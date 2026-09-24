library(sf)
library(terra)
library(tidyverse)
library(data.table)
library(exactextractr)
library(here)

# =========================================================================
# CONFIGURATIE & CONFIG-FUNCTIES
# =========================================================================
totale_tijd_start <- Sys.time()
epsilon <- 0.00001 # Klein getal voor 'aanwezigheid' (indices > 4)

# Percentage-toewijzing voor de eerste 4 EENH-kolommen
get_eenh_pct <- function(n, index) {
  case_when(
    n == 1 & index == 1 ~ 1.00,
    n == 2 & index == 1 ~ 0.70, n == 2 & index == 2 ~ 0.30,
    n == 3 & index == 1 ~ 0.60, n == 3 & index == 2 ~ 0.20, n == 3 & index == 3 ~ 0.20,
    n >= 4 & index == 1 ~ 0.60, n >= 4 & index == 2 ~ 0.20, n >= 4 & index == 3 ~ 0.10, n >= 4 & index == 4 ~ 0.10,
    TRUE ~ 0.00
  )
}

# =========================================================================
# STAP 1: INPUT BESTANDEN LADEN (VOLLEDIG VLAANDEREN - GEEN CLIPPING)
# =========================================================================
message("=== STAP 1: GEOMETRIEËN EN REFERENCE-RASTER LADEN ===")

# 1. Master grid (10m x 10m Vlaanderen)
master_grid_vlaanderen <- rast(here("data/input/Raster_Vlaanderen/Vlaanderen_MasterGrid_10m.tif"))

# 2. De originele BWK Basiskaart (BWK_2025.shp)
message("-> BWK basiskaart voor héél Vlaanderen inladen en valideren...")
bwk_sf <- st_read(here("data/references/Shapefiles/BWK_2025.shp"), quiet = TRUE) %>% 
  st_transform(31370) %>% 
  st_make_valid()

# =========================================================================
# STAP 2: HERSTRUCTUREREN NAAR LONG FORMAT (8 EENH EN 5 HAB CODES)
# =========================================================================
message("=== STAP 2: HET HERSTRUCTUREREN VAN EENH (1-8) EN HAB (1-5) CODES ===")

bwk_sf <- bwk_sf %>% mutate(poly_id = row_number())

# 2A. Bereken HAB oppervlakte fracties (PHAB1 - PHAB4, PHAB5 bestaat niet in shape -> 0/epsilon)
for(i in 1:5) {
  phab_col       <- paste0("PHAB", i)
  target_opp_col <- paste0("OPP_HAB", i)
  bwk_sf[[target_opp_col]] <- if(phab_col %in% names(bwk_sf)) {
    (replace_na(as.numeric(bwk_sf[[phab_col]]), 0) / 100)
  } else 0
}

# 2B. Bereken EENH oppervlakte fracties (EENH1 - EENH4 gebruiken verdeelsleutel)
exist_eenh_first4 <- intersect(c("EENH1", "EENH2", "EENH3", "EENH4"), names(bwk_sf))
if(length(exist_eenh_first4) > 0) {
  eenh_matrix <- as.matrix(st_drop_geometry(bwk_sf[, exist_eenh_first4]))
  bwk_sf$n_filled <- rowSums(!is.na(eenh_matrix) & eenh_matrix != "" & eenh_matrix != " ")
  for(i in 1:4) bwk_sf[[paste0("OPP_EENH", i)]] <- get_eenh_pct(bwk_sf$n_filled, i)
}

# Zoek alle beschikbare EENH (1-8) en HAB (1-5) kolommen in het shapefile
exist_eenh_all <- intersect(paste0("EENH", 1:8), names(bwk_sf))
exist_hab_all  <- intersect(paste0("HAB", 1:5), names(bwk_sf))

# Pivot EENH codes (1 t/m 8) naar long format
df_eenh_finaal <- bwk_sf %>% 
  st_drop_geometry() %>% 
  select(poly_id, any_of(exist_eenh_all)) %>% 
  pivot_longer(cols = any_of(exist_eenh_all), names_to = "source_col", values_to = "CODE") %>% 
  mutate(idx = as.numeric(sub("EENH", "", source_col))) %>% 
  left_join(
    bwk_sf %>% st_drop_geometry() %>% 
      select(poly_id, any_of(paste0("OPP_EENH", 1:4))) %>% 
      pivot_longer(cols = any_of(paste0("OPP_EENH", 1:4)), names_to = "opp_col", values_to = "OPP_VAL") %>% 
      mutate(idx = as.numeric(sub("OPP_EENH", "", opp_col))) %>% 
      select(poly_id, idx, OPP_VAL), 
    by = c("poly_id", "idx")
  ) %>% 
  # EENH 1-4 krijgen de berekende oppervlaktes, EENH 5-8 krijgen epsilon
  mutate(BWK_FRAC = case_when(
    idx <= 4 ~ coalesce(as.numeric(OPP_VAL), 0.00), 
    idx >= 5 ~ epsilon, 
    TRUE ~ 0.00
  ))

# Pivot HAB codes (1 t/m 5) naar long format
df_hab_finaal <- bwk_sf %>% 
  st_drop_geometry() %>% 
  select(poly_id, any_of(exist_hab_all)) %>% 
  pivot_longer(cols = any_of(exist_hab_all), names_to = "source_col", values_to = "CODE") %>% 
  mutate(idx = as.numeric(sub("HAB", "", source_col))) %>% 
  left_join(
    bwk_sf %>% st_drop_geometry() %>% 
      select(poly_id, any_of(paste0("OPP_HAB", 1:4))) %>% 
      pivot_longer(cols = any_of(paste0("OPP_HAB", 1:4)), names_to = "opp_col", values_to = "OPP_VAL") %>% 
      mutate(idx = as.numeric(sub("OPP_HAB", "", opp_col))) %>% 
      select(poly_id, idx, OPP_VAL), 
    by = c("poly_id", "idx")
  ) %>% 
  # HAB 1-4 krijgen PHAB% / 100, HAB 5 krijgt epsilon
  mutate(BWK_FRAC = case_when(
    idx <= 4 ~ coalesce(as.numeric(OPP_VAL), 0.00), 
    idx >= 5 ~ epsilon, 
    TRUE ~ 0.00
  ))

# Combineer, filter ongeldige codes en neem de maximale fractie per polygoon & code
df_bwk_clean <- bind_rows(df_eenh_finaal, df_hab_finaal) %>% 
  filter(!is.na(CODE) & CODE != "" & CODE != " " & CODE != "NA") %>% 
  mutate(CODE = tolower(trimws(as.character(CODE)))) %>% 
  group_by(poly_id, CODE) %>% 
  summarise(BWK_FRAC = max(BWK_FRAC, na.rm = TRUE), .groups = "drop")

# Koppel opgeschoonde attributes terug aan sf object
vec_bwk_sf <- bwk_sf %>% 
  select(poly_id) %>% 
  inner_join(df_bwk_clean, by = "poly_id") %>% 
  st_make_valid() %>% 
  filter(!st_is_empty(.))

vlakken_index <- st_dimension(vec_bwk_sf) == 2
vec_bwk_sf   <- vec_bwk_sf[vlakken_index, ]
vec_bwk_sf   <- st_cast(vec_bwk_sf, "MULTIPOLYGON")
vec_bwk_sf   <- vec_bwk_sf %>% filter(as.numeric(st_area(.)) > 0.001)

# Schoon het RAM-geheugen op voor de zware exact_extract loop
rm(bwk_sf, df_eenh_finaal, df_hab_finaal, df_bwk_clean)
gc()

# =========================================================================
# STAP 3: RASTER EXTRACTIE WASSTRAAT (EXACT_EXTRACT VLAANDEREN-BREED)
# =========================================================================
message("=== STAP 3: START VERRASTERING EN DEKKINGSFRACTIE BEREKENING ===")

unieke_codes <- sort(unique(vec_bwk_sf$CODE))
n_codes      <- length(unieke_codes)

lijst_code_tabellen <- list()

for(i in 1:n_codes) {
  h_code <- unieke_codes[i]
  sub_sf <- vec_bwk_sf[vec_bwk_sf$CODE == h_code, ]
  
  if(nrow(sub_sf) > 0) {
    extractie <- exact_extract(master_grid_vlaanderen[[1]], sub_sf, include_cell = TRUE, progress = FALSE)
    df_extract <- if(is.data.frame(extractie)) as.data.table(extractie)[, polygon_id := 1] else rbindlist(extractie, idcol = "polygon_id")
    
    if(nrow(df_extract) > 0) {
      df_extract[, poly_id := sub_sf$poly_id[polygon_id]]
      df_extract[, BWK_FRAC_poly := sub_sf$BWK_FRAC[polygon_id]]
      df_extract[, BWK_FRAC := coverage_fraction * BWK_FRAC_poly]
      
      df_cel_som <- df_extract[, .(BWK_FRAC = sum(BWK_FRAC, na.rm = TRUE)), by = .(cell)]
      if(nrow(df_cel_som) > 0) {
        setnames(df_cel_som, "cell", "cel_id")
        df_cel_som[, CODE := h_code]
        lijst_code_tabellen[[h_code]] <- df_cel_som[BWK_FRAC > 0]
      }
    }
  }
}

# =========================================================================
# STAP 4: SAMENVOEGEN EN WEGSCHRIJVEN ALS RDS
# =========================================================================
message("=== STAP 4: FINALE TABEL OPBOUWEN EN EXPORTEREN ===")

tabel_bwk_finaal <- rbindlist(lijst_code_tabellen)
tabel_bwk_finaal[, BWK_FRAC := pmin(BWK_FRAC, 1.00)]

export_pad <- here("data/input/Raster_Vlaanderen/Test_BWK_TidyTabel_Smal_Vlaanderen_2025.rds")
saveRDS(tabel_bwk_finaal, export_pad)

totale_duur <- difftime(Sys.time(), totale_tijd_start, units="mins")
message(paste("\n🏁 FINISH! Volledige Vlaanderen BWK rastertabel opgeslagen als:", export_pad, "in", round(totale_duur, 2), "minuten."))
