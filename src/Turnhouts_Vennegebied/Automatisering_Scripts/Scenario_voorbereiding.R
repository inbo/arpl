library(sf)
library(terra)
library(tidyverse)
library(data.table)
library(exactextractr)
library(here)

# Start de totale klok voor alle scenario's samen
totale_tijd_start <- Sys.time()
epsilon <- 0.00001

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
# CONFIGURATIE: GEEF HIER AL JE SCENARIO'S OP
# =========================================================================
# Drop hier alle bestanden die je achter elkaar wilt verrasteren.
scenario_bestanden <- c(
  "Scenario/streefbeelden_incl_allocatie_bosbehoudss_ss31fix_tvg_vrij_20260107_r2.gpkg",
  "Scenario/streefbeelden_incl_allocatie_bosbehoudss_ss31fix_kij_vrij_20260107_r2.gpkg"
)

# =========================================================================
# DEEL I: EÉNMALIGE STARTBLOK (BUITEN DE LOOP - REKENTIJD BESPAREN)
# =========================================================================
message("=== STAP 1: GLOBALE GEOMETRIEËN EN REFERENCE-RASTERS LADEN ===")

master_grid_vlaanderen <- rast(here("data/input/Raster_Vlaanderen/Vlaanderen_MasterGrid_10m.tif"))
target_crs             <- terra::crs(master_grid_vlaanderen)

area_shape_proj        <- project(vect(here("data/input/Turnhouts_Vennegebied.shp")), target_crs)
area_shape_sf          <- st_as_sf(area_shape_proj)

message("-> De gigantische BWK-basiskaart éénmalig inladen en valideren...")
bwk_sf_global  <- st_read(here("data/references/Shapefiles/BwkHab.shp"), quiet = TRUE) %>% 
  st_transform(31370) %>% 
  st_make_valid()

message("-> Regionale 50km-buffer basistabel éénmalig inladen...")
tabel_bwk_lokaal_global <- readRDS(here("data/input/Raster_Vlaanderen/BWK_Tabel_Turnhout_Plus_50km.rds"))
setDT(tabel_bwk_lokaal_global)

# =========================================================================
# DEEL II: DE MULTI-SCENARIO WASSTRAAT (FOR-LOOP)
# =========================================================================
message("\n=== STAP 2: START MULTI-SCENARIO VERRASTERINGSSCHEMA ===")

for(scen_path in scenario_bestanden) {
  
  # Check of het bestand fysiek bestaat om crashes te voorkomen
  if(!file.exists(scen_path)) {
    message(paste("⚠️ Waarschuwing: Bestand", scen_path, "niet gevonden. Overslaan..."))
    next
  }
  
  # Automatische naamontleding voor de export (.rds)
  bestandsnaam <- basename(scen_path)
  basis_naam   <- sub("\\.(gpkg|shp)$", "", bestandsnaam)
  
  message(paste("\n▶ START VERRASTERING VOOR:", toupper(basis_naam)))
  scen_tijd_start <- Sys.time()
  
  # 1. Scenario inladen en direct filteren op gevulde informatie
  scen_sf <- st_read(scen_path, quiet = TRUE) %>% 
    st_transform(31370) %>% 
    st_make_valid() %>% 
    filter(!is.na(nsb1) & nsb1 != "" & nsb1 != " " & nsb1 != "geen vegetatie NSB")
  
  if(nrow(scen_sf) == 0) {
    message("⚠️ Dit scenario bevat geen geldige vegetatiepolygonen. Overslaan...")
    next
  }
  
  # 2. Vector overlay met de reeds ingeladen globale bwk_sf
  sf_overlap  <- scen_sf
  sf_rest_bwk <- st_difference(bwk_sf_global, st_union(sf_overlap))
  sf_rest_bwk <- st_intersection(sf_rest_bwk, area_shape_sf)
  
  # 3. Herstructureren naar Long Format (Deel A: Scenario)
  sf_overlap <- sf_overlap %>% mutate(scen_poly_id = row_number())
  
  df_scen_codes <- sf_overlap %>% st_drop_geometry() %>% 
    select(scen_poly_id, matches("^nsb[1-8]$")) %>% 
    pivot_longer(cols = matches("^nsb[1-8]$"), names_to = "source_col", values_to = "CODE") %>% 
    mutate(idx = as.numeric(sub("nsb", "", source_col)))
  
  df_scen_opp <- sf_overlap %>% st_drop_geometry() %>% 
    select(scen_poly_id, matches("^pnsb[1-8]$")) %>% 
    pivot_longer(cols = matches("^pnsb[1-8]$"), names_to = "opp_col", values_to = "PCT_VAL") %>% 
    mutate(idx = as.numeric(sub("pnsb", "", opp_col))) %>% 
    select(scen_poly_id, idx, PCT_VAL)
  
  df_overlap_long <- df_scen_codes %>% 
    left_join(df_scen_opp, by = c("scen_poly_id", "idx")) %>% 
    mutate(BWK_FRAC = coalesce(as.numeric(PCT_VAL), 0) / 100) %>% 
    filter(!is.na(CODE) & CODE != "" & CODE != " " & CODE != "NA") %>% 
    mutate(CODE = tolower(trimws(as.character(CODE)))) %>% 
    group_by(scen_poly_id, CODE) %>% 
    summarise(BWK_FRAC = sum(BWK_FRAC, na.rm = TRUE), .groups = "drop")
  
  vec_scen_finaal_sf <- sf_overlap %>% select(scen_poly_id) %>% inner_join(df_overlap_long, by = "scen_poly_id") %>% 
    mutate(poly_id = scen_poly_id) %>% select(poly_id, CODE, BWK_FRAC)
  
  # 3. Herstructureren naar Long Format (Deel B: Rest-BWK)
  sf_rest_bwk <- sf_rest_bwk %>% mutate(rest_poly_id = row_number())
  
  for(i in 1:4) {
    phab_col <- paste0("PHAB", i)
    target_opp_col <- paste0("OPP_HAB", i)
    sf_rest_bwk[[target_opp_col]] <- if(phab_col %in% names(sf_rest_bwk)) (replace_na(as.numeric(sf_rest_bwk[[phab_col]]), 0) / 100) else 0
  }
  
  exist_eenh <- intersect(c("EENH1", "EENH2", "EENH3", "EENH4"), names(sf_rest_bwk))
  if(length(exist_eenh) > 0) {
    eenh_matrix <- as.matrix(st_drop_geometry(sf_rest_bwk[, exist_eenh]))
    sf_rest_bwk$n_filled <- rowSums(!is.na(eenh_matrix) & eenh_matrix != "" & eenh_matrix != " ")
    for(i in 1:4) sf_rest_bwk[[paste0("OPP_EENH", i)]] <- get_eenh_pct(sf_rest_bwk$n_filled, i)
  }
  
  exist_hab <- intersect(paste0("HAB", 1:5), names(sf_rest_bwk))
  
  df_eenh_finaal <- sf_rest_bwk %>% st_drop_geometry() %>% select(rest_poly_id, any_of(exist_eenh)) %>% 
    pivot_longer(cols = any_of(exist_eenh), names_to = "source_col", values_to = "CODE") %>% mutate(idx = as.numeric(sub("EENH", "", source_col))) %>% 
    left_join(sf_rest_bwk %>% st_drop_geometry() %>% select(rest_poly_id, any_of(paste0("OPP_EENH", 1:4))) %>% pivot_longer(cols = any_of(paste0("OPP_EENH", 1:4)), names_to = "opp_col", values_to = "OPP_VAL") %>% mutate(idx = as.numeric(sub("OPP_EENH", "", opp_col))) %>% select(rest_poly_id, idx, OPP_VAL), by = c("rest_poly_id", "idx")) %>% 
    mutate(BWK_FRAC = case_when(idx <= 4 ~ coalesce(as.numeric(OPP_VAL), 0.00), idx >= 5 ~ epsilon, TRUE ~ 0.00))
  
  df_hab_finaal <- sf_rest_bwk %>% st_drop_geometry() %>% select(rest_poly_id, any_of(exist_hab)) %>% 
    pivot_longer(cols = any_of(exist_hab), names_to = "source_col", values_to = "CODE") %>% mutate(idx = as.numeric(sub("HAB", "", source_col))) %>% 
    left_join(sf_rest_bwk %>% st_drop_geometry() %>% select(rest_poly_id, any_of(paste0("OPP_HAB", 1:4))) %>% pivot_longer(cols = any_of(paste0("OPP_HAB", 1:4)), names_to = "opp_col", values_to = "OPP_VAL") %>% mutate(idx = as.numeric(sub("OPP_HAB", "", opp_col))) %>% select(rest_poly_id, idx, OPP_VAL), by = c("rest_poly_id", "idx")) %>% 
    mutate(BWK_FRAC = case_when(idx <= 4 ~ coalesce(as.numeric(OPP_VAL), 0.00), idx >= 5 ~ epsilon, TRUE ~ 0.00))
  
  df_rest_clean <- bind_rows(df_eenh_finaal, df_hab_finaal) %>% 
    filter(!is.na(CODE) & CODE != "" & CODE != " " & CODE != "NA") %>% 
    mutate(CODE = tolower(trimws(as.character(CODE)))) %>% 
    group_by(rest_poly_id, CODE) %>% summarise(BWK_FRAC = max(BWK_FRAC, na.rm = TRUE), .groups = "drop")
  
  max_scen_id <- max(sf_overlap$scen_poly_id, default = 0)
  vec_bwk_finaal_sf <- sf_rest_bwk %>% select(rest_poly_id) %>% inner_join(df_rest_clean, by = "rest_poly_id") %>% 
    mutate(poly_id = rest_poly_id + max_scen_id) %>% select(poly_id, CODE, BWK_FRAC)
  
  # Samenvoegen tot gecombineerd target
  vec_sf_long <- bind_rows(vec_scen_finaal_sf, vec_bwk_finaal_sf)
  
  # 4. Kogelvrije geometrie-reiniging per scenario run
  vec_sf_long <- vec_sf_long %>% 
    st_make_valid() %>% 
    filter(!st_is_empty(.))
  
  vlakken_index <- st_dimension(vec_sf_long) == 2
  vec_sf_long   <- vec_sf_long[vlakken_index, ]
  vec_sf_long   <- st_cast(vec_sf_long, "MULTIPOLYGON")
  
  vec_sf_long <- vec_sf_long %>% filter(as.numeric(st_area(.)) > 0.001)
  
  df_long_clean <- st_drop_geometry(vec_sf_long)
  unieke_codes  <- sort(unique(df_long_clean$CODE))
  n_codes       <- length(unieke_codes)
  
  # 5. Extractie wasstraat
  lijst_code_tabellen <- list()
  for(i in 1:n_codes) {
    h_code <- unieke_codes[i]
    sub_sf <- vec_sf_long[vec_sf_long$CODE == h_code, ]
    
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
  
  tabel_TV_exact <- rbindlist(lijst_code_tabellen)
  tabel_TV_exact[, BWK_FRAC := pmin(BWK_FRAC, 1.00)]
  
  # 6. Integratie met de reeds ingeladen globale 50km tabel
  target_cellen_tv <- unique(tabel_TV_exact$cel_id)
  tabel_bwk_gezuiverd <- tabel_bwk_lokaal_global[!cel_id %in% target_cellen_tv]
  tabel_finaal_scenario <- rbindlist(list(tabel_bwk_gezuiverd, tabel_TV_exact), use.names = TRUE)
  
  # 7. Unieke bestandsnaam wegschrijven op basis van de gpkg-naam
  export_naam <- paste0("TV_Scenario_", basis_naam, ".rds")
  saveRDS(tabel_finaal_scenario, export_naam)
  
  scen_duur <- difftime(Sys.time(), scen_tijd_start, units="mins")
  message(paste("✓ Scenario", basis_naam, "is succesvol verrasterd en opgeslagen als:", export_naam, "in", round(scen_duur, 2), "minuten."))
  
  # Tussentijdse schoonmaak van tijdelijke scenario-objecten
  rm(scen_sf, sf_overlap, sf_rest_bwk, df_scen_codes, df_scen_opp, df_overlap_long, 
     vec_scen_finaal_sf, vec_bwk_finaal_sf, vec_sf_long, df_long_clean, tabel_TV_exact, 
     tabel_bwk_gezuiverd, tabel_finaal_scenario, lijst_code_tabellen)
  gc()
}

totale_duur <- difftime(Sys.time(), totale_tijd_start, units="mins")
message(paste("\n🏁 FINISH! Alle scenario's zijn verwerkt. Totale run-tijd:", round(totale_duur, 2), "minuten."))