library(here)
# CRUCIALE FIX: Dwing R Markdown om te werken vanaf de hoofdmap (arpl/)
# Hierdoor werken relatieve paden met hier() overal hetzelfde.

library(knitr)
library(tidyverse)
library(sf)
library(terra)
library(readxl)
library(tidyterra)
library(leaflet)
library(kableExtra)
library(data.table)

conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::select)
conflicted::conflicts_prefer(dplyr::first)
conflicted::conflicts_prefer(terra::intersect)
conflicted::conflicts_prefer(terra::any)

# Exacte hectare-berekening behouden (cellSize voor 100% precisie)
calc_ha_exact <- function(r) {
  if(is.null(r)) return(0)
  if(all(is.na(terra::values(r, mat=FALSE)))) return(0)
  area_raster <- r * terra::cellSize(r, unit = "ha")
  val <- terra::global(area_raster, "sum", na.rm = TRUE)[[1]]
  return(as.numeric(val))
}

# Geoptimaliseerde get_stats helper
get_stats <- function(cat_id, label) {
  cid <- cat_id
  target_mask <- r_status == cid
  f_pix <- terra::freq(target_mask)
  n_pix <- if(nrow(f_pix) > 0) sum(f_pix$count[f_pix$value == 1], na.rm=TRUE) else 0
  area_ha <- (n_pix * 100) / 10000
  
  cl_src <- if(cid %in% c(1, 2)) cl_id_max else cl_id_opp
  cl_zone <- terra::mask(cl_src, target_mask)
  f_cl <- terra::freq(cl_zone)
  n_cl <- if(!is.null(f_cl) && nrow(f_cl) > 0) nrow(f_cl) else 0
  
  return(data.frame(Type = label, Clusters = n_cl, Oppervlakte_ha = round(area_ha, 2)))
}

maak_kaart_laag <- function(res_list) {
  if(is.null(res_list)) return(NULL)
  kernen <- terra::ifel(!is.na(res_list$kern), 1, NA)
  if(!is.null(res_list$bouw)) {
    bouw <- terra::ifel(!is.na(res_list$bouw), 2, NA)
    finaal <- terra::cover(kernen, bouw)
  } else {
    finaal = kernen
  }
  return(finaal)
}

cluster_filter_compleet <- function(masker, opp_laag, drempel_m2, dist_m, werkelijk = FALSE) {
  # 1. Controleer of het raster leeg is
  if (terra::global(is.na(masker), "sum")[[1]] == terra::ncell(masker)) {
    return(list(raster = masker * NA, clusters = masker * NA))
  }
  
  # ----------------------------------------------------------------------------
  # STAP 1: VEILIGE & SNELLER NETWERKVORMING ZONDER RAM-CRASH
  # ----------------------------------------------------------------------------
  if (dist_m > 0) {
    # Maak binaire kaart (1 = biotoop, NA = rest)
    r_binair <- terra::ifel(!is.na(masker) & masker > 0, 1, NA)
    
    # Buffer de BINAIR kaart (dit kost vrijwel geen geheugen!)
    # Volle afstand dist_m zorgt dat plukjes binnen dist_m gegarandeerd samensmelten
    r_buffered <- terra::buffer(r_binair, width = dist_m / 2)
    
    # Maak unieke Netwerk-ID's op de gebufferde zones
    cl_network <- terra::patches(r_buffered, directions = 4, zeroAsNA = TRUE)
    
    # Snijd de Netwerk-ID's direct terug naar waar de originele biotooppixels liggen
    cl_biotoop_only <- terra::mask(cl_network, masker)
  } else {
    cl_network <- terra::patches(masker, directions = 8, zeroAsNA = TRUE)
    cl_biotoop_only <- cl_network
  }
  
  # ----------------------------------------------------------------------------
  # STAP 2: OPPERVLAKTE-OPTELSOM PER GEKOPPELD NETWERK
  # ----------------------------------------------------------------------------
  if(werkelijk) {
    stats_df <- terra::zonal(opp_laag, cl_biotoop_only, fun = "sum", na.rm = TRUE)
    colnames(stats_df) <- c("ID", "Waarde")
    stats_df$Area_m2 <- stats_df$Waarde * 100 
  } else {
    f <- terra::freq(cl_biotoop_only)
    stats_df <- data.frame(ID = f$value, Waarde = f$count)
    stats_df$Area_m2 <- stats_df$Waarde * 100 
  }
  
  stats_df <- stats_df[!is.na(stats_df$ID), ]
  if(nrow(stats_df) == 0) return(list(raster = masker * NA, clusters = masker * NA))
  
  # ----------------------------------------------------------------------------
  # STAP 3: FILTEREN OP TOTALE NETWERK-OPPERVLAKTE >= DREMPEL
  # ----------------------------------------------------------------------------
  voldoet_ids <- stats_df$ID[stats_df$Area_m2 >= drempel_m2]
  if(length(voldoet_ids) == 0) return(list(raster = masker * NA, clusters = masker * NA))
  
  # Behaal alleen de winnende netwerk-ID's
  masker_binair <- cl_biotoop_only %in% voldoet_ids
  final_network_mask <- terra::ifel(masker_binair == 1, 1, NA)
  
  r_finaal  <- terra::mask(masker, final_network_mask)
  cl_finaal <- terra::mask(cl_biotoop_only, r_finaal) 
  
  return(list(raster = r_finaal, clusters = cl_finaal))
}

terraOptions(
  memfrac = 0.8,        # Dwing terra om tot max. 80% van het RAM-geheugen te gebruiken
  tempdir = tempdir(),  # Geef toestemming voor automatische disk-swapping bij zware rasters
  verbose = FALSE
)

df <- read_excel(here::here("data/input/Excel_files/Soorten_bwk_afstanden.xlsx"))
soort <- "venwitsnuitlibel"

# Scenario pad en naam bepalen

# --- DYNAMISCHE SCENARIO PARAMETER CHECK ---
if (!exists("params") || is.null(params$scenario_rds_path)) {
  scenario_rds_path <- "data/input/Scenario_rds/TV_Scenario_BWK_2025.rds"
} else {
  scenario_rds_path <- params$scenario_rds_path
}

p_raw <- gsub("^([.][.]/)+", "", scenario_rds_path)
scenario_path <- here::here(p_raw)


if (!file.exists(scenario_path)) {
  stop(paste("❌ FOUT: Scenario RDS bestand NIET gevonden op:", scenario_path))
}

scen_volledig <- basename(scenario_path)
scenario_naam <- gsub("^TV_Scenario_|^Scenario_|.rds$", "", scen_volledig)

message(paste("Verwerken van soort:", soort, "binnen scenario:", scenario_naam))

resultaat <- df %>%
  filter(tolower(trimws(Soort)) == soort) %>%
  select(Type, MinOpp_ha, AfstandBiotopen_m, Dispersiecap_m)

# Variabelen definiëren
buffer_m       <- resultaat$Dispersiecap_m[1]
straal_water_m <- 500  # Harde afsnij-afstand van land rond het water (500m)

print(resultaat)
rm(df, resultaat)

area_shape  <- vect(here("data/input/Turnhouts_Vennegebied.shp"))
master_grid <- rast(here("data/input/Raster_Vlaanderen/Vlaanderen_MasterGrid_10m.tif"))[[1]]

df_namen_sleutel <- read_csv(here("data/input/Excel_files/BWK_Laag_Namen_2025.csv"), show_col_types = FALSE)
gouden_namenlijst <- tolower(trimws(df_namen_sleutel$Laagnaam))

area_shape_proj <- project(area_shape, crs(master_grid))
area_buffer_fix <- buffer(area_shape_proj, width = buffer_m)

message("-> Vertaalraster voor globale/lokale cellen opbouwen via snelle MASK methode...")
id_raster_TV <- crop(master_grid, area_buffer_fix, snap = "near")

globale_id_raster <- master_grid
values(globale_id_raster) <- 1:ncell(globale_id_raster)

id_raster_TV_globale_values <- crop(globale_id_raster, area_buffer_fix, snap = "near")
id_raster_TV_masked <- mask(id_raster_TV_globale_values, area_buffer_fix)

message("-> Vertaaltabel bliksemsnel opbouwen via C++ dataframe extractie...")

df_extractie <- as.data.frame(id_raster_TV_masked, cells = TRUE)
vertaal_df <- as.data.table(df_extractie)
setnames(vertaal_df, c(1, 2), c("lokale_id", "globale_id"))

vertaal_df <- vertaal_df[!is.na(globale_id)]
studiegebied_globale_ids <- unique(vertaal_df$globale_id)

values(id_raster_TV) <- NA
template_TV <- terra::rasterize(area_buffer_fix, id_raster_TV, field = 1, background = 0)

cat("Gecorrigeerd aantal pixels in template_TV: ", sum(terra::values(template_TV) == 1, na.rm=TRUE), "\n")

grens_web <- sf::st_as_sf(terra::project(area_shape, "EPSG:4326"))

rm(globale_id_raster, id_raster_TV_globale_values, id_raster_TV_masked, df_extractie)
gc()

# 1. LAAD DE BRON-CROSSWALK/DICTIONARY IN
df_nieuw <- read_csv(here("data/input/Excel_files/Resultaten_Totaal_Samengevoegd.csv"), show_col_types = FALSE)

resultaten_gegroepeerd <- df_nieuw %>%
  mutate(Soort_clean = tolower(trimws(Soort))) %>%
  filter(Soort_clean == soort) %>%
  group_by(Type) %>%
  nest(Data = c(Code, Match))

# Lees de scenario-RDS in
tabel_vlaanderen <- readRDS(scenario_path)
setDT(tabel_vlaanderen)
tabel_vlaanderen[, CODE := tolower(trimws(CODE))]

lijst_matches      <- list()
lijst_oppervlaktes <- list()

for(i in 1:nrow(resultaten_gegroepeerd)) {
  h_type <- resultaten_gegroepeerd$Type[i]
  h_data <- resultaten_gegroepeerd$Data[[i]]
  
  exact_codes <- tolower(trimws(h_data$Code[h_data$Match == "exact"]))
  bevat_codes <- tolower(trimws(h_data$Code[h_data$Match == "bevat"]))
  
  if(length(bevat_codes) > 0) {
    bevat_codes_escaped <- gsub("([\\.\\^\\$\\*\\+\\?\\(\\)\\[\\{\\\\\\|])", "\\\\\\1", bevat_codes)
    bevat_codes_anchored <- paste0("^", bevat_codes_escaped)
    regex_term <- paste0(bevat_codes_anchored, collapse = "|")
    
    tabel_gefilterd <- tabel_vlaanderen[CODE %in% exact_codes | grepl(regex_term, CODE)]
  } else {
    tabel_gefilterd <- tabel_vlaanderen[CODE %in% exact_codes]
  }
  
  tabel_TV <- tabel_gefilterd[cel_id %in% studiegebied_globale_ids]
  tabel_TV_unique <- unique(tabel_TV, by = c("cel_id", "CODE"))

  tabel_cel_som <- tabel_TV_unique[, .(Oppervlakte = pmin(sum(BWK_FRAC, na.rm = TRUE), 1.0)), by = .(cel_id)]
  
  r_match_type <- id_raster_TV * NA
  r_opp_type   <- id_raster_TV * NA
  
  if(nrow(tabel_cel_som) > 0) {
    tabel_cel_som[, Oppervlakte := pmin(Oppervlakte, 1)]
    tabel_cel_som[, Match := ifelse(Oppervlakte >= 0.01, 1, 0)]
    
    setnames(tabel_cel_som, "cel_id", "globale_id")
    tabel_finaal_mapping <- merge(tabel_cel_som, vertaal_df, by = "globale_id", all.x = TRUE)
    tabel_finaal_mapping <- tabel_finaal_mapping[!is.na(lokale_id)]
    
    if(nrow(tabel_finaal_mapping) > 0) {
      tabel_matches_clean <- tabel_finaal_mapping[Match == 1]
      if(nrow(tabel_matches_clean) > 0) {
        r_match_type[tabel_matches_clean$lokale_id] <- tabel_matches_clean$Match
      }
      r_opp_type[tabel_finaal_mapping$lokale_id] <- tabel_finaal_mapping$Oppervlakte
      
      lijst_matches[[h_type]]      <- r_match_type
      lijst_oppervlaktes[[h_type]] <- r_opp_type
    }
  }
}

waterbiotoop_bwk1_max        <- lijst_matches[["waterbiotoop_bwk1"]]
waterbiotoop_bwk1_opp        <- lijst_oppervlaktes[["waterbiotoop_bwk1"]]
waterbiotoop_notbwk_max      <- lijst_matches[["waterbiotoop_notbwk"]]
waterbiotoop_notbwk_opp      <- lijst_oppervlaktes[["waterbiotoop_notbwk"]]
landbiotoop_bwk_max          <- lijst_matches[["landbiotoop_bwk"]]
landbiotoop_bwk_opp          <- lijst_oppervlaktes[["landbiotoop_bwk"]]
landbiotoop_bwk_ruim_max     <- lijst_matches[["landbiotoop_bwk_ruim"]]
landbiotoop_bwk_ruim_opp     <- lijst_oppervlaktes[["landbiotoop_bwk_ruim"]]
landbiotoop_bwk_ruim_bos_max <- lijst_matches[["landbiotoop_bwk_ruim_bos"]]
landbiotoop_bwk_ruim_bos_opp <- lijst_oppervlaktes[["landbiotoop_bwk_ruim_bos"]]

rm(tabel_vlaanderen, vertaal_df, lijst_matches, lijst_oppervlaktes)
gc()

# ==============================================================================
# BEREIDING WATERBIOTOOP VENWITSNUITLIBEL (Minimaal 100 m² - PARALLEL)
# ==============================================================================
message("-> Waterbiotoop filteren (inclusief uitsluitingen)...")

water_max <- (waterbiotoop_bwk1_max == 1) & (waterbiotoop_notbwk_max == 0 | is.na(waterbiotoop_notbwk_max))
water_max <- terra::ifel(water_max == 1, 1, NA)

water_opp_raw <- terra::mask(waterbiotoop_bwk1_opp, water_max)
water_opp     <- terra::ifel(!is.na(water_opp_raw) & water_opp_raw > 0, water_opp_raw, NA)

message("-> Waterbiotoop patches groeperen (50m fuzzy) en filteren op min. 100 m²...")

# Spoor A: MAX
waterbiotoop_clusters_max <- cluster_filter_compleet(
  masker     = water_max,
  opp_laag   = water_opp,
  drempel_m2 = 100,   
  dist_m     = 50,       
  werkelijk  = FALSE      
)

# Spoor B: OPP (Parallel & Autonoom)
r_binair_water_opp <- terra::ifel(!is.na(water_opp) & water_opp > 0, 1, NA)
waterbiotoop_clusters_opp <- cluster_filter_compleet(
  masker     = r_binair_water_opp,
  opp_laag   = water_opp,
  drempel_m2 = 100,   
  dist_m     = 50,       
  werkelijk  = TRUE      
)

rm(water_max, water_opp_raw, r_binair_water_opp)
gc()

# ==============================================================================
# VOORBEREIDING LANDBIOTOOP NABIJ (Minimaal 5 ha - PARALLEL)
# ==============================================================================
message("-> Landbiotoop nabij patches groeperen (20m) en filteren op min. 5 ha (50.000 m²)...")

# Spoor A: MAX
land_nabij_clusters_max <- cluster_filter_compleet(
  masker     = landbiotoop_bwk_max,
  opp_laag   = landbiotoop_bwk_opp,
  drempel_m2 = 50000,   
  dist_m     = 20,       
  werkelijk  = FALSE      
)
venwitsnuitlibel_land_nabij_max <- land_nabij_clusters_max$raster

# Spoor B: OPP (Parallel & Autonoom)
r_binair_land_nabij_opp <- terra::ifel(!is.na(landbiotoop_bwk_opp) & landbiotoop_bwk_opp > 0, 1, NA)
land_nabij_clusters_opp <- cluster_filter_compleet(
  masker     = r_binair_land_nabij_opp,
  opp_laag   = landbiotoop_bwk_opp,
  drempel_m2 = 50000,   
  dist_m     = 20,       
  werkelijk  = TRUE      
)
venwitsnuitlibel_land_nabij_opp <- land_nabij_clusters_opp$raster

rm(land_nabij_clusters_max, land_nabij_clusters_opp, r_binair_land_nabij_opp)
gc()

# ==============================================================================
# BEREIDING RUIME OMGEVING & ZUIVERE PERCENTAGE-CHECK (>60% BOS - PARALLEL)
# ==============================================================================
message("-> Landbiotoop ruime omgeving patches groeperen (50m) en filteren op min. 20 ha...")

# --- SPOOR A: MAX ---
land_ruim_clusters_max <- cluster_filter_compleet(
  masker     = landbiotoop_bwk_ruim_max,
  opp_laag   = landbiotoop_bwk_ruim_opp,
  drempel_m2 = 200000,   
  dist_m     = 50,       
  werkelijk  = FALSE      
)
r_patch_ids_max <- land_ruim_clusters_max$clusters
venwitsnuitlibel_land_ruim_max <- template_TV * NA

if (!all(is.na(terra::values(r_patch_ids_max, mat=FALSE)))) {
  area_ha_raster_max <- r_patch_ids_max * terra::cellSize(r_patch_ids_max, unit = "ha")
  totale_opp_per_patch_max <- terra::zonal(area_ha_raster_max, r_patch_ids_max, fun = "sum", na.rm = TRUE)
  colnames(totale_opp_per_patch_max) <- c("Patch_ID", "Totaal_ha")
  
  bos_opp_raster_max <- terra::mask(landbiotoop_bwk_ruim_bos_max, r_patch_ids_max)
  bos_ha_raster_max  <- bos_opp_raster_max * terra::cellSize(bos_opp_raster_max, unit = "ha")
  
  bos_opp_per_patch_max <- terra::zonal(bos_ha_raster_max, r_patch_ids_max, fun = "sum", na.rm = TRUE)
  colnames(bos_opp_per_patch_max) <- c("Patch_ID", "Bos_ha")
  
  patch_stats_max <- merge(totale_opp_per_patch_max, bos_opp_per_patch_max, by = "Patch_ID", all.x = TRUE)
  patch_stats_max$Bos_ha[is.na(patch_stats_max$Bos_ha)] <- 0
  patch_stats_max$Percentage_Bos <- (patch_stats_max$Bos_ha / patch_stats_max$Totaal_ha) * 100
  
  goedgekeurde_patch_ids_max <- patch_stats_max$Patch_ID[patch_stats_max$Percentage_Bos >= 60]
  
  if (length(goedgekeurde_patch_ids_max) > 0) {
    venwitsnuitlibel_land_ruim_max <- terra::ifel(r_patch_ids_max %in% goedgekeurde_patch_ids_max, 1, NA)
  }
  rm(area_ha_raster_max, totale_opp_per_patch_max, bos_opp_raster_max, bos_ha_raster_max, bos_opp_per_patch_max, patch_stats_max)
}

# --- SPOOR B: OPP (Parallel & Autonoom) ---
r_binair_ruim_opp_init <- terra::ifel(!is.na(landbiotoop_bwk_ruim_opp) & landbiotoop_bwk_ruim_opp > 0, 1, NA)

land_ruim_clusters_opp <- cluster_filter_compleet(
  masker     = r_binair_ruim_opp_init,
  opp_laag   = landbiotoop_bwk_ruim_opp,
  drempel_m2 = 200000,   
  dist_m     = 50,       
  werkelijk  = TRUE      
)
r_patch_ids_opp <- land_ruim_clusters_opp$clusters
venwitsnuitlibel_land_ruim_opp <- template_TV * NA

if (!all(is.na(terra::values(r_patch_ids_opp, mat=FALSE)))) {
  totale_opp_per_patch_opp <- terra::zonal(landbiotoop_bwk_ruim_opp, r_patch_ids_opp, fun = "sum", na.rm = TRUE)
  colnames(totale_opp_per_patch_opp) <- c("Patch_ID", "Totaal_ha")
  totale_opp_per_patch_opp$Totaal_ha <- totale_opp_per_patch_opp$Totaal_ha * 0.01
  
  bos_in_ruim_opp <- terra::mask(landbiotoop_bwk_ruim_bos_opp, r_patch_ids_opp)
  bos_opp_per_patch_opp <- terra::zonal(bos_in_ruim_opp, r_patch_ids_opp, fun = "sum", na.rm = TRUE)
  colnames(bos_opp_per_patch_opp) <- c("Patch_ID", "Bos_ha")
  bos_opp_per_patch_opp$Bos_ha <- bos_opp_per_patch_opp$Bos_ha * 0.01
  
  patch_stats_opp <- merge(totale_opp_per_patch_opp, bos_opp_per_patch_opp, by = "Patch_ID", all.x = TRUE)
  patch_stats_opp$Bos_ha[is.na(patch_stats_opp$Bos_ha)] <- 0
  patch_stats_opp$Percentage_Bos <- (patch_stats_opp$Bos_ha / patch_stats_opp$Totaal_ha) * 100
  
  goedgekeurde_patch_ids_opp <- patch_stats_opp$Patch_ID[patch_stats_opp$Percentage_Bos >= 60]
  
  if (length(goedgekeurde_patch_ids_opp) > 0) {
    venwitsnuitlibel_land_ruim_opp <- terra::mask(land_ruim_clusters_opp$raster, r_patch_ids_opp %in% goedgekeurde_patch_ids_opp)
  }
  rm(totale_opp_per_patch_opp, bos_in_ruim_opp, bos_opp_per_patch_opp, patch_stats_opp)
}

rm(land_ruim_clusters_max, land_ruim_clusters_opp, r_patch_ids_max, r_patch_ids_opp, r_binair_ruim_opp_init,
   landbiotoop_bwk_ruim_max, landbiotoop_bwk_ruim_opp, landbiotoop_bwk_ruim_bos_max, landbiotoop_bwk_ruim_bos_opp)
gc()

message("-> Ruimtelijke interacties uitvoeren + Water-First 500m vector-afsnijding...")

straal_water_m <- 500

r_water1_max_clean <- waterbiotoop_clusters_max$raster
r_water1_opp_clean <- waterbiotoop_clusters_opp$raster

# ==============================================================================
# VERWERKING SPOOR 1: MAXIMALE POTENTIE
# ==============================================================================
water_buffer_10m_max <- terra::buffer(r_water1_max_clean, width = 10)
land_buffer_10m_max  <- terra::buffer(venwitsnuitlibel_land_nabij_max, width = 10)

waterbiotoop2_max <- terra::mask(r_water1_max_clean, land_buffer_10m_max)
landbiotoop2_max  <- terra::mask(venwitsnuitlibel_land_nabij_max, water_buffer_10m_max)

ruim_land_buffer_500m_max <- terra::buffer(venwitsnuitlibel_land_ruim_max, width = 500)

waterbiotoop_finaal_max_ruw <- terra::mask(waterbiotoop2_max, ruim_land_buffer_500m_max)
landbiotoop_finaal_max_ruw  <- terra::mask(landbiotoop2_max, ruim_land_buffer_500m_max)

# --- WATER-FIRST VECTORIELE AFSNIJDING OP 500M ---
if (!all(is.na(terra::values(waterbiotoop_finaal_max_ruw, mat = FALSE)))) {
  poly_water_max  <- terra::as.polygons(waterbiotoop_finaal_max_ruw, aggregate = TRUE)
  poly_buffer_max <- terra::buffer(poly_water_max, width = straal_water_m)
  landbiotoop_finaal_max <- terra::mask(landbiotoop_finaal_max_ruw, poly_buffer_max)
} else {
  landbiotoop_finaal_max <- landbiotoop_finaal_max_ruw
}
waterbiotoop_finaal_max <- waterbiotoop_finaal_max_ruw

venwitsnuitlibel_leefgebied_max <- terra::cover(waterbiotoop_finaal_max, landbiotoop_finaal_max)
venwitsnuitlibel_leefgebied_max <- terra::crop(venwitsnuitlibel_leefgebied_max, template_TV)
leefgebied_max    <- venwitsnuitlibel_leefgebied_max
finaal_max_binair <- terra::ifel(!is.na(venwitsnuitlibel_leefgebied_max) & venwitsnuitlibel_leefgebied_max > 0, 1, NA)
cl_id_max         <- terra::patches(finaal_max_binair, directions = 8, zeroAsNA = TRUE)


# ==============================================================================
# VERWERKING SPOOR 2: WERKELIJKE OPPERVLAKTE (PARALLEL)
# ==============================================================================
water_buffer_10m_opp <- terra::buffer(!is.na(r_water1_opp_clean) & r_water1_opp_clean > 0, width = 10)
land_buffer_10m_opp  <- terra::buffer(!is.na(venwitsnuitlibel_land_nabij_opp) & venwitsnuitlibel_land_nabij_opp > 0, width = 10)

waterbiotoop2_opp <- terra::mask(r_water1_opp_clean, land_buffer_10m_opp)
landbiotoop2_opp  <- terra::mask(venwitsnuitlibel_land_nabij_opp, water_buffer_10m_opp)

ruim_land_buffer_500m_opp <- terra::buffer(!is.na(venwitsnuitlibel_land_ruim_opp) & venwitsnuitlibel_land_ruim_opp > 0, width = 500)

waterbiotoop_finaal_opp_ruw <- terra::mask(waterbiotoop2_opp, ruim_land_buffer_500m_opp)
landbiotoop_finaal_opp_ruw  <- terra::mask(landbiotoop2_opp, ruim_land_buffer_500m_opp)

# --- WATER-FIRST VECTORIELE AFSNIJDING OP 500M ---
if (!all(is.na(terra::values(waterbiotoop_finaal_opp_ruw, mat = FALSE)))) {
  r_bin_water_opp <- terra::ifel(!is.na(waterbiotoop_finaal_opp_ruw) & waterbiotoop_finaal_opp_ruw > 0, 1, NA)
  poly_water_opp  <- terra::as.polygons(r_bin_water_opp, aggregate = TRUE)
  poly_buffer_opp <- terra::buffer(poly_water_opp, width = straal_water_m)
  landbiotoop_finaal_opp <- terra::mask(landbiotoop_finaal_opp_ruw, poly_buffer_opp)
} else {
  landbiotoop_finaal_opp <- landbiotoop_finaal_opp_ruw
}
waterbiotoop_finaal_opp <- waterbiotoop_finaal_opp_ruw

w_clean <- terra::ifel(is.na(waterbiotoop_finaal_opp), 0, waterbiotoop_finaal_opp)
l_clean <- terra::ifel(is.na(landbiotoop_finaal_opp), 0, landbiotoop_finaal_opp)
som_raw <- w_clean + l_clean
som_cl  <- terra::clamp(som_raw, upper = 1.0)

venwitsnuitlibel_leefgebied_opp <- terra::ifel(som_cl > 0, som_cl, NA)
venwitsnuitlibel_leefgebied_opp <- terra::crop(venwitsnuitlibel_leefgebied_opp, template_TV)

leefgebied_opp    <- venwitsnuitlibel_leefgebied_opp
finaal_opp_binair <- terra::ifel(!is.na(leefgebied_opp) & leefgebied_opp > 0, 1, NA)
cl_id_opp         <- terra::patches(finaal_opp_binair, directions = 8, zeroAsNA = TRUE)

rm(water_buffer_10m_max, land_buffer_10m_max, waterbiotoop2_max, landbiotoop2_max, ruim_land_buffer_500m_max,
   water_buffer_10m_opp, land_buffer_10m_opp, waterbiotoop2_opp, landbiotoop2_opp, ruim_land_buffer_500m_opp,
   w_clean, l_clean, som_raw, som_cl, r_water1_max_clean, r_water1_opp_clean,
   venwitsnuitlibel_land_nabij_max, venwitsnuitlibel_land_nabij_opp, venwitsnuitlibel_land_ruim_max, venwitsnuitlibel_land_ruim_opp)
gc()

message("-> Venwitsnuitlibel modellering voltooid met water-first afsnijding op 500m!")


# ==============================================================================
# SCHONE EXPORT BIOTOOP EN ANALYTISCH ID-RASTER (VOOR SCRIPT 2 / ARPL)
# ==============================================================================
base_dir <- here::here("data/output/Turnhouts_Vennegebied/Rasters_Soorten", scenario_naam)

folders <- list(
  potentie  = file.path(base_dir, "01_Maximale_Potentie"),
  werkelijk = file.path(base_dir, "02_Werkelijke_Oppervlaktes"),
  id_raster = file.path(base_dir, "00_ID_Rasters")
)
purrr::walk(folders, ~if (!dir.exists(.x)) dir.create(.x, showWarnings = FALSE, recursive = TRUE))

# 1. Bepaal Maximale Potentie Raster
potentie_export_rast <- if (exists("grutto_kaart_A") && !is.null(grutto_kaart_A) && !all(is.na(terra::values(grutto_kaart_A, mat=FALSE)))) {
  terra::ifel(grutto_kaart_A > 0, 1, NA)
} else if (exists("final_max") && !all(is.na(terra::values(final_max, mat=FALSE)))) {
  terra::ifel(!is.na(final_max) & final_max > 0, 1, NA)
} else {
  terra::rast(template_TV, vals = NA)
}

# 2. Bepaal Werkelijke Oppervlakte Raster
werkelijk_export_rast <- if (exists("resB_strikt") && !is.null(resB_strikt) && (!all(is.na(terra::values(resB_strikt$kern, mat=FALSE))) || !all(is.na(terra::values(resB_strikt$bouw, mat=FALSE))))) {
  r_net_totaal_opp <- terra::cover(resB_strikt$kern, resB_strikt$bouw)
  terra::ifel(!is.na(r_net_totaal_opp) & r_net_totaal_opp > 0, 1, NA)
} else if (exists("final_opp") && !all(is.na(terra::values(final_opp, mat=FALSE)))) {
  terra::ifel(!is.na(final_opp) & final_opp > 0, 1, NA)
} else {
  terra::rast(template_TV, vals = NA)
}

# 3. Bepaal Analytisch Metacluster ID-raster
if (exists("resB_strikt") && !is.null(resB_strikt) && (!all(is.na(terra::values(resB_strikt$kern, mat=FALSE))) || !all(is.na(terra::values(resB_strikt$bouw, mat=FALSE))))) {
  r_net_totaal <- terra::cover(resB_strikt$kern, resB_strikt$bouw)
  r_buf        <- terra::buffer(r_net_totaal, width = 100)
  id_export_rast <- terra::mask(terra::patches(r_buf, directions = 8, zeroAsNA = TRUE), r_net_totaal)
} else if (exists("cl_max") && !all(is.na(terra::values(cl_max, mat=FALSE)))) {
  id_export_rast <- cl_max
} else if (exists("cl_opp") && !all(is.na(terra::values(cl_opp, mat=FALSE)))) {
  id_export_rast <- cl_opp
} else {
  id_export_rast <- terra::rast(template_TV, vals = NA)
}

# --- EXPORT LUS VOOR VISUELE BINAIR RASTERS ---
export_config <- list(
  list(rast = potentie_export_rast,  naam = "Maximale_Potentie",        folder = folders$potentie),
  list(rast = werkelijk_export_rast, naam = "Werkelijke_Oppervlaktes", folder = folders$werkelijk)
)

for(item in export_config) {
  suffix <- if (exists("zoek_sleutel") && grepl("_wv$", zoek_sleutel)) "_wv.tif" else ".tif"
  file_path <- file.path(item$folder, paste0("Habitat_", item$naam, "_", soort, suffix))
  
  export_rast <- terra::deepcopy(item$rast)
  export_rast <- terra::extend(export_rast, master_grid, fill = NA)
  
  terra::writeRaster(
    export_rast, 
    filename = file_path, 
    overwrite = TRUE, 
    gdal = c("COMPRESS=LZW"), 
    datatype = "INT1U", 
    NAflag = 255
  )
  message(paste("    [OK] Geëxporteerd:", basename(file_path)))
}

# --- EXPORT VOOR ANALYTISCH ID-RASTER ---
suffix <- if (exists("zoek_sleutel") && grepl("_wv$", zoek_sleutel)) "_wv.tif" else ".tif"
file_path_id <- file.path(folders$id_raster, paste0("ID_Netwerken_", soort, suffix))

export_id_rast <- terra::extend(id_export_rast, master_grid, fill = NA)

terra::writeRaster(
  export_id_rast, 
  filename = file_path_id, 
  overwrite = TRUE, 
  gdal = c("COMPRESS=LZW"), 
  datatype = "INT4U", 
  NAflag = 0
)
message(paste("    [OK] Analytisch ID-raster geëxporteerd:", basename(file_path_id)))

suppressWarnings(
  rm(potentie_export_rast, werkelijk_export_rast, id_export_rast, export_config, export_rast, export_id_rast)
)
gc()

message(paste("🏁 SCENARIO EXPORT VOLLEDIG AFGEROND VOOR:", toupper(soort)))

