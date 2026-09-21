library(here)
library(knitr)
library(tidyverse)
library(sf)
library(terra)
library(readxl)
library(tidyterra)
library(data.table)

conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::select)
conflicted::conflicts_prefer(dplyr::first)
conflicted::conflicts_prefer(terra::intersect)
conflicted::conflicts_prefer(terra::any)

calc_ha_exact <- function(r) {
  if(is.null(r)) return(0)
  if(all(is.na(terra::values(r, mat=FALSE)))) return(0)
  area_raster <- r * terra::cellSize(r, unit = "ha")
  val <- terra::global(area_raster, "sum", na.rm = TRUE)[[1]]
  return(as.numeric(val))
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
  if (terra::global(is.na(masker), "sum")[[1]] == terra::ncell(masker)) {
    return(list(raster = masker * NA, clusters = masker * NA))
  }
  
  if (dist_m > 0) {
    r_binair <- terra::ifel(!is.na(masker) & masker > 0, 1, NA)
    r_buffered <- terra::buffer(r_binair, width = dist_m / 2)
    cl_network <- terra::patches(r_buffered, directions = 4, zeroAsNA = TRUE)
    cl_biotoop_only <- terra::mask(cl_network, masker)
  } else {
    cl_network <- terra::patches(masker, directions = 8, zeroAsNA = TRUE)
    cl_biotoop_only <- cl_network
  }
  
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
  
  voldoet_ids <- stats_df$ID[stats_df$Area_m2 >= drempel_m2]
  if(length(voldoet_ids) == 0) return(list(raster = masker * NA, clusters = masker * NA))
  
  masker_binair <- cl_biotoop_only %in% voldoet_ids
  final_network_mask <- terra::ifel(masker_binair == 1, 1, NA)
  
  r_finaal  <- terra::mask(masker, final_network_mask)
  cl_finaal <- terra::mask(cl_biotoop_only, r_finaal) 
  
  return(list(raster = r_finaal, clusters = cl_finaal))
}

terraOptions(
  memfrac = 0.8,
  tempdir = tempdir(),
  verbose = FALSE
)

soort <- "venglazenmaker"

# --- DYNAMISCHE SCENARIO PARAMETER CHECK ---
if (exists("SCENARIO_RDS_PAD") && !is.null(SCENARIO_RDS_PAD)) {
  scenario_rds_path <- SCENARIO_RDS_PAD
} else if (exists("params") && !is.null(params$scenario_rds_path)) {
  scenario_rds_path <- params$scenario_rds_path
} else {
  scenario_rds_path <- "data/input/Scenario_rds/KH_Scenario_BWK_2025.rds"
}

p_raw <- gsub("^([.][.]/)+", "", scenario_rds_path)
scenario_path <- here::here(p_raw)

if (!file.exists(scenario_path)) {
  stop(paste("❌ FOUT: Scenario RDS bestand NIET gevonden op:", scenario_path))
}

scen_volledig <- basename(scenario_path)
scenario_naam <- gsub("^KH_Scenario_|^Scenario_|.rds$", "", scen_volledig)

message(paste("Verwerken van soort:", soort, "binnen scenario:", scenario_naam))

df <- read_excel(here::here("data/input/Excel_files/Soorten_bwk_afstanden.xlsx"))
resultaat <- df %>%
  filter(tolower(trimws(Soort)) == soort) %>%
  select(Type, MinOpp_ha, AfstandBiotopen_m, Dispersiecap_m)

buffer_m       <- resultaat$Dispersiecap_m[1]
straal_water_m <- 500  # Harde afsnij-afstand van land rond het water (500m)

rm(df, resultaat)

area_shape  <- vect(here("data/input/Kalmthoutse_Heide.shp"))
master_grid <- rast(here("data/input/Raster_Vlaanderen/Vlaanderen_MasterGrid_10m.tif"))[[1]]

df_namen_sleutel <- read_csv(here("data/input/Excel_files/BWK_Laag_Namen_2025.csv"), show_col_types = FALSE)
gouden_namenlijst <- tolower(trimws(df_namen_sleutel$Laagnaam))

area_shape_proj <- project(area_shape, crs(master_grid))
area_buffer_fix <- buffer(area_shape_proj, width = buffer_m)

message("-> Vertaalraster voor globale/lokale cellen opbouwen via snelle MASK methode...")
id_raster_KH <- crop(master_grid, area_buffer_fix, snap = "near")

globale_id_raster <- master_grid
globale_id_raster <- terra::init(globale_id_raster, fun = "cell")

id_raster_KH_globale_values <- crop(globale_id_raster, area_buffer_fix, snap = "near")
id_raster_KH_masked <- mask(id_raster_KH_globale_values, area_buffer_fix)

message("-> Vertaaltabel bliksemsnel opbouwen via C++ dataframe extractie...")

df_extractie <- as.data.frame(id_raster_KH_masked, cells = TRUE)
vertaal_df <- as.data.table(df_extractie)
setnames(vertaal_df, c(1, 2), c("lokale_id", "globale_id"))

vertaal_df <- vertaal_df[!is.na(globale_id)]
studiegebied_globale_ids <- unique(vertaal_df$globale_id)

values(id_raster_KH) <- NA
template_KH <- terra::rasterize(area_buffer_fix, id_raster_KH, field = 1, background = 0)

grens_web <- sf::st_as_sf(terra::project(area_shape, "EPSG:4326"))

rm(globale_id_raster, id_raster_KH_globale_values, id_raster_KH_masked, df_extractie)
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
    bevat_codes_escaped  <- gsub("([\\.\\^\\$\\*\\+\\?\\(\\)\\[\\{\\\\\\|])", "\\\\\\1", bevat_codes)
    bevat_codes_anchored <- paste0("^", bevat_codes_escaped)
    regex_term           <- paste0(bevat_codes_anchored, collapse = "|")
    
    tabel_gefilterd <- tabel_vlaanderen[CODE %in% exact_codes | grepl(regex_term, CODE)]
  } else {
    tabel_gefilterd <- tabel_vlaanderen[CODE %in% exact_codes]
  }
  
  tabel_KH <- tabel_gefilterd[cel_id %in% studiegebied_globale_ids]
  tabel_KH_unique <- unique(tabel_KH, by = c("cel_id", "CODE"))
  
  tabel_cel_som <- tabel_KH_unique[, .(Oppervlakte = pmin(sum(BWK_FRAC, na.rm = TRUE), 1.0)), by = .(cel_id)]
  
  r_match_type <- id_raster_KH * NA
  r_opp_type   <- id_raster_KH * NA
  
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
landbiotoop_bwk1_max         <- lijst_matches[["landbiotoop_bwk1"]]
landbiotoop_bwk1_opp         <- lijst_oppervlaktes[["landbiotoop_bwk1"]]
landbiotoop_notbwk_max       <- lijst_matches[["landbiotoop_notbwk"]]
landbiotoop_notbwk_opp       <- lijst_oppervlaktes[["landbiotoop_notbwk"]]
landbiotoop_nabij_bos_max    <- lijst_matches[["landbiotoop_nabij_bos"]]
landbiotoop_nabij_bos_opp    <- lijst_oppervlaktes[["landbiotoop_nabij_bos"]]
landbiotoop_bwk_ruim_max     <- lijst_matches[["landbiotoop_bwk_ruim"]]
landbiotoop_bwk_ruim_opp     <- lijst_oppervlaktes[["landbiotoop_bwk_ruim"]]
landbiotoop_bos_ruim_max     <- lijst_matches[["landbiotoop_bos_ruim"]]
landbiotoop_bos_ruim_opp     <- lijst_oppervlaktes[["landbiotoop_bos_ruim"]]

rm(tabel_vlaanderen, vertaal_df, lijst_matches, lijst_oppervlaktes, df_nieuw, resultaten_gegroepeerd)
gc()

# ==============================================================================
# BEREIDING WATERBIOTOOP VENGLAZENMAKER (Minimaal 400 m² - PARALLEL)
# ==============================================================================
message("-> Waterbiotoop filteren (combinatie van geschikte en uitgesloten BWK-codes)...")

venglazenmaker_water_max <- (waterbiotoop_bwk1_max == 1) & (waterbiotoop_notbwk_max == 0 | is.na(waterbiotoop_notbwk_max))
venglazenmaker_water_max <- terra::ifel(venglazenmaker_water_max == 1, 1, NA)

venglazenmaker_water_opp_raw <- terra::mask(waterbiotoop_bwk1_opp, venglazenmaker_water_max)
venglazenmaker_water_opp     <- terra::ifel(!is.na(venglazenmaker_water_opp_raw) & venglazenmaker_water_opp_raw > 0, venglazenmaker_water_opp_raw, NA)

# Spoor A: MAX
waterbiotoop_clusters_max <- cluster_filter_compleet(
  masker     = venglazenmaker_water_max,
  opp_laag   = venglazenmaker_water_opp,
  drempel_m2 = 400,   
  dist_m     = 50,        
  werkelijk  = FALSE      
)

# Spoor B: OPP
r_binair_water_opp <- terra::ifel(!is.na(venglazenmaker_water_opp) & venglazenmaker_water_opp > 0, 1, NA)
waterbiotoop_clusters_opp <- cluster_filter_compleet(
  masker     = r_binair_water_opp,
  opp_laag   = venglazenmaker_water_opp,
  drempel_m2 = 400,   
  dist_m     = 50,        
  werkelijk  = TRUE       
)

gc()

# ==============================================================================
# BEREIDING LANDBIOTOOP NABIJ (>40% BOS)
# ==============================================================================
message("-> Landbiotoop nabij filteren en uitsluiten van ongeschikte codes...")

land_nabij_filter_max <- (landbiotoop_bwk1_max == 1) & (landbiotoop_notbwk_max == 0 | is.na(landbiotoop_notbwk_max))
land_nabij_filter_max <- terra::ifel(land_nabij_filter_max == 1, 1, NA)

land_nabij_filter_opp_raw <- terra::mask(landbiotoop_bwk1_opp, land_nabij_filter_max)
land_nabij_filter_opp     <- terra::ifel(!is.na(land_nabij_filter_opp_raw) & land_nabij_filter_opp_raw > 0, land_nabij_filter_opp_raw, NA)

# SPOOR A: MAX
land_nabij_clusters_max <- cluster_filter_compleet(
  masker     = land_nabij_filter_max,
  opp_laag   = land_nabij_filter_opp,
  drempel_m2 = 100000,    
  dist_m     = 50,        
  werkelijk  = FALSE      
)
r_nabij_patch_ids_max <- land_nabij_clusters_max$clusters
venglazenmaker_land_nabij_max <- template_KH * NA
goedgekeurde_patch_ids_nabij_max <- c()

if (!all(is.na(suppressWarnings(terra::minmax(r_nabij_patch_ids_max))))) {
  r_binair_patch_max <- terra::ifel(!is.na(r_nabij_patch_ids_max), 1, NA)
  totale_opp_per_patch_nabij_max <- terra::zonal(r_binair_patch_max, r_nabij_patch_ids_max, fun = "sum", na.rm = TRUE)
  colnames(totale_opp_per_patch_nabij_max) <- c("Patch_ID", "Totaal_ha")
  totale_opp_per_patch_nabij_max$Totaal_ha <- totale_opp_per_patch_nabij_max$Totaal_ha * 0.01
  
  bos_binair_max <- terra::ifel(!is.na(landbiotoop_nabij_bos_max) & landbiotoop_nabij_bos_max > 0, 1, NA)
  bos_opp_raster_nabij_max <- terra::mask(bos_binair_max, r_nabij_patch_ids_max)
  
  bos_opp_per_patch_nabij_max <- terra::zonal(bos_opp_raster_nabij_max, r_nabij_patch_ids_max, fun = "sum", na.rm = TRUE)
  colnames(bos_opp_per_patch_nabij_max) <- c("Patch_ID", "Bos_ha")
  bos_opp_per_patch_nabij_max$Bos_ha <- bos_opp_per_patch_nabij_max$Bos_ha * 0.01
  
  patch_stats_nabij_max <- merge(totale_opp_per_patch_nabij_max, bos_opp_per_patch_nabij_max, by = "Patch_ID", all.x = TRUE)
  patch_stats_nabij_max$Bos_ha[is.na(patch_stats_nabij_max$Bos_ha)] <- 0
  patch_stats_nabij_max$Percentage_Bos <- round((patch_stats_nabij_max$Bos_ha / patch_stats_nabij_max$Totaal_ha) * 100, 2)
  
  goedgekeurde_patch_ids_nabij_max <- patch_stats_nabij_max$Patch_ID[patch_stats_nabij_max$Percentage_Bos >= 40]
  
  if (length(goedgekeurde_patch_ids_nabij_max) > 0) {
    m_match_nabij_max <- terra::match(r_nabij_patch_ids_max, goedgekeurde_patch_ids_nabij_max)
    venglazenmaker_land_nabij_max <- terra::ifel(!is.na(m_match_nabij_max), 1, NA)
  }
  rm(r_binair_patch_max, totale_opp_per_patch_nabij_max, bos_binair_max, bos_opp_raster_nabij_max, bos_opp_per_patch_nabij_max)
}

# SPOOR B: OPP
r_binair_land_nabij_opp <- terra::ifel(!is.na(land_nabij_filter_opp) & land_nabij_filter_opp > 0, 1, NA)
land_nabij_clusters_opp <- cluster_filter_compleet(
  masker     = r_binair_land_nabij_opp,
  opp_laag   = land_nabij_filter_opp,
  drempel_m2 = 100000,    
  dist_m     = 50,        
  werkelijk  = TRUE       
)
r_nabij_patch_ids_opp <- land_nabij_clusters_opp$clusters
venglazenmaker_land_nabij_opp <- template_KH * NA
goedgekeurde_patch_ids_nabij_opp <- c()

if (!all(is.na(suppressWarnings(terra::minmax(r_nabij_patch_ids_opp))))) {
  totale_opp_per_patch_nabij_opp <- terra::zonal(land_nabij_filter_opp, r_nabij_patch_ids_opp, fun = "sum", na.rm = TRUE)
  colnames(totale_opp_per_patch_nabij_opp) <- c("Patch_ID", "Totaal_ha")
  totale_opp_per_patch_nabij_opp$Totaal_ha <- totale_opp_per_patch_nabij_opp$Totaal_ha * 0.01
  
  bos_in_zone_opp <- terra::mask(landbiotoop_nabij_bos_opp, r_nabij_patch_ids_opp)
  bos_opp_per_patch_nabij_opp <- terra::zonal(bos_in_zone_opp, r_nabij_patch_ids_opp, fun = "sum", na.rm = TRUE)
  colnames(bos_opp_per_patch_nabij_opp) <- c("Patch_ID", "Bos_ha")
  bos_opp_per_patch_nabij_opp$Bos_ha <- bos_opp_per_patch_nabij_opp$Bos_ha * 0.01
  
  patch_stats_nabij_opp <- merge(totale_opp_per_patch_nabij_opp, bos_opp_per_patch_nabij_opp, by = "Patch_ID", all.x = TRUE)
  patch_stats_nabij_opp$Bos_ha[is.na(patch_stats_nabij_opp$Bos_ha)] <- 0
  patch_stats_nabij_opp$Percentage_Bos <- round((patch_stats_nabij_opp$Bos_ha / patch_stats_nabij_opp$Totaal_ha) * 100, 2)
  
  goedgekeurde_patch_ids_nabij_opp <- patch_stats_nabij_opp$Patch_ID[patch_stats_nabij_opp$Percentage_Bos >= 40]
  
  if (length(goedgekeurde_patch_ids_nabij_opp) > 0) {
    m_match_nabij_opp <- terra::match(r_nabij_patch_ids_opp, goedgekeurde_patch_ids_nabij_opp)
    winnende_masker_opp <- terra::ifel(!is.na(m_match_nabij_opp), 1, NA)
    venglazenmaker_land_nabij_opp <- terra::mask(land_nabij_filter_opp, winnende_masker_opp)
  }
  rm(totale_opp_per_patch_nabij_opp, bos_in_zone_opp, bos_opp_per_patch_nabij_opp)
}

gc()

message("-> Landbiotoop ruime omgeving patches groeperen (50m) en filteren op min. 30 ha (PARALLEL)...")

# SPOOR A: MAX
land_ruim_clusters_max <- cluster_filter_compleet(
  masker     = landbiotoop_bwk_ruim_max,
  opp_laag   = landbiotoop_bwk_ruim_opp,
  drempel_m2 = 300000,    
  dist_m     = 50,        
  werkelijk  = FALSE      
)
r_patch_ids_max <- land_ruim_clusters_max$clusters
venglazenmaker_land_ruim_max <- template_KH * NA

if (!all(is.na(suppressWarnings(terra::minmax(r_patch_ids_max))))) {
  area_ha_raster_max <- r_patch_ids_max * terra::cellSize(r_patch_ids_max, unit = "ha")
  totale_opp_per_patch_max <- terra::zonal(area_ha_raster_max, r_patch_ids_max, fun = "sum", na.rm = TRUE)
  colnames(totale_opp_per_patch_max) <- c("Patch_ID", "Totaal_ha")
  
  bos_opp_raster_max <- terra::mask(landbiotoop_bos_ruim_max, r_patch_ids_max)
  bos_ha_raster_max  <- bos_opp_raster_max * terra::cellSize(bos_opp_raster_max, unit = "ha")
  
  bos_opp_per_patch_max <- terra::zonal(bos_ha_raster_max, r_patch_ids_max, fun = "sum", na.rm = TRUE)
  colnames(bos_opp_per_patch_max) <- c("Patch_ID", "Bos_ha")
  
  patch_stats_max <- merge(totale_opp_per_patch_max, bos_opp_per_patch_max, by = "Patch_ID", all.x = TRUE)
  patch_stats_max$Bos_ha[is.na(patch_stats_max$Bos_ha)] <- 0
  patch_stats_max$Percentage_Bos <- (patch_stats_max$Bos_ha / patch_stats_max$Totaal_ha) * 100
  
  goedgekeurde_patch_ids_max <- patch_stats_max$Patch_ID[patch_stats_max$Percentage_Bos >= 40]
  
  if (length(goedgekeurde_patch_ids_max) > 0) {
    m_match_max <- terra::match(r_patch_ids_max, goedgekeurde_patch_ids_max)
    venglazenmaker_land_ruim_max <- terra::ifel(!is.na(m_match_max), 1, NA)
  }
}

# SPOOR B: OPP
r_binair_ruim_opp_init <- terra::ifel(!is.na(landbiotoop_bwk_ruim_opp) & landbiotoop_bwk_ruim_opp > 0, 1, NA)

land_ruim_clusters_opp <- cluster_filter_compleet(
  masker     = r_binair_ruim_opp_init,
  opp_laag   = landbiotoop_bwk_ruim_opp,
  drempel_m2 = 300000,    
  dist_m     = 50,        
  werkelijk  = TRUE       
)
r_patch_ids_opp <- land_ruim_clusters_opp$clusters
venglazenmaker_land_ruim_opp <- template_KH * NA

if (!all(is.na(suppressWarnings(terra::minmax(r_patch_ids_opp))))) {
  totale_opp_per_patch_opp <- terra::zonal(landbiotoop_bwk_ruim_opp, r_patch_ids_opp, fun = "sum", na.rm = TRUE)
  colnames(totale_opp_per_patch_opp) <- c("Patch_ID", "Totaal_ha")
  totale_opp_per_patch_opp$Totaal_ha <- totale_opp_per_patch_opp$Totaal_ha * 0.01
  
  bos_in_ruim_opp <- terra::mask(landbiotoop_bos_ruim_opp, r_patch_ids_opp)
  bos_opp_per_patch_opp <- terra::zonal(bos_in_ruim_opp, r_patch_ids_opp, fun = "sum", na.rm = TRUE)
  colnames(bos_opp_per_patch_opp) <- c("Patch_ID", "Bos_ha")
  bos_opp_per_patch_opp$Bos_ha <- bos_opp_per_patch_opp$Bos_ha * 0.01
  
  patch_stats_opp <- merge(totale_opp_per_patch_opp, bos_opp_per_patch_opp, by = "Patch_ID", all.x = TRUE)
  patch_stats_opp$Bos_ha[is.na(patch_stats_opp$Bos_ha)] <- 0
  patch_stats_opp$Percentage_Bos <- (patch_stats_opp$Bos_ha / patch_stats_opp$Totaal_ha) * 100
  
  goedgekeurde_patch_ids_opp <- patch_stats_opp$Patch_ID[patch_stats_opp$Percentage_Bos >= 40]
  
  if (length(goedgekeurde_patch_ids_opp) > 0) {
    m_match_ruim <- terra::match(r_patch_ids_opp, goedgekeurde_patch_ids_opp)
    winnende_masker_ruim <- terra::ifel(!is.na(m_match_ruim), 1, NA)
    venglazenmaker_land_ruim_opp <- terra::mask(landbiotoop_bwk_ruim_opp, winnende_masker_ruim)
  }
}

gc()

message("-> Ruimtelijke interacties uitvoeren + Water-First 500m vector-afsnijding...")

straal_water_m <- 500

r_water1_max_clean <- waterbiotoop_clusters_max$raster
r_water1_opp_clean <- waterbiotoop_clusters_opp$raster

# SPOOR 1: MAX
water_buf20_max_bin <- terra::buffer(!is.na(r_water1_max_clean) & r_water1_max_clean > 0, width = 20)

if (exists("r_nabij_patch_ids_max") && !all(is.na(suppressWarnings(terra::minmax(r_nabij_patch_ids_max))))) {
  r_land_in_water_buf_max     <- terra::mask(r_nabij_patch_ids_max, water_buf20_max_bin)
  winnende_land_patch_ids_max <- unique(na.omit(terra::values(r_land_in_water_buf_max, mat = FALSE)))
  actieve_land_ids_max        <- intersect(winnende_land_patch_ids_max, goedgekeurde_patch_ids_nabij_max)
} else { actieve_land_ids_max <- c() }

if (length(actieve_land_ids_max) > 0) {
  m_land_gekoppeld_max <- terra::match(r_nabij_patch_ids_max, actieve_land_ids_max)
  venglazenmaker_landbiotoop_max_ruw <- terra::ifel(!is.na(m_land_gekoppeld_max), 1, NA)
} else { venglazenmaker_landbiotoop_max_ruw <- template_KH * NA }

land_buf20_max_bin              <- terra::buffer(!is.na(venglazenmaker_land_nabij_max) & venglazenmaker_land_nabij_max > 0, width = 20)
venglazenmaker_waterbiotoop_max <- terra::mask(r_water1_max_clean, land_buf20_max_bin)

# --- WATER-FIRST AFSNIJDING OP 500M ---
if (!all(is.na(suppressWarnings(terra::minmax(venglazenmaker_waterbiotoop_max))))) {
  poly_water_max  <- terra::as.polygons(venglazenmaker_waterbiotoop_max, aggregate = TRUE)
  poly_buffer_max <- terra::buffer(poly_water_max, width = straal_water_m)
  venglazenmaker_landbiotoop_max <- terra::mask(venglazenmaker_landbiotoop_max_ruw, poly_buffer_max)
} else {
  venglazenmaker_landbiotoop_max <- venglazenmaker_landbiotoop_max_ruw
}

waterbiotoop_finaal_max <- venglazenmaker_waterbiotoop_max
landbiotoop_finaal_max  <- venglazenmaker_landbiotoop_max

venglazenmaker_leefgebied_max <- terra::cover(waterbiotoop_finaal_max, landbiotoop_finaal_max)
venglazenmaker_leefgebied_max <- terra::crop(venglazenmaker_leefgebied_max, template_KH)

# SPOOR 2: OPP
r_water1_opp_bin    <- terra::ifel(!is.na(r_water1_opp_clean) & r_water1_opp_clean > 0, 1, NA)
water_buf20_opp_bin <- terra::buffer(r_water1_opp_bin, width = 20)

if (exists("r_nabij_patch_ids_opp") && !all(is.na(suppressWarnings(terra::minmax(r_nabij_patch_ids_opp))))) {
  r_land_in_water_buf_opp     <- terra::mask(r_nabij_patch_ids_opp, water_buf20_opp_bin)
  winnende_land_patch_ids_opp <- unique(na.omit(terra::values(r_land_in_water_buf_opp, mat = FALSE)))
  actieve_land_ids_opp        <- intersect(winnende_land_patch_ids_opp, goedgekeurde_patch_ids_nabij_opp)
} else { actieve_land_ids_opp <- c() }

if (length(actieve_land_ids_opp) > 0) {
  m_land_gekoppeld_opp <- terra::match(r_nabij_patch_ids_opp, actieve_land_ids_opp)
  winnende_masker_opp  <- terra::ifel(!is.na(m_land_gekoppeld_opp), 1, NA)
  venglazenmaker_landbiotoop_opp_ruw <- terra::mask(land_nabij_filter_opp, winnende_masker_opp)
} else { venglazenmaker_landbiotoop_opp_ruw <- template_KH * NA }

venglazenmaker_land_nabij_opp_bin <- terra::ifel(!is.na(venglazenmaker_land_nabij_opp) & venglazenmaker_land_nabij_opp > 0, 1, NA)
land_buf20_opp_bin                <- terra::buffer(venglazenmaker_land_nabij_opp_bin, width = 20)
venglazenmaker_waterbiotoop_opp   <- terra::mask(r_water1_opp_clean, land_buf20_opp_bin)

# --- WATER-FIRST AFSNIJDING OP 500M ---
if (!all(is.na(suppressWarnings(terra::minmax(venglazenmaker_waterbiotoop_opp))))) {
  r_bin_water_opp <- terra::ifel(!is.na(venglazenmaker_waterbiotoop_opp) & venglazenmaker_waterbiotoop_opp > 0, 1, NA)
  poly_water_opp  <- terra::as.polygons(r_bin_water_opp, aggregate = TRUE)
  poly_buffer_opp <- terra::buffer(poly_water_opp, width = straal_water_m)
  venglazenmaker_landbiotoop_opp <- terra::mask(venglazenmaker_landbiotoop_opp_ruw, poly_buffer_opp)
} else {
  venglazenmaker_landbiotoop_opp <- venglazenmaker_landbiotoop_opp_ruw
}

waterbiotoop_finaal_opp <- venglazenmaker_waterbiotoop_opp
landbiotoop_finaal_opp  <- venglazenmaker_landbiotoop_opp

opp_stack <- c(waterbiotoop_finaal_opp, landbiotoop_finaal_opp)
venglazenmaker_leefgebied_opp <- terra::app(opp_stack, fun = function(x) {
  if (all(is.na(x))) return(NA)
  return(min(sum(x, na.rm = TRUE), 1.0))
})

venglazenmaker_leefgebied_opp <- terra::crop(venglazenmaker_leefgebied_opp, template_KH)

# DEFINITIEVE MODELUITGANGEN
final_max <- venglazenmaker_leefgebied_max
final_opp <- venglazenmaker_leefgebied_opp

if (!all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  cl_opp <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE)
} else {
  cl_opp <- template_KH * NA
}

if (!all(is.na(suppressWarnings(terra::minmax(final_max))))) {
  cl_max <- terra::patches(final_max, directions = 8, zeroAsNA = TRUE)
} else {
  cl_max <- template_KH * NA
}

rm(water_buffer_10m_max, land_buffer_10m_max, waterbiotoop2_max, landbiotoop2_max, ruim_land_buffer_500m_max,
   water_buffer_10m_opp, land_buffer_10m_opp, waterbiotoop2_opp, landbiotoop2_opp, ruim_land_buffer_500m_opp,
   w_clean, l_clean, som_raw, som_cl, r_water1_max_clean, r_water1_opp_clean,
   venwitsnuitlibel_land_nabij_max, venwitsnuitlibel_land_nabij_opp, venwitsnuitlibel_land_ruim_max, venwitsnuitlibel_land_ruim_opp)
gc()

# ==============================================================================
# SCHONE EXPORT BIOTOOP EN ANALYTISCH ID-RASTER (VOOR SCRIPT 2 / ARPL)
# ==============================================================================
base_dir <- here::here("data/output/Kalmthoutse_Heide/Rasters_Soorten", scenario_naam)

folders <- list(
  potentie  = file.path(base_dir, "01_Maximale_Potentie"),
  werkelijk = file.path(base_dir, "02_Werkelijke_Oppervlaktes"),
  id_raster = file.path(base_dir, "00_ID_Rasters")
)
purrr::walk(folders, ~if (!dir.exists(.x)) dir.create(.x, showWarnings = FALSE, recursive = TRUE))

# 1. Bepaal Maximale Potentie Raster
potentie_export_rast <- if (exists("final_max") && !is.null(final_max) && !all(is.na(suppressWarnings(terra::minmax(final_max))))) {
  terra::ifel(!is.na(final_max) & final_max > 0, 1, NA)
} else {
  terra::rast(template_KH, vals = NA)
}

# 2. Bepaal Werkelijke Oppervlakte Raster
werkelijk_export_rast <- if (exists("final_opp") && !is.null(final_opp) && !all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  terra::ifel(!is.na(final_opp) & final_opp > 0, 1, NA)
} else {
  terra::rast(template_KH, vals = NA)
}

# 3. Bepaal Analytisch Metacluster ID-raster (EXCLUSIEF OP BASIS VAN WERKELIJKE OPPERVLAKTES)
if (exists("cl_opp") && !is.null(cl_opp) && !all(is.na(suppressWarnings(terra::minmax(cl_opp))))) {
  id_export_rast <- cl_opp
} else if (exists("final_opp") && !is.null(final_opp) && !all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  id_export_rast <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE)
} else {
  id_export_rast <- template_KH * NA
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
  message(paste("    [OK] Geëxporteerd naar scenariomap:", basename(file_path)))
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
