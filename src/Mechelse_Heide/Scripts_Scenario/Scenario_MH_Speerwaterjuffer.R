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

soort <- "speerwaterjuffer"

# --- DYNAMISCHE SCENARIO PARAMETER CHECK ---
if (exists("SCENARIO_RDS_PAD") && !is.null(SCENARIO_RDS_PAD)) {
  scenario_rds_path <- SCENARIO_RDS_PAD
} else if (exists("params") && !is.null(params$scenario_rds_path)) {
  scenario_rds_path <- params$scenario_rds_path
} else {
  scenario_rds_path <- "data/input/Scenario_rds/MH_Scenario_BWK_2025.rds"
}

p_raw <- gsub("^([.][.]/)+", "", scenario_rds_path)
scenario_path <- here::here(p_raw)

if (!file.exists(scenario_path)) {
  stop(paste("❌ FOUT: Scenario RDS bestand NIET gevonden op:", scenario_path))
}

scen_volledig <- basename(scenario_path)
scenario_naam <- gsub("^MH_Scenario_|^Scenario_|.rds$", "", scen_volledig)

message(paste("Verwerken van soort:", soort, "binnen scenario:", scenario_naam))

df <- read_excel(here::here("data/input/Excel_files/Soorten_bwk_afstanden.xlsx"))
resultaat <- df %>%
  filter(tolower(trimws(Soort)) == soort) %>%
  select(Type, MinOpp_ha, AfstandBiotopen_m, Dispersiecap_m)

buffer_m       <- resultaat$Dispersiecap_m[1]
straal_water_m <- 500  # Harde afsnij-afstand van land rond het water (500m)

rm(df, resultaat)

area_shape  <- vect(here("data/input/Mechelse_Heide.shp"))
master_grid <- rast(here("data/input/Raster_Vlaanderen/Vlaanderen_MasterGrid_10m.tif"))[[1]]

df_namen_sleutel <- read_csv(here("data/input/Excel_files/BWK_Laag_Namen_2025.csv"), show_col_types = FALSE)
gouden_namenlijst <- tolower(trimws(df_namen_sleutel$Laagnaam))

area_shape_proj <- project(area_shape, crs(master_grid))
area_buffer_fix <- buffer(area_shape_proj, width = buffer_m)

message("-> Vertaalraster voor globale/lokale cellen opbouwen via snelle MASK methode...")
id_raster_MH <- crop(master_grid, area_buffer_fix, snap = "near")

globale_id_raster <- master_grid
globale_id_raster <- terra::init(globale_id_raster, fun = "cell")

id_raster_MH_globale_values <- crop(globale_id_raster, area_buffer_fix, snap = "near")
id_raster_MH_masked <- mask(id_raster_MH_globale_values, area_buffer_fix)

message("-> Vertaaltabel bliksemsnel opbouwen via C++ dataframe extractie...")

df_extractie <- as.data.frame(id_raster_MH_masked, cells = TRUE)
vertaal_df <- as.data.table(df_extractie)
setnames(vertaal_df, c(1, 2), c("lokale_id", "globale_id"))

vertaal_df <- vertaal_df[!is.na(globale_id)]
studiegebied_globale_ids <- unique(vertaal_df$globale_id)

values(id_raster_MH) <- NA
template_MH <- terra::rasterize(area_buffer_fix, id_raster_MH, field = 1, background = 0)

rm(globale_id_raster, id_raster_MH_globale_values, id_raster_MH_masked, df_extractie)
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
  
  tabel_MH <- tabel_gefilterd[cel_id %in% studiegebied_globale_ids]
  tabel_MH_unique <- unique(tabel_MH, by = c("cel_id", "CODE"))
  
  tabel_cel_som <- tabel_MH_unique[, .(Oppervlakte = pmin(sum(BWK_FRAC, na.rm = TRUE), 1.0)), by = .(cel_id)]
  
  r_match_type <- id_raster_MH * NA
  r_opp_type   <- id_raster_MH * NA
  
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

water1_max         <- lijst_matches[["waterbiotoop_bwk1"]]
water1_opp         <- lijst_oppervlaktes[["waterbiotoop_bwk1"]]
not_water_max      <- lijst_matches[["waterbiotoop_notbwk"]]
not_water_opp      <- lijst_oppervlaktes[["waterbiotoop_notbwk"]]
land_nabij_max     <- lijst_matches[["landbiotoop_nabij_bwk"]]
land_nabij_opp     <- lijst_oppervlaktes[["landbiotoop_nabij_bwk"]]
not_land_nabij_max <- lijst_matches[["landbiotoop_nabij_notbwk"]]
not_land_nabij_opp <- lijst_oppervlaktes[["landbiotoop_nabij_notbwk"]]
land_nabij_bos_max <- lijst_matches[["landbiotoop_nabij_bos_bwk"]]
land_nabij_bos_opp <- lijst_oppervlaktes[["landbiotoop_nabij_bos_bwk"]]
land_ruim_max      <- lijst_matches[["landbiotoop_bwk_ruim"]]
land_ruim_opp      <- lijst_oppervlaktes[["landbiotoop_bwk_ruim"]]
land_ruim_bos_max  <- lijst_matches[["landbiotoop_bwk_ruim_bos"]]
land_ruim_bos_opp  <- lijst_oppervlaktes[["landbiotoop_bwk_ruim_bos"]]

rm(vertaal_df, lijst_matches, lijst_oppervlaktes, tabel_vlaanderen, resultaten_gegroepeerd, df_nieuw)
gc()

# ==============================================================================
# STAP 1: WATERBIOTOOP BEREKENEN
# ==============================================================================
waterbiotoop_max <- (water1_max == 1) & (not_water_max == 0 | is.na(not_water_max))
waterbiotoop_max <- terra::ifel(waterbiotoop_max == 1, 1, NA)

waterbiotoop_opp_raw <- terra::mask(water1_opp, waterbiotoop_max)
waterbiotoop_opp     <- terra::ifel(!is.na(waterbiotoop_opp_raw) & waterbiotoop_opp_raw > 0, waterbiotoop_opp_raw, NA)

waterbiotoop_clusters_max <- cluster_filter_compleet(
  masker     = waterbiotoop_max,
  opp_laag   = waterbiotoop_opp,
  drempel_m2 = 200,   
  dist_m     = 10,        
  werkelijk  = FALSE      
)

r_binair_water_opp <- terra::ifel(!is.na(waterbiotoop_opp) & waterbiotoop_opp > 0, 1, NA)
waterbiotoop_clusters_opp <- cluster_filter_compleet(
  masker     = r_binair_water_opp,
  opp_laag   = waterbiotoop_opp,
  drempel_m2 = 200,   
  dist_m     = 10,        
  werkelijk  = TRUE       
)

rm(r_binair_water_opp, waterbiotoop_opp_raw)
gc()

# ==============================================================================
# STAP 2: LANDBIOTOOP NABIJE OMGEVING BEREKENEN
# ==============================================================================
land_nabij_filter_max <- (land_nabij_max == 1) & (not_land_nabij_max == 0 | is.na(not_land_nabij_max))
land_nabij_filter_max <- terra::ifel(land_nabij_filter_max == 1, 1, NA)

land_nabij_filter_opp_raw <- terra::mask(land_nabij_opp, land_nabij_filter_max)
land_nabij_filter_opp     <- terra::ifel(!is.na(land_nabij_filter_opp_raw) & land_nabij_filter_opp_raw > 0, land_nabij_filter_opp_raw, NA)

land_nabij_clusters_max <- cluster_filter_compleet(
  masker     = land_nabij_filter_max,
  opp_laag   = land_nabij_filter_opp,
  drempel_m2 = 100000,    
  dist_m     = 20,        
  werkelijk  = FALSE      
)

r_binair_land_nabij_opp <- terra::ifel(!is.na(land_nabij_filter_opp) & land_nabij_filter_opp > 0, 1, NA)
land_nabij_clusters_opp <- cluster_filter_compleet(
  masker     = r_binair_land_nabij_opp,
  opp_laag   = land_nabij_filter_opp,
  drempel_m2 = 100000,    
  dist_m     = 20,        
  werkelijk  = TRUE       
)

rm(r_binair_land_nabij_opp, land_nabij_filter_opp_raw)
gc()

# ==============================================================================
# STAP 3: 60% BOS-MOZAÏEK CHECK OP HET NABIJE LANDBIOTOOP (PARALLEL)
# ==============================================================================
message("-> Nabij landbiotoop filteren op basis van netto 60% bos-eis binnenshuis...")

# --- SPOOR A: MAX ---
cl_land_nabij_max <- land_nabij_clusters_max$clusters
land_nabij_finaal_max <- template_MH * NA

if (!all(is.na(suppressWarnings(terra::minmax(cl_land_nabij_max))))) {
  l_stats_max <- terra::freq(cl_land_nabij_max)
  df_l_stats_max <- data.frame(ID = l_stats_max$value, Land_ha = l_stats_max$count * 0.01)
  
  w_50m <- terra::focalMat(cl_land_nabij_max, 50, type = 'circle')
  w_50m[w_50m > 0] <- 1
  cl_land_buf_max <- terra::focal(cl_land_nabij_max, w = w_50m, fun = "max", na.rm = TRUE)
  
  r_bos_clean_max <- land_nabij_bos_max
  r_bos_clean_max[is.na(r_bos_clean_max)] <- 0
  
  bos_in_zone_max <- terra::mask(r_bos_clean_max, cl_land_buf_max)
  zonal_bos_max   <- terra::zonal(bos_in_zone_max * 0.01, cl_land_buf_max, fun = "sum", na.rm = TRUE)
  colnames(zonal_bos_max) <- c("ID", "Bos_ha")
  
  df_criteria_max <- merge(df_l_stats_max, zonal_bos_max, by = "ID", all.x = TRUE)
  df_criteria_max[is.na(df_criteria_max)] <- 0
  df_criteria_max$Percentage_Bos <- (df_criteria_max$Bos_ha / df_criteria_max$Land_ha) * 100
  
  goedgekeurde_land_ids_max <- df_criteria_max$ID[df_criteria_max$Land_ha >= 10 & df_criteria_max$Percentage_Bos >= 60]
  
  if (length(goedgekeurde_land_ids_max) > 0) {
    masker_land_bin_max <- cl_land_nabij_max %in% goedgekeurde_land_ids_max
    land_nabij_finaal_max <- terra::mask(land_nabij_clusters_max$raster, terra::ifel(masker_land_bin_max == 1, 1, NA))
  }
  rm(l_stats_max, df_l_stats_max, cl_land_buf_max, bos_in_zone_max, zonal_bos_max, df_criteria_max, goedgekeurde_land_ids_max, r_bos_clean_max)
}

# --- SPOOR B: OPP ---
cl_land_nabij_opp <- land_nabij_clusters_opp$clusters
land_nabij_finaal_opp <- template_MH * NA

if (!all(is.na(suppressWarnings(terra::minmax(cl_land_nabij_opp))))) {
  r_binair_land_cl_opp <- terra::ifel(!is.na(cl_land_nabij_opp) & cl_land_nabij_opp > 0, 1, NA)
  
  l_stats_opp <- terra::zonal(land_nabij_filter_opp, cl_land_nabij_opp, fun = "sum", na.rm = TRUE)
  colnames(l_stats_opp) <- c("ID", "Land_ha")
  l_stats_opp$Land_ha <- l_stats_opp$Land_ha * 0.01
  
  if (!exists("w_50m")) {
    w_50m <- terra::focalMat(cl_land_nabij_opp, 50, type = 'circle')
    w_50m[w_50m > 0] <- 1
  }
  cl_land_buf_opp <- terra::focal(cl_land_nabij_opp, w = w_50m, fun = "max", na.rm = TRUE)
  
  r_bos_clean_opp <- land_nabij_bos_opp
  r_bos_clean_opp[is.na(r_bos_clean_opp)] <- 0
  
  bos_in_zone_opp <- terra::mask(r_bos_clean_opp, cl_land_buf_opp)
  zonal_bos_opp   <- terra::zonal(bos_in_zone_opp, cl_land_buf_opp, fun = "sum", na.rm = TRUE)
  colnames(zonal_bos_opp) <- c("ID", "Bos_Fractie_Sum")
  zonal_bos_opp$Bos_ha    <- zonal_bos_opp$Bos_Fractie_Sum * 0.01
  
  df_criteria_opp <- merge(l_stats_opp, zonal_bos_opp, by = "ID", all.x = TRUE)
  df_criteria_opp[is.na(df_criteria_opp)] <- 0
  df_criteria_opp$Percentage_Bos <- (df_criteria_opp$Bos_ha / df_criteria_opp$Land_ha) * 100
  
  goedgekeurde_land_ids_opp <- df_criteria_opp$ID[df_criteria_opp$Land_ha >= 10 & df_criteria_opp$Percentage_Bos >= 60]
  
  if (length(goedgekeurde_land_ids_opp) > 0) {
    masker_land_bin_opp <- cl_land_nabij_opp %in% goedgekeurde_land_ids_opp
    land_nabij_finaal_opp <- terra::mask(land_nabij_clusters_opp$raster, terra::ifel(masker_land_bin_opp == 1, 1, NA))
  }
  rm(l_stats_opp, cl_land_buf_opp, bos_in_zone_opp, zonal_bos_opp, df_criteria_opp, goedgekeurde_land_ids_opp, r_bos_clean_opp, r_binair_land_cl_opp)
}

if (exists("w_50m")) rm(w_50m)
rm(land_nabij_bos_max, land_nabij_bos_opp)
gc()

# ==============================================================================
# STAP 4: RUIM LANDBIOTOOP 70% BOS CHECK
# ==============================================================================
message("-> Ruim landbiotoop filteren op basis van de 70% bos-samenstellingscheck (PARALLEL)...")

# --- SPOOR A: MAX ---
r_ruim_land_src_max <- land_ruim_max
r_ruim_land_src_max[is.na(r_ruim_land_src_max) | r_ruim_land_src_max == 0] <- NA

r_ruim_bos_src_max <- land_ruim_bos_max
r_ruim_bos_src_max[is.na(r_ruim_bos_src_max) | r_ruim_bos_src_max == 0] <- NA

land_ruim_clusters_list_max <- cluster_filter_compleet(
  masker     = r_ruim_land_src_max, 
  opp_laag   = land_ruim_max, 
  drempel_m2 = 500000, 
  dist_m     = 50, 
  werkelijk  = FALSE
)
cl_id_ruim_max <- land_ruim_clusters_list_max$clusters

speerwaterjuffer_land_ruim_max <- template_MH * NA

if (!all(is.na(suppressWarnings(terra::minmax(cl_id_ruim_max))))) {
  stats_ruim_totaal_max <- terra::zonal(r_ruim_land_src_max * 0.01, cl_id_ruim_max, fun = "sum", na.rm = TRUE)
  colnames(stats_ruim_totaal_max) <- c("Cluster_ID", "Totaal_Ruim_ha")
  
  bos_in_ruim_max <- terra::mask(r_ruim_bos_src_max, cl_id_ruim_max)
  stats_bos_in_ruim_max <- terra::zonal(bos_in_ruim_max * 0.01, cl_id_ruim_max, fun = "sum", na.rm = TRUE)
  colnames(stats_bos_in_ruim_max) <- c("Cluster_ID", "Bos_ha")
  
  df_samenstelling_ruim_max <- merge(stats_ruim_totaal_max, stats_bos_in_ruim_max, by = "Cluster_ID", all.x = TRUE)
  df_samenstelling_ruim_max$Bos_ha[is.na(df_samenstelling_ruim_max$Bos_ha)] <- 0
  df_samenstelling_ruim_max$Percentage_Bos <- (df_samenstelling_ruim_max$Bos_ha / df_samenstelling_ruim_max$Totaal_Ruim_ha) * 100
  
  valide_ruim_ids_max <- df_samenstelling_ruim_max$Cluster_ID[df_samenstelling_ruim_max$Percentage_Bos >= 70]
  
  if(length(valide_ruim_ids_max) > 0) {
    masker_ruim_bin_max <- cl_id_ruim_max %in% valide_ruim_ids_max
    speerwaterjuffer_land_ruim_max <- terra::mask(r_ruim_land_src_max, terra::ifel(masker_ruim_bin_max == 1, 1, NA))
  }
  rm(stats_ruim_totaal_max, bos_in_ruim_max, stats_bos_in_ruim_max, df_samenstelling_ruim_max, valide_ruim_ids_max)
}

# --- SPOOR B: OPP ---
r_binair_ruim_opp_init <- terra::ifel(!is.na(land_ruim_opp) & land_ruim_opp > 0, 1, NA)

land_ruim_clusters_list_opp <- cluster_filter_compleet(
  masker     = r_binair_ruim_opp_init, 
  opp_laag   = land_ruim_opp, 
  drempel_m2 = 500000, 
  dist_m     = 50, 
  werkelijk  = TRUE
)
cl_id_ruim_opp <- land_ruim_clusters_list_opp$clusters

speerwaterjuffer_land_ruim_opp <- template_MH * NA

if (!all(is.na(suppressWarnings(terra::minmax(cl_id_ruim_opp))))) {
  stats_ruim_totaal_opp <- terra::zonal(land_ruim_opp, cl_id_ruim_opp, fun = "sum", na.rm = TRUE)
  colnames(stats_ruim_totaal_opp) <- c("Cluster_ID", "Totaal_Ruim_ha")
  stats_ruim_totaal_opp$Totaal_Ruim_ha <- stats_ruim_totaal_opp$Totaal_Ruim_ha * 0.01
  
  bos_in_ruim_opp <- terra::mask(land_ruim_bos_opp, cl_id_ruim_opp)
  stats_bos_in_ruim_opp <- terra::zonal(bos_in_ruim_opp, cl_id_ruim_opp, fun = "sum", na.rm = TRUE)
  colnames(stats_bos_in_ruim_opp) <- c("Cluster_ID", "Bos_ha")
  stats_bos_in_ruim_opp$Bos_ha <- stats_bos_in_ruim_opp$Bos_ha * 0.01
  
  df_samenstelling_ruim_opp <- merge(stats_ruim_totaal_opp, stats_bos_in_ruim_opp, by = "Cluster_ID", all.x = TRUE)
  df_samenstelling_ruim_opp$Bos_ha[is.na(df_samenstelling_ruim_opp$Bos_ha)] <- 0
  df_samenstelling_ruim_opp$Percentage_Bos <- (df_samenstelling_ruim_opp$Bos_ha / df_samenstelling_ruim_opp$Totaal_Ruim_ha) * 100
  
  valide_ruim_ids_opp <- df_samenstelling_ruim_opp$Cluster_ID[df_samenstelling_ruim_opp$Percentage_Bos >= 70]
  
  if(length(valide_ruim_ids_opp) > 0) {
    masker_ruim_bin_opp <- cl_id_ruim_opp %in% valide_ruim_ids_opp
    speerwaterjuffer_land_ruim_opp <- terra::mask(land_ruim_opp, terra::ifel(masker_ruim_bin_opp == 1, 1, NA))
  }
  rm(stats_ruim_totaal_opp, bos_in_ruim_opp, stats_bos_in_ruim_opp, df_samenstelling_ruim_opp, valide_ruim_ids_opp)
}

rm(r_ruim_land_src_max, r_ruim_bos_src_max, land_ruim_clusters_list_max, cl_id_ruim_max, 
   r_binair_ruim_opp_init, land_ruim_clusters_list_opp, cl_id_ruim_opp,
   land_ruim_max, land_ruim_opp, land_ruim_bos_max, land_ruim_bos_opp)
gc()

message("-> Ruimtelijke interacties uitvoeren + Water-First 500m vector-afsnijding...")

r_water1_max_clean <- waterbiotoop_clusters_max$raster
r_water1_opp_clean <- waterbiotoop_clusters_opp$raster

# ==============================================================================
# STAP 5: RUIMTELIJKE INTERACTIES EN WATER-FIRST AFSNIJDING (500M)
# ==============================================================================

# --- SPOOR A: MAX ---
water_buffer10_max <- terra::buffer(r_water1_max_clean, width = 10)
land_buffer10_max  <- terra::buffer(land_nabij_finaal_max, width = 10)

waterbiotoop2_max <- terra::mask(r_water1_max_clean, land_buffer10_max)
landbiotoop2_max  <- terra::mask(land_nabij_finaal_max, water_buffer10_max)

ruim_buffer500_max <- terra::buffer(speerwaterjuffer_land_ruim_max, width = 500)

waterbiotoop_finaal_max_ruw <- terra::mask(waterbiotoop2_max, ruim_buffer500_max)
landbiotoop_finaal_max_ruw  <- terra::mask(landbiotoop2_max, ruim_buffer500_max)

if (!all(is.na(suppressWarnings(terra::minmax(waterbiotoop_finaal_max_ruw))))) {
  poly_water_max  <- terra::as.polygons(waterbiotoop_finaal_max_ruw, aggregate = TRUE)
  poly_buffer_max <- terra::buffer(poly_water_max, width = straal_water_m)
  landbiotoop_finaal_max <- terra::mask(landbiotoop_finaal_max_ruw, poly_buffer_max)
} else {
  landbiotoop_finaal_max <- landbiotoop_finaal_max_ruw
}
waterbiotoop_finaal_max <- waterbiotoop_finaal_max_ruw

speerwaterjuffer_leefgebied_max <- terra::cover(waterbiotoop_finaal_max, landbiotoop_finaal_max)

# --- SPOOR B: OPP ---
water_buffer10_opp <- terra::buffer(!is.na(r_water1_opp_clean) & r_water1_opp_clean > 0, width = 10)
land_buffer10_opp  <- terra::buffer(!is.na(land_nabij_finaal_opp) & land_nabij_finaal_opp > 0, width = 10)

waterbiotoop2_opp <- terra::mask(r_water1_opp_clean, land_buffer10_opp)
landbiotoop2_opp  <- terra::mask(land_nabij_finaal_opp, water_buffer10_opp)

ruim_buffer500_opp <- terra::buffer(!is.na(speerwaterjuffer_land_ruim_opp) & speerwaterjuffer_land_ruim_opp > 0, width = 500)

waterbiotoop_finaal_opp_ruw <- terra::mask(waterbiotoop2_opp, ruim_buffer500_opp)
landbiotoop_finaal_opp_ruw  <- terra::mask(landbiotoop2_opp, ruim_buffer500_opp)

if (!all(is.na(suppressWarnings(terra::minmax(waterbiotoop_finaal_opp_ruw))))) {
  r_bin_water_opp <- terra::ifel(!is.na(waterbiotoop_finaal_opp_ruw) & waterbiotoop_finaal_opp_ruw > 0, 1, NA)
  poly_water_opp  <- terra::as.polygons(r_bin_water_opp, aggregate = TRUE)
  poly_buffer_opp <- terra::buffer(poly_water_opp, width = straal_water_m)
  landbiotoop_finaal_opp <- terra::mask(landbiotoop_finaal_opp_ruw, poly_buffer_opp)
} else {
  landbiotoop_finaal_opp <- landbiotoop_finaal_opp_ruw
}
waterbiotoop_finaal_opp <- waterbiotoop_finaal_opp_ruw

speerwaterjuffer_leefgebied_opp <- terra::cover(waterbiotoop_finaal_opp, landbiotoop_finaal_opp)

# DEFINITIEVE MODELUITGANGEN
final_max <- speerwaterjuffer_leefgebied_max %>% terra::crop(template_MH)
final_opp <- speerwaterjuffer_leefgebied_opp %>% terra::crop(template_MH)

if (!all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  cl_opp <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE)
} else {
  cl_opp <- template_MH * NA
}

if (!all(is.na(suppressWarnings(terra::minmax(final_max))))) {
  cl_max <- terra::patches(final_max, directions = 8, zeroAsNA = TRUE)
} else {
  cl_max <- template_MH * NA
}

rm(water_buffer10_max, land_buffer10_max, waterbiotoop2_max, landbiotoop2_max, ruim_buffer500_max,
   water_buffer10_opp, land_buffer10_opp, waterbiotoop2_opp, landbiotoop2_opp, ruim_buffer500_opp,
   r_water1_max_clean, r_water1_opp_clean, land_nabij_finaal_max, land_nabij_finaal_opp,
   speerwaterjuffer_land_ruim_max, speerwaterjuffer_land_ruim_opp)
gc()

# ==============================================================================
# SCHONE EXPORT BIOTOOP EN ANALYTISCH ID-RASTER (VOOR SCRIPT 2 / ARPL)
# ==============================================================================
base_dir <- here::here("data/output/Mechelse_Heide/Rasters_Soorten", scenario_naam)

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
  terra::rast(template_MH, vals = NA)
}

# 2. Bepaal Werkelijke Oppervlakte Raster
werkelijk_export_rast <- if (exists("final_opp") && !is.null(final_opp) && !all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  terra::ifel(!is.na(final_opp) & final_opp > 0, 1, NA)
} else {
  terra::rast(template_MH, vals = NA)
}

# 3. Bepaal Analytisch Metacluster ID-raster (EXCLUSIEF OP BASIS VAN WERKELIJKE OPPERVLAKTES)
if (exists("cl_opp") && !is.null(cl_opp) && !all(is.na(suppressWarnings(terra::minmax(cl_opp))))) {
  id_export_rast <- cl_opp
} else if (exists("final_opp") && !is.null(final_opp) && !all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  id_export_rast <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE)
} else {
  id_export_rast <- template_MH * NA
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
