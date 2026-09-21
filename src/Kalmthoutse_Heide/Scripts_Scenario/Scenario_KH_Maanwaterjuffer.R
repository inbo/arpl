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

soort <- "maanwaterjuffer"

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

buffer_m <- resultaat$Dispersiecap_m[1]

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
    min_drempel <- if(grepl("waterbiotoop_bwk1", h_type, ignore.case = TRUE)) 0.20 else 0.01
    
    tabel_cel_som[, Match := ifelse(Oppervlakte >= min_drempel, 1, 0)]
    
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

waterbiotoop_bwk1_max     <- lijst_matches[["waterbiotoop_bwk1"]]
waterbiotoop_bwk1_opp     <- lijst_oppervlaktes[["waterbiotoop_bwk1"]]
waterbiotoop_notbwk_max   <- lijst_matches[["waterbiotoop_notbwk"]]
waterbiotoop_notbwk_opp   <- lijst_oppervlaktes[["waterbiotoop_notbwk"]]
landbiotoop_bwk1_max      <- lijst_matches[["landbiotoop_bwk1"]]
landbiotoop_bwk1_opp      <- lijst_oppervlaktes[["landbiotoop_bwk1"]]
landbiotoop_notbwk_max    <- lijst_matches[["landbiotoop_notbwk"]]
landbiotoop_notbwk_opp    <- lijst_oppervlaktes[["landbiotoop_notbwk"]]
landbiotoop_bos_nabij_max <- lijst_matches[["landbiotoop_nabij_bos"]]
landbiotoop_bos_nabij_opp <- lijst_oppervlaktes[["landbiotoop_nabij_bos"]]
landbiotoop_bwk_ruim_max  <- lijst_matches[["landbiotoop_bwk_ruim"]]
landbiotoop_bwk_ruim_opp  <- lijst_oppervlaktes[["landbiotoop_bwk_ruim"]]
landbiotoop_bos_ruim_max  <- lijst_matches[["landbiotoop_ruim_bos"]]
landbiotoop_bos_ruim_opp  <- lijst_oppervlaktes[["landbiotoop_ruim_bos"]]

rm(tabel_vlaanderen, vertaal_df, lijst_matches, lijst_oppervlaktes, df_nieuw, resultaten_gegroepeerd)
gc()

message("-> Waterbiotoop berekenen op basis van expert-tabel...")

waterbiotoop_bwk_basis_max <- terra::ifel(!is.na(waterbiotoop_bwk1_max) & is.na(waterbiotoop_notbwk_max), 1, NA)
waterbiotoop_bwk_basis_opp <- terra::ifel(!is.na(waterbiotoop_bwk1_opp) & waterbiotoop_bwk1_opp > 0 & (is.na(waterbiotoop_notbwk_opp) | waterbiotoop_notbwk_opp == 0), 1, NA)

waterbiotoop1_max_list <- cluster_filter_compleet(
  masker     = waterbiotoop_bwk_basis_max, 
  opp_laag   = waterbiotoop_bwk1_max, 
  drempel_m2 = 200, 
  dist_m     = 10, 
  werkelijk  = FALSE
)
waterbiotoop1_max <- waterbiotoop1_max_list$raster

r_binair_water_opp <- terra::ifel(!is.na(waterbiotoop_bwk_basis_opp) & waterbiotoop_bwk_basis_opp > 0, 1, NA)

waterbiotoop1_opp_list <- cluster_filter_compleet(
  masker     = r_binair_water_opp, 
  opp_laag   = waterbiotoop_bwk1_opp, 
  drempel_m2 = 200, 
  dist_m     = 10, 
  werkelijk  = TRUE
)
waterbiotoop1_opp <- waterbiotoop1_opp_list$raster

rm(r_binair_water_opp, waterbiotoop_bwk_basis_max, waterbiotoop_bwk_basis_opp, waterbiotoop1_max_list, waterbiotoop1_opp_list)
gc()

message("-> Landbiotoop ruimere omgeving berekenen (Minimaal 50 ha & >50% bos/struweel)...")

fix_bos_raster <- function(r_bos, template) {
  if (is.null(r_bos) || !inherits(r_bos, "SpatRaster")) {
    return(terra::rast(template, vals = NA))
  }
  return(terra::ifel(!is.na(r_bos) & r_bos > 0, 1, NA))
}

r_bos_ruim_max_clean <- fix_bos_raster(landbiotoop_bos_ruim_max, template_KH)
r_bos_ruim_opp_clean <- fix_bos_raster(landbiotoop_bos_ruim_opp, template_KH)

# SPOOR 1: MAXIMALE POTENTIE
landbiotoop_ruim_max_list <- cluster_filter_compleet(
  masker     = landbiotoop_bwk_ruim_max, 
  opp_laag   = landbiotoop_bwk_ruim_max, 
  drempel_m2 = 500000, # 50 ha
  dist_m     = 50, 
  werkelijk  = FALSE
)
r_ruim_max_src <- landbiotoop_ruim_max_list$raster

if(!all(is.na(suppressWarnings(terra::minmax(r_ruim_max_src))))) {
  cl_id_ruim_max <- terra::patches(r_ruim_max_src, directions = 8, zeroAsNA = TRUE)
  
  stats_ruim_totaal_max <- terra::zonal(r_ruim_max_src * 0.01, cl_id_ruim_max, fun = "sum", na.rm = TRUE)
  colnames(stats_ruim_totaal_max) <- c("Cluster_ID", "Totaal_Ruim_ha")
  
  bos_in_ruim_max <- terra::mask(r_bos_ruim_max_clean, cl_id_ruim_max)
  stats_bos_in_ruim_max <- terra::zonal(bos_in_ruim_max * 0.01, cl_id_ruim_max, fun = "sum", na.rm = TRUE)
  colnames(stats_bos_in_ruim_max) <- c("Cluster_ID", "Bos_ha")
  
  df_ruim_samenstelling_max <- merge(stats_ruim_totaal_max, stats_bos_in_ruim_max, by = "Cluster_ID", all.x = TRUE)
  df_ruim_samenstelling_max$Bos_ha[is.na(df_ruim_samenstelling_max$Bos_ha)] <- 0
  df_ruim_samenstelling_max$Percentage_Bos <- (df_ruim_samenstelling_max$Bos_ha / df_ruim_samenstelling_max$Totaal_Ruim_ha) * 100
  
  valide_ruim_ids_max <- df_ruim_samenstelling_max$Cluster_ID[df_ruim_samenstelling_max$Percentage_Bos > 50]
  
  if(length(valide_ruim_ids_max) > 0) {
    masker_ruim_max <- cl_id_ruim_max %in% valide_ruim_ids_max
    landbiotoop_ruim_max <- terra::mask(r_ruim_max_src, terra::ifel(masker_ruim_max, 1, NA))
  } else {
    landbiotoop_ruim_max <- template_KH * NA
  }
} else {
  landbiotoop_ruim_max <- template_KH * NA
}

# SPOOR 2: WERKELIJKE OPPERVLAKTE
r_binair_ruim_opp <- terra::ifel(!is.na(landbiotoop_bwk_ruim_opp) & landbiotoop_bwk_ruim_opp > 0, 1, NA)

landbiotoop_ruim_opp_list <- cluster_filter_compleet(
  masker     = r_binair_ruim_opp, 
  opp_laag   = landbiotoop_bwk_ruim_opp, 
  drempel_m2 = 500000, # 50 ha
  dist_m     = 50, 
  werkelijk  = TRUE
)
r_ruim_opp_src <- landbiotoop_ruim_opp_list$raster

if(!all(is.na(suppressWarnings(terra::minmax(r_ruim_opp_src))))) {
  cl_id_ruim_opp <- terra::patches(r_ruim_opp_src, directions = 8, zeroAsNA = TRUE)
  
  stats_ruim_totaal_opp <- terra::zonal(r_ruim_opp_src * 0.01, cl_id_ruim_opp, fun = "sum", na.rm = TRUE)
  colnames(stats_ruim_totaal_opp) <- c("Cluster_ID", "Totaal_Ruim_ha")
  
  bos_in_ruim_opp <- terra::mask(r_bos_ruim_opp_clean, cl_id_ruim_opp)
  stats_bos_in_ruim_opp <- terra::zonal(bos_in_ruim_opp * 0.01, cl_id_ruim_opp, fun = "sum", na.rm = TRUE)
  colnames(stats_bos_in_ruim_opp) <- c("Cluster_ID", "Bos_ha")
  
  df_ruim_samenstelling_opp <- merge(stats_ruim_totaal_opp, stats_bos_in_ruim_opp, by = "Cluster_ID", all.x = TRUE)
  df_ruim_samenstelling_opp$Bos_ha[is.na(df_ruim_samenstelling_opp$Bos_ha)] <- 0
  df_ruim_samenstelling_opp$Percentage_Bos <- (df_ruim_samenstelling_opp$Bos_ha / df_ruim_samenstelling_opp$Totaal_Ruim_ha) * 100
  
  valide_ruim_ids_opp <- df_ruim_samenstelling_opp$Cluster_ID[df_ruim_samenstelling_opp$Percentage_Bos > 50]
  
  if(length(valide_ruim_ids_opp) > 0) {
    masker_ruim_opp <- cl_id_ruim_opp %in% valide_ruim_ids_opp
    landbiotoop_ruim_opp <- terra::mask(r_ruim_opp_src, terra::ifel(masker_ruim_opp, 1, NA))
  } else {
    landbiotoop_ruim_opp <- template_KH * NA
  }
} else {
  landbiotoop_ruim_opp <- template_KH * NA
}

rm(r_binair_ruim_opp, landbiotoop_ruim_max_list, landbiotoop_ruim_opp_list, r_ruim_max_src, r_ruim_opp_src)
gc()

message("-> WATER-FIRST FILTERING: Vectoriële buffer voorkomt foute 0-waarden...")

straal_water_m <- 500

# ==============================================================================
# SPOOR 1: MAXIMALE POTENTIE
# ==============================================================================

if (exists("waterbiotoop1_max") && !all(is.na(suppressWarnings(terra::minmax(waterbiotoop1_max))))) {
  r_water_max_echt <- terra::ifel(!is.na(waterbiotoop1_max) & waterbiotoop1_max > 0, 1, NA)
} else {
  r_water_max_echt <- template_KH * NA
}

if (!all(is.na(suppressWarnings(terra::minmax(r_water_max_echt))))) {
  poly_water_max <- terra::as.polygons(r_water_max_echt, aggregate = TRUE)
  poly_buffer_max <- terra::buffer(poly_water_max, width = straal_water_m)
  
  r_land_basis_max <- terra::ifel(!is.na(landbiotoop_bwk1_max) & landbiotoop_bwk1_max > 0 & is.na(landbiotoop_notbwk_max), 1, NA)
  land_afgesneden_500m_max <- terra::mask(r_land_basis_max, poly_buffer_max)
  
  res_land_10ha_max <- cluster_filter_compleet(
    masker     = land_afgesneden_500m_max, 
    opp_laag   = land_afgesneden_500m_max, 
    drempel_m2 = 100000, # 10 ha
    dist_m     = 20, 
    werkelijk  = FALSE
  )
  land_10ha_max <- res_land_10ha_max$raster
  cl_id_10ha_max <- res_land_10ha_max$clusters
  
  if (!all(is.na(suppressWarnings(terra::minmax(land_10ha_max))))) {
    r_bos_nabij_max_clean <- fix_bos_raster(landbiotoop_bos_nabij_max, template_KH)
    stats_totaal_max <- terra::zonal(land_10ha_max * 0.01, cl_id_10ha_max, fun = "sum", na.rm = TRUE)
    colnames(stats_totaal_max) <- c("Cluster_ID", "Totaal_ha")
    bos_in_cl_max <- terra::mask(r_bos_nabij_max_clean * 0.01, cl_id_10ha_max)
    stats_bos_max <- terra::zonal(bos_in_cl_max, cl_id_10ha_max, fun = "sum", na.rm = TRUE)
    colnames(stats_bos_max) <- c("Cluster_ID", "Bos_ha")
    df_bos_max <- merge(stats_totaal_max, stats_bos_max, by = "Cluster_ID", all.x = TRUE)
    df_bos_max$Bos_ha[is.na(df_bos_max$Bos_ha)] <- 0
    df_bos_max$Pct_Bos <- (df_bos_max$Bos_ha / df_bos_max$Totaal_ha) * 100
    valide_bos_ids_max <- df_bos_max$Cluster_ID[df_bos_max$Pct_Bos >= 40]
    
    land_finaal_max <- if(length(valide_bos_ids_max) > 0) terra::mask(land_10ha_max, terra::ifel(cl_id_10ha_max %in% valide_bos_ids_max, 1, NA)) else template_KH * NA
  } else { land_finaal_max <- template_KH * NA }
  
} else {
  land_finaal_max <- template_KH * NA
}

water_finaal_max  <- r_water_max_echt
leefgebied_max    <- terra::cover(water_finaal_max, land_finaal_max)

# ==============================================================================
# SPOOR 2: WERKELIJKE OPPERVLAKTE
# ==============================================================================

if (exists("waterbiotoop1_opp") && !all(is.na(suppressWarnings(terra::minmax(waterbiotoop1_opp))))) {
  r_water_opp_echt <- terra::ifel(!is.na(waterbiotoop1_opp) & waterbiotoop1_opp > 0, waterbiotoop1_opp, NA)
  r_water_opp_bin  <- terra::ifel(!is.na(r_water_opp_echt) & r_water_opp_echt > 0, 1, NA)
} else {
  r_water_opp_echt <- template_KH * NA
  r_water_opp_bin  <- template_KH * NA
}

if (!all(is.na(suppressWarnings(terra::minmax(r_water_opp_bin))))) {
  poly_water_opp <- terra::as.polygons(r_water_opp_bin, aggregate = TRUE)
  poly_buffer_opp <- terra::buffer(poly_water_opp, width = straal_water_m)
  
  r_land_opp_basis <- terra::ifel(!is.na(landbiotoop_bwk1_opp) & landbiotoop_bwk1_opp > 0 & (is.na(landbiotoop_notbwk_opp) | landbiotoop_notbwk_opp == 0), landbiotoop_bwk1_opp, NA)
  land_afgesneden_500m_opp <- terra::mask(r_land_opp_basis, poly_buffer_opp)
  r_bin_land_afgesneden <- terra::ifel(!is.na(land_afgesneden_500m_opp) & land_afgesneden_500m_opp > 0, 1, NA)
  
  res_land_10ha_opp <- cluster_filter_compleet(
    masker     = r_bin_land_afgesneden, 
    opp_laag   = land_afgesneden_500m_opp, 
    drempel_m2 = 100000, # 10 ha
    dist_m     = 20, 
    werkelijk  = TRUE
  )
  land_10ha_opp <- res_land_10ha_opp$raster
  cl_id_10ha_opp <- res_land_10ha_opp$clusters
  
  if (!all(is.na(suppressWarnings(terra::minmax(land_10ha_opp))))) {
    r_bos_nabij_opp_clean <- fix_bos_raster(landbiotoop_bos_nabij_opp, template_KH)
    stats_totaal_opp <- terra::zonal(land_10ha_opp, cl_id_10ha_opp, fun = "sum", na.rm = TRUE)
    colnames(stats_totaal_opp) <- c("Cluster_ID", "Totaal_m2")
    bos_in_cl_opp <- terra::mask(r_bos_nabij_opp_clean, cl_id_10ha_opp)
    stats_bos_opp <- terra::zonal(bos_in_cl_opp, cl_id_10ha_opp, fun = "sum", na.rm = TRUE)
    colnames(stats_bos_opp) <- c("Cluster_ID", "Bos_m2")
    df_bos_opp <- merge(stats_totaal_opp, stats_bos_opp, by = "Cluster_ID", all.x = TRUE)
    df_bos_opp$Bos_m2[is.na(df_bos_opp$Bos_m2)] <- 0
    df_bos_opp$Pct_Bos <- (df_bos_opp$Bos_m2 / df_bos_opp$Totaal_m2) * 100
    valide_bos_ids_opp <- df_bos_opp$Cluster_ID[df_bos_opp$Pct_Bos >= 40]
    
    land_finaal_opp <- if(length(valide_bos_ids_opp) > 0) terra::mask(land_10ha_opp, terra::ifel(cl_id_10ha_opp %in% valide_bos_ids_opp, 1, NA)) else template_KH * NA
  } else { land_finaal_opp <- template_KH * NA }
  
} else {
  land_finaal_opp <- template_KH * NA
}

water_finaal_opp  <- r_water_opp_echt
leefgebied_opp    <- terra::cover(water_finaal_opp, land_finaal_opp)

# DEFINITIEVE MODELUITGANGEN
final_max <- leefgebied_max
final_opp <- leefgebied_opp

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
  id_export_rast <- terra::rast(template_KH, vals = NA)
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
