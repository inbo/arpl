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
conflicted::conflicts_prefer(terra::intersect)
conflicted::conflicts_prefer(terra::any)

calc_ha_exact <- function(r) {
  if(is.null(r)) return(0)
  if(all(is.na(terra::values(r, mat=FALSE)))) return(0)
  area_raster <- r * terra::cellSize(r, unit = "ha")
  val <- terra::global(area_raster, "sum", na.rm = TRUE)[[1]]
  return(as.numeric(val))
}

# ==============================================================================
# HOOFDFUNCTIE: GEHEUGENZUINIG CLUSTEREN & FILTEREN
# ==============================================================================
cluster_filter_compleet <- function(masker, opp_laag, drempel_m2, dist_m, werkelijk = FALSE, straal_m = 0) {
  if (terra::global(is.na(masker), "sum")[[1]] == terra::ncell(masker)) {
    return(list(raster = masker * NA, clusters = masker * NA))
  }
  
  if (exists("id_raster_TV")) {
    masker   <- terra::crop(masker, id_raster_TV, snap = "out")
    opp_laag <- terra::crop(opp_laag, id_raster_TV, snap = "out")
  }
  
  if (straal_m > 0) {
    r_binair <- terra::ifel(!is.na(masker) & masker > 0, 1, NA)
    if (terra::global(is.na(r_binair), "sum")[[1]] == terra::ncell(r_binair)) {
      return(list(raster = masker * NA, clusters = masker * NA))
    }
    
    v_biotoop <- terra::as.polygons(r_binair, aggregate = TRUE)
    v_eroded  <- terra::buffer(v_biotoop, width = -straal_m)
    
    if (length(v_eroded) == 0 || terra::geomtype(v_eroded) == "none") {
      return(list(raster = masker * NA, clusters = masker * NA))
    }
    
    v_dilated <- terra::buffer(v_eroded, width = straal_m)
    r_zuiver  <- terra::rasterize(v_dilated, r_binair, field = 1, background = NA)
    r_masker_werkelijk <- terra::mask(r_zuiver, masker)
    rm(v_biotoop, v_eroded, v_dilated, r_zuiver)
  } else {
    r_masker_werkelijk <- masker
  }
  
  if (terra::global(is.na(r_masker_werkelijk), "sum")[[1]] == terra::ncell(r_masker_werkelijk)) {
    return(list(raster = masker * NA, clusters = masker * NA))
  }
  
  if (dist_m > 0) {
    r_binair_z <- terra::ifel(!is.na(r_masker_werkelijk) & r_masker_werkelijk > 0, 1, NA)
    r_buffered <- terra::buffer(r_binair_z, width = dist_m / 2)
    cl_network <- terra::patches(r_buffered, directions = 4, zeroAsNA = TRUE)
    cl_biotoop_only <- terra::mask(cl_network, r_masker_werkelijk)
  } else {
    cl_network <- terra::patches(r_masker_werkelijk, directions = 8, zeroAsNA = TRUE)
    cl_biotoop_only <- cl_network
  }
  
  if(werkelijk) {
    opp_laag_sub <- terra::crop(opp_laag, cl_biotoop_only)
    stats_df     <- terra::zonal(opp_laag_sub, cl_biotoop_only, fun = "sum", na.rm = TRUE)
    colnames(stats_df) <- c("ID", "Waarde")
    stats_df$Area_m2  <- stats_df$Waarde * 100 
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
  
  r_finaal  <- terra::mask(r_masker_werkelijk, final_network_mask)
  cl_finaal <- terra::mask(cl_biotoop_only, r_finaal) 
  
  r_finaal  <- terra::deepcopy(r_finaal)
  cl_finaal <- terra::deepcopy(cl_finaal)
  
  return(list(raster = r_finaal, clusters = cl_finaal))
}

terraOptions(
  memfrac = 0.8,
  tempdir = tempdir(),
  verbose = FALSE
)

df <- read_excel(here::here("data/input/Excel_files/Soorten_bwk_afstanden.xlsx"))
soort <- "beekrombout"

# --- DYNAMISCHE SCENARIO PARAMETER CHECK ---
if (exists("SCENARIO_RDS_PAD") && !is.null(SCENARIO_RDS_PAD)) {
  scenario_rds_path <- SCENARIO_RDS_PAD
} else if (exists("params") && !is.null(params$scenario_rds_path)) {
  scenario_rds_path <- params$scenario_rds_path
} else {
  scenario_rds_path <- "data/input/Scenario_rds/TV_Scenario_BWK_2025.rds"
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

oppervlakte_ha <- resultaat$MinOpp_ha
afstand_m      <- resultaat$AfstandBiotopen_m
buffer_m       <- resultaat$Dispersiecap_m[1]
straal_water_m <- 500  # Afsnij-afstand van land rond de beken (500m)

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

waterlopen1_bwk_max  <- lijst_matches[["waterlopen1_bwk"]]
waterlopen1_bwk_opp  <- lijst_oppervlaktes[["waterlopen1_bwk"]]
landbiotoop_bos_max  <- lijst_matches[["landbiotoop_bos_bwk"]]
landbiotoop_bos_opp  <- lijst_oppervlaktes[["landbiotoop_bos_bwk"]]
landbiotoop_rest_max <- lijst_matches[["landbiotoop_rest_bwk"]]
landbiotoop_rest_opp <- lijst_oppervlaktes[["landbiotoop_rest_bwk"]]

rm(tabel_vlaanderen, vertaal_df, lijst_matches, lijst_oppervlaktes)
gc()

# ==============================================================================
# HYDROGRAFISCHE DATA (HUET-ZONES) INLADEN EN SYNCHRONISEREN
# ==============================================================================
message("-> Huet-zone waterlopen inlezen...")
r_huetzon_vlaanderen <- rast(here("data/input/Raster_Vlaanderen/vlaanderen_huetzon_10m.tif"))

r_huetzon_local <- r_huetzon_vlaanderen %>% 
  terra::crop(area_buffer_fix) %>% 
  terra::resample(template_TV, method = "near")

legende_huet  <- terra::cats(r_huetzon_local)[[1]]
target_labels <- c("Bg", "BgK", "Bk1", "Bk2", "BkK1", "BkK2")

valide_ids <- legende_huet$value[legende_huet$Label %in% target_labels]
message(paste("-> Gekoppelde Huet-zone ID's op basis van labels:", paste(valide_ids, collapse = ", ")))

r_waterlopen_bin <- r_huetzon_local %in% valide_ids
r_waterlopen_bin <- terra::ifel(r_waterlopen_bin == 1, 1, NA)

# ==============================================================================
# PARALLELLE FUSIE: SPOOR A (MAX) EN SPOOR B (OPP)
# ==============================================================================
message("-> Strikt gescheiden parallelle fusie uitvoeren...")

# SPOOR A: MAXIMALE POTENTIE
waterlopen_max <- terra::cover(waterlopen1_bwk_max, r_waterlopen_bin)
names(waterlopen_max) <- "Match_Max"

# SPOOR B: WERKELIJKE OPPERVLAKTE
r_waterlopen_opp <- terra::ifel(r_waterlopen_bin == 1, 0.01, NA)
waterlopen_opp <- terra::cover(waterlopen1_bwk_opp, r_waterlopen_opp)
names(waterlopen_opp) <- "Oppervlakte_Real"

rm(r_huetzon_vlaanderen, r_huetzon_local, r_waterlopen_bin, r_waterlopen_opp)
gc()

# ==============================================================================
# RUIMTELIJKE CLUSTERING & OPPERVLAKTEFILTER WATERLOPEN
# ==============================================================================
message("-> Start ruimtelijke clustering waterlopen...")

w_matrix <- matrix(1, nrow = 3, ncol = 3) 

# SPOOR A: MAX
r_bridge_max <- terra::focal(waterlopen_max, w = w_matrix, fun = "max", na.rm = TRUE)
cl_id_max    <- terra::patches(r_bridge_max, directions = 8, zeroAsNA = TRUE)

habitat_pixels_max <- waterlopen_max * 0.01
stats_max <- terra::zonal(habitat_pixels_max, cl_id_max, fun = "sum", na.rm = TRUE)
colnames(stats_max) <- c("Cluster_ID", "Potentie_ha")

valide_cluster_ids_max <- stats_max$Cluster_ID[stats_max$Potentie_ha >= 0.4]

if(length(valide_cluster_ids_max) > 0) {
  masker_max <- terra::ifel(cl_id_max %in% valide_cluster_ids_max, 1, NA)
  waterlopen_max <- terra::mask(waterlopen_max, masker_max)
} else {
  waterlopen_max <- terra::rast(template_TV, vals = NA)
}
names(waterlopen_max) <- "Leefgebied_Max"

# SPOOR B: OPP
r_binair_water_opp <- terra::ifel(!is.na(waterlopen_opp) & waterlopen_opp > 0, 1, NA)

r_bridge_opp <- terra::focal(r_binair_water_opp, w = w_matrix, fun = "max", na.rm = TRUE)
cl_id_opp    <- terra::patches(r_bridge_opp, directions = 8, zeroAsNA = TRUE)

stats_opp <- terra::zonal(waterlopen_opp, cl_id_opp, fun = "sum", na.rm = TRUE)
colnames(stats_opp) <- c("Cluster_ID", "Werkelijke_ha")

valide_cluster_ids_opp <- stats_opp$Cluster_ID[stats_opp$Werkelijke_ha >= 0.4]

if(length(valide_cluster_ids_opp) > 0) {
  masker_opp <- cl_id_opp %in% valide_cluster_ids_opp
  waterlopen_opp <- terra::mask(waterlopen_opp, terra::ifel(masker_opp, 1, NA))
} else {
  waterlopen_opp <- terra::rast(template_TV, vals = NA)
}
names(waterlopen_opp) <- "Leefgebied_Opp"

# ==============================================================================
# LANDBIOTOOP CLUSTEREN & FILTEREN
# ==============================================================================
message("-> Clusteren en filteren van landbiotoop-lagen...")

bos_res_max <- cluster_filter_compleet(landbiotoop_bos_max, landbiotoop_bos_opp, dist_m = 50, drempel_m2 = 5 * 10000, werkelijk = FALSE)
bos_res_opp <- cluster_filter_compleet(landbiotoop_bos_max, landbiotoop_bos_opp, dist_m = 50, drempel_m2 = 5 * 10000, werkelijk = TRUE)

bos_cluster_max <- bos_res_max$raster
bos_cluster_opp <- bos_res_opp$raster

rest_res_max <- cluster_filter_compleet(landbiotoop_rest_max, opp_laag = landbiotoop_rest_opp, drempel_m2 = 500000, dist_m = 50, werkelijk = FALSE)
rest_res_opp <- cluster_filter_compleet(landbiotoop_rest_opp, opp_laag = landbiotoop_rest_opp, drempel_m2 = 500000, dist_m = 50, werkelijk = TRUE)

rest_cluster_max <- rest_res_max$raster
rest_cluster_opp <- rest_res_opp$raster

filter_bos_in_landmatrice <- function(r_bos_clusters, r_rest_clusters) {
  if (is.null(r_bos_clusters) || is.null(r_rest_clusters)) return(template_TV * NA)
  if (all(is.na(terra::values(r_bos_clusters, mat = FALSE)))) return(template_TV * NA)
  if (all(is.na(terra::values(r_rest_clusters, mat = FALSE)))) return(template_TV * NA)
  
  w_matrix <- matrix(1, nrow = 3, ncol = 3)
  rest_zone <- terra::focal(!is.na(r_rest_clusters) & r_rest_clusters > 0, w = w_matrix, fun = "max", na.rm = TRUE)
  
  r_bos_bin <- terra::ifel(!is.na(r_bos_clusters) & r_bos_clusters > 0, 1, NA)
  cl_bos_id <- terra::patches(r_bos_bin, directions = 8, zeroAsNA = TRUE)
  
  overlap_tabel <- terra::zonal(rest_zone, cl_bos_id, fun = "max", na.rm = TRUE)
  colnames(overlap_tabel) <- c("Bos_Cluster_ID", "Rest_Aanwezig")
  
  geldige_bos_ids <- overlap_tabel$Bos_Cluster_ID[!is.na(overlap_tabel$Rest_Aanwezig) & overlap_tabel$Rest_Aanwezig > 0]
  
  if (length(geldige_bos_ids) == 0) {
    return(template_TV * NA)
  }
  
  r_bos_goedgekeurd <- terra::mask(r_bos_clusters, cl_bos_id %in% geldige_bos_ids, maskvalues = FALSE)
  return(r_bos_goedgekeurd)
}

bos_goedgekeurd_max <- filter_bos_in_landmatrice(bos_cluster_max, rest_cluster_max)
bos_goedgekeurd_opp <- filter_bos_in_landmatrice(bos_cluster_opp, rest_cluster_opp)

# ==============================================================================
# AFSTANDSFILTERS LAND EN WATER
# ==============================================================================
message("-> Wederzijdse afstandsfilters tussen water en boskernen...")

if (terra::hasValues(waterlopen_max) && terra::hasValues(bos_goedgekeurd_max) &&
    terra::global(waterlopen_max, "notNA")$notNA > 0 && 
    terra::global(bos_goedgekeurd_max, "notNA")$notNA > 0) {
  
  r_omgekeerd_bos <- terra::ifel(is.na(bos_goedgekeurd_max), 1, NA)
  dist_naar_bos_max <- terra::distance(r_omgekeerd_bos, target = 1, maxdist = 505)
  zone_bos_500m_max <- terra::ifel(dist_naar_bos_max <= 500, 1, NA)
  
  beekrombout_waterbiotoop_max <- terra::mask(waterlopen_max, zone_bos_500m_max)
  
  if (!all(is.na(terra::values(beekrombout_waterbiotoop_max, mat = FALSE)))) {
    poly_water_max  <- terra::as.polygons(beekrombout_waterbiotoop_max, aggregate = TRUE)
    poly_buffer_max <- terra::buffer(poly_water_max, width = straal_water_m)
    beekrombout_landbiotoop_max <- terra::mask(bos_goedgekeurd_max, poly_buffer_max)
  } else {
    beekrombout_landbiotoop_max <- template_TV * NA
  }
} else {
  beekrombout_waterbiotoop_max <- template_TV * NA
  beekrombout_landbiotoop_max  <- template_TV * NA
}

if (terra::hasValues(waterlopen_opp) && terra::hasValues(bos_goedgekeurd_opp) &&
    terra::global(waterlopen_opp, "notNA")$notNA > 0 && 
    terra::global(bos_goedgekeurd_opp, "notNA")$notNA > 0) {
  
  r_bos_bin <- terra::ifel(!is.na(bos_goedgekeurd_opp) & bos_goedgekeurd_opp > 0, 1, NA)
  r_omgekeerd_bos_opp <- terra::ifel(is.na(r_bos_bin), 1, NA)
  dist_naar_bos_opp <- terra::distance(r_omgekeerd_bos_opp, target = 1, maxdist = 505)
  zone_bos_500m_opp <- terra::ifel(dist_naar_bos_opp <= 500, 1, NA)
  
  beekrombout_waterbiotoop_opp <- terra::mask(waterlopen_opp, zone_bos_500m_opp)
  
  if (!all(is.na(terra::values(beekrombout_waterbiotoop_opp, mat = FALSE)))) {
    r_bin_water_opp <- terra::ifel(!is.na(beekrombout_waterbiotoop_opp) & beekrombout_waterbiotoop_opp > 0, 1, NA)
    poly_water_opp  <- terra::as.polygons(r_bin_water_opp, aggregate = TRUE)
    poly_buffer_opp <- terra::buffer(poly_water_opp, width = straal_water_m)
    beekrombout_landbiotoop_opp <- terra::mask(bos_goedgekeurd_opp, poly_buffer_opp)
  } else {
    beekrombout_landbiotoop_opp <- template_TV * NA
  }
} else {
  beekrombout_waterbiotoop_opp <- template_TV * NA
  beekrombout_landbiotoop_opp  <- template_TV * NA
}

gc()

# ==============================================================================
# DEFINITIEVE SAMENVOEGING LEEFGEBIEDEN
# ==============================================================================
final_max <- terra::cover(beekrombout_waterbiotoop_max, beekrombout_landbiotoop_max)
final_opp <- terra::cover(beekrombout_waterbiotoop_opp, beekrombout_landbiotoop_opp)

if (!all(is.na(terra::values(final_max, mat=FALSE)))) {
  cl_max <- terra::patches(final_max, directions = 8, zeroAsNA = TRUE)
} else {
  cl_max <- template_TV * NA
}

if (!all(is.na(terra::values(final_opp, mat=FALSE)))) {
  cl_opp <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE)
} else {
  cl_opp <- template_TV * NA
}

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
potentie_export_rast <- if (exists("final_max") && !all(is.na(terra::values(final_max, mat=FALSE)))) {
  terra::ifel(!is.na(final_max) & final_max > 0, 1, NA)
} else {
  terra::rast(template_TV, vals = NA)
}

# 2. Bepaal Werkelijke Oppervlakte Raster
werkelijk_export_rast <- if (exists("final_opp") && !all(is.na(terra::values(final_opp, mat=FALSE)))) {
  terra::ifel(!is.na(final_opp) & final_opp > 0, 1, NA)
} else {
  terra::rast(template_TV, vals = NA)
}

# 3. Bepaal Analytisch Metacluster ID-raster
if (exists("cl_opp") && !is.null(cl_opp) && !all(is.na(terra::values(cl_opp, mat=FALSE)))) {
  id_export_rast <- cl_opp
} else if (exists("final_opp") && !is.null(final_opp) && !all(is.na(terra::values(final_opp, mat=FALSE)))) {
  id_export_rast <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE)
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
