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

soort           <- "boomkikker"

# --- DYNAMISCHE SCENARIO PARAMETER CHECK ---
if (exists("SCENARIO_RDS_PAD") && !is.null(SCENARIO_RDS_PAD)) {
  scenario_rds_path <- SCENARIO_RDS_PAD
} else if (exists("params") && !is.null(params$scenario_rds_path)) {
  scenario_rds_path <- params$scenario_rds_path
} else {
  scenario_rds_path <- "data/input/Scenario_rds/DM_Scenario_BWK_2025.rds"
}

p_raw <- gsub("^([.][.]/)+", "", scenario_rds_path)
scenario_path <- here::here(p_raw)

if (!file.exists(scenario_path)) {
  stop(paste("❌ FOUT: Scenario RDS bestand NIET gevonden op:", scenario_path))
}

scen_volledig <- basename(scenario_path)
scenario_naam <- gsub("^DM_Scenario_|^Scenario_|.rds$", "", scen_volledig)

message(paste("Verwerken van soort:", soort, "binnen scenario:", scenario_naam))

# Boomkikker expert-parameters
opp_water_ha   <- 0.01   # MinOpp waterbiotoop (0,01 ha = 100 m2)
opp_land_ha    <- 5.0    # MinOpp landbiotoop per eiland (5 ha)
buffer_m       <- 4000   # Dispersiecapaciteit (4000 m)
min_totale_ha  <- 160    # Minimale totale netwerkoppervlakte (160 ha)
min_plassen    <- 10     # Minimaal aantal geschikte plassen in het netwerk
afstand_m      <- 0      # Afstand land tot water (0 m)

area_shape  <- vect(here("data/input/De_Maten.shp"))
master_grid <- rast(here("data/input/Raster_Vlaanderen/Vlaanderen_MasterGrid_10m.tif"))[[1]]

df_namen_sleutel <- read_csv(here("data/input/Excel_files/BWK_Laag_Namen_2025.csv"), show_col_types = FALSE)
gouden_namenlijst <- tolower(trimws(df_namen_sleutel$Laagnaam))

area_shape_proj <- project(area_shape, crs(master_grid))
area_buffer_fix <- buffer(area_shape_proj, width = buffer_m)

message("-> Vertaalraster voor globale/lokale cellen opbouwen via snelle MASK methode...")
id_raster_DM <- crop(master_grid, area_buffer_fix, snap = "near")

globale_id_raster <- master_grid
globale_id_raster <- terra::init(globale_id_raster, fun = "cell")

id_raster_DM_globale_values <- crop(globale_id_raster, area_buffer_fix, snap = "near")
id_raster_DM_masked <- mask(id_raster_DM_globale_values, area_buffer_fix)

message("-> Vertaaltabel bliksemsnel opbouwen via C++ dataframe extractie...")
df_extractie <- as.data.frame(id_raster_DM_masked, cells = TRUE)
vertaal_df <- as.data.table(df_extractie)
setnames(vertaal_df, c(1, 2), c("lokale_id", "globale_id"))

vertaal_df <- vertaal_df[!is.na(globale_id)]
studiegebied_globale_ids <- unique(vertaal_df$globale_id)

values(id_raster_DM) <- NA
template_DM <- terra::rasterize(area_buffer_fix, id_raster_DM, field = 1, background = 0)

rm(globale_id_raster, id_raster_DM_globale_values, id_raster_DM_masked, df_extractie)
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
  
  tabel_DM <- tabel_gefilterd[cel_id %in% studiegebied_globale_ids]
  tabel_DM_unique <- unique(tabel_DM, by = c("cel_id", "CODE"))

  tabel_cel_som <- tabel_DM_unique[, .(Oppervlakte = pmin(sum(BWK_FRAC, na.rm = TRUE), 1.0)), by = .(cel_id)]
  
  r_match_type <- id_raster_DM * NA
  r_opp_type   <- id_raster_DM * NA
  
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

waterbiotoop_bwk_max <- lijst_matches[["waterbiotoop_bwk"]]
waterbiotoop_bwk_opp <- lijst_oppervlaktes[["waterbiotoop_bwk"]]
landbiotoop_hms_max  <- lijst_matches[["landbiotoop_hms"]]
landbiotoop_hms_opp  <- lijst_oppervlaktes[["landbiotoop_hms"]]
landbiotoop_k_max    <- lijst_matches[["landbiotoop_k"]]
landbiotoop_k_opp    <- lijst_oppervlaktes[["landbiotoop_k"]]
akkers_max           <- lijst_matches[["akkers"]]
akkers_opp           <- lijst_oppervlaktes[["akkers"]]

raster_simpel_final  <- id_raster_DM

rm(tabel_vlaanderen, vertaal_df, lijst_matches, lijst_oppervlaktes)
gc()

message("-> Watervlakken shapefile inlezen en verwerken...")

watervlakken_v    <- vect(here("data/input/ASCI Files/watervlakken2024.shp"))
watervlakken_proj <- project(watervlakken_v, crs(template_DM))

area_buffer_clean  <- aggregate(makeValid(area_buffer_fix))
watervlakken_clean <- makeValid(watervlakken_proj)
watervlakken_DM    <- crop(watervlakken_clean, area_buffer_clean)

# Fuzzy clustering voor kleine waters
watervlakken_buffered       <- buffer(watervlakken_DM, width = 5)
watervlakken_fuzzy_clusters <- aggregate(watervlakken_buffered, by = NULL)
watervlakken_fuzzy_clusters$fuzzy_area_m2 <- terra::expanse(watervlakken_fuzzy_clusters)

kleine_clusters_v <- watervlakken_fuzzy_clusters[watervlakken_fuzzy_clusters$fuzzy_area_m2 <= 4000, ]

if (!is.null(kleine_clusters_v) && nrow(kleine_clusters_v) > 0) {
  kleine_watervlakken_finaal <- crop(watervlakken_DM, kleine_clusters_v)
  
  if (!is.null(kleine_watervlakken_finaal) && nrow(kleine_watervlakken_finaal) > 0) {
    r_boomkikker_kleinwater <- terra::rasterize(kleine_watervlakken_finaal, template_DM, field = 1, background = NA)
    
    r_boomkikker_waterbiotoop1     <- terra::cover(waterbiotoop_bwk_max, r_boomkikker_kleinwater)
    r_boomkikker_waterbiotoop1_opp <- terra::cover(waterbiotoop_bwk_opp, r_boomkikker_kleinwater)
  } else {
    r_boomkikker_waterbiotoop1     <- waterbiotoop_bwk_max
    r_boomkikker_waterbiotoop1_opp <- waterbiotoop_bwk_opp
  }
} else {
  r_boomkikker_waterbiotoop1     <- waterbiotoop_bwk_max
  r_boomkikker_waterbiotoop1_opp <- waterbiotoop_bwk_opp
}

suppressWarnings(
  rm(watervlakken_v, watervlakken_buffered, watervlakken_fuzzy_clusters, kleine_clusters_v, kleine_watervlakken_finaal, r_boomkikker_kleinwater)
)
gc()

message("-> Landbiotoop opbouwen en akkers uitsluiten...")

r_hms_zuiver_max <- landbiotoop_hms_max
r_hms_zuiver_max[is.na(r_hms_zuiver_max) | r_hms_zuiver_max == 0] <- NA

r_k_zuiver_max <- landbiotoop_k_max
r_k_zuiver_max[is.na(r_k_zuiver_max) | r_k_zuiver_max == 0] <- NA

dist_k   <- terra::distance(r_k_zuiver_max)
dist_hms <- terra::distance(r_hms_zuiver_max)

landbiotoop_bwk1_max <- r_hms_zuiver_max == 1 & dist_k <= 500
landbiotoop_bwk2_max <- r_k_zuiver_max == 1   & dist_hms <= 500

landbiotoop_bwk3_max <- landbiotoop_bwk1_max | landbiotoop_bwk2_max
landbiotoop_bwk3_max <- terra::ifel(landbiotoop_bwk3_max == 1, 1, NA)

landbiotoop_bwk1_opp <- terra::mask(landbiotoop_hms_opp, landbiotoop_bwk1_max)
landbiotoop_bwk2_opp <- terra::mask(landbiotoop_k_opp, landbiotoop_bwk2_max)
landbiotoop_bwk3_opp <- terra::cover(landbiotoop_bwk1_opp, landbiotoop_bwk2_opp)

r_akkers_zuiver <- akkers_max
r_akkers_zuiver[is.na(r_akkers_zuiver) | r_akkers_zuiver == 0] <- NA

landbiotoop_na_max <- terra::mask(landbiotoop_bwk3_max, r_akkers_zuiver, inverse = TRUE)
landbiotoop_na_opp <- terra::mask(landbiotoop_bwk3_opp, r_akkers_zuiver, inverse = TRUE)

suppressWarnings(
  rm(dist_k, dist_hms, r_hms_zuiver_max, r_k_zuiver_max, r_akkers_zuiver,
     landbiotoop_bwk1_max, landbiotoop_bwk2_max, landbiotoop_bwk3_max,
     landbiotoop_bwk1_opp, landbiotoop_bwk2_opp, landbiotoop_bwk3_opp)
)
gc()

message("-> Landbiotoop filteren op minimale eiland-oppervlakte (>= 5 ha) voor MAX en OPP parallel...")

drempel_land_m2 <- opp_land_ha * 10000 # 5 ha = 50.000 m²

# --- SPOOR MAX ---
if (!all(is.na(suppressWarnings(terra::minmax(landbiotoop_na_max))))) {
  land_patches_max <- terra::patches(landbiotoop_na_max, directions = 8, zeroAsNA = TRUE)
  land_freq_max    <- as.data.frame(terra::freq(land_patches_max))
  ids_land_ok_max  <- land_freq_max$value[(land_freq_max$count * 100) >= drempel_land_m2]
  
  if (length(ids_land_ok_max) > 0) {
    boomkikker_land_gefilterd_max <- terra::mask(landbiotoop_na_max, land_patches_max %in% ids_land_ok_max)
  } else {
    boomkikker_land_gefilterd_max <- template_DM * NA
  }
} else {
  boomkikker_land_gefilterd_max <- template_DM * NA
}

# --- SPOOR OPP ---
if (!all(is.na(suppressWarnings(terra::minmax(landbiotoop_na_opp))))) {
  r_binair_land_opp <- terra::ifel(!is.na(landbiotoop_na_opp) & landbiotoop_na_opp > 0, 1, NA)
  land_patches_opp  <- terra::patches(r_binair_land_opp, directions = 8, zeroAsNA = TRUE)
  
  zonal_land_opp <- terra::zonal(landbiotoop_na_opp, land_patches_opp, fun = "sum", na.rm = TRUE)
  colnames(zonal_land_opp) <- c("ID", "Werk_ha")
  ids_land_ok_opp <- zonal_land_opp$ID[(zonal_land_opp$Werk_ha * 100) >= drempel_land_m2]
  
  if (length(ids_land_ok_opp) > 0) {
    boomkikker_land_gefilterd_opp <- terra::mask(landbiotoop_na_opp, land_patches_opp %in% ids_land_ok_opp)
  } else {
    boomkikker_land_gefilterd_opp <- template_DM * NA
  }
} else {
  boomkikker_land_gefilterd_opp <- template_DM * NA
}

message("-> Ruimtelijke koppeling land- en waterbiotoop (AfstandBiotopen = 0 m) voor MAX en OPP parallel...")

# --- SPOOR MAX ---
if (!all(is.na(suppressWarnings(terra::minmax(boomkikker_land_gefilterd_max)))) && 
    !all(is.na(suppressWarnings(terra::minmax(r_boomkikker_waterbiotoop1))))) {
  
  r_land_mask_max  <- boomkikker_land_gefilterd_max
  r_water_mask_max <- r_boomkikker_waterbiotoop1
  
  bk_land_gekoppeld_max  <- terra::mask(boomkikker_land_gefilterd_max, r_water_mask_max)
  bk_water_gekoppeld_max <- terra::mask(r_boomkikker_waterbiotoop1, r_land_mask_max)
  
  boomkikker_leefgebied1_max <- bk_land_gekoppeld_max | bk_water_gekoppeld_max
  boomkikker_leefgebied1_max <- terra::ifel(boomkikker_leefgebied1_max == 1, 1, NA)
} else {
  boomkikker_leefgebied1_max <- template_DM * NA
}

# --- SPOOR OPP ---
if (!all(is.na(suppressWarnings(terra::minmax(boomkikker_land_gefilterd_opp)))) && 
    !all(is.na(suppressWarnings(terra::minmax(r_boomkikker_waterbiotoop1_opp))))) {
  
  r_land_mask_opp  <- terra::ifel(!is.na(boomkikker_land_gefilterd_opp) & boomkikker_land_gefilterd_opp > 0, 1, NA)
  r_water_mask_opp <- terra::ifel(!is.na(r_boomkikker_waterbiotoop1_opp) & r_boomkikker_waterbiotoop1_opp > 0, 1, NA)
  
  bk_land_gekoppeld_opp  <- terra::mask(boomkikker_land_gefilterd_opp, r_water_mask_opp)
  bk_water_gekoppeld_opp <- terra::mask(r_boomkikker_waterbiotoop1_opp, r_land_mask_opp)
  
  water_opp_clean <- terra::ifel(is.na(bk_water_gekoppeld_opp), 0, bk_water_gekoppeld_opp)
  land_opp_clean  <- terra::ifel(is.na(bk_land_gekoppeld_opp), 0, bk_land_gekoppeld_opp)

  som_opp <- terra::clamp(water_opp_clean + land_opp_clean, upper = 1.0)
  boomkikker_leefgebied1_opp <- terra::ifel(som_opp > 0, som_opp, NA)
} else {
  boomkikker_leefgebied1_opp <- template_DM * NA
}

message("-> Netwerkvorming binnen 4000m en controle op >= 160 ha EN >= 10 plassen...")

drempel_netwerk_m2 <- min_totale_ha * 10000 # 160 ha = 1.600.000 m²

# --- SPOOR MAX ---
res_netwerk_max <- cluster_filter_compleet(
  masker     = boomkikker_leefgebied1_max,
  opp_laag   = boomkikker_leefgebied1_max,
  drempel_m2 = drempel_netwerk_m2,
  dist_m     = buffer_m,
  werkelijk  = FALSE
)

cl_netwerk_max <- res_netwerk_max$clusters
r_netwerk_max  <- res_netwerk_max$raster

if (!all(is.na(suppressWarnings(terra::minmax(cl_netwerk_max))))) {
  netwerk_poly_max <- terra::as.polygons(cl_netwerk_max, dissolve = TRUE) %>% terra::makeValid()
  watervlakken_in_netwerk_max <- terra::crop(watervlakken_clean, netwerk_poly_max)
  
  if (!is.null(watervlakken_in_netwerk_max) && nrow(watervlakken_in_netwerk_max) > 0) {
    rel_max <- matrix(terra::is.related(netwerk_poly_max, watervlakken_in_netwerk_max, "intersects"), 
                      nrow = nrow(netwerk_poly_max), ncol = nrow(watervlakken_in_netwerk_max))
    netwerk_poly_max$Aantal_Plassen <- rowSums(rel_max)
    id_col_max <- names(netwerk_poly_max)[1]
    valid_ids_max <- as.numeric(netwerk_poly_max[[id_col_max]][netwerk_poly_max$Aantal_Plassen >= min_plassen, 1])
  } else {
    valid_ids_max <- c()
  }
  
  if (length(na.omit(valid_ids_max)) > 0) {
    boomkikker_leefgebied_max <- terra::mask(r_netwerk_max, cl_netwerk_max %in% valid_ids_max)
    cl_max                    <- terra::mask(cl_netwerk_max, cl_netwerk_max %in% valid_ids_max)
  } else {
    boomkikker_leefgebied_max <- template_DM * NA
    cl_max                    <- template_DM * NA
  }
} else {
  boomkikker_leefgebied_max <- template_DM * NA
  cl_max                    <- template_DM * NA
}

# --- SPOOR OPP ---
r_binair_leef1_opp <- terra::ifel(!is.na(boomkikker_leefgebied1_opp) & boomkikker_leefgebied1_opp > 0, 1, NA)

res_netwerk_opp <- cluster_filter_compleet(
  masker     = r_binair_leef1_opp,
  opp_laag   = boomkikker_leefgebied1_opp,
  drempel_m2 = drempel_netwerk_m2,
  dist_m     = buffer_m,
  werkelijk  = TRUE
)

cl_netwerk_opp <- res_netwerk_opp$clusters
r_netwerk_opp  <- res_netwerk_opp$raster

if (!all(is.na(suppressWarnings(terra::minmax(cl_netwerk_opp))))) {
  netwerk_poly_opp <- terra::as.polygons(cl_netwerk_opp, dissolve = TRUE) %>% terra::makeValid()
  watervlakken_in_netwerk_opp <- terra::crop(watervlakken_clean, netwerk_poly_opp)
  
  if (!is.null(watervlakken_in_netwerk_opp) && nrow(watervlakken_in_netwerk_opp) > 0) {
    rel_opp <- matrix(terra::is.related(netwerk_poly_opp, watervlakken_in_netwerk_opp, "intersects"), 
                      nrow = nrow(netwerk_poly_opp), ncol = nrow(watervlakken_in_netwerk_opp))
    netwerk_poly_opp$Aantal_Plassen <- rowSums(rel_opp)
    id_col_opp <- names(netwerk_poly_opp)[1]
    valid_ids_opp <- as.numeric(netwerk_poly_opp[[id_col_opp]][netwerk_poly_opp$Aantal_Plassen >= min_plassen, 1])
  } else {
    valid_ids_opp <- c()
  }
  
  if (length(na.omit(valid_ids_opp)) > 0) {
    boomkikker_leefgebied_opp <- terra::mask(boomkikker_leefgebied1_opp, cl_netwerk_opp %in% valid_ids_opp)
    cl_opp                    <- terra::mask(cl_netwerk_opp, cl_netwerk_opp %in% valid_ids_opp)
  } else {
    boomkikker_leefgebied_opp <- template_DM * NA
    cl_opp                    <- template_DM * NA
  }
} else {
  boomkikker_leefgebied_opp <- template_DM * NA
  cl_opp                    <- template_DM * NA
}

# ==============================================================================
# DEFINITIEVE TOEWIJSING EXPORT VARIABELEN (CRUCIALE FIX)
# ==============================================================================
final_max <- boomkikker_leefgebied_max
final_opp <- boomkikker_leefgebied_opp

if (is.null(cl_max) || all(is.na(terra::values(cl_max, mat=FALSE)))) {
  if (!is.null(final_max) && !all(is.na(suppressWarnings(terra::minmax(final_max))))) {
    cl_max <- terra::patches(final_max, directions = 8, zeroAsNA = TRUE)
  } else {
    cl_max <- template_DM * NA
  }
}

if (is.null(cl_opp) || all(is.na(terra::values(cl_opp, mat=FALSE)))) {
  if (!is.null(final_opp) && !all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
    cl_opp <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE)
  } else {
    cl_opp <- template_DM * NA
  }
}

# ==============================================================================
# SCHONE EXPORT BIOTOOP EN ANALYTISCH ID-RASTER (VOOR SCRIPT 2 / ARPL)
# ==============================================================================
base_dir <- here::here("data/output/De_Maten/Rasters_Soorten", scenario_naam)

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
  terra::rast(template_DM, vals = NA)
}

# 2. Bepaal Werkelijke Oppervlakte Raster
werkelijk_export_rast <- if (exists("final_opp") && !is.null(final_opp) && !all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  terra::ifel(!is.na(final_opp) & final_opp > 0, 1, NA)
} else {
  terra::rast(template_DM, vals = NA)
}

# 3. Bepaal Analytisch Metacluster ID-raster
# 3. Bepaal Analytisch Metacluster ID-raster (ENKEL OP BASIS VAN WERKELIJKE OPPERVLAKTES)
if (exists("cl_opp") && !is.null(cl_opp) && !all(is.na(suppressWarnings(terra::minmax(cl_opp))))) {
  id_export_rast <- cl_opp
} else if (exists("final_opp") && !is.null(final_opp) && !all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  id_export_rast <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE)
} else {
  id_export_rast <- terra::rast(template_DM, vals = NA)
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
