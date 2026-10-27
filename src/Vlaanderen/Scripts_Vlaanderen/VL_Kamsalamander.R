# ==============================================================================
# LEEFGEBIEDENGRID VOOR KAMSALAMANDER - HEEL VLAANDEREN (WERKELIJKE OPPERVLAKTE)
# Ecologische analyse & netwerkeisen voor amfibieën (INBO)
# Author: Bert Van Hecke / Update 2026
# ==============================================================================

# 1. BIBLIOTHEKEN LADEN & CONFLICTEN BEHEREN -----------------------------------
library(here)
library(knitr)
library(tidyverse)
library(sf)
library(terra)
library(readxl)
library(tidyterra)
library(data.table)
library(arrow)
library(conflicted)

conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::select)
conflicted::conflicts_prefer(dplyr::first)
conflicted::conflicts_prefer(terra::intersect)
conflicted::conflicts_prefer(terra::any)

# Instellingen voor terra (RAM-geheugen beheer)
terraOptions(
  memfrac = 0.8,       # Max. 80% van het RAM-geheugen gebruiken
  todisk = TRUE,
  tempdir = tempdir(), 
  verbose = FALSE
)

# 2. HELPER FUNCTIES -----------------------------------------------------------
calc_ha_exact <- function(r) {
  if (is.null(r) || all(is.na(suppressWarnings(terra::minmax(r))))) return(0)
  area_raster <- r * terra::cellSize(r, unit = "ha")
  val <- terra::global(area_raster, "sum", na.rm = TRUE)[[1]]
  return(as.numeric(val))
}

cluster_filter_compleet <- function(masker, opp_laag, drempel_m2, dist_m, werkelijk = FALSE) {
  if (all(is.na(suppressWarnings(terra::minmax(masker))))) {
    return(list(raster = masker * NA, clusters = masker * NA))
  }
  
  # 1. Binaire kaart maken
  r_binair <- terra::ifel(!is.na(masker) & masker > 0, 1, NA)
  
  # 2. Netwerkvorming via tijdelijke schijfbestanden
  if (dist_m > 0) {
    stral_cellen <- ceiling((dist_m / 2) / 10)
    f_matrix     <- matrix(1, nrow = (2 * stral_cellen + 1), ncol = (2 * stral_cellen + 1))
    
    tmp_focal  <- tempfile(fileext = ".tif")
    r_buffered <- terra::focal(r_binair, w = f_matrix, fun = "max", na.rm = TRUE, 
                               filename = tmp_focal, overwrite = TRUE)
    r_buffered <- terra::ifel(r_buffered > 0, 1, NA)
    
    tmp_patch  <- tempfile(fileext = ".tif")
    cl_network <- terra::patches(r_buffered, directions = 4, zeroAsNA = TRUE,
                                 filename = tmp_patch, overwrite = TRUE)
    
    cl_biotoop_only <- terra::mask(cl_network, masker)
    
    unlink(c(tmp_focal, tmp_patch))
    rm(r_buffered, cl_network, r_binair)
  } else {
    tmp_patch       <- tempfile(fileext = ".tif")
    cl_biotoop_only <- terra::patches(masker, directions = 8, zeroAsNA = TRUE,
                                      filename = tmp_patch, overwrite = TRUE)
  }
  
  # 3. Slimme extractie: Enkel actieve cellen verwerken
  df_cl <- terra::as.data.frame(cl_biotoop_only, cells = TRUE)
  if (nrow(df_cl) == 0) return(list(raster = masker * NA, clusters = masker * NA))
  colnames(df_cl) <- c("cell", "ID")
  
  if (werkelijk) {
    opp_waarden <- terra::extract(opp_laag, df_cl$cell)[[1]]
    dt_calc <- data.table(ID = df_cl$ID, Waarde = opp_waarden)
  } else {
    dt_calc <- data.table(ID = df_cl$ID, Waarde = 1)
  }
  
  stats_dt <- dt_calc[!is.na(ID) & !is.na(Waarde), .(Area_m2 = sum(Waarde, na.rm = TRUE) * 100), by = ID]
  voldoet_ids <- stats_dt[Area_m2 >= drempel_m2, ID]
  
  if (length(voldoet_ids) == 0) return(list(raster = masker * NA, clusters = masker * NA))
  
  # 4. Resultaat filteren
  cl_finaal <- terra::ifel(cl_biotoop_only %in% voldoet_ids, cl_biotoop_only, NA)
  r_finaal  <- terra::mask(masker, cl_finaal)
  
  rm(df_cl, dt_calc, stats_dt, cl_biotoop_only)
  gc()
  
  return(list(raster = r_finaal, clusters = cl_finaal))
}

# 3. SOORT PARAMETERS OPHALEN --------------------------------------------------
soort         <- "kamsalamander"

opp_water_ha  <- 0.01   # MinOpp waterbiotoop (0,01 ha = 100 m2)
opp_land_ha   <- 5.0    # MinOpp landbiotoop per eiland (5 ha)
buffer_m      <- 1000   # Dispersiecapaciteit (1000 m)
min_totale_ha <- 14     # Minimale totale netwerkoppervlakte (14 ha)
min_plassen   <- 5      # Minimaal aantal geschikte plassen in het netwerk
afstand_m     <- 0      # Afstand land tot water (0 m)

# 4. GEBIED VOORBEREIDEN (MASTERGRID VLAANDEREN) -------------------------------
message("-> Mastergrid Vlaanderen inlezen...")
master_grid <- rast(here("data/input/Raster_Vlaanderen/Vlaanderen_MasterGrid_10m.tif"))[[1]]

template_Vlaanderen <- master_grid
values(template_Vlaanderen) <- NA

# 5. BIOTOOPFILTERING VIA PARQUET (BWK TYPETABELLEN) ---------------------------
message("-> BWK Parquet dataset koppelen via Arrow...")

df_nieuw <- read_csv(
  here("data/input/Excel_files/Resultaten_Totaal_Samengevoegd.csv"), 
  show_col_types = FALSE,
  progress = FALSE
)

resultaten_gegroepeerd <- df_nieuw %>%
  mutate(Soort_clean = tolower(trimws(Soort))) %>%
  filter(Soort_clean == soort) %>%
  group_by(Type) %>%
  nest(Data = c(Code, Match))

ds_vlaanderen <- open_dataset(here("data/input/Raster_Vlaanderen/BWK_TidyTabel_Smal_Vlaanderen_2025.parquet"))

lijst_oppervlaktes <- list()

for(i in 1:nrow(resultaten_gegroepeerd)) {
  h_type <- resultaten_gegroepeerd$Type[i]
  h_data <- resultaten_gegroepeerd$Data[[i]]
  
  exact_codes <- tolower(trimws(h_data$Code[h_data$Match == "exact"]))
  bevat_codes <- tolower(trimws(h_data$Code[h_data$Match == "bevat"]))
  
  if (length(bevat_codes) > 0) {
    bevat_codes_escaped  <- gsub("([\\.\\^\\$\\*\\+\\?\\(\\)\\[\\{\\\\\\|])", "\\\\\\1", bevat_codes)
    bevat_codes_anchored <- paste0("^", bevat_codes_escaped)
    regex_term           <- paste0(bevat_codes_anchored, collapse = "|")
    
    tabel_sub <- ds_vlaanderen %>%
      filter(tolower(CODE) %in% exact_codes | str_detect(tolower(CODE), regex_term)) %>%
      collect() %>%
      as.data.table()
  } else {
    tabel_sub <- ds_vlaanderen %>%
      filter(tolower(CODE) %in% exact_codes) %>%
      collect() %>%
      as.data.table()
  }
  
  if (nrow(tabel_sub) > 0) {
    tabel_sub[, CODE := tolower(trimws(CODE))]
    tabel_unique  <- unique(tabel_sub, by = c("cel_id", "CODE"))
    tabel_cel_som <- tabel_unique[, .(Oppervlakte = pmin(sum(BWK_FRAC, na.rm = TRUE), 1.0)), by = .(cel_id)]
    
    r_opp_type <- template_Vlaanderen
    r_opp_type[tabel_cel_som$cel_id] <- tabel_cel_som$Oppervlakte
    lijst_oppervlaktes[[h_type]]     <- r_opp_type
  } else {
    lijst_oppervlaktes[[h_type]]     <- template_Vlaanderen
  }
}

waterbiotoop_bwk_opp <- if("waterbiotoop_bwk" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["waterbiotoop_bwk"]] else template_Vlaanderen
landbiotoop_bwk_opp  <- if("landbiotoop_bwk" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["landbiotoop_bwk"]] else template_Vlaanderen

rm(ds_vlaanderen, resultaten_gegroepeerd, df_nieuw, lijst_oppervlaktes)
gc()

# 6. WATERVLAKKEN SHAPEFILE INLEZEN EN VERWERKEN -------------------------------
message("-> Watervlakken shapefile inlezen voor heel Vlaanderen...")

watervlakken_path <- here("data/input/ASCI Files/watervlakken2024.shp")

if (file.exists(watervlakken_path)) {
  watervlakken_v     <- vect(watervlakken_path)
  watervlakken_proj  <- project(watervlakken_v, crs(template_Vlaanderen))
  watervlakken_clean <- makeValid(watervlakken_proj)
  
  watervlakken_buffered       <- buffer(watervlakken_clean, width = 5)
  watervlakken_fuzzy_clusters <- aggregate(watervlakken_buffered, by = NULL)
  watervlakken_fuzzy_clusters$fuzzy_area_m2 <- terra::expanse(watervlakken_fuzzy_clusters)
  
  kleine_clusters_v <- watervlakken_fuzzy_clusters[watervlakken_fuzzy_clusters$fuzzy_area_m2 <= 4000, ]
  
  if (!is.null(kleine_clusters_v) && nrow(kleine_clusters_v) > 0) {
    kleine_watervlakken_finaal <- crop(watervlakken_clean, kleine_clusters_v)
    
    if (!is.null(kleine_watervlakken_finaal) && nrow(kleine_watervlakken_finaal) > 0) {
      r_ks_kleinwater <- terra::rasterize(kleine_watervlakken_finaal, template_Vlaanderen, field = 1, background = NA)
      r_kamsalamander_waterbiotoop1_opp <- terra::cover(waterbiotoop_bwk_opp, r_ks_kleinwater)
    } else {
      r_kamsalamander_waterbiotoop1_opp <- waterbiotoop_bwk_opp
    }
  } else {
    r_kamsalamander_waterbiotoop1_opp <- waterbiotoop_bwk_opp
  }
} else {
  r_kamsalamander_waterbiotoop1_opp <- waterbiotoop_bwk_opp
  watervlakken_clean                <- NULL
}

suppressWarnings(
  rm(watervlakken_v, watervlakken_buffered, watervlakken_fuzzy_clusters, kleine_clusters_v, kleine_watervlakken_finaal, r_ks_kleinwater)
)
gc()

# 7. LANDBIOTOOP FILTEREN OP MINIMALE EILANDOPPERVLAKTE (>= 5 HA) --------------
message("-> Landbiotoop filteren op minimale eiland-oppervlakte (>= 5 ha)...")

drempel_land_m2 <- opp_land_ha * 10000 # 5 ha = 50.000 m²

if (!is.null(landbiotoop_bwk_opp) && !all(is.na(suppressWarnings(terra::minmax(landbiotoop_bwk_opp))))) {
  r_binair_land_opp <- terra::ifel(!is.na(landbiotoop_bwk_opp) & landbiotoop_bwk_opp > 0, 1, NA)
  
  tmp_patch_land   <- tempfile(fileext = ".tif")
  land_patches_opp <- terra::patches(r_binair_land_opp, directions = 8, zeroAsNA = TRUE,
                                     filename = tmp_patch_land, overwrite = TRUE)
  
  zonal_land_opp <- terra::zonal(landbiotoop_bwk_opp, land_patches_opp, fun = "sum", na.rm = TRUE)
  colnames(zonal_land_opp) <- c("ID", "Werk_ha")
  ids_land_ok_opp <- zonal_land_opp$ID[(zonal_land_opp$Werk_ha * 100) >= drempel_land_m2]
  
  unlink(tmp_patch_land)
  
  if (length(ids_land_ok_opp) > 0) {
    kamsalamander_land_gefilterd_opp <- terra::mask(landbiotoop_bwk_opp, land_patches_opp %in% ids_land_ok_opp)
  } else {
    kamsalamander_land_gefilterd_opp <- template_Vlaanderen
  }
} else {
  kamsalamander_land_gefilterd_opp <- template_Vlaanderen
}

# 8. RUIMTELIJKE KOPPELING LAND- EN WATERBIOTOOP ------------------------------
message("-> Ruimtelijke koppeling land- en waterbiotoop...")

if (!is.null(kamsalamander_land_gefilterd_opp) && !all(is.na(suppressWarnings(terra::minmax(kamsalamander_land_gefilterd_opp)))) && 
    !is.null(r_kamsalamander_waterbiotoop1_opp) && !all(is.na(suppressWarnings(terra::minmax(r_kamsalamander_waterbiotoop1_opp))))) {
  
  r_land_mask_opp  <- terra::ifel(!is.na(kamsalamander_land_gefilterd_opp) & kamsalamander_land_gefilterd_opp > 0, 1, NA)
  r_water_mask_opp <- terra::ifel(!is.na(r_kamsalamander_waterbiotoop1_opp) & r_kamsalamander_waterbiotoop1_opp > 0, 1, NA)
  
  ks_land_gekoppeld_opp  <- terra::mask(kamsalamander_land_gefilterd_opp, r_water_mask_opp)
  ks_water_gekoppeld_opp <- terra::mask(r_kamsalamander_waterbiotoop1_opp, r_land_mask_opp)
  
  water_opp_clean <- terra::ifel(is.na(ks_water_gekoppeld_opp), 0, ks_water_gekoppeld_opp)
  land_opp_clean  <- terra::ifel(is.na(ks_land_gekoppeld_opp), 0, ks_land_gekoppeld_opp)
  
  som_opp <- terra::clamp(water_opp_clean + land_opp_clean, upper = 1.0)
  kamsalamander_leefgebied1_opp <- terra::ifel(som_opp > 0, som_opp, NA)
} else {
  kamsalamander_leefgebied1_opp <- template_Vlaanderen
}

# 9. NETWERKVORMING BINNEN 1000M EN PLASSEN-DREMPEL (>= 14 HA EN >= 5 PLASSEN)
message("-> Netwerkvorming binnen 1000m en controle op >= 14 ha EN >= 5 plassen...")

drempel_netwerk_m2 <- min_totale_ha * 10000 # 14 ha = 140.000 m²
r_binair_leef1_opp <- terra::ifel(!is.na(kamsalamander_leefgebied1_opp) & kamsalamander_leefgebied1_opp > 0, 1, NA)

res_netwerk_opp <- cluster_filter_compleet(
  masker     = r_binair_leef1_opp,
  opp_laag   = kamsalamander_leefgebied1_opp,
  drempel_m2 = drempel_netwerk_m2,
  dist_m     = buffer_m,
  werkelijk  = TRUE
)

cl_netwerk_opp <- res_netwerk_opp$clusters
r_netwerk_opp  <- res_netwerk_opp$raster

if (!all(is.na(suppressWarnings(terra::minmax(cl_netwerk_opp))))) {
  netwerk_poly_opp <- terra::as.polygons(cl_netwerk_opp, dissolve = TRUE) %>% terra::makeValid()
  watervlakken_in_netwerk_opp <- if(!is.null(watervlakken_clean)) terra::crop(watervlakken_clean, netwerk_poly_opp) else NULL
  
  if (!is.null(watervlakken_in_netwerk_opp) && nrow(watervlakken_in_netwerk_opp) > 0) {
    rel_opp <- matrix(terra::is.related(netwerk_poly_opp, watervlakken_in_netwerk_opp, "intersects"), 
                      nrow = nrow(netwerk_poly_opp), ncol = nrow(watervlakken_in_netwerk_opp))
    netwerk_poly_opp$Aantal_Plassen <- rowSums(rel_opp)
    id_col_opp    <- names(netwerk_poly_opp)[1]
    valid_ids_opp <- as.numeric(netwerk_poly_opp[[id_col_opp]][netwerk_poly_opp$Aantal_Plassen >= min_plassen, 1])
  } else {
    valid_ids_opp <- c()
  }
  
  if (length(na.omit(valid_ids_opp)) > 0) {
    final_opp <- terra::mask(kamsalamander_leefgebied1_opp, cl_netwerk_opp %in% valid_ids_opp)
    cl_opp    <- terra::mask(cl_netwerk_opp, cl_netwerk_opp %in% valid_ids_opp)
  } else {
    final_opp <- template_Vlaanderen
    cl_opp    <- template_Vlaanderen
  }
} else {
  final_opp <- template_Vlaanderen
  cl_opp    <- template_Vlaanderen
}

rm(kamsalamander_leefgebied1_opp, r_binair_leef1_opp, res_netwerk_opp, watervlakken_clean)
gc()

# 10. EXPORT NAAR DISK (ENKEL 00_ID_Rasters EN 01_Werkelijke_Oppervlaktes) -----
message("-> Geformatteerde export starten...")

base_dir <- here::here("data/output/Vlaanderen/Rasters_Soorten")

folders <- list(
  id_raster = file.path(base_dir, "00_ID_Rasters"),
  werkelijk = file.path(base_dir, "01_Werkelijke_Oppervlaktes")
)
purrr::walk(folders, ~if (!dir.exists(.x)) dir.create(.x, showWarnings = FALSE, recursive = TRUE))

# 1. Export Werkelijke Oppervlakte Raster
file_path_opp <- file.path(folders$werkelijk, paste0("Habitat_Werkelijke_Oppervlaktes_", soort, ".tif"))
terra::writeRaster(
  final_opp, 
  filename = file_path_opp, 
  overwrite = TRUE, 
  gdal = c("COMPRESS=LZW"), 
  datatype = "FLT4S", 
  NAflag = -9999
)
message(paste("    [OK] Geëxporteerd:", basename(file_path_opp)))

# 2. Export Analytisch Netwerk ID Raster
file_path_id <- file.path(folders$id_raster, paste0("ID_Netwerken_", soort, ".tif"))
terra::writeRaster(
  cl_opp, 
  filename = file_path_id, 
  overwrite = TRUE, 
  gdal = c("COMPRESS=LZW"), 
  datatype = "INT4U", 
  NAflag = 0
)
message(paste("    [OK] Analytisch ID-raster geëxporteerd:", basename(file_path_id)))

message("🏁 LEEFGEBIEDENGRID KAMSALAMANDER SUCCESVOL GEËXPORTEERD VOOR HEEL VLAANDEREN!")
