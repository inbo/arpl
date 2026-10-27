# ==============================================================================
# LEEFGEBIEDENGRID VOOR HOOGVEENGLANSLIBEL - HEEL VLAANDEREN (WERKELIJKE OPPERVLAKTE)
# Ecologische analyse water, bosrijk landbiotoop & ruime omgevingsfactor (INBO)
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
  if (is.null(r)) return(0)
  if (all(is.na(terra::values(r, mat = FALSE)))) return(0)
  area_raster <- r * terra::cellSize(r, unit = "ha")
  val <- terra::global(area_raster, "sum", na.rm = TRUE)[[1]]
  return(as.numeric(val))
}

cluster_filter_compleet <- function(masker, opp_laag, drempel_m2, dist_m, werkelijk = FALSE) {
  if (terra::global(is.na(masker), "sum")[[1]] == terra::ncell(masker)) {
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

# 3. SOORT INFORMATIE OPHALEN --------------------------------------------------
soort <- "hoogveenglanslibel"

buffer_m        <- 5000   # Dispersiecapaciteit voor smear/buffer (5000m)
straal_water_m  <- 500    # Harde afsnij-afstand rond het water (500m)

drempel_water_m2 <- 100    # Min. 0.01 ha water (100 m²)
dist_water_m     <- 50     # 50m overbrugging voor water

drempel_land_m2  <- 50000  # Min. 5 ha land
dist_land_m      <- 10     # 10m overbrugging voor aanpalend land

drempel_ruim_m2  <- 500000 # Min. 50 ha ruim land
dist_ruim_m      <- 50     # 50m overbrugging voor ruim land

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

water1_opp    <- if("waterbiotoop_bwk1" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["waterbiotoop_bwk1"]] else template_Vlaanderen
not_water_opp <- if("waterbiotoop_notbwk" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["waterbiotoop_notbwk"]] else template_Vlaanderen
land1_opp     <- if("landbiotoop_bwk1" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["landbiotoop_bwk1"]] else template_Vlaanderen
not_land_opp  <- if("landbiotoop_notbwk" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["landbiotoop_notbwk"]] else template_Vlaanderen
bos_opp       <- if("landbiotoop_bos" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["landbiotoop_bos"]] else template_Vlaanderen
land_ruim_opp <- if("landbiotoop_bwk_ruim" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["landbiotoop_bwk_ruim"]] else template_Vlaanderen

rm(ds_vlaanderen, resultaten_gegroepeerd, df_nieuw, lijst_oppervlaktes)
gc()

# 6. WATERBIOTOOP EN LANDBIOTOOP NETWERKEN CLUSTEREN --------------------------
message("-> Waterbiotoop zuiveren en clusteren (0.01 ha, 50m netwerk)...")

r_not_water_zero <- terra::ifel(is.na(not_water_opp), 0, not_water_opp)
r_water_rest     <- terra::clamp(water1_opp - r_not_water_zero, lower = 0, upper = 1)
waterbiotoop_opp <- terra::ifel(r_water_rest > 0, r_water_rest, NA)

rm(r_not_water_zero, r_water_rest, water1_opp, not_water_opp)
gc()

r_binair_water_opp        <- terra::ifel(!is.na(waterbiotoop_opp) & waterbiotoop_opp > 0, 1, NA)
waterbiotoop_clusters_opp <- cluster_filter_compleet(
  masker     = r_binair_water_opp,
  opp_laag   = waterbiotoop_opp,
  drempel_m2 = drempel_water_m2, # 100 m²
  dist_m     = dist_water_m,    # 50m
  werkelijk  = TRUE
)

rm(r_binair_water_opp)
gc()

message("-> Landbiotoop zuiveren en clusteren (5 ha, 10m netwerk)...")

r_not_land_zero <- terra::ifel(is.na(not_land_opp), 0, not_land_opp)
r_land_rest     <- terra::clamp(land1_opp - r_not_land_zero, lower = 0, upper = 1)
landbiotoop_opp <- terra::ifel(r_land_rest > 0, r_land_rest, NA)

rm(r_not_land_zero, r_land_rest, land1_opp, not_land_opp)
gc()

r_binair_land_opp        <- terra::ifel(!is.na(landbiotoop_opp) & landbiotoop_opp > 0, 1, NA)
landbiotoop_clusters_opp <- cluster_filter_compleet(
  masker     = r_binair_land_opp,
  opp_laag   = landbiotoop_opp,
  drempel_m2 = drempel_land_m2, # 5 ha
  dist_m     = dist_land_m,     # 10m
  werkelijk  = TRUE
)

rm(r_binair_land_opp)
gc()

message("-> Ruim landbiotoop clusteren (50 ha, 50m netwerk)...")

r_binair_ruim_opp        <- terra::ifel(!is.na(land_ruim_opp) & land_ruim_opp > 0, 1, NA)
ruim_land_clusters_opp <- cluster_filter_compleet(
  masker     = r_binair_ruim_opp,
  opp_laag   = land_ruim_opp,
  drempel_m2 = drempel_ruim_m2, # 50 ha
  dist_m     = dist_ruim_m,     # 50m
  werkelijk  = TRUE
)

rm(r_binair_ruim_opp, land_ruim_opp)
gc()

# 7. BOSPATCHES EN BOSPERCENTAGES PER CLUSTER (MIN. 60% BOS) --------------------
message("-> Samenstellingsanalyse landbiotoopclusters (min. 60% bos)...")

r_land_zuiver_opp <- landbiotoop_clusters_opp$raster
r_binair_land_opp <- terra::ifel(!is.na(r_land_zuiver_opp) & r_land_zuiver_opp > 0, 1, NA)

if (!all(is.na(suppressWarnings(terra::minmax(r_binair_land_opp))))) {
  tmp_patch_l    <- tempfile(fileext = ".tif")
  cl_id_land_opp <- terra::patches(r_binair_land_opp, directions = 8, zeroAsNA = TRUE, filename = tmp_patch_l, overwrite = TRUE)
  
  totale_opp_per_patch_opp <- terra::zonal(landbiotoop_opp, cl_id_land_opp, fun = "sum", na.rm = TRUE)
  colnames(totale_opp_per_patch_opp) <- c("Cluster_ID", "Totaal_Land_ha_sum")
  totale_opp_per_patch_opp$Totaal_Land_ha <- totale_opp_per_patch_opp$Totaal_Land_ha_sum * 0.01
  
  bos_in_land_opp       <- terra::mask(bos_opp, cl_id_land_opp)
  stats_bos_in_land_opp <- terra::zonal(bos_in_land_opp, cl_id_land_opp, fun = "sum", na.rm = TRUE)
  colnames(stats_bos_in_land_opp) <- c("Cluster_ID", "Bos_ha_sum")
  stats_bos_in_land_opp$Bos_ha    <- stats_bos_in_land_opp$Bos_ha_sum * 0.01
  
  df_samenstelling_land_opp <- merge(totale_opp_per_patch_opp, stats_bos_in_land_opp, by = "Cluster_ID", all.x = TRUE)
  df_samenstelling_land_opp$Bos_ha[is.na(df_samenstelling_land_opp$Bos_ha)] <- 0
  df_samenstelling_land_opp$Percentage_Bos   <- round((df_samenstelling_land_opp$Bos_ha / df_samenstelling_land_opp$Totaal_Land_ha) * 100, 2)
  
  valide_land_ids_opp <- df_samenstelling_land_opp$Cluster_ID[df_samenstelling_land_opp$Percentage_Bos >= 60]
  
  if(length(valide_land_ids_opp) > 0) {
    winnende_masker_land_opp <- cl_id_land_opp %in% valide_land_ids_opp
    land_nabij_opp           <- terra::mask(landbiotoop_opp, winnende_masker_land_opp)
  } else {
    land_nabij_opp           <- template_Vlaanderen
  }
  unlink(tmp_patch_l)
} else {
  land_nabij_opp             <- template_Vlaanderen
}

rm(r_land_zuiver_opp, r_binair_land_opp, bos_in_land_opp, stats_bos_in_land_opp, landbiotoop_opp, bos_opp, landbiotoop_clusters_opp)
gc()

# 8. WATER-FIRST LANDBIOTOOP FILTERING (500M AFSNIJDING) ----------------------
message("-> Landbiotoop (min. 60% bos) hard afsnijden op ", straal_water_m, "m van geschikt water...")

if (exists("waterbiotoop_clusters_opp") && !is.null(waterbiotoop_clusters_opp$raster) &&
    !all(is.na(suppressWarnings(terra::minmax(waterbiotoop_clusters_opp$raster))))) {
  
  r_water_src_opp <- waterbiotoop_clusters_opp$raster
  r_bin_water_opp <- terra::ifel(!is.na(r_water_src_opp) & r_water_src_opp > 0, 1, NA)
  poly_water_opp  <- terra::as.polygons(r_bin_water_opp, aggregate = TRUE)
  poly_buffer_opp <- terra::buffer(poly_water_opp, width = straal_water_m)
  
  if (exists("land_nabij_opp") && !all(is.na(suppressWarnings(terra::minmax(land_nabij_opp))))) {
    landbiotoop2_opp  <- terra::mask(land_nabij_opp, poly_buffer_opp)
    waterbiotoop2_opp <- r_water_src_opp
  } else {
    landbiotoop2_opp  <- template_Vlaanderen
    waterbiotoop2_opp <- template_Vlaanderen
  }
} else {
  landbiotoop2_opp  <- template_Vlaanderen
  waterbiotoop2_opp <- template_Vlaanderen
}

rm(waterbiotoop_clusters_opp)
gc()

# 9. RUIME OMGEVING EN FINALE SAMENSMELTING ------------------------------------
message("-> Ruime omgevingsfactor (50 ha complex binnen 500m) toepassen op water...")

r_ruim_land_src_opp   <- ruim_land_clusters_opp$raster
r_binair_ruim_opp     <- terra::ifel(!is.na(r_ruim_land_src_opp) & r_ruim_land_src_opp > 0, 1, NA)

if (!all(is.na(suppressWarnings(terra::minmax(r_binair_ruim_opp))))) {
  tmp_dist_ruim         <- tempfile(fileext = ".tif")
  dist_to_ruim_land_opp <- terra::distance(r_binair_ruim_opp, filename = tmp_dist_ruim, overwrite = TRUE)
  
  m_water_ruim_opp                    <- dist_to_ruim_land_opp <= 500
  hoogveenglanslibel_waterbiotoop_opp <- terra::mask(waterbiotoop2_opp, terra::ifel(m_water_ruim_opp, 1, NA))
  unlink(tmp_dist_ruim)
} else {
  hoogveenglanslibel_waterbiotoop_opp <- template_Vlaanderen
}

hoogveenglanslibel_landbiotoop_opp  <- landbiotoop2_opp

w_clean <- terra::ifel(is.na(hoogveenglanslibel_waterbiotoop_opp), 0, hoogveenglanslibel_waterbiotoop_opp)
l_clean <- terra::ifel(is.na(hoogveenglanslibel_landbiotoop_opp), 0, hoogveenglanslibel_landbiotoop_opp)

som_fracties <- w_clean + l_clean
final_opp_raw <- terra::clamp(som_fracties, upper = 1.0)
final_opp     <- terra::ifel(final_opp_raw > 0, final_opp_raw, NA)

if (!all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  tmp_patch_cl <- tempfile(fileext = ".tif")
  cl_opp       <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE, filename = tmp_patch_cl, overwrite = TRUE)
} else {
  cl_opp <- template_Vlaanderen
}

rm(ruim_land_clusters_opp, r_ruim_land_src_opp, r_binair_ruim_opp, waterbiotoop2_opp, landbiotoop2_opp,
   hoogveenglanslibel_waterbiotoop_opp, hoogveenglanslibel_landbiotoop_opp, w_clean, l_clean, som_fracties, final_opp_raw)
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

message("🏁 LEEFGEBIEDENGRID HOOGVEENGLANSLIBEL SUCCESVOL GEËXPORTEERD VOOR HEEL VLAANDEREN!")
