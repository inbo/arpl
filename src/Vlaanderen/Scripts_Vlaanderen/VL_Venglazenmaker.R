# ==============================================================================
# LEEFGEBIEDENGRID VOOR VENGLAZENMAKER - HEEL VLAANDEREN (WERKELIJKE OPPERVLAKTE)
# Ecologische analyse water-first landkoppeling & bospercentage-eisen (INBO)
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
soort          <- "venglazenmaker"
straal_water_m <- 500 # Harde afsnij-afstand van land rond het water (500m)

df <- read_excel(here::here("data/input/Excel_files/Soorten_bwk_afstanden.xlsx"))
resultaat <- df %>%
  filter(tolower(trimws(Soort)) == soort) %>%
  select(Type, MinOpp_ha, AfstandBiotopen_m, Dispersiecap_m)

buffer_m <- resultaat$Dispersiecap_m[1]

print(resultaat)
rm(df, resultaat)

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

waterbiotoop_bwk1_opp     <- if("waterbiotoop_bwk1" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["waterbiotoop_bwk1"]] else template_Vlaanderen
waterbiotoop_notbwk_opp   <- if("waterbiotoop_notbwk" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["waterbiotoop_notbwk"]] else template_Vlaanderen
landbiotoop_bwk1_opp      <- if("landbiotoop_bwk1" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["landbiotoop_bwk1"]] else template_Vlaanderen
landbiotoop_notbwk_opp    <- if("landbiotoop_notbwk" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["landbiotoop_notbwk"]] else template_Vlaanderen
landbiotoop_nabij_bos_opp <- if("landbiotoop_nabij_bos" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["landbiotoop_nabij_bos"]] else template_Vlaanderen
landbiotoop_bwk_ruim_opp  <- if("landbiotoop_bwk_ruim" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["landbiotoop_bwk_ruim"]] else template_Vlaanderen
landbiotoop_bos_ruim_opp  <- if("landbiotoop_bos_ruim" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["landbiotoop_bos_ruim"]] else template_Vlaanderen

rm(ds_vlaanderen, resultaten_gegroepeerd, df_nieuw, lijst_oppervlaktes)
gc()

# 6. WATERBIOTOOP BEREKENEN & CLUSTEREN (MIN. 0.04 HA / 400 M²) ----------------
message("-> Waterbiotoop filteren en clusteren (0.04 ha, 50m netwerk)...")

r_not_water_zero <- terra::ifel(is.na(waterbiotoop_notbwk_opp), 0, waterbiotoop_notbwk_opp)
r_water_rest     <- terra::clamp(waterbiotoop_bwk1_opp - r_not_water_zero, lower = 0, upper = 1)
venglazenmaker_water_opp <- terra::ifel(r_water_rest > 0, r_water_rest, NA)

rm(r_not_water_zero, r_water_rest, waterbiotoop_bwk1_opp, waterbiotoop_notbwk_opp)
gc()

r_binair_water_opp        <- terra::ifel(!is.na(venglazenmaker_water_opp) & venglazenmaker_water_opp > 0, 1, NA)
waterbiotoop_clusters_opp <- cluster_filter_compleet(
  masker     = r_binair_water_opp,
  opp_laag   = venglazenmaker_water_opp,
  drempel_m2 = 400, # 0.04 ha
  dist_m     = 50, 
  werkelijk  = TRUE
)
r_water1_opp_clean <- waterbiotoop_clusters_opp$raster

rm(r_binair_water_opp, venglazenmaker_water_opp, waterbiotoop_clusters_opp)
gc()

# 7. LANDBIOTOOP NABIJE OMGEVING BEREKENEN (MIN. 10 HA & >=40% BOS) ------------
message("-> Landbiotoop nabije omgeving filteren (10 ha, 50m netwerk, >=40% bos)...")

r_not_land_zero       <- terra::ifel(is.na(landbiotoop_notbwk_opp), 0, landbiotoop_notbwk_opp)
r_land_rest           <- terra::clamp(landbiotoop_bwk1_opp - r_not_land_zero, lower = 0, upper = 1)
land_nabij_filter_opp <- terra::ifel(r_land_rest > 0, r_land_rest, NA)

rm(r_not_land_zero, r_land_rest, landbiotoop_bwk1_opp, landbiotoop_notbwk_opp)
gc()

r_binair_land_nabij_opp <- terra::ifel(!is.na(land_nabij_filter_opp) & land_nabij_filter_opp > 0, 1, NA)

land_nabij_clusters_opp <- cluster_filter_compleet(
  masker     = r_binair_land_nabij_opp,
  opp_laag   = land_nabij_filter_opp,
  drempel_m2 = 100000, # 10 ha
  dist_m     = 50, 
  werkelijk  = TRUE
)

r_nabij_patch_ids_opp          <- land_nabij_clusters_opp$clusters
venglazenmaker_land_nabij_opp <- template_Vlaanderen
goedgekeurde_patch_ids_nabij_opp <- c()

if (!all(is.na(suppressWarnings(terra::minmax(r_nabij_patch_ids_opp))))) {
  totale_opp_per_patch_nabij_opp <- terra::zonal(land_nabij_filter_opp, r_nabij_patch_ids_opp, fun = "sum", na.rm = TRUE)
  colnames(totale_opp_per_patch_nabij_opp) <- c("Patch_ID", "Totaal_ha_sum")
  totale_opp_per_patch_nabij_opp$Totaal_ha <- totale_opp_per_patch_nabij_opp$Totaal_ha_sum * 0.01
  
  bos_in_zone_opp             <- terra::mask(landbiotoop_nabij_bos_opp, r_nabij_patch_ids_opp)
  bos_opp_per_patch_nabij_opp <- terra::zonal(bos_in_zone_opp, r_nabij_patch_ids_opp, fun = "sum", na.rm = TRUE)
  colnames(bos_opp_per_patch_nabij_opp) <- c("Patch_ID", "Bos_ha_sum")
  bos_opp_per_patch_nabij_opp$Bos_ha    <- bos_opp_per_patch_nabij_opp$Bos_ha_sum * 0.01
  
  dt_tot_nabij <- as.data.table(totale_opp_per_patch_nabij_opp)
  dt_bos_nabij <- as.data.table(bos_opp_per_patch_nabij_opp)
  
  patch_stats_nabij_opp <- dt_bos_nabij[dt_tot_nabij, on = "Patch_ID"]
  patch_stats_nabij_opp[is.na(Bos_ha), Bos_ha := 0]
  patch_stats_nabij_opp[, Percentage_Bos := round((Bos_ha / Totaal_ha) * 100, 2)]
  
  goedgekeurde_patch_ids_nabij_opp <- patch_stats_nabij_opp[Percentage_Bos >= 40, Patch_ID]
  
  if (length(goedgekeurde_patch_ids_nabij_opp) > 0) {
    venglazenmaker_land_nabij_opp <- terra::mask(land_nabij_filter_opp, r_nabij_patch_ids_opp %in% goedgekeurde_patch_ids_nabij_opp)
  }
  rm(totale_opp_per_patch_nabij_opp, bos_in_zone_opp, bos_opp_per_patch_nabij_opp, dt_tot_nabij, dt_bos_nabij, patch_stats_nabij_opp)
}

rm(r_binair_land_nabij_opp, land_nabij_clusters_opp, land_nabij_bos_opp)
gc()

# 8. RUIM LANDBIOTOOP BEREKENEN (MIN. 30 HA & >=40% BOS) ------------------------
message("-> Ruim landbiotoop filteren (30 ha, 50m netwerk, >=40% bos)...")

r_binair_ruim_opp_init  <- terra::ifel(!is.na(landbiotoop_bwk_ruim_opp) & landbiotoop_bwk_ruim_opp > 0, 1, NA)
land_ruim_clusters_opp  <- cluster_filter_compleet(
  masker     = r_binair_ruim_opp_init,
  opp_laag   = landbiotoop_bwk_ruim_opp,
  drempel_m2 = 300000, # 30 ha
  dist_m     = 50, 
  werkelijk  = TRUE
)
r_patch_ids_opp              <- land_ruim_clusters_opp$clusters
venglazenmaker_land_ruim_opp <- template_Vlaanderen

if (!all(is.na(suppressWarnings(terra::minmax(r_patch_ids_opp))))) {
  totale_opp_per_patch_opp <- terra::zonal(landbiotoop_bwk_ruim_opp, r_patch_ids_opp, fun = "sum", na.rm = TRUE)
  colnames(totale_opp_per_patch_opp) <- c("Patch_ID", "Totaal_ha_sum")
  totale_opp_per_patch_opp$Totaal_ha <- totale_opp_per_patch_opp$Totaal_ha_sum * 0.01
  
  bos_in_ruim_opp       <- terra::mask(landbiotoop_bos_ruim_opp, r_patch_ids_opp)
  bos_opp_per_patch_opp <- terra::zonal(bos_in_ruim_opp, r_patch_ids_opp, fun = "sum", na.rm = TRUE)
  colnames(bos_opp_per_patch_opp) <- c("Patch_ID", "Bos_ha_sum")
  bos_opp_per_patch_opp$Bos_ha    <- bos_opp_per_patch_opp$Bos_ha_sum * 0.01
  
  dt_tot_ruim <- as.data.table(totale_opp_per_patch_opp)
  dt_bos_ruim <- as.data.table(bos_opp_per_patch_opp)
  
  patch_stats_opp <- dt_bos_ruim[dt_tot_ruim, on = "Patch_ID"]
  patch_stats_opp[is.na(Bos_ha), Bos_ha := 0]
  patch_stats_opp[, Percentage_Bos := (Bos_ha / Totaal_ha) * 100]
  
  goedgekeurde_patch_ids_opp <- patch_stats_opp[Percentage_Bos >= 40, Patch_ID]
  
  if (length(goedgekeurde_patch_ids_opp) > 0) {
    venglazenmaker_land_ruim_opp <- terra::mask(landbiotoop_bwk_ruim_opp, r_patch_ids_opp %in% goedgekeurde_patch_ids_opp)
  }
  rm(totale_opp_per_patch_opp, bos_in_ruim_opp, bos_opp_per_patch_opp, dt_tot_ruim, dt_bos_ruim, patch_stats_opp)
}

rm(r_binair_ruim_opp_init, land_ruim_clusters_opp, r_patch_ids_opp, landbiotoop_bwk_ruim_opp, landbiotoop_bos_ruim_opp)
gc()

# 9. RUIMTELIJKE INTERACTIES & WATER-FIRST AFSNIJDING (500M) --------------------
message("-> Ruimtelijke interacties (20m water/land koppeling, 500m omgevingsbuffer) & Water-First afsnijding op 500m...")

r_water1_opp_bin    <- terra::ifel(!is.na(r_water1_opp_clean) & r_water1_opp_clean > 0, 1, NA)

tmp_buf_w20         <- tempfile(fileext = ".tif")
water_buf20_opp_bin <- terra::buffer(r_water1_opp_bin, width = 20, filename = tmp_buf_w20, overwrite = TRUE)

if (exists("r_nabij_patch_ids_opp") && !all(is.na(suppressWarnings(terra::minmax(r_nabij_patch_ids_opp))))) {
  r_land_in_water_buf_opp     <- terra::mask(r_nabij_patch_ids_opp, water_buf20_opp_bin)
  df_w_land                   <- terra::as.data.frame(r_land_in_water_buf_opp, na.rm = TRUE)
  winnende_land_patch_ids_opp <- unique(df_w_land[[1]])
  actieve_land_ids_opp        <- intersect(winnende_land_patch_ids_opp, goedgekeurde_patch_ids_nabij_opp)
} else { actieve_land_ids_opp <- c() }

unlink(tmp_buf_w20)

if (length(actieve_land_ids_opp) > 0) {
  venglazenmaker_landbiotoop_opp_ruw <- terra::mask(land_nabij_filter_opp, r_nabij_patch_ids_opp %in% actieve_land_ids_opp)
} else { 
  venglazenmaker_landbiotoop_opp_ruw <- template_Vlaanderen 
}

venglazenmaker_land_nabij_opp_bin <- terra::ifel(!is.na(venglazenmaker_land_nabij_opp) & venglazenmaker_land_nabij_opp > 0, 1, NA)

tmp_buf_l20        <- tempfile(fileext = ".tif")
land_buf20_opp_bin <- terra::buffer(venglazenmaker_land_nabij_opp_bin, width = 20, filename = tmp_buf_l20, overwrite = TRUE)

venglazenmaker_waterbiotoop_opp <- terra::mask(r_water1_opp_clean, land_buf20_opp_bin)

unlink(tmp_buf_l20)

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

w_clean <- terra::ifel(is.na(waterbiotoop_finaal_opp), 0, waterbiotoop_finaal_opp)
l_clean <- terra::ifel(is.na(landbiotoop_finaal_opp), 0, landbiotoop_finaal_opp)

final_opp_som <- w_clean + l_clean
final_opp_raw <- terra::clamp(final_opp_som, upper = 1.0)
final_opp     <- terra::ifel(final_opp_raw > 0, final_opp_raw, NA)

if (!all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  tmp_patch_cl <- tempfile(fileext = ".tif")
  cl_opp       <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE, filename = tmp_patch_cl, overwrite = TRUE)
} else {
  cl_opp <- template_Vlaanderen
}

rm(r_water1_opp_clean, venglazenmaker_land_nabij_opp, land_nabij_filter_opp, venglazenmaker_landbiotoop_opp_ruw, 
   witsnuit_water_finaal_opp, w_clean, l_clean, final_opp_som, final_opp_raw)
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

message("🏁 LEEFGEBIEDENGRID VENGLAZENMAKER SUCCESVOL GEËXPORTEERD VOOR HEEL VLAANDEREN!")
