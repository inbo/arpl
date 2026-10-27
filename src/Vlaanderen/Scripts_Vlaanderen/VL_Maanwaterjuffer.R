# ==============================================================================
# LEEFGEBIEDENGRID VOOR MAANWATERJUFFER - HEEL VLAANDEREN (WERKELIJKE OPPERVLAKTE)
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

fix_bos_raster <- function(r_bos, template) {
  if (is.null(r_bos) || !inherits(r_bos, "SpatRaster")) {
    return(terra::rast(template, vals = NA))
  }
  return(terra::ifel(!is.na(r_bos) & r_bos > 0, 1, NA))
}

# 3. SOORT INFORMATIE OPHALEN --------------------------------------------------
soort          <- "maanwaterjuffer"
straal_water_m <- 500

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
    
    min_drempel  <- if(grepl("waterbiotoop_bwk1", h_type, ignore.case = TRUE)) 0.20 else 0.01
    tabel_cel_som <- tabel_unique[, .(Oppervlakte = pmin(sum(BWK_FRAC, na.rm = TRUE), 1.0)), by = .(cel_id)][Oppervlakte >= min_drempel]
    
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
landbiotoop_bos_nabij_opp <- if("landbiotoop_nabij_bos" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["landbiotoop_nabij_bos"]] else template_Vlaanderen
landbiotoop_bwk_ruim_opp  <- if("landbiotoop_bwk_ruim" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["landbiotoop_bwk_ruim"]] else template_Vlaanderen
landbiotoop_bos_ruim_opp  <- if("landbiotoop_ruim_bos" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["landbiotoop_ruim_bos"]] else template_Vlaanderen

rm(ds_vlaanderen, resultaten_gegroepeerd, df_nieuw, lijst_oppervlaktes)
gc()

# 6. WATERBIOTOOP BEREKENEN & CLUSTEREN (MIN. 0.02 HA) ------------------------
message("-> Waterbiotoop filteren en clusteren (0.02 ha, 10m netwerk)...")

r_not_water_zero <- terra::ifel(is.na(waterbiotoop_notbwk_opp), 0, waterbiotoop_notbwk_opp)
r_water_rest     <- terra::clamp(waterbiotoop_bwk1_opp - r_not_water_zero, lower = 0, upper = 1)
waterbiotoop_bwk_basis_opp <- terra::ifel(r_water_rest > 0, r_water_rest, NA)

rm(r_not_water_zero, r_water_rest, waterbiotoop_bwk1_opp, waterbiotoop_notbwk_opp)
gc()

r_binair_water_opp     <- terra::ifel(!is.na(waterbiotoop_bwk_basis_opp) & waterbiotoop_bwk_basis_opp > 0, 1, NA)
waterbiotoop1_opp_list <- cluster_filter_compleet(
  masker     = r_binair_water_opp, 
  opp_laag   = waterbiotoop_bwk_basis_opp, 
  drempel_m2 = 200, # 0.02 ha
  dist_m     = 10, 
  werkelijk  = TRUE
)
waterbiotoop1_opp <- waterbiotoop1_opp_list$raster

rm(r_binair_water_opp, waterbiotoop_bwk_basis_opp, waterbiotoop1_opp_list)
gc()

# 7. RUIME OMGEVING LANDBIOTOOP BEREKENEN (MIN. 50 HA & >50% BOS/STRUWEEL) ------
message("-> Landbiotoop ruimere omgeving berekenen (min. 50 ha & >50% bos/struweel)...")

r_bos_ruim_opp_clean <- fix_bos_raster(landbiotoop_bos_ruim_opp, template_Vlaanderen)

r_binair_ruim_opp         <- terra::ifel(!is.na(landbiotoop_bwk_ruim_opp) & landbiotoop_bwk_ruim_opp > 0, 1, NA)
landbiotoop_ruim_opp_list <- cluster_filter_compleet(
  masker     = r_binair_ruim_opp, 
  opp_laag   = landbiotoop_bwk_ruim_opp, 
  drempel_m2 = 500000, # 50 ha
  dist_m     = 50, 
  werkelijk  = TRUE
)
r_ruim_opp_src <- landbiotoop_ruim_opp_list$raster

if(!all(is.na(suppressWarnings(terra::minmax(r_ruim_opp_src))))) {
  tmp_patch_r    <- tempfile(fileext = ".tif")
  cl_id_ruim_opp <- terra::patches(r_ruim_opp_src, directions = 8, zeroAsNA = TRUE, filename = tmp_patch_r, overwrite = TRUE)
  
  stats_ruim_totaal_opp <- terra::zonal(r_ruim_opp_src, cl_id_ruim_opp, fun = "sum", na.rm = TRUE)
  colnames(stats_ruim_totaal_opp) <- c("Cluster_ID", "Totaal_Ruim_ha_sum")
  stats_ruim_totaal_opp$Totaal_Ruim_ha <- stats_ruim_totaal_opp$Totaal_Ruim_ha_sum * 0.01
  
  bos_in_ruim_opp       <- terra::mask(r_bos_ruim_opp_clean, cl_id_ruim_opp)
  stats_bos_in_ruim_opp <- terra::zonal(bos_in_ruim_opp, cl_id_ruim_opp, fun = "sum", na.rm = TRUE)
  colnames(stats_bos_in_ruim_opp) <- c("Cluster_ID", "Bos_ha_sum")
  stats_bos_in_ruim_opp$Bos_ha    <- stats_bos_in_ruim_opp$Bos_ha_sum * 0.01
  
  df_ruim_samenstelling_opp <- merge(stats_ruim_totaal_opp, stats_bos_in_ruim_opp, by = "Cluster_ID", all.x = TRUE)
  df_ruim_samenstelling_opp$Bos_ha[is.na(df_ruim_samenstelling_opp$Bos_ha)] <- 0
  df_ruim_samenstelling_opp$Percentage_Bos   <- (df_ruim_samenstelling_opp$Bos_ha / df_ruim_samenstelling_opp$Totaal_Ruim_ha) * 100
  
  valide_ruim_ids_opp <- df_ruim_samenstelling_opp$Cluster_ID[df_ruim_samenstelling_opp$Percentage_Bos > 50]
  
  if(length(valide_ruim_ids_opp) > 0) {
    masker_ruim_opp     <- cl_id_ruim_opp %in% valide_ruim_ids_opp
    landbiotoop_ruim_opp <- terra::mask(r_ruim_opp_src, terra::ifel(masker_ruim_opp, 1, NA))
  } else {
    landbiotoop_ruim_opp <- template_Vlaanderen
  }
  unlink(tmp_patch_r)
} else {
  landbiotoop_ruim_opp   <- template_Vlaanderen
}

rm(r_binair_ruim_opp, landbiotoop_ruim_opp_list, r_ruim_opp_src, landbiotoop_bwk_ruim_opp, landbiotoop_bos_ruim_opp, r_bos_ruim_opp_clean)
gc()

# 8. WATER-FIRST FILTERING & LANDBIOTOOP CLUSTERING (10 HA & >=40% BOS) --------
message("-> Water-first filtering: landbiotoop afsnijden op 500m van geschikt water...")

if (exists("waterbiotoop1_opp") && !all(is.na(suppressWarnings(terra::minmax(waterbiotoop1_opp))))) {
  r_water_opp_echt <- terra::ifel(!is.na(waterbiotoop1_opp) & waterbiotoop1_opp > 0, waterbiotoop1_opp, NA)
  r_water_opp_bin  <- terra::ifel(!is.na(r_water_opp_echt) & r_water_opp_echt > 0, 1, NA)
} else {
  r_water_opp_echt <- template_Vlaanderen
  r_water_opp_bin  <- template_Vlaanderen
}

if (!all(is.na(suppressWarnings(terra::minmax(r_water_opp_bin))))) {
  poly_water_opp  <- terra::as.polygons(r_water_opp_bin, aggregate = TRUE)
  poly_buffer_opp <- terra::buffer(poly_water_opp, width = straal_water_m)
  
  r_not_land_zero  <- terra::ifel(is.na(landbiotoop_notbwk_opp), 0, landbiotoop_notbwk_opp)
  r_land_rest      <- terra::clamp(landbiotoop_bwk1_opp - r_not_land_zero, lower = 0, upper = 1)
  r_land_opp_basis <- terra::ifel(r_land_rest > 0, r_land_rest, NA)
  
  land_afgesneden_500m_opp <- terra::mask(r_land_opp_basis, poly_buffer_opp)
  r_bin_land_afgesneden    <- terra::ifel(!is.na(land_afgesneden_500m_opp) & land_afgesneden_500m_opp > 0, 1, NA)
  
  res_land_10ha_opp <- cluster_filter_compleet(
    masker     = r_bin_land_afgesneden, 
    opp_laag   = land_afgesneden_500m_opp, 
    drempel_m2 = 100000, # 10 ha
    dist_m     = 20, 
    werkelijk  = TRUE
  )
  land_10ha_opp  <- res_land_10ha_opp$raster
  cl_id_10ha_opp <- res_land_10ha_opp$clusters
  
  if (!all(is.na(suppressWarnings(terra::minmax(land_10ha_opp))))) {
    r_bos_nabij_opp_clean <- fix_bos_raster(landbiotoop_bos_nabij_opp, template_Vlaanderen)
    stats_totaal_opp      <- terra::zonal(land_10ha_opp, cl_id_10ha_opp, fun = "sum", na.rm = TRUE)
    colnames(stats_totaal_opp) <- c("Cluster_ID", "Totaal_m2")
    
    bos_in_cl_opp         <- terra::mask(r_bos_nabij_opp_clean, cl_id_10ha_opp)
    stats_bos_opp         <- terra::zonal(bos_in_cl_opp, cl_id_10ha_opp, fun = "sum", na.rm = TRUE)
    colnames(stats_bos_opp)    <- c("Cluster_ID", "Bos_m2")
    
    df_bos_opp            <- merge(stats_totaal_opp, stats_bos_opp, by = "Cluster_ID", all.x = TRUE)
    df_bos_opp$Bos_m2[is.na(df_bos_opp$Bos_m2)] <- 0
    df_bos_opp$Pct_Bos    <- (df_bos_opp$Bos_m2 / df_bos_opp$Totaal_m2) * 100
    valide_bos_ids_opp    <- df_bos_opp$Cluster_ID[df_bos_opp$Pct_Bos >= 40]
    
    land_finaal_opp       <- if(length(valide_bos_ids_opp) > 0) terra::mask(land_10ha_opp, terra::ifel(cl_id_10ha_opp %in% valide_bos_ids_opp, 1, NA)) else template_Vlaanderen
  } else { land_finaal_opp <- template_Vlaanderen }
  
} else {
  land_finaal_opp <- template_Vlaanderen
}

water_finaal_opp <- r_water_opp_echt

w_clean <- terra::ifel(is.na(water_finaal_opp), 0, water_finaal_opp)
l_clean <- terra::ifel(is.na(land_finaal_opp), 0, land_finaal_opp)

som_fracties  <- w_clean + l_clean
final_opp_raw <- terra::clamp(som_fracties, upper = 1.0)
final_opp     <- terra::ifel(final_opp_raw > 0, final_opp_raw, NA)

if (!all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  tmp_patch_cl <- tempfile(fileext = ".tif")
  cl_opp       <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE, filename = tmp_patch_cl, overwrite = TRUE)
} else {
  cl_opp <- template_Vlaanderen
}

rm(waterbiotoop1_opp, landbiotoop_bwk1_opp, landbiotoop_notbwk_opp, landbiotoop_bos_nabij_opp, w_clean, l_clean, som_fracties, final_opp_raw)
gc()

# 9. EXPORT NAAR DISK (ENKEL 00_ID_Rasters EN 01_Werkelijke_Oppervlaktes) -----
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

message("🏁 LEEFGEBIEDENGRID MAANWATERJUFFER SUCCESVOL GEËXPORTEERD VOOR HEEL VLAANDEREN!")
