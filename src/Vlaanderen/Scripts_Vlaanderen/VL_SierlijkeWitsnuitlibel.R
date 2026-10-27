# ==============================================================================
# LEEFGEBIEDENGRID VOOR SIERLIJKE WITSNUITLIBEL - HEEL VLAANDEREN (WERKELIJKE OPPERVLAKTE)
# Ecologische analyse water-first landkoppeling, netwerkconnectiviteit & struweelaandeel (INBO)
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
soort          <- "sierlijkewitsnuitlibel"
straal_water_m <- 800  # Harde afsnij-afstand van land rond het water (800m)

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

voortplanting_bwk_opp <- if("voortplanting_bwk" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["voortplanting_bwk"]] else template_Vlaanderen
foerageer0_notbwk_opp <- if("foerageer0" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["foerageer0"]] else template_Vlaanderen
struweel1_opp         <- if("struweel1" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["struweel1"]] else template_Vlaanderen

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
      r_kleinwater                 <- terra::rasterize(kleine_watervlakken_finaal, template_Vlaanderen, field = 1, background = NA)
      r_witsnuit_waterbiotoop1_opp <- terra::cover(voortplanting_bwk_opp, r_kleinwater)
    } else {
      r_witsnuit_waterbiotoop1_opp <- voortplanting_bwk_opp
    }
  } else {
    r_witsnuit_waterbiotoop1_opp <- voortplanting_bwk_opp
  }
} else {
  r_witsnuit_waterbiotoop1_opp <- voortplanting_bwk_opp
  watervlakken_clean           <- NULL
}

suppressWarnings(
  rm(watervlakken_v, watervlakken_buffered, watervlakken_fuzzy_clusters, kleine_clusters_v, kleine_watervlakken_finaal, r_kleinwater)
)
gc()

# 7. STAP 2: VOORTPLANTINGSBIOTOOP FILTEREN OP CONNECTIVITEIT (>= 4 PLASSEN) -----
message("-> Stap 2: Voortplantingsbiotoop filteren op netwerkconnectiviteit (4 plassen binnen 500m)...")

r_binair_water_opp   <- terra::ifel(!is.na(r_witsnuit_waterbiotoop1_opp) & r_witsnuit_waterbiotoop1_opp > 0, 1, NA)

tmp_patch_w          <- tempfile(fileext = ".tif")
r_water_clusters_opp <- terra::patches(r_binair_water_opp, directions = 8, zeroAsNA = TRUE, filename = tmp_patch_w, overwrite = TRUE)

witsnuit_water_finaal_opp <- template_Vlaanderen

if (!all(is.na(suppressWarnings(terra::minmax(r_water_clusters_opp))))) {
  r_poel_kernen_opp <- r_water_clusters_opp * NA
  df_cells_opp      <- terra::as.data.frame(r_water_clusters_opp, cells = TRUE)
  colnames(df_cells_opp) <- c("cell", "ID")
  
  dt_kernen_opp <- as.data.table(df_cells_opp)[, .(cell = cell[1]), by = ID]
  r_poel_kernen_opp[dt_kernen_opp$cell] <- dt_kernen_opp$ID
  
  xy_kernen_opp <- terra::xyFromCell(r_poel_kernen_opp, dt_kernen_opp$cell)
  n_plassen_opp <- nrow(xy_kernen_opp)
  
  if (n_plassen_opp >= 4) {
    dist_matrix_opp    <- stats::dist(xy_kernen_opp, method = "euclidean") %>% as.matrix()
    cluster_counts_opp <- rowSums(dist_matrix_opp <= 500)
    
    df_scores_opp        <- data.table(ID = dt_kernen_opp$ID, Buren_Count = cluster_counts_opp)
    goedgekeurde_ids_opp <- df_scores_opp[Buren_Count >= 4, ID]
    
    if (length(goedgekeurde_ids_opp) > 0) {
      masker_water_bin_opp      <- r_water_clusters_opp %in% goedgekeurde_ids_opp
      witsnuit_water_finaal_opp <- terra::mask(r_witsnuit_waterbiotoop1_opp, terra::ifel(masker_water_bin_opp, 1, NA))
    }
    rm(dist_matrix_opp, df_scores_opp, goedgekeurde_ids_opp, cluster_counts_opp)
  }
  rm(r_poel_kernen_opp, df_cells_opp, dt_kernen_opp, xy_kernen_opp)
}

unlink(tmp_patch_w)
rm(r_water_clusters_opp, r_binair_water_opp, r_witsnuit_waterbiotoop1_opp)
gc()

# 8. STAP 3: FOERAGEERHABITAT FILTEREN EN STRUWEEL-AANDEEL (>= 50 HA & >= 70% STRUWEEL)
message("-> Stap 3: Foerageergebieden filteren op grootte (>= 50 ha) en struweel-aandeel (>= 70%)...")

witsnuit_foerageer_finaal_opp <- template_Vlaanderen

if (!all(is.na(suppressWarnings(terra::minmax(struweel1_opp))))) {
  
  r_binair_struweel <- terra::ifel(!is.na(struweel1_opp) & struweel1_opp > 0, 1, NA)
  
  tmp_buf_foer       <- tempfile(fileext = ".tif")
  foerageer_buffered <- terra::buffer(r_binair_struweel, width = 25, filename = tmp_buf_foer, overwrite = TRUE)
  
  tmp_patch_foer     <- tempfile(fileext = ".tif")
  cl_foerageer       <- terra::patches(foerageer_buffered, directions = 8, zeroAsNA = TRUE, filename = tmp_patch_foer, overwrite = TRUE) %>% 
    terra::mask(r_binair_struweel)
  
  unlink(c(tmp_buf_foer, tmp_patch_foer))
  
  df_cl_foer <- terra::as.data.frame(cl_foerageer, cells = TRUE)
  
  if (nrow(df_cl_foer) > 0) {
    colnames(df_cl_foer) <- c("cell", "ID")
    dt_foer_cl          <- as.data.table(df_cl_foer)
    
    # 1. Bereken totale foerageeroppervlakte per cluster
    dt_stats            <- dt_foer_cl[, .(Foerageer_ha = .N * 0.01), by = ID]
    
    # 2. Bereken struweeloppervlakte in 50m zone rond cluster
    w_50m <- terra::focalMat(template_Vlaanderen, 50, type = 'circle')
    w_50m[w_50m > 0] <- 1
    
    tmp_focal_foer   <- tempfile(fileext = ".tif")
    cl_foerageer_buf <- terra::focal(cl_foerageer, w = w_50m, fun = "max", na.rm = TRUE, filename = tmp_focal_foer, overwrite = TRUE)
    
    struweel_in_zone_opp <- terra::mask(struweel1_opp, cl_foerageer_buf)
    zonal_struweel_opp   <- terra::zonal(struweel_in_zone_opp, cl_foerageer_buf, fun = "sum", na.rm = TRUE)
    
    unlink(tmp_focal_foer)
    
    dt_struweel <- as.data.table(zonal_struweel_opp)
    colnames(dt_struweel) <- c("ID", "Struweel_Fractie_Sum")
    dt_struweel[, Struweel_ha := Struweel_Fractie_Sum * 0.01]
    
    dt_criteria <- merge(dt_stats, dt_struweel, by = "ID", all.x = TRUE)
    dt_criteria[is.na(Struweel_ha), Struweel_ha := 0]
    dt_criteria[, Percentage_Struweel := (Struweel_ha / Foerageer_ha) * 100]
    
    goedgekeurde_foer_ids <- dt_criteria[Foerageer_ha >= 50 & Percentage_Struweel >= 70, ID]
    
    if (length(goedgekeurde_foer_ids) > 0) {
      masker_foer_bin               <- cl_foerageer %in% goedgekeurde_foer_ids
      witsnuit_foerageer_finaal_opp <- terra::mask(struweel1_opp, terra::ifel(masker_foer_bin, 1, NA))
    }
    rm(dt_foer_cl, dt_stats, dt_struweel, dt_criteria, goedgekeurde_foer_ids)
  }
  rm(cl_foerageer, r_binair_struweel)
}

rm(foerageer0_notbwk_opp, struweel1_opp)
gc()

# 9. STAP 4: RUIMTELIJK RELATIE WATER - FOERAGEERGEBIED (800M) -----------------
message("-> Stap 4: Ruimtelijke koppelingsafstand tussen water en foerageergebied (800m)...")

witsnuit_water_relatie_opp <- template_Vlaanderen
witsnuit_foer_relatie_opp  <- template_Vlaanderen

if (!all(is.na(suppressWarnings(terra::minmax(witsnuit_water_finaal_opp))))) && 
  !all(is.na(suppressWarnings(terra::minmax(witsnuit_foerageer_finaal_opp))))) {
    
    r_binair_water_finaal_opp <- terra::ifel(!is.na(witsnuit_water_finaal_opp) & witsnuit_water_finaal_opp > 0, 1, NA)
    r_binair_foer_finaal_opp  <- terra::ifel(!is.na(witsnuit_foerageer_finaal_opp) & witsnuit_foerageer_finaal_opp > 0, 1, NA)
    
    tmp_dist_w800 <- tempfile(fileext = ".tif")
    tmp_dist_f800 <- tempfile(fileext = ".tif")
    
    water_buffer_800m_opp <- terra::distance(r_binair_water_finaal_opp, filename = tmp_dist_w800, overwrite = TRUE) <= 800
    foer_buffer_800m_opp  <- terra::distance(r_binair_foer_finaal_opp, filename = tmp_dist_f800, overwrite = TRUE) <= 800
    
    witsnuit_water_relatie_opp <- terra::mask(witsnuit_water_finaal_opp, foer_buffer_800m_opp, maskvalues = FALSE)
    witsnuit_foer_relatie_opp  <- terra::mask(witsnuit_foerageer_finaal_opp, water_buffer_800m_opp, maskvalues = FALSE)
    
    unlink(c(tmp_dist_w800, tmp_dist_f800))
    rm(r_binair_water_finaal_opp, r_binair_foer_finaal_opp, water_buffer_800m_opp, foer_buffer_800m_opp)
  }

rm(witsnuit_water_finaal_opp, witsnuit_foerageer_finaal_opp)
gc()

# 10. STAP 5: FINAAL LEEFGEBIED COMPLEX (35 HA) & WATER-FIRST AFSNIJDING -------
message("-> Stap 5: Finaal leefgebied complex filteren op 35 ha met 800m water-first afsnijding...")

w_clean <- terra::ifel(is.na(witsnuit_water_relatie_opp), 0, witsnuit_water_relatie_opp)
f_clean <- terra::ifel(is.na(witsnuit_foer_relatie_opp), 0, witsnuit_foer_relatie_opp)

som_fracties                 <- w_clean + f_clean
som_gemedieerd               <- terra::clamp(som_fracties, upper = 1.0)
witsnuit_leefgebied1_opp_raw <- terra::ifel(som_gemedieerd > 0, som_gemedieerd, NA)

r_binair_leef1_opp <- terra::ifel(!is.na(witsnuit_leefgebied1_opp_raw) & witsnuit_leefgebied1_opp_raw > 0, 1, NA)

if (!all(is.na(suppressWarnings(terra::minmax(r_binair_leef1_opp))))) {
  
  tmp_buf_leef             <- tempfile(fileext = ".tif")
  leefgebied_buffered_opp  <- terra::buffer(r_binair_leef1_opp, width = 25, filename = tmp_buf_leef, overwrite = TRUE)
  
  tmp_patch_leef           <- tempfile(fileext = ".tif")
  r_leefgebied_clusters_opp<- terra::patches(leefgebied_buffered_opp, directions = 8, filename = tmp_patch_leef, overwrite = TRUE)
  r_scherpe_clusters_opp   <- terra::mask(r_leefgebied_clusters_opp, r_binair_leef1_opp)
  
  unlink(c(tmp_buf_leef, tmp_patch_leef))
  
  cluster_stats_opp <- terra::zonal(witsnuit_leefgebied1_opp_raw, r_scherpe_clusters_opp, fun = "sum", na.rm = TRUE)
  colnames(cluster_stats_opp) <- c("value", "ha_exact")
  
  dt_cl_stats          <- as.data.table(cluster_stats_opp)
  goedgekeurde_ids_opp <- dt_cl_stats[ha_exact >= 35, value]
  
  if (length(goedgekeurde_ids_opp) > 0) {
    masker_opp_leef     <- r_scherpe_clusters_opp %in% goedgekeurde_ids_opp
    leef_35ha_mask      <- terra::ifel(masker_opp_leef, 1, NA)
    
    water_finaal_opp    <- terra::mask(witsnuit_water_relatie_opp, leef_35ha_mask)
    land_finaal_opp_ruw <- terra::mask(witsnuit_foer_relatie_opp, leef_35ha_mask)
  } else {
    water_finaal_opp    <- template_Vlaanderen
    land_finaal_opp_ruw <- template_Vlaanderen
  }
} else {
  water_finaal_opp    <- template_Vlaanderen
  land_finaal_opp_ruw <- template_Vlaanderen
}

if (!all(is.na(suppressWarnings(terra::minmax(water_finaal_opp))))) {
  r_bin_water_opp <- terra::ifel(!is.na(water_finaal_opp) & water_finaal_opp > 0, 1, NA)
  poly_water_opp  <- terra::as.polygons(r_bin_water_opp, aggregate = TRUE)
  poly_buffer_opp <- terra::buffer(poly_water_opp, width = straal_water_m)
  land_finaal_opp <- terra::mask(land_finaal_opp_ruw, poly_buffer_opp)
} else {
  land_finaal_opp <- land_finaal_opp_ruw
}

w_f_clean <- terra::ifel(is.na(water_finaal_opp), 0, water_finaal_opp)
l_f_clean <- terra::ifel(is.na(land_finaal_opp), 0, land_finaal_opp)

final_opp_som <- w_f_clean + l_f_clean
final_opp_raw <- terra::clamp(final_opp_som, upper = 1.0)
final_opp     <- terra::ifel(final_opp_raw > 0, final_opp_raw, NA)

if (!all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  tmp_patch_cl <- tempfile(fileext = ".tif")
  cl_opp       <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE, filename = tmp_patch_cl, overwrite = TRUE)
} else {
  cl_opp <- template_Vlaanderen
}

rm(witsnuit_water_relatie_opp, witsnuit_foer_relatie_opp, w_clean, f_clean, som_fracties, som_gemedieerd, 
   witsnuit_leefgebied1_opp_raw, r_binair_leef1_opp, w_f_clean, l_f_clean, final_opp_som, final_opp_raw)
gc()

# 11. EXPORT NAAR DISK (ENKEL 00_ID_Rasters EN 01_Werkelijke_Oppervlaktes) -----
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

message("🏁 LEEFGEBIEDENGRID SIERLIJKE WITSNUITLIBEL SUCCESVOL GEËXPORTEERD VOOR HEEL VLAANDEREN!")
