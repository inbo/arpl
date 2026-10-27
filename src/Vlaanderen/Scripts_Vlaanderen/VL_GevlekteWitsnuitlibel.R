# ==============================================================================
# LEEFGEBIEDENGRID VOOR GEVLEKTE WITSNUITLIBEL - HEEL VLAANDEREN (WERKELIJKE OPPERVLAKTE)
# Ecologische analyse water-first landkoppeling & struweelaandeel (INBO)
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

# 3. SOORT PARAMETERS OPHALEN --------------------------------------------------
soort                <- "gevlektewitsnuitlibel"
buffer_m             <- 20000  # Dispersiecapaciteit voor smear/buffer
straal_water_m       <- 200    # Afsnij-afstand van land rond het water (200m)

drempel_water_min_m2 <- 200    # Min. 200 m² (0.02 ha)
drempel_water_max_m2 <- 10000  # Max. 10.000 m² (1.0 ha)
dist_water_m          <- 10     # 10m overbrugging voor water

drempel_land_m2      <- 50000   # Min. 5 ha land
dist_land_m          <- 10     # 10m overbrugging voor aanpalend land

drempel_ruim_m2      <- 300000  # Min. 30 ha ruim land
dist_ruim_m          <- 50     # 50m overbrugging voor ruim land

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
land_nabij_opp        <- if("land_nabij" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["land_nabij"]] else template_Vlaanderen
land_ruim_opp         <- if("land_ruim" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["land_ruim"]] else template_Vlaanderen
struweel1_opp         <- if("struweel1" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["struweel1"]] else template_Vlaanderen

rm(ds_vlaanderen, resultaten_gegroepeerd, df_nieuw, lijst_oppervlaktes)
gc()

# 6. WATERVLAKKEN SHAPEFILE INLEZEN EN BUNDELEN MET BWK-WATER ------------------
message("-> Watervlakken shapefile inlezen voor heel Vlaanderen...")

watervlakken_v    <- vect(here("data/input/ASCI Files/watervlakken2024.shp"))
watervlakken_proj <- project(watervlakken_v, crs(template_Vlaanderen))
watervlakken_clean <- makeValid(watervlakken_proj)

# Fuzzy clustering op kleine wateren (<= 200m²)
watervlakken_buffered       <- buffer(watervlakken_clean, width = 5)
watervlakken_fuzzy_clusters <- aggregate(watervlakken_buffered, by = NULL)
watervlakken_fuzzy_clusters$fuzzy_area_m2 <- terra::expanse(watervlakken_fuzzy_clusters)

kleine_clusters_v <- watervlakken_fuzzy_clusters[watervlakken_fuzzy_clusters$fuzzy_area_m2 <= 200, ]

if (!is.null(kleine_clusters_v) && nrow(kleine_clusters_v) > 0) {
  kleine_watervlakken_finaal <- crop(watervlakken_clean, kleine_clusters_v)
  if (!is.null(kleine_watervlakken_finaal) && nrow(kleine_watervlakken_finaal) > 0) {
    r_witsnuit_kleinwater <- terra::rasterize(kleine_watervlakken_finaal, template_Vlaanderen, field = 1, background = NA)
    r_witsnuit_waterbiotoop1_opp <- terra::cover(voortplanting_bwk_opp, r_witsnuit_kleinwater)
  } else {
    r_witsnuit_waterbiotoop1_opp <- voortplanting_bwk_opp
  }
} else {
  r_witsnuit_waterbiotoop1_opp <- voortplanting_bwk_opp
}

suppressWarnings(
  rm(watervlakken_v, watervlakken_proj, watervlakken_clean, watervlakken_buffered, watervlakken_fuzzy_clusters, kleine_clusters_v, kleine_watervlakken_finaal, r_witsnuit_kleinwater, voortplanting_bwk_opp)
)
gc()

# 7. WATER FILTEREN: 0.02 HA <= OPP <= 1 HA + DICHTHEIDSEIS (3 BUREN) ----------
message("-> Water filteren: 0.02 ha <= Opp <= 1 ha (10m netwerk) + dichtheidseis...")

if (!is.null(r_witsnuit_waterbiotoop1_opp) && !all(is.na(suppressWarnings(terra::minmax(r_witsnuit_waterbiotoop1_opp))))) {
  r_binair_water_opp <- terra::ifel(!is.na(r_witsnuit_waterbiotoop1_opp) & r_witsnuit_waterbiotoop1_opp > 0, 1, NA)
  
  res_water_opp <- cluster_filter_compleet(
    masker     = r_binair_water_opp, 
    opp_laag   = r_witsnuit_waterbiotoop1_opp, 
    drempel_m2 = drempel_water_min_m2, 
    dist_m     = dist_water_m, 
    werkelijk  = TRUE
  )
  r_water_cand_opp  <- res_water_opp$raster
  cl_water_cand_opp <- res_water_opp$clusters
  
  if (!all(is.na(suppressWarnings(terra::minmax(cl_water_cand_opp))))) {
    stats_water_opp <- terra::zonal(r_water_cand_opp, cl_water_cand_opp, fun = "sum", na.rm = TRUE)
    colnames(stats_water_opp) <- c("ID", "Waarde")
    stats_water_opp$Area_m2 <- stats_water_opp$Waarde * 100
    
    valide_opp_ids_opp <- stats_water_opp$ID[stats_water_opp$Area_m2 <= drempel_water_max_m2]
    
    if (length(valide_opp_ids_opp) > 0) {
      r_water_cand_opp <- terra::mask(r_water_cand_opp, terra::ifel(cl_water_cand_opp %in% valide_opp_ids_opp, 1, NA))
    } else { r_water_cand_opp <- template_Vlaanderen }
  } else { r_water_cand_opp <- template_Vlaanderen }
  
  if (!all(is.na(suppressWarnings(terra::minmax(r_water_cand_opp))))) {
    tmp_patch_w      <- tempfile(fileext = ".tif")
    r_water_cl_opp   <- terra::patches(r_water_cand_opp, directions = 8, zeroAsNA = TRUE, filename = tmp_patch_w, overwrite = TRUE)
    
    if (!all(is.na(suppressWarnings(terra::minmax(r_water_cl_opp))))) {
      watervlakken_poly <- terra::as.polygons(r_water_cl_opp, dissolve = TRUE)
      
      if (!is.null(watervlakken_poly) && nrow(watervlakken_poly) > 0) {
        centroids        <- terra::centroids(watervlakken_poly)
        centroids_buffer <- terra::buffer(centroids, width = 200)
        intersect_matrix <- terra::relate(centroids_buffer, centroids, "intersects")
        cluster_counts   <- rowSums(intersect_matrix)
        
        watervlakken_poly$unieke_buren_count <- cluster_counts
        goedgekeurde_clusters_poly <- watervlakken_poly[watervlakken_poly$unieke_buren_count >= 4, ]
        
        if (nrow(goedgekeurde_clusters_poly) > 0) {
          gevlektewitsnuitlibel_voortplanting_opp <- terra::mask(r_water_cand_opp, terra::rasterize(goedgekeurde_clusters_poly, template_Vlaanderen, field = 1))
        } else {
          gevlektewitsnuitlibel_voortplanting_opp <- template_Vlaanderen
        }
      } else {
        gevlektewitsnuitlibel_voortplanting_opp <- template_Vlaanderen
      }
      unlink(tmp_patch_w)
    } else {
      gevlektewitsnuitlibel_voortplanting_opp <- template_Vlaanderen
    }
  } else {
    gevlektewitsnuitlibel_voortplanting_opp <- template_Vlaanderen
  }
} else {
  gevlektewitsnuitlibel_voortplanting_opp <- template_Vlaanderen
}

rm(r_witsnuit_waterbiotoop1_opp)
gc()

# 8. WATER-FIRST LANDBIOTOOP FILTERING -----------------------------------------
message("-> Landbiotoop filteren (10m & 50m netwerken, 30% bos) & afsnijden op 200m van water...")

# --- STAP 8.1: RUIM LANDBIOTOOP FILTEREN (30 HA) ---
if (!all(is.na(suppressWarnings(terra::minmax(land_ruim_opp))))) {
  r_bin_ruim_opp <- terra::ifel(!is.na(land_ruim_opp) & land_ruim_opp > 0, 1, NA)
  res_ruim_opp   <- cluster_filter_compleet(
    masker     = r_bin_ruim_opp,
    opp_laag   = land_ruim_opp,
    drempel_m2 = drempel_ruim_m2, 
    dist_m     = dist_ruim_m,      
    werkelijk  = TRUE
  )
  land_ruim_ok_opp <- res_ruim_opp$raster
} else { land_ruim_ok_opp <- template_Vlaanderen }

# --- STAP 8.2: LAND NABIJDHEIDFILTER & 30% BOSSTRUCTUUR (5 HA) ---
if (!all(is.na(suppressWarnings(terra::minmax(gevlektewitsnuitlibel_voortplanting_opp))))) {
  
  r_water_bin_opp <- terra::ifel(!is.na(gevlektewitsnuitlibel_voortplanting_opp) & gevlektewitsnuitlibel_voortplanting_opp > 0, 1, NA)
  poly_water_opp  <- terra::as.polygons(r_water_bin_opp, aggregate = TRUE)
  poly_buffer_opp <- terra::buffer(poly_water_opp, width = straal_water_m)
  
  if (!all(is.na(suppressWarnings(terra::minmax(land_nabij_opp))))) {
    
    land_afgesneden_opp   <- terra::mask(land_nabij_opp, poly_buffer_opp)
    r_bin_land_afgesneden <- terra::ifel(!is.na(land_afgesneden_opp) & land_afgesneden_opp > 0, 1, NA)
    
    res_land_5ha_opp <- cluster_filter_compleet(
      masker     = r_bin_land_afgesneden, 
      opp_laag   = land_afgesneden_opp, 
      drempel_m2 = drempel_land_m2, 
      dist_m     = dist_land_m,      
      werkelijk  = TRUE
    )
    land_5ha_opp  <- res_land_5ha_opp$raster
    cl_id_5ha_opp <- res_land_5ha_opp$clusters
    
    if (!all(is.na(suppressWarnings(terra::minmax(land_5ha_opp))))) {
      r_struweel_opp_clean <- fix_bos_raster(struweel1_opp, template_Vlaanderen)
      stats_totaal_opp     <- terra::zonal(land_5ha_opp, cl_id_5ha_opp, fun = "sum", na.rm = TRUE)
      colnames(stats_totaal_opp) <- c("Cluster_ID", "Totaal_m2")
      struweel_in_cl_opp   <- terra::mask(r_struweel_opp_clean, cl_id_5ha_opp)
      stats_struweel_opp   <- terra::zonal(struweel_in_cl_opp, cl_id_5ha_opp, fun = "sum", na.rm = TRUE)
      colnames(stats_struweel_opp) <- c("Cluster_ID", "Struweel_m2")
      
      df_check_opp <- merge(stats_totaal_opp, stats_struweel_opp, by = "Cluster_ID", all.x = TRUE)
      df_check_opp$Struweel_m2[is.na(df_check_opp$Struweel_m2)] <- 0
      df_check_opp$Pct_Bos <- (df_check_opp$Struweel_m2 / df_check_opp$Totaal_m2) * 100
      valide_ids_opp <- df_check_opp$Cluster_ID[df_check_opp$Pct_Bos >= 30] 
      
      land_ok_opp <- if(length(valide_ids_opp) > 0) terra::mask(land_5ha_opp, terra::ifel(cl_id_5ha_opp %in% valide_ids_opp, 1, NA)) else template_Vlaanderen
    } else { land_ok_opp <- template_Vlaanderen }
    
    if (!all(is.na(suppressWarnings(terra::minmax(land_ok_opp)))) && !all(is.na(suppressWarnings(terra::minmax(land_ruim_ok_opp))))) {
      gevlektewitsnuitlibel_land_finaal_opp <- terra::mask(land_ok_opp, terra::ifel(!is.na(land_ruim_ok_opp) & land_ruim_ok_opp > 0, 1, NA))
    } else { gevlektewitsnuitlibel_land_finaal_opp <- land_ok_opp }
    
  } else { gevlektewitsnuitlibel_land_finaal_opp <- template_Vlaanderen }
} else { gevlektewitsnuitlibel_land_finaal_opp <- template_Vlaanderen }

water_finaal_opp <- if(exists("gevlektewitsnuitlibel_voortplanting_opp")) gevlektewitsnuitlibel_voortplanting_opp else template_Vlaanderen
final_opp        <- terra::cover(water_finaal_opp, gevlektewitsnuitlibel_land_finaal_opp)

if (!all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  tmp_patch_cl <- tempfile(fileext = ".tif")
  cl_opp       <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE, filename = tmp_patch_cl, overwrite = TRUE)
} else {
  cl_opp <- template_Vlaanderen
}

rm(land_ruim_opp, land_nabij_opp, struweel1_opp, gevlektewitsnuitlibel_voortplanting_opp, gevlektewitsnuitlibel_land_finaal_opp)
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

message("🏁 LEEFGEBIEDENGRID GEVLEKTE WITSNUITLIBEL SUCCESVOL GEËXPORTEERD VOOR HEEL VLAANDEREN!")
