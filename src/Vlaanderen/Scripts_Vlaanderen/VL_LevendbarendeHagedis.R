# ==============================================================================
# LEEFGEBIEDENGRID VOOR LEVENDBARENDE HAGEDIS - HEEL VLAANDEREN (WERKELIJKE OPPERVLAKTE)
# Ecologische analyse vlek- en netwerkdrempels op werkelijke fracties (INBO)
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

# 3. SOORT INFORMATIE OPHALEN --------------------------------------------------
soort <- "levendbarendehagedis"

oppervlakte_ha <- 5    # Minimale vlek-oppervlakte (MinOpp)
min_totale_ha  <- 10   # Minimale netwerk-oppervlakte (Min totale Opp)
afstand_m      <- 0    # Afstand tussen biotopen
buffer_m       <- 100  # Dispersiecapaciteit

drempel_vlek_m2    <- oppervlakte_ha * 10000  # 5 ha = 50.000 m²
drempel_netwerk_m2 <- min_totale_ha * 10000   # 10 ha = 100.000 m²

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

if(nrow(resultaten_gegroepeerd) > 0) {
  for(i in 1:nrow(resultaten_gegroepeerd)) {
    h_type <- resultaten_gegroepeerd$Type[i]
    h_data <- resultaten_gegroepeerd$Data[[i]]
    
    exact_codes <- tolower(trimws(h_data$Code[h_data$Match == "exact"]))
    bevat_codes <- tolower(trimws(h_data$Code[h_data$Match == "bevat"]))
    
    if(length(bevat_codes) > 0) {
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
    
    if(nrow(tabel_sub) > 0) {
      tabel_sub[, CODE := tolower(trimws(CODE))]
      tabel_unique  <- unique(tabel_sub, by = c("cel_id", "CODE"))
      tabel_cel_som <- tabel_unique[, .(Oppervlakte = pmin(sum(BWK_FRAC, na.rm = TRUE), 1.0)), by = .(cel_id)]
      
      r_opp_type <- template_Vlaanderen
      r_opp_type[tabel_cel_som$cel_id] <- tabel_cel_som$Oppervlakte
      lijst_oppervlaktes[[h_type]]     <- r_opp_type
    }
  }
}

lijst_oppervlaktes <- purrr::compact(lijst_oppervlaktes)

if (length(lijst_oppervlaktes) == 0) {
  message("⚠️ Geen enkele doel-BWK code voor '", soort, "' aangetroffen in de database. Leeg raster wordt geëxporteerd.")
  opp_bron <- template_Vlaanderen
} else {
  opp_bron <- lijst_oppervlaktes[[1]]
}

rm(ds_vlaanderen, resultaten_gegroepeerd, df_nieuw, lijst_oppervlaktes)
gc()

# 6. STAP 1: VOORFILTER - LOSSE VLEKKEN < 5 HA FILTEREN ------------------------
message("-> Stap 1: Voorfilter losse biotopen < 5 ha op werkelijke oppervlaktes...")

if (!all(is.na(suppressWarnings(terra::minmax(opp_bron))))) {
  r_binair_bron <- terra::ifel(!is.na(opp_bron) & opp_bron > 0, 1, NA)
  
  tmp_patch_vlek <- tempfile(fileext = ".tif")
  vlek_patches   <- terra::patches(r_binair_bron, directions = 8, zeroAsNA = TRUE, 
                                   filename = tmp_patch_vlek, overwrite = TRUE)
  
  df_vlek <- terra::as.data.frame(vlek_patches, cells = TRUE)
  if (nrow(df_vlek) > 0) {
    colnames(df_vlek) <- c("cell", "ID")
    opp_val_vlek     <- terra::extract(opp_bron, df_vlek$cell)[[1]]
    dt_vlek          <- data.table(ID = df_vlek$ID, Waarde = opp_val_vlek)
    
    stats_vlek            <- dt_vlek[!is.na(ID) & !is.na(Waarde), .(Area_m2 = sum(Waarde, na.rm = TRUE) * 100), by = ID]
    gevalideerde_vlek_ids <- stats_vlek[Area_m2 >= drempel_vlek_m2, ID]
    
    if (length(gevalideerde_vlek_ids) > 0) {
      opp_gefilterd <- terra::mask(opp_bron, vlek_patches %in% gevalideerde_vlek_ids)
    } else {
      opp_gefilterd <- template_Vlaanderen
    }
    rm(df_vlek, dt_vlek, stats_vlek)
  } else {
    opp_gefilterd <- template_Vlaanderen
  }
  unlink(tmp_patch_vlek)
  rm(r_binair_bron)
} else {
  opp_gefilterd <- template_Vlaanderen
}

gc()

# 7. STAP 2: NETWERKCLUSTERING - VERBIND VLEKKEN BINNEN 100M EN CHECK >= 10 HA -
message("-> Stap 2: Netwerkclustering uitvoeren op 100m (min. 10 ha netwerk)...")

r_binair_gefilterd <- terra::ifel(!is.na(opp_gefilterd) & opp_gefilterd > 0, 1, NA)

res_opp <- cluster_filter_compleet(
  masker     = r_binair_gefilterd, 
  opp_laag   = opp_gefilterd, 
  drempel_m2 = drempel_netwerk_m2, 
  dist_m     = buffer_m, 
  werkelijk  = TRUE
)

final_opp <- res_opp$raster
cl_opp    <- res_opp$clusters

if (is.null(cl_opp) || all(is.na(terra::values(cl_opp, mat = FALSE)))) {
  if (!is.null(final_opp) && !all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
    tmp_patch <- tempfile(fileext = ".tif")
    cl_opp    <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE, filename = tmp_patch, overwrite = TRUE)
  } else {
    cl_opp <- template_Vlaanderen
  }
}

rm(opp_bron, opp_gefilterd, r_binair_gefilterd, res_opp)
gc()

# 8. EXPORT NAAR DISK (ENKEL 00_ID_Rasters EN 01_Werkelijke_Oppervlaktes) -----
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

message("🏁 LEEFGEBIEDENGRID LEVENDBARENDE HAGEDIS SUCCESVOL GEËXPORTEERD VOOR HEEL VLAANDEREN!")
