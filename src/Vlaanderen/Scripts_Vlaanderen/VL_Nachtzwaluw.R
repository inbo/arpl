# ==============================================================================
# LEEFGEBIEDENGRID VOOR NACHTZWALUW - HEEL VLAANDEREN (WERKELIJKE OPPERVLAKTE)
# Ecologische analyse bosranden, open heide/kapvlakten & ruimtelijke clustering (INBO)
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

bepaal_bosrand <- function(r_voortplanting, r_bos_cluster) {
  v_bin <- terra::ifel(!is.na(r_voortplanting) & r_voortplanting > 0, 1, 0)
  b_bin <- terra::ifel(!is.na(r_bos_cluster) & r_bos_cluster > 0, 1, 0)
  
  tmp_focal  <- tempfile(fileext = ".tif")
  bos_buffer <- terra::focal(b_bin, w = 3, fun = "max", na.rm = TRUE, filename = tmp_focal, overwrite = TRUE)
  bosrand    <- terra::ifel(v_bin == 1 & bos_buffer == 1, 1, NA)
  
  unlink(tmp_focal)
  return(bosrand)
}

# 3. SOORT INFORMATIE OPHALEN --------------------------------------------------
soort <- "nachtzwaluw"

df <- read_excel(here("data/input/Excel_files/Soorten_bwk_afstanden.xlsx"))
resultaat <- df %>%
  filter(tolower(trimws(Soort)) == soort) %>%
  select(Type, MinOpp_ha, AfstandBiotopen_m, Dispersiecap_m)

oppervlakte_ha <- resultaat$MinOpp_ha[1]
afstand_m      <- resultaat$AfstandBiotopen_m[1]
buffer_m       <- resultaat$Dispersiecap_m[1]

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

voortplanting_opp <- if("voortplanting" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["voortplanting"]] else template_Vlaanderen
bos_opp          <- if("bos" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["bos"]] else template_Vlaanderen

rm(ds_vlaanderen, resultaten_gegroepeerd, df_nieuw, lijst_oppervlaktes)
gc()

# 6. STAP 1: BOSCLUSTERING & 5 HA FILTER ----------------------------------------
message("-> Stap 1: Bosclustering en 5 ha filter toepassen op werkelijke oppervlaktes...")

drempel_bos_m2      <- 5 * 10000 
cluster_afstand_bos <- 50    

r_binair_bos_opp <- terra::ifel(!is.na(bos_opp) & bos_opp > 0, 1, NA)

bos_opp_gefilterd_lijst <- cluster_filter_compleet(
  masker     = r_binair_bos_opp, 
  opp_laag   = bos_opp, 
  drempel_m2 = drempel_bos_m2, 
  dist_m     = cluster_afstand_bos, 
  werkelijk  = TRUE
)
bos_cluster_5ha_opp <- bos_opp_gefilterd_lijst$raster

rm(r_binair_bos_opp, bos_opp_gefilterd_lijst, bos_opp)
gc()

# 7. STAP 2: BOSRAND DETECTIE & FILTER (MIN. 1500M / 150 PIXELS) ---------------
message("-> Stap 2: Bosranden detecteren en filteren op minstens 1500m...")

min_rand_pixels <- 150 

rand_opp <- bepaal_bosrand(voortplanting_opp, bos_cluster_5ha_opp)

r_binair_vp_opp <- terra::ifel(!is.na(voortplanting_opp) & voortplanting_opp > 0, 1, NA)

tmp_buf_vp <- tempfile(fileext = ".tif")
r_buf_vp   <- terra::buffer(r_binair_vp_opp, width = 25, filename = tmp_buf_vp, overwrite = TRUE)

tmp_patch_vp         <- tempfile(fileext = ".tif")
voortplanting_cl_opp <- terra::patches(r_buf_vp, directions = 8, zeroAsNA = TRUE, filename = tmp_patch_vp, overwrite = TRUE) %>% 
  terra::mask(r_binair_vp_opp)

unlink(c(tmp_buf_vp, tmp_patch_vp))

rand_stats_opp           <- terra::zonal(rand_opp, voortplanting_cl_opp, fun = "sum", na.rm = TRUE)
colnames(rand_stats_opp) <- c("ID", "RandPixels")

voldoende_rand_ids_opp <- rand_stats_opp$ID[rand_stats_opp$RandPixels >= min_rand_pixels]

if (length(voldoende_rand_ids_opp) > 0) {
  leefgebied1_opp <- terra::mask(voortplanting_opp, voortplanting_cl_opp %in% voldoende_rand_ids_opp)
} else {
  leefgebied1_opp <- template_Vlaanderen
}

rm(rand_opp, voortplanting_cl_opp, rand_stats_opp, r_binair_vp_opp, bos_cluster_5ha_opp, voortplanting_opp)
gc()

# 8. STAP 3: FINALE CLUSTERING (25 HA) -----------------------------------------
message("-> Stap 3: Finale clustering op 25 ha (50m netwerk)...")

drempel_leefgebied_m2      <- 25 * 10000 
cluster_afstand_leefgebied <- 50       

r_binair_leef1_opp <- terra::ifel(!is.na(leefgebied1_opp) & leefgebied1_opp > 0, 1, NA)

leefgebied_opp_lijst <- cluster_filter_compleet(
  masker     = r_binair_leef1_opp, 
  opp_laag   = leefgebied1_opp, 
  drempel_m2 = drempel_leefgebied_m2, 
  dist_m     = cluster_afstand_leefgebied, 
  werkelijk  = TRUE
)

final_opp <- leefgebied_opp_lijst$raster
cl_opp    <- leefgebied_opp_lijst$clusters

if (is.null(cl_opp) || all(is.na(terra::values(cl_opp, mat = FALSE)))) {
  if (!is.null(final_opp) && !all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
    tmp_patch <- tempfile(fileext = ".tif")
    cl_opp    <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE, filename = tmp_patch, overwrite = TRUE)
  } else {
    cl_opp <- template_Vlaanderen
  }
}

rm(leefgebied_opp_lijst, r_binair_leef1_opp, leefgebied1_opp)
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

message("🏁 LEEFGEBIEDENGRID NACHTZWALUW SUCCESVOL GEËXPORTEERD VOOR HEEL VLAANDEREN!")
