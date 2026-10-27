# ==============================================================================
# LEEFGEBIEDENGRID VOOR PAAPJE - HEEL VLAANDEREN (WERKELIJKE OPPERVLAKTE)
# Ecologische analyse open graslanden, bosrandverstoring & ruimtelijke clustering (INBO)
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
df <- read_excel(here("data/input/Excel_files/Soorten_bwk_afstanden.xlsx"))
soort <- "paapje"

resultaat <- df %>%
  filter(tolower(trimws(Soort)) == soort) %>%
  select(Type, MinOpp_ha, AfstandBiotopen_m, Dispersiecap_m)

oppervlakte_ha <- resultaat$MinOpp_ha[1]
afstand_m      <- resultaat$AfstandBiotopen_m[1]

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

bwk_opp    <- if("bwk" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["bwk"]] else template_Vlaanderen
bwkbos_opp <- if("bwkbos" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["bwkbos"]] else template_Vlaanderen

rm(ds_vlaanderen, resultaten_gegroepeerd, df_nieuw, lijst_oppervlaktes)
gc()

# 6. STAP 1: BIOTOOP EN BOSLAGEN COMBINEREN EN FILTEREN (GEFILTERD OP >= 20 ARE)
message("-> Stap 1: Boslagen verwerken en filteren op >= 20 are...")

groenkaart_raw <- terra::rast(here("data/input/ASCI Files/Groenkaart_2021.tif"))
if(!terra::compareGeom(groenkaart_raw, template_Vlaanderen, stopOnError = FALSE)) {
  groenkaart_aligned <- terra::resample(groenkaart_raw, template_Vlaanderen, method = "near")
} else {
  groenkaart_aligned <- groenkaart_raw
}

groenkaart_bos <- terra::ifel(groenkaart_aligned == 1, 1.0, NA)

bwk_binair          <- terra::ifel(!is.na(bwk_opp) & bwk_opp > 0, 1, NA)
bwkbos_gecorrigeerd <- terra::mask(bwkbos_opp, bwk_binair, inverse = TRUE)

b1_clean <- terra::ifel(is.na(bwkbos_gecorrigeerd), 0, bwkbos_gecorrigeerd)
b2_clean <- terra::ifel(is.na(groenkaart_bos), 0, groenkaart_bos)

paapje_bos_ruw_opp <- terra::clamp(b1_clean + b2_clean, upper = 1.0)
paapje_bos_ruw_bin <- terra::ifel(paapje_bos_ruw_opp > 0, 1, NA)

# OPPERVLAKTEFILTER VOOR BOS: 20 ARE = 20 CELLEN VAN 10X10M
drempel_bos_cellen <- 20 

tmp_patch_bos <- tempfile(fileext = ".tif")
bos_patches   <- terra::patches(paapje_bos_ruw_bin, directions = 8, zeroAsNA = TRUE, filename = tmp_patch_bos, overwrite = TRUE)

df_bos_cl <- terra::as.data.frame(bos_patches, cells = TRUE)
if (nrow(df_bos_cl) > 0) {
  colnames(df_bos_cl) <- c("cell", "ID")
  stats_bos_count     <- data.table(ID = df_bos_cl$ID)[, .(count = .N), by = ID]
  grote_bos_ids       <- stats_bos_count[count >= drempel_bos_cellen, ID]
  paapje_bos_bin      <- terra::ifel(bos_patches %in% grote_bos_ids, 1, NA)
} else {
  paapje_bos_bin      <- template_Vlaanderen
}

unlink(tmp_patch_bos)
rm(groenkaart_raw, groenkaart_aligned, groenkaart_bos, bwkbos_gecorrigeerd, paapje_bos_ruw_opp, paapje_bos_ruw_bin, 
   bos_patches, df_bos_cl, b1_clean, b2_clean, bwk_binair, bwkbos_opp)
gc()

# 7. STAP 2: 100 METER BOSRAND-BUFFER TOEPASSEN ---------------------------------
message("-> Stap 2: 100 meter bosrand-buffer berekenen en toepassen op biotoop...")

if (!all(is.na(suppressWarnings(terra::minmax(paapje_bos_bin))))) {
  tmp_dist_bos            <- tempfile(fileext = ".tif")
  afstand_tot_bos         <- terra::distance(paapje_bos_bin, filename = tmp_dist_bos, overwrite = TRUE)
  geschikt_buiten_bosrand <- terra::ifel(afstand_tot_bos >= 100, 1, NA)
  
  paapje_leefgebied1_opp  <- terra::mask(bwk_opp, geschikt_buiten_bosrand)
  unlink(tmp_dist_bos)
} else {
  paapje_leefgebied1_opp  <- bwk_opp
}

rm(paapje_bos_bin, bwk_opp)
gc()

# 8. STAP 3: CLUSTERING EN MINIMALE OPPERVLAKTE (MIN. 5 HA) ---------------------
message("-> Stap 3: Clustering en minimale oppervlaktefilter (min. ", oppervlakte_ha, " ha)...")

drempel_m2         <- oppervlakte_ha * 10000 
zoekafstand_m      <- 50 
r_binair_leef1_opp <- terra::ifel(!is.na(paapje_leefgebied1_opp) & paapje_leefgebied1_opp > 0, 1, NA)

cluster_output_opp <- cluster_filter_compleet(
  masker     = r_binair_leef1_opp, 
  opp_laag   = paapje_leefgebied1_opp, 
  drempel_m2 = drempel_m2, 
  dist_m     = zoekafstand_m, 
  werkelijk  = TRUE
)

final_opp <- cluster_output_opp$raster
cl_opp    <- cluster_output_opp$clusters

if (is.null(cl_opp) || all(is.na(terra::values(cl_opp, mat = FALSE)))) {
  if (!is.null(final_opp) && !all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
    tmp_patch <- tempfile(fileext = ".tif")
    cl_opp    <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE, filename = tmp_patch, overwrite = TRUE)
  } else {
    cl_opp <- template_Vlaanderen
  }
}

rm(r_binair_leef1_opp, paapje_leefgebied1_opp, cluster_output_opp)
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

message("🏁 LEEFGEBIEDENGRID PAAPJE SUCCESVOL GEËXPORTEERD VOOR HEEL VLAANDEREN!")
