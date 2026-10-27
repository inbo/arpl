# ==============================================================================
# LEEFGEBIEDENGRID VOOR VLIEGEND HERT - HEEL VLAANDEREN (WERKELIJKE OPPERVLAKTE)
# Ecologische analyse loofbossen, bosleeftijd, hellingen & bosrand-dieptebuffers (INBO)
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
soort <- "vliegendhert"

df <- read_excel(here::here("data/input/Excel_files/Soorten_bwk_afstanden.xlsx"))
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

naaldbos_opp <- if("naaldbos" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["naaldbos"]] else template_Vlaanderen

rm(ds_vlaanderen, resultaten_gegroepeerd, df_nieuw, lijst_oppervlaktes)
gc()

# 6. EXTERNE GIS-LAGEN INLEZEN EN SYNCHRONISEREN ------------------------------
message("-> Externe GIS-lagen inlezen voor het vliegend hert...")

load_and_sync <- function(pad) {
  if (file.exists(pad)) {
    r_raw <- terra::rast(pad)
    if(!terra::compareGeom(r_raw, template_Vlaanderen, stopOnError = FALSE)) {
      r_sync <- terra::resample(r_raw, template_Vlaanderen, method = "near")
    } else {
      r_sync <- r_raw
    }
    return(r_sync)
  } else {
    message("⚠️ Waarschuwing: Bestand niet gevonden op pad: ", pad)
    return(template_Vlaanderen)
  }
}

r_groenkaart     <- load_and_sync(here("data/input/ASCI Files/Groenkaart_2021.tif"))
r_groen_hoog_bin <- terra::ifel(r_groenkaart == 1, 1, NA)
r_groen_open_bin <- terra::ifel(r_groenkaart %in% c(2, 3, 4), 1, NA)

r_bosleeftijd_raw <- load_and_sync(here("data/input/Raster_Vlaanderen/vlaanderen_bosleeftijd_10m.tif"))
r_oud_bos_bin     <- terra::ifel(r_bosleeftijd_raw %in% c(1, 2), 1, NA)

r_helling_raw     <- load_and_sync(here("data/input/Raster_Vlaanderen/vliegend_hert_4gr_DTMVLII_10m_comp.tif"))
r_helling_bin     <- terra::ifel(r_helling_raw == 1, 1, NA)

rm(r_groenkaart, r_bosleeftijd_raw, r_helling_raw)
gc()

# 7. STAP 1: GESCHIKT LOOFBOS AFBAKENEN (MAX. 40% NAALDHOUT) -------------------
message("-> Geschikt loofbos afbakenen met 40% naaldhout-tolerantie...")

r_naaldbos_opp_clean <- terra::ifel(is.na(naaldbos_opp), 0, naaldbos_opp)
r_naaldbos_dominant  <- r_naaldbos_opp_clean > 0.40

r_bossen1_bin <- terra::ifel(r_groen_hoog_bin == 1 & !r_naaldbos_dominant, 1, NA)

# Afstandstoetsing oud bos (<= 1000m)
if (!all(is.na(suppressWarnings(terra::minmax(r_oud_bos_bin))))) {
  tmp_dist_ob      <- tempfile(fileext = ".tif")
  dist_oud_bos_opp <- terra::distance(r_oud_bos_bin, filename = tmp_dist_ob, overwrite = TRUE)
  r_bossen2_bin    <- terra::mask(r_bossen1_bin, dist_oud_bos_opp <= 1000)
  unlink(tmp_dist_ob)
} else {
  r_bossen2_bin    <- r_bossen1_bin
}

# Afstandstoetsing hellingen (<= 1000m)
if (!all(is.na(suppressWarnings(terra::minmax(r_helling_bin))))) {
  tmp_dist_hel       <- tempfile(fileext = ".tif")
  dist_helling_opp   <- terra::distance(r_helling_bin, filename = tmp_dist_hel, overwrite = TRUE)
  r_bossen_gefilterd <- terra::mask(r_bossen2_bin, dist_helling_opp <= 1000)
  unlink(tmp_dist_hel)
} else {
  r_bossen_gefilterd <- r_bossen2_bin
}

# Bosclustering op 10 ha (50m netwerk)
cl_bos_10ha_opp <- cluster_filter_compleet(
  masker     = r_bossen_gefilterd,
  opp_laag   = r_bossen_gefilterd,
  drempel_m2 = 10 * 10000, 
  dist_m     = 50,
  werkelijk  = FALSE
)
vh_bossen3_mask <- cl_bos_10ha_opp$raster

r_loof_frac    <- terra::clamp(1.0 - r_naaldbos_opp_clean, lower = 0.0, upper = 1.0)
vh_bossen3_opp <- terra::mask(r_loof_frac, vh_bossen3_mask)

rm(r_naaldbos_opp_clean, r_naaldbos_dominant, r_bossen1_bin, r_bossen2_bin, r_bossen_gefilterd, 
   cl_bos_10ha_opp, r_loof_frac, r_oud_bos_bin, r_helling_bin)
gc()

# 8. STAP 2: OPEN GROENGEBIEDEN CLUSTEREN & BOSRANDEN BEPALEN ------------------
message("-> Open groengebieden clusteren (min. 5 ha) en bosranden identificeren...")

cl_open_5ha_opp <- cluster_filter_compleet(
  masker     = r_groen_open_bin,
  opp_laag   = r_groen_open_bin,
  drempel_m2 = 5 * 10000, 
  dist_m     = 50,
  werkelijk  = FALSE
)
vh_open1_mask <- cl_open_5ha_opp$raster

vh_bosrand2_opp <- template_Vlaanderen

if (!all(is.na(suppressWarnings(terra::minmax(vh_bossen3_mask))))) && !all(is.na(suppressWarnings(terra::minmax(vh_open1_mask))))) {
  r_bin_open1        <- terra::ifel(!is.na(vh_open1_mask) & vh_open1_mask > 0, 1, NA)
  tmp_dist_open      <- tempfile(fileext = ".tif")
  dist_tot_open_opp  <- terra::distance(r_bin_open1, filename = tmp_dist_open, overwrite = TRUE)
  
  r_bosrand_rauw_opp <- terra::mask(vh_bossen3_mask, dist_tot_open_opp <= 15)
  unlink(tmp_dist_open)
  rm(r_bin_open1)
  
  cl_bosrand_opp <- cluster_filter_compleet(
    masker     = r_bosrand_rauw_opp,
    opp_laag   = r_bosrand_rauw_opp,
    drempel_m2 = 0.25 * 10000, # 0.25 ha (100m randlengte)
    dist_m     = 10,           
    werkelijk  = FALSE
  )
  vh_bosrand2_mask <- cl_bosrand_opp$raster
  vh_bosrand2_opp  <- terra::mask(vh_bossen3_opp, vh_bosrand2_mask)
  
  rm(r_bosrand_rauw_opp, cl_bosrand_opp)
}

rm(cl_open_5ha_opp, r_groen_open_bin, vh_open1_mask)
gc()

# 9. STAP 3: BOSDIEPTE-BUFFER VAN 40M VANAF DE BOSRAND -----------------------
message("-> Bosdiepte-buffer van 40m berekenen vanaf de kwalitatieve bosrand...")

if (!all(is.na(suppressWarnings(terra::minmax(vh_bosrand2_opp))))) {
  r_bin_rand2       <- terra::ifel(!is.na(vh_bosrand2_opp) & vh_bosrand2_opp > 0, 1, NA)
  tmp_dist_rand     <- tempfile(fileext = ".tif")
  dist_tot_rand_opp <- terra::distance(r_bin_rand2, filename = tmp_dist_rand, overwrite = TRUE)
  
  vh_bosrand3_mask  <- terra::mask(vh_bossen3_mask, dist_tot_rand_opp <= 40)
  unlink(tmp_dist_rand)
  rm(r_bin_rand2)
  
  b1_clean <- terra::ifel(is.na(vh_bosrand2_opp), 0, vh_bosrand2_opp)
  b2_clean <- terra::ifel(is.na(terra::mask(vh_bossen3_opp, vh_bosrand3_mask)), 0, terra::mask(vh_bossen3_opp, vh_bosrand3_mask))
  
  bruto_opp_raw           <- terra::clamp(b1_clean + b2_clean, upper = 1.0)
  vh_leefgebied_bruto_opp <- terra::ifel(bruto_opp_raw > 0, bruto_opp_raw, NA)
  
  rm(b1_clean, b2_clean, bruto_opp_raw, vh_bosrand3_mask)
} else {
  vh_leefgebied_bruto_opp <- template_Vlaanderen
}

r_binair_bruto_opp <- terra::ifel(!is.na(vh_leefgebied_bruto_opp) & vh_leefgebied_bruto_opp > 0, 1, NA)

cl_finaal_opp <- cluster_filter_compleet(
  masker     = r_binair_bruto_opp,
  opp_laag   = vh_leefgebied_bruto_opp,
  drempel_m2 = 5 * 10000, # 5 ha
  dist_m     = 50,
  werkelijk  = TRUE
)

final_opp <- cl_finaal_opp$raster
cl_opp    <- cl_finaal_opp$clusters

if (is.null(cl_opp) || all(is.na(terra::values(cl_opp, mat = FALSE)))) {
  if (!is.null(final_opp) && !all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
    tmp_patch <- tempfile(fileext = ".tif")
    cl_opp    <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE, filename = tmp_patch, overwrite = TRUE)
  } else {
    cl_opp <- template_Vlaanderen
  }
}

rm(r_binair_bruto_opp, vh_leefgebied_bruto_opp, cl_finaal_opp, vh_bossen3_opp, vh_bossen3_mask, vh_bosrand2_opp)
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

message("🏁 LEEFGEBIEDENGRID VLIEGEND HERT SUCCESVOL GEËXPORTEERD VOOR HEEL VLAANDEREN!")
