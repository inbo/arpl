# ==============================================================================
# LEEFGEBIEDENGRID VOOR WULP - HEEL VLAANDEREN (WERKELIJKE OPPERVLAKTE)
# Ecologische analyse valleigraslanden, verbossingsfilter & 100ha clustering (INBO)
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
soort <- "wulp"

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

biotoop_grasland_nat_opp <- if("biotoop_grasland_nat" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["biotoop_grasland_nat"]] else template_Vlaanderen
biotoop_overig_opp       <- if("biotoop_overig" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["biotoop_overig"]] else template_Vlaanderen

message("-> Akkers (b%) isoleren uit Parquet...")
tabel_akkers <- ds_vlaanderen %>%
  filter(str_detect(tolower(CODE), "^b")) %>%
  collect() %>%
  as.data.table()

r_akkers_opp <- template_Vlaanderen
if (nrow(tabel_akkers) > 0) {
  tabel_akkers[, CODE := tolower(trimws(CODE))]
  tabel_akkers_unique <- unique(tabel_akkers, by = c("cel_id", "CODE"))
  tabel_akkers_som    <- tabel_akkers_unique[, .(Oppervlakte = pmin(sum(BWK_FRAC, na.rm = TRUE), 1.0)), by = .(cel_id)]
  r_akkers_opp[tabel_akkers_som$cel_id] <- tabel_akkers_som$Oppervlakte
}

rm(resultaten_gegroepeerd, df_nieuw, lijst_oppervlaktes, tabel_akkers)
gc()

# 6. EXTERNE DATASETS INLEZEN EN SYNCHRONISEREN -------------------------------
message("-> Externe omgevingslagen en groenkaart inlezen...")

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

r_ovstr      <- load_and_sync(here("data/input/Raster_Vlaanderen/vlaanderen_ovstrg_10m.tif"))
r_drainage   <- load_and_sync(here("data/input/Raster_Vlaanderen/vlaanderen_drainage_n2khab_10m.tif"))
r_profiel    <- load_and_sync(here("data/input/Raster_Vlaanderen/Vlaanderen_profiel_10m.tif"))
r_groenkaart <- load_and_sync(here("data/input/ASCI Files/Groenkaart_2021.tif"))

message("-> Verbossingsfilter berekenen: hooggroencomplexen groter dan 20 are uitsluiten...")
r_hooggroen_bin <- terra::ifel(r_groenkaart == 1, 1, NA)

tmp_patch_groen    <- tempfile(fileext = ".tif")
verbossing_clusters <- terra::patches(r_hooggroen_bin, directions = 8, zeroAsNA = TRUE, filename = tmp_patch_groen, overwrite = TRUE)

df_groen_cl <- terra::as.data.frame(verbossing_clusters, cells = TRUE)
r_verbossing_masker <- template_Vlaanderen

if (nrow(df_groen_cl) > 0) {
  colnames(df_groen_cl) <- c("cell", "ID")
  stats_groen_count     <- data.table(ID = df_groen_cl$ID)[, .(count = .N), by = ID]
  
  # Groen-patches van minstens 20 are (20 cellen van 10x10m = 2000 m2)
  grote_groen_ids       <- stats_groen_count[count >= 20, ID]
  if (length(grote_groen_ids) > 0) {
    r_verbossing_masker <- terra::ifel(verbossing_clusters %in% grote_groen_ids, 1, NA)
  }
}

unlink(tmp_patch_groen)
rm(r_hooggroen_bin, verbossing_clusters, df_groen_cl, r_groenkaart)
gc()

# 7. STAP 2: VALLEIFILTER & AKKER-EXCLUSIE --------------------------------------
message("-> Valleicriteria configureren op basis van biologische klassen...")

gewenste_drainage_letters <- c("d", "e", "f", "g", "h", "i", "e-f", "h-i", "e-i")
gewenste_profiel_types    <- c("p", "p+x")
gewenste_overstroming     <- c(1, 2)

drain_cats       <- terra::cats(r_drainage)[[1]]
ids_drainage_nat <- drain_cats$value[drain_cats$Label %in% gewenste_drainage_letters]

prof_cats        <- terra::cats(r_profiel)[[1]]
ids_profiel_val  <- prof_cats$value[prof_cats$Label %in% gewenste_profiel_types]

drainage_nat_ok  <- r_drainage %in% ids_drainage_nat
profiel_vallei_ok <- r_profiel %in% ids_profiel_val
overstroming_ok  <- r_ovstr %in% gewenste_overstroming

vallei_indicator <- overstroming_ok | drainage_nat_ok | profiel_vallei_ok
vallei_masker    <- terra::ifel(vallei_indicator, 1, NA)

rm(r_drainage, r_profiel, r_ovstr, drain_cats, prof_cats, ids_drainage_nat, ids_profiel_val, 
   drainage_nat_ok, profiel_vallei_ok, overstroming_ok, vallei_indicator, ds_vlaanderen)

r_akkers_masker             <- terra::ifel(!is.na(r_akkers_opp) & r_akkers_opp > 0, 1, NA)
r_grasland_gecorrigeerd_opp <- terra::mask(biotoop_grasland_nat_opp, r_akkers_masker, inverse = TRUE)
r_grasland_vallei_opp       <- terra::mask(r_grasland_gecorrigeerd_opp, vallei_masker)

g_clean <- terra::ifel(is.na(r_grasland_vallei_opp), 0, r_grasland_vallei_opp)
o_clean <- terra::ifel(is.na(biotoop_overig_opp), 0, biotoop_overig_opp)

som_raw <- g_clean + o_clean
som_cl  <- terra::clamp(som_raw, upper = 1.0)

r_biotoop_finaal_opp_raw <- terra::ifel(som_cl > 0, som_cl, NA)
r_biotoop_finaal_opp     <- terra::mask(r_biotoop_finaal_opp_raw, r_verbossing_masker, inverse = TRUE)

rm(biotoop_grasland_nat_opp, biotoop_overig_opp, r_akkers_opp, r_akkers_masker,
   r_grasland_gecorrigeerd_opp, r_grasland_vallei_opp, vallei_masker, r_verbossing_masker, 
   g_clean, o_clean, som_raw, som_cl, r_biotoop_finaal_opp_raw)
gc()

# 8. STAP 3: RUIMTELIJKE CLUSTERING (100M FUZZY, MIN. 100 HA) -------------------
message("-> Landschappelijke clustering (100m fuzzy, min. 100 ha)...")

drempel_100ha_m2     <- 100 * 10000
r_binair_biotoop_opp <- terra::ifel(!is.na(r_biotoop_finaal_opp) & r_biotoop_finaal_opp > 0, 1, NA)

wulp_clusters_opp <- cluster_filter_compleet(
  masker     = r_binair_biotoop_opp,
  opp_laag   = r_biotoop_finaal_opp,
  drempel_m2 = drempel_100ha_m2,
  dist_m     = 100,     
  werkelijk  = TRUE    
)

final_opp <- wulp_clusters_opp$raster
cl_opp    <- wulp_clusters_opp$clusters

if (is.null(cl_opp) || all(is.na(terra::values(cl_opp, mat = FALSE)))) {
  if (!is.null(final_opp) && !all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
    tmp_patch <- tempfile(fileext = ".tif")
    cl_opp    <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE, filename = tmp_patch, overwrite = TRUE)
  } else {
    cl_opp <- template_Vlaanderen
  }
}

rm(r_binair_biotoop_opp, r_biotoop_finaal_opp, wulp_clusters_opp)
gc()

# 9. EXPORT NAAR DISK (ENKEL 00_ID_Rasters EN 01_Werkelijke_Oppervlaktes) ------
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

message("🏁 LEEFGEBIEDENGRID WULP SUCCESVOL GEËXPORTEERD VOOR HEEL VLAANDEREN!")
