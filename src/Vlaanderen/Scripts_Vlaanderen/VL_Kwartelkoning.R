# ==============================================================================
# LEEFGEBIEDENGRID VOOR KWARTELKONING - HEEL VLAANDEREN (WERKELIJKE OPPERVLAKTE)
# Ecologische analyse valleigronden, hydrografische profilering & graslandaandeel (INBO)
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

cluster_filter_compleet <- function(masker, opp_laag = NULL, drempel_m2, dist_m, werkelijk = FALSE) {
  if (is.null(masker) || terra::global(is.na(masker), "sum")[[1]] == terra::ncell(masker)) {
    return(list(raster = masker * NA, clusters = masker * NA))
  }
  
  r_binair <- terra::ifel(!is.na(masker) & masker > 0, 1, NA)
  if (terra::global(!is.na(r_binair), "sum")[[1]] == 0) {
    return(list(raster = masker * NA, clusters = masker * NA))
  }
  
  tmp_patch <- tempfile(fileext = ".tif")
  
  if (dist_m > 0) {
    stral_cellen <- ceiling((dist_m / 2) / 10)
    f_matrix     <- matrix(1, nrow = (2 * stral_cellen + 1), ncol = (2 * stral_cellen + 1))
    
    tmp_focal  <- tempfile(fileext = ".tif")
    r_buffered <- terra::focal(r_binair, w = f_matrix, fun = "max", na.rm = TRUE, 
                               filename = tmp_focal, overwrite = TRUE)
    r_buffered <- terra::ifel(r_buffered > 0, 1, NA)
    
    cl_network <- terra::patches(r_buffered, directions = 4, zeroAsNA = TRUE, 
                                 filename = tmp_patch, overwrite = TRUE)
    cl_biotoop_only <- terra::mask(cl_network, masker)
    
    unlink(tmp_focal)
  } else {
    cl_biotoop_only <- terra::patches(r_binair, directions = 8, zeroAsNA = TRUE, 
                                      filename = tmp_patch, overwrite = TRUE)
  }
  
  df_cl <- terra::as.data.frame(cl_biotoop_only, cells = TRUE)
  if (nrow(df_cl) == 0) {
    unlink(tmp_patch)
    return(list(raster = masker * NA, clusters = masker * NA))
  }
  colnames(df_cl) <- c("cell", "ID")
  
  if (werkelijk && !is.null(opp_laag)) {
    opp_waarden <- terra::extract(opp_laag, df_cl$cell)[[1]]
    dt_calc     <- data.table(ID = df_cl$ID, Waarde = opp_waarden)
  } else {
    dt_calc     <- data.table(ID = df_cl$ID, Waarde = 1)
  }
  
  stats_dt    <- dt_calc[!is.na(ID) & !is.na(Waarde), .(Area_m2 = sum(Waarde, na.rm = TRUE) * 100), by = ID]
  voldoet_ids <- stats_dt[Area_m2 >= drempel_m2, ID]
  
  unlink(tmp_patch)
  rm(df_cl, dt_calc, stats_dt)
  gc()
  
  if (length(voldoet_ids) == 0) return(list(raster = masker * NA, clusters = masker * NA))
  
  cl_finaal <- terra::ifel(cl_biotoop_only %in% voldoet_ids, cl_biotoop_only, NA)
  r_finaal  <- terra::mask(masker, cl_finaal)
  
  return(list(raster = r_finaal, clusters = cl_finaal))
}

# 3. SOORT INFORMATIE OPHALEN --------------------------------------------------
df <- read_excel(here("data/input/Excel_files/Soorten_bwk_afstanden.xlsx"))
soort <- "kwartelkoning"

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

grasland_opp <- if("grasland" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["grasland"]] else template_Vlaanderen
moeras_opp   <- if("moeras" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["moeras"]] else template_Vlaanderen

rm(resultaten_gegroepeerd, df_nieuw, lijst_oppervlaktes)
gc()

# 6. EXTERNE STRONGLAYERS (GROENKAART EN TERRILS) ------------------------------
message("-> Groenkaart synchroniseren voor heel Vlaanderen...")
hooggroen_raw <- rast(here("data/input/ASCI Files/Groenkaart_2021.tif"))
if(!terra::compareGeom(hooggroen_raw, template_Vlaanderen, stopOnError = FALSE)) {
  hooggroen_sync <- terra::resample(hooggroen_raw, template_Vlaanderen, method = "near")
} else {
  hooggroen_sync <- hooggroen_raw
}

tf_groen_sync <- tempfile(fileext = ".tif")
hooggroen_sync <- terra::writeRaster(hooggroen_sync, filename = tf_groen_sync, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
rm(hooggroen_raw)

message("-> Terrils (kg%) isoleren uit Parquet...")
tabel_terril <- ds_vlaanderen %>%
  filter(str_detect(tolower(CODE), "^kg")) %>%
  collect() %>%
  as.data.table()

r_terril_opp <- template_Vlaanderen
if (nrow(tabel_terril) > 0) {
  tabel_terril[, CODE := tolower(trimws(CODE))]
  tabel_terril_unique <- unique(tabel_terril, by = c("cel_id", "CODE"))
  tabel_terril_som    <- tabel_terril_unique[, .(Oppervlakte = pmin(sum(BWK_FRAC, na.rm = TRUE), 1.0)), by = .(cel_id)]
  r_terril_opp[tabel_terril_som$cel_id] <- tabel_terril_som$Oppervlakte
}
rm(tabel_terril)

# 7. STAP 1: ONGESCHIKTE ZONES UITSLUITEN --------------------------------------
message("-> Stap 1: Initieel biotoop opbouwen en ongeschikte zones uitsluiten...")

g_clean <- terra::ifel(is.na(grasland_opp), 0, grasland_opp)
m_clean <- terra::ifel(is.na(moeras_opp), 0, moeras_opp)

biotoop_initieel_frac_raw <- terra::clamp(g_clean + m_clean, upper = 1.0)
biotoop_initieel_frac     <- terra::ifel(biotoop_initieel_frac_raw > 0, biotoop_initieel_frac_raw, NA)

ongeschikt_masker     <- (hooggroen_sync == 1) | (!is.na(r_terril_opp) & r_terril_opp > 0)
biotoop_fase1_frac    <- terra::mask(biotoop_initieel_frac, ongeschikt_masker, inverse = TRUE)

rm(g_clean, m_clean, biotoop_initieel_frac_raw, biotoop_initieel_frac, ongeschikt_masker, r_terril_opp)
gc()

# 8. STAP 2: VALLEIGRONDEN EN HYDROGRAFIE STUURKAART ---------------------------
message("-> Stap 2: Profiel-, drainage- en overstromingsstuurkaarten koppelen...")

raw_profiel      <- rast(here("data/input/Raster_Vlaanderen/Vlaanderen_profiel_10m.tif"))
raw_drainage     <- rast(here("data/input/Raster_Vlaanderen/vlaanderen_drainage_n2khab_10m.tif"))
raw_overstroming <- rast(here("data/input/Raster_Vlaanderen/vlaanderen_ovstrg_10m.tif")) 

if(!terra::compareGeom(raw_profiel, template_Vlaanderen, stopOnError = FALSE)) {
  profiel_sync      <- terra::resample(raw_profiel, template_Vlaanderen, method = "near")
  drainage_sync     <- terra::resample(raw_drainage, template_Vlaanderen, method = "near")
  overstroming_sync <- terra::resample(raw_overstroming, template_Vlaanderen, method = "near")
} else {
  profiel_sync      <- raw_profiel
  drainage_sync     <- raw_drainage
  overstroming_sync <- raw_overstroming
}

gewenste_profielen  <- c("p", "p+x")
legenda_profiel     <- terra::cats(profiel_sync)[[1]]
profiel_nummers     <- legenda_profiel$value[tolower(trimws(legenda_profiel$Label)) %in% tolower(gewenste_profielen)]

gewenste_drainages  <- c("d", "e", "f", "g", "h", "i", "e-f", "h-i", "e-i") 
legenda_drainage    <- terra::cats(drainage_sync)[[1]]
drainage_nummers    <- legenda_drainage$value[tolower(trimws(legenda_drainage$Label)) %in% tolower(gewenste_drainages)]

overstroming_nummers <- c(1, 2)

mask_profiel      <- profiel_sync %in% profiel_nummers
mask_drainage     <- drainage_sync %in% drainage_nummers
mask_overstroming <- overstroming_sync %in% overstroming_nummers

vallei_geschikt        <- (mask_profiel == 1) | (mask_drainage == 1) | (mask_overstroming == 1)
vallei_geschikt_binair <- terra::ifel(vallei_geschikt == 1, 1, NA)

biotoop_fase2_frac     <- terra::mask(biotoop_fase1_frac, vallei_geschikt_binair)

rm(raw_profiel, raw_drainage, raw_overstroming, profiel_sync, drainage_sync, overstroming_sync, 
   mask_profiel, mask_drainage, mask_overstroming, vallei_geschikt, vallei_geschikt_binair, biotoop_fase1_frac)
gc()

# 9. STAP 3: MAJORITY FILTER (20M STRUKTUURVOORKOMEN) --------------------------
message("-> Stap 3: Majority/focal-structuurfilter toepassen op 20m straal...")

f_matrix               <- terra::focalMat(template_Vlaanderen, d = 20, type = "circle")
f_matrix[f_matrix > 0] <- 1
drempel                <- ceiling(sum(f_matrix) / 2)

binaire_aanwezigheid   <- terra::ifel(!is.na(biotoop_fase2_frac) & biotoop_fase2_frac > 0, 1, 0)

temp_focal_file        <- tempfile(fileext = ".tif")
buurt_som              <- terra::focal(binaire_aanwezigheid, w = f_matrix, fun = "sum", na.rm = TRUE, filename = temp_focal_file, overwrite = TRUE)

vorm_fase3             <- terra::ifel(buurt_som >= drempel, 1, NA)
biotoop_fase3_frac     <- terra::mask(biotoop_fase2_frac, vorm_fase3)

unlink(temp_focal_file)
rm(binaire_aanwezigheid, buurt_som, vorm_fase3, biotoop_fase2_frac)
gc()

# 10. STAP 4: MINIMAAL 3 HA FILTER (50M FUZZY NETWERK) ------------------------
message("-> Stap 4: Minimaal 3 ha clusterdrempel toepassen...")

r_binair_fase3_frac <- terra::ifel(!is.na(biotoop_fase3_frac) & biotoop_fase3_frac > 0, 1, NA)

fase4_lijst_frac    <- cluster_filter_compleet(
  masker     = r_binair_fase3_frac, 
  opp_laag   = biotoop_fase3_frac, 
  drempel_m2 = 30000, # 3 ha
  dist_m     = 50, 
  werkelijk  = TRUE
)

biotoop_fase4_frac <- fase4_lijst_frac$raster

rm(r_binair_fase3_frac, biotoop_fase3_frac, fase4_lijst_frac)
gc()

# 11. STAP 5: HOOGGROEN BUFFERFILTER (100M AFSTANDSTOETS) --------------------
message("-> Stap 5: Hooggroen-storing (20 are cluster + 100m buffer) uitsluiten...")

r_groen_bin <- terra::ifel(!is.na(hooggroen_sync) & hooggroen_sync > 0, 1, NA)

w_matrix        <- matrix(1, nrow = 5, ncol = 5)
temp_groen_file <- tempfile(fileext = ".tif")

groen_som <- terra::focal(
  r_groen_bin, 
  w = w_matrix, 
  fun = "sum", 
  na.rm = TRUE,
  filename = temp_groen_file,
  overwrite = TRUE
)

r_groen_20are            <- terra::ifel(groen_som >= 12 & r_groen_bin == 1, 1, NA)
tmp_buf_g                <- tempfile(fileext = ".tif")
buf_groen                <- terra::buffer(r_groen_20are, width = 100, filename = tmp_buf_g, overwrite = TRUE)
masker_groen_100m_cluster <- terra::ifel(!is.na(buf_groen) & buf_groen > 0, 1, NA)

biotoop_fase5_frac       <- terra::mask(biotoop_fase4_frac, masker_groen_100m_cluster, inverse = TRUE)

unlink(c(temp_groen_file, tmp_buf_g, tf_groen_sync))
rm(r_groen_bin, groen_som, r_groen_20are, buf_groen, masker_groen_100m_cluster, hooggroen_sync, biotoop_fase4_frac)
gc()

# 12. STAP 6 & 7: FUNCTIONEEL LEEFGEBIED (10 HA) & GRASLANDFRACTIE-CONTROLE (>= 75%)
message("-> Stap 6 & 7: Functioneel leefgebied (10 ha) & graslandaandeel (>= 75%) toetsen...")

r_binair_fase5_frac <- terra::ifel(!is.na(biotoop_fase5_frac) & biotoop_fase5_frac > 0, 1, NA)

fase6_lijst_frac <- cluster_filter_compleet(
  masker     = r_binair_fase5_frac, 
  opp_laag   = biotoop_fase5_frac, 
  drempel_m2 = 100000, # 10 ha
  dist_m     = 50, 
  werkelijk  = TRUE
)

biotoop_fase6_frac <- fase6_lijst_frac$raster
cl_id_fase6_frac   <- fase6_lijst_frac$clusters

if (!all(is.na(suppressWarnings(terra::minmax(cl_id_fase6_frac))))) {
  cl_id_fase6_frac_gefilterd <- terra::mask(cl_id_fase6_frac, biotoop_fase6_frac)
  
  tot_biotoop_frac  <- terra::zonal(biotoop_fase5_frac, cl_id_fase6_frac_gefilterd, fun = "sum", na.rm = TRUE)
  tot_grasland_frac <- terra::zonal(grasland_opp, cl_id_fase6_frac_gefilterd, fun = "sum", na.rm = TRUE)
  colnames(tot_biotoop_frac)  <- c("ID", "Tot_Biotoop")
  colnames(tot_grasland_frac) <- c("ID", "Tot_Grasland")
  
  dt_bio  <- as.data.table(tot_biotoop_frac)
  dt_gras <- as.data.table(tot_grasland_frac)
  
  dt_check <- dt_gras[dt_bio, on = "ID"]
  dt_check[is.na(Tot_Grasland), Tot_Grasland := 0]
  dt_check[, Grasland_Pct := (Tot_Grasland / Tot_Biotoop) * 100]
  
  ids_ok_frac <- dt_check[Grasland_Pct >= 75, ID]
  
  if(length(ids_ok_frac) > 0) {
    final_opp <- terra::mask(biotoop_fase6_frac, cl_id_fase6_frac_gefilterd %in% ids_ok_frac)
    cl_opp    <- terra::mask(cl_id_fase6_frac_gefilterd, cl_id_fase6_frac_gefilterd %in% ids_ok_frac)
  } else {
    final_opp <- template_Vlaanderen
    cl_opp    <- template_Vlaanderen
  }
} else {
  final_opp <- template_Vlaanderen
  cl_opp    <- template_Vlaanderen
}

rm(ds_vlaanderen, biotoop_fase5_frac, biotoop_fase6_frac, cl_id_fase6_frac, grasland_opp, moeras_opp)
gc()

# 13. EXPORT NAAR DISK (ENKEL 00_ID_Rasters EN 01_Werkelijke_Oppervlaktes) -----
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

message("🏁 LEEFGEBIEDENGRID KWARTELKONING SUCCESVOL GEËXPORTEERD VOOR HEEL VLAANDEREN!")
