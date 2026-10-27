# ==============================================================================
# LEEFGEBIEDENGRID VOOR ROERDOMP - HEEL VLAANDEREN (WERKELIJKE OPPERVLAKTE)
# Ecologische analyse rietmoerassen, ecotonen & metapopulatie-eisen (INBO)
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
soort <- "roerdomp"

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

openwater_opp <- if("openwater" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["openwater"]] else template_Vlaanderen
rietland_opp  <- if("rietland" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["rietland"]] else template_Vlaanderen
zegge_opp     <- if("zegge" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["zegge"]] else template_Vlaanderen
moeras_opp    <- if("moeras" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["moeras"]] else template_Vlaanderen
grasland_opp  <- if("grasland" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["grasland"]] else template_Vlaanderen
struweel_opp  <- if("struweel" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["struweel"]] else template_Vlaanderen

rm(resultaten_gegroepeerd, df_nieuw, lijst_oppervlaktes)
gc()

# 6. EXTERNE DATASETS (GROENKAART EN ESTUARIUM 1130) ---------------------------
message("-> Groenkaart synchroniseren voor heel Vlaanderen...")
hooggroen_raw <- rast(here("data/input/ASCI Files/Groenkaart_2021.tif"))
if(!terra::compareGeom(hooggroen_raw, template_Vlaanderen, stopOnError = FALSE)) {
  hooggroen_aligned <- terra::resample(hooggroen_raw, template_Vlaanderen, method = "near")
} else {
  hooggroen_aligned <- hooggroen_raw
}

r_bos_bin <- terra::ifel(hooggroen_aligned == 1, 1, NA)

message("-> Estuarium code (1130) isoleren uit Parquet...")
tabel_estuarium <- ds_vlaanderen %>%
  filter(tolower(CODE) == "1130") %>%
  collect() %>%
  as.data.table()

r_estuarium <- template_Vlaanderen
if (nrow(tabel_estuarium) > 0) {
  tabel_estuarium[, CODE := tolower(trimws(CODE))]
  tabel_est_unique <- unique(tabel_estuarium, by = c("cel_id", "CODE"))
  r_estuarium[tabel_est_unique$cel_id] <- 1
}

stoorlaag <- terra::cover(r_bos_bin, r_estuarium)

rm(hooggroen_raw, hooggroen_aligned, r_bos_bin, r_estuarium, tabel_estuarium)
gc()

# 7. STAP 1: BASISHABITAT BUNDELEN EN STORING UITSLUITEN -----------------------
message("-> Stap 1: Basishabitat bundelen op werkelijke oppervlaktes en stoorlagen uitsluiten...")

opp_lagen <- list(openwater_opp, rietland_opp, zegge_opp, moeras_opp, grasland_opp, struweel_opp)

biotoop_opp_som <- template_Vlaanderen
values(biotoop_opp_som) <- 0

for(i in 1:length(opp_lagen)) {
  clean_opp <- terra::ifel(is.na(opp_lagen[[i]]), 0, opp_lagen[[i]])
  biotoop_opp_som <- biotoop_opp_som + clean_opp
}

biotoop_opp_raw          <- terra::clamp(biotoop_opp_som, upper = 1.0)
biotoop_opp              <- terra::ifel(biotoop_opp_raw > 0, biotoop_opp_raw, NA)

roerdomp_leefgebied1_opp <- terra::mask(biotoop_opp, stoorlaag, inverse = TRUE)
rietland_clean_opp       <- terra::mask(rietland_opp, stoorlaag, inverse = TRUE)
openwater_clean_opp      <- terra::mask(openwater_opp, stoorlaag, inverse = TRUE)

rm(biotoop_opp_som, biotoop_opp_raw, biotoop_opp, stoorlaag, clean_opp, opp_lagen, 
   zegge_opp, moeras_opp, grasland_opp, struweel_opp)
gc()

# 8. STAP 2: OPEN WATER PARTIJEN FILTEREN OP MINIMALE GROOTTE (10 HA) ----------
message("-> Stap 2: Open water partijen filteren op minstens 10 ha...")

drempel_10ha_m2     <- 10 * 10000
r_binair_openw_opp  <- terra::ifel(!is.na(openwater_clean_opp) & openwater_clean_opp > 0, 1, NA)

water_clusters_opp  <- cluster_filter_compleet(
  masker     = r_binair_openw_opp,
  opp_laag   = openwater_clean_opp,
  drempel_m2 = drempel_10ha_m2,
  dist_m     = 50,
  werkelijk  = TRUE
)
r_grote_plassen_opp <- water_clusters_opp$raster

if (!all(is.na(suppressWarnings(terra::minmax(r_grote_plassen_opp))))) {
  r_bin_grote_plassen <- terra::ifel(!is.na(r_grote_plassen_opp) & r_grote_plassen_opp > 0, 1, NA)
  
  tmp_dist_w          <- tempfile(fileext = ".tif")
  dist_tot_water_opp  <- terra::distance(r_bin_grote_plassen, filename = tmp_dist_w, overwrite = TRUE)
  mal_nabij_water_opp <- terra::ifel(dist_tot_water_opp <= 100, 1, NA)
  
  roerdomp_leefgebied2_opp <- terra::mask(roerdomp_leefgebied1_opp, mal_nabij_water_opp)
  
  unlink(tmp_dist_w)
  rm(r_bin_grote_plassen, dist_tot_water_opp, mal_nabij_water_opp)
} else {
  roerdomp_leefgebied2_opp <- template_Vlaanderen
}

rm(water_clusters_opp, r_grote_plassen_opp, r_binair_openw_opp, roerdomp_leefgebied1_opp)
gc()

# 9. STAP 3: RIETKRAGEN / ECOTONEN ANALYSEREN (MIN. 500M / 0.5 HA WATERKANT) ----
message("-> Stap 3: Rietkragen/ecotonen analyseren (min. 0.5 ha waterkant in rietbuffer)...")

drempel_waterkant_m2 <- 5000 # 0.5 ha = 5000 m2

r_binair_riet_opp <- terra::ifel(!is.na(rietland_clean_opp) & rietland_clean_opp > 0, 1, NA)
r_binair_open_opp <- terra::ifel(!is.na(openwater_clean_opp) & openwater_clean_opp > 0, 1, NA)

tmp_buf_riet     <- tempfile(fileext = ".tif")
riet_buffer_opp  <- terra::buffer(r_binair_riet_opp, width = 10, filename = tmp_buf_riet, overwrite = TRUE)
r_waterkant_opp  <- terra::mask(openwater_clean_opp, riet_buffer_opp)

waterkant_gefilterd_opp <- cluster_filter_compleet(
  masker     = r_binair_open_opp & !is.na(r_waterkant_opp),
  opp_laag   = r_waterkant_opp,
  drempel_m2 = drempel_waterkant_m2,
  dist_m     = 50,
  werkelijk  = TRUE
)
r_robuuste_waterkant_opp <- waterkant_gefilterd_opp$raster

unlink(tmp_buf_riet)

if (!all(is.na(suppressWarnings(terra::minmax(r_robuuste_waterkant_opp))))) {
  r_bin_rob_water           <- terra::ifel(!is.na(r_robuuste_waterkant_opp) & r_robuuste_waterkant_opp > 0, 1, NA)
  
  tmp_dist_wk               <- tempfile(fileext = ".tif")
  dist_tot_waterkant_opp    <- terra::distance(r_bin_rob_water, filename = tmp_dist_wk, overwrite = TRUE)
  mal_nabij_waterkant_opp   <- terra::ifel(dist_tot_waterkant_opp <= 100, 1, NA)
  
  roerdomp_leefgebied3_opp  <- terra::mask(roerdomp_leefgebied2_opp, mal_nabij_waterkant_opp)
  
  unlink(tmp_dist_wk)
  rm(r_bin_rob_water, dist_tot_waterkant_opp, mal_nabij_waterkant_opp)
} else {
  roerdomp_leefgebied3_opp  <- template_Vlaanderen
}

rm(riet_buffer_opp, r_waterkant_opp, waterkant_gefilterd_opp, r_robuuste_waterkant_opp, 
   r_binair_riet_opp, r_binair_open_opp, roerdomp_leefgebied2_opp, openwater_clean_opp)
gc()

# 10. STAP 4: NESTKERNEN FILTEREN (MINIMAAL 0.5 HA ZUIVER RIETLAND) -------------
message("-> Stap 4: Nestkernen filteren (min. 0.5 ha zuiver rietland)...")

drempel_05ha_m2     <- 0.5 * 10000
r_binair_riet2_opp  <- terra::ifel(!is.na(rietland_clean_opp) & rietland_clean_opp > 0, 1, NA)

riet_kernen_opp     <- cluster_filter_compleet(
  masker     = r_binair_riet2_opp,
  opp_laag   = rietland_clean_opp,
  drempel_m2 = drempel_05ha_m2,
  dist_m     = 50,
  werkelijk  = TRUE
)
r_riet_kernen_opp   <- riet_kernen_opp$raster

if (!all(is.na(suppressWarnings(terra::minmax(r_riet_kernen_opp))))) {
  r_bin_rietkern          <- terra::ifel(!is.na(r_riet_kernen_opp) & r_riet_kernen_opp > 0, 1, NA)
  
  tmp_dist_rk             <- tempfile(fileext = ".tif")
  dist_tot_rietkern_opp   <- terra::distance(r_bin_rietkern, filename = tmp_dist_rk, overwrite = TRUE)
  mal_nabij_rietkern_opp  <- terra::ifel(dist_tot_rietkern_opp <= 100, 1, NA)
  
  roerdomp_leefgebied4_opp <- terra::mask(roerdomp_leefgebied3_opp, mal_nabij_rietkern_opp)
  
  unlink(tmp_dist_rk)
  rm(r_bin_rietkern, dist_tot_rietkern_opp, mal_nabij_rietkern_opp)
} else {
  roerdomp_leefgebied4_opp <- template_Vlaanderen
}

rm(riet_kernen_opp, r_riet_kernen_opp, r_binair_riet2_opp, roerdomp_leefgebied3_opp, rietland_clean_opp)
gc()

# 11. STAP 5: GROTE METAPOPULATIE-CLUSTERING ON 25 HECTARE ---------------------
message("-> Stap 5: Finale metapopulatie-clustering op 25 ha (500m netwerk)...")

drempel_25ha_m2    <- 25 * 10000
r_binair_leef4_opp <- terra::ifel(!is.na(roerdomp_leefgebied4_opp) & roerdomp_leefgebied4_opp > 0, 1, NA)

finaal_roerdomp_opp <- cluster_filter_compleet(
  masker     = r_binair_leef4_opp,
  opp_laag   = roerdomp_leefgebied4_opp,
  drempel_m2 = drempel_25ha_m2,
  dist_m     = 500,
  werkelijk  = TRUE
)

final_opp <- finaal_roerdomp_opp$raster
cl_opp    <- finaal_roerdomp_opp$clusters

if (is.null(cl_opp) || all(is.na(terra::values(cl_opp, mat = FALSE)))) {
  if (!is.null(final_opp) && !all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
    tmp_patch <- tempfile(fileext = ".tif")
    cl_opp    <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE, filename = tmp_patch, overwrite = TRUE)
  } else {
    cl_opp <- template_Vlaanderen
  }
}

rm(finaal_roerdomp_opp, roerdomp_leefgebied4_opp, r_binair_leef4_opp, ds_vlaanderen)
gc()

# 12. EXPORT NAAR DISK (ENKEL 00_ID_Rasters EN 01_Werkelijke_Oppervlaktes) -----
message("-> Geformatteerde export starten...")

base_dir <- here::here("data/output/Vlaanderen/Rasters_Soorten")

folders <- list(
  id_raster = file.path(base_dir, "00_ID_Rasters"),
  werkelijk = file.path(base_dir, "01_Werkelijke_Oppervlaktes")
)
purrr::walk(folders, ~if (!dir.exists(.x)) dir.create(.x, showWarnings = FALSE, recursive = TRUE))

# 1. Export Werkelijke Oppervlakte Raster
file_path_opp <- file.path(folders$werkelijk, paste0("Habitat_Werkelijke_Oppervlaktes_", soort, "_wv.tif"))
terra::writeRaster(
  final_opp, 
  filename = file_path_opp, 
  overwrite = TRUE, 
  gdal = c("COMPRESS=LZW"), 
  datatype = "FLT4S", 
  NAflag = -9999
)
message(paste("    [OK] Geëxporteerd naar scenariomap:", basename(file_path_opp)))

# 2. Export Analytisch Netwerk ID Raster
file_path_id <- file.path(folders$id_raster, paste0("ID_Netwerken_", soort, "_wv.tif"))
terra::writeRaster(
  cl_opp, 
  filename = file_path_id, 
  overwrite = TRUE, 
  gdal = c("COMPRESS=LZW"), 
  datatype = "INT4U", 
  NAflag = 0
)
message(paste("    [OK] Analytisch ID-raster geëxporteerd:", basename(file_path_id)))

message("🏁 LEEFGEBIEDENGRID ROERDOMP SUCCESVOL GEËXPORTEERD VOOR HEEL VLAANDEREN!")
