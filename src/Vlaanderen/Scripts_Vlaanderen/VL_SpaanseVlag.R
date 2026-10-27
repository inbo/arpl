# ==============================================================================
# LEEFGEBIEDENGRID VOOR SPAANSE VLAG - HEEL VLAANDEREN (WERKELIJKE OPPERVLAKTE)
# Ecologische analyse voortplantingshabitat x bosrand-contactzones (INBO)
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
soort <- "spaansevlag"

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

voortplanting_bwk_opp <- if("voortplanting_bwk" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["voortplanting_bwk"]] else template_Vlaanderen
bos_bwk_opp          <- if("bos_bwk" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["bos_bwk"]] else template_Vlaanderen

rm(ds_vlaanderen, resultaten_gegroepeerd, df_nieuw, lijst_oppervlaktes)
gc()

# 6. STAP 1: VOORTPLANTINGSGEBIEDEN AFBAKENEN (RUPSEN) -------------------------
message("-> Stap 1: Voortplantingsgebieden rupsen verwerken op werkelijke oppervlaktes (min. 1 ha)...")

r_binair_vp_opp              <- terra::ifel(!is.na(voortplanting_bwk_opp) & voortplanting_bwk_opp > 0, 1, NA)
voortplanting_geclusterd_opp <- cluster_filter_compleet(
  masker     = r_binair_vp_opp,
  opp_laag   = voortplanting_bwk_opp,
  drempel_m2 = 10000, # 1 ha
  dist_m     = 50,
  werkelijk  = TRUE
)
spaansevlag_voortplanting2_opp <- voortplanting_geclusterd_opp$raster

rm(r_binair_vp_opp, voortplanting_geclusterd_opp, voortplanting_bwk_opp)
gc()

# 7. STAP 2: BOSSEN AFBAKENEN (ADULTE VLINDERS) --------------------------------
message("-> Stap 2: Boshabitats verwerken op werkelijke oppervlaktes (min. 1 ha)...")

r_binair_bos_opp     <- terra::ifel(!is.na(bos_bwk_opp) & bos_bwk_opp > 0, 1, NA)
bos_geclusterd_opp   <- cluster_filter_compleet(
  masker     = r_binair_bos_opp,
  opp_laag   = bos_bwk_opp,
  drempel_m2 = 10000, # 1 ha
  dist_m     = 50,
  werkelijk  = TRUE
)
spaansevlag_bos_opp  <- bos_geclusterd_opp$raster

rm(r_binair_bos_opp, bos_geclusterd_opp, bos_bwk_opp)
gc()

# 8. STAP 3: SCHERPE BOSRANDEN (CONTACTZONES) BEPALEN --------------------------
message("-> Stap 3: Contactzones scherp identificeren binnen 10m afstand...")

if (!all(is.na(suppressWarnings(terra::minmax(spaansevlag_bos_opp))))) {
  r_bin_bos           <- terra::ifel(!is.na(spaansevlag_bos_opp) & spaansevlag_bos_opp > 0, 1, NA)
  tmp_dist_bos        <- tempfile(fileext = ".tif")
  dist_tot_bos_opp    <- terra::distance(r_bin_bos, filename = tmp_dist_bos, overwrite = TRUE)
  spaansevlag_bosrand_vpt_opp <- terra::mask(spaansevlag_voortplanting2_opp, dist_tot_bos_opp <= 10)
  unlink(tmp_dist_bos)
  rm(r_bin_bos)
} else {
  spaansevlag_bosrand_vpt_opp <- template_Vlaanderen
}

if (!all(is.na(suppressWarnings(terra::minmax(spaansevlag_voortplanting2_opp))))) {
  r_bin_vp            <- terra::ifel(!is.na(spaansevlag_voortplanting2_opp) & spaansevlag_voortplanting2_opp > 0, 1, NA)
  tmp_dist_vp         <- tempfile(fileext = ".tif")
  dist_tot_vpt_opp    <- terra::distance(r_bin_vp, filename = tmp_dist_vp, overwrite = TRUE)
  spaansevlag_bosrand_bos_opp <- terra::mask(spaansevlag_bos_opp, dist_tot_vpt_opp <= 10)
  unlink(tmp_dist_vp)
  rm(r_bin_vp)
} else {
  spaansevlag_bosrand_bos_opp <- template_Vlaanderen
}

v1_clean <- terra::ifel(is.na(spaansevlag_bosrand_vpt_opp), 0, spaansevlag_bosrand_vpt_opp)
v2_clean <- terra::ifel(is.na(spaansevlag_bosrand_bos_opp), 0, spaansevlag_bosrand_bos_opp)

rand1_opp_raw            <- terra::clamp(v1_clean + v2_clean, upper = 1.0)
spaansevlag_bosrand1_opp <- terra::ifel(rand1_opp_raw > 0, rand1_opp_raw, NA)

rm(spaansevlag_bosrand_vpt_opp, spaansevlag_bosrand_bos_opp, v1_clean, v2_clean, rand1_opp_raw)
gc()

# 9. STAP 4: ECHTE AANEENGESLOTEN LENGTE FILTEREN (>= 100M / 0.25 HA) -----------
message("-> Stap 4: Bosranden filteren op aaneengesloten lengte (>= 100m)...")

if (!all(is.na(suppressWarnings(terra::minmax(spaansevlag_bosrand1_opp))))) {
  r_binair_rand1_opp          <- terra::ifel(!is.na(spaansevlag_bosrand1_opp) & spaansevlag_bosrand1_opp > 0, 1, NA)
  tmp_patch_rand              <- tempfile(fileext = ".tif")
  r_scherpe_rand_clusters_opp <- terra::patches(r_binair_rand1_opp, directions = 8, zeroAsNA = TRUE, filename = tmp_patch_rand, overwrite = TRUE)
  
  df_rand_cl <- terra::as.data.frame(r_scherpe_rand_clusters_opp, cells = TRUE)
  
  if (nrow(df_rand_cl) > 0) {
    colnames(df_rand_cl) <- c("cell", "ID")
    opp_rand_vals        <- terra::extract(spaansevlag_bosrand1_opp, df_rand_cl$cell)[[1]]
    dt_rand              <- data.table(ID = df_rand_cl$ID, Waarde = opp_rand_vals)
    
    stats_rand           <- dt_rand[!is.na(ID) & !is.na(Waarde), .(Area_ha = sum(Waarde, na.rm = TRUE) * 0.01), by = ID]
    goedgekeurde_rand_ids_opp <- stats_rand[Area_ha >= 0.25, ID] # 0.25 ha
    
    if (length(goedgekeurde_rand_ids_opp) > 0) {
      masker_rand_opp          <- r_scherpe_rand_clusters_opp %in% goedgekeurde_rand_ids_opp
      spaansevlag_bosrand2_opp <- terra::mask(spaansevlag_bosrand1_opp, terra::ifel(masker_rand_opp, 1, NA))
    } else {
      spaansevlag_bosrand2_opp <- template_Vlaanderen
    }
    rm(df_rand_cl, dt_rand, stats_rand)
  } else {
    spaansevlag_bosrand2_opp <- template_Vlaanderen
  }
  unlink(tmp_patch_rand)
  rm(r_binair_rand1_opp, r_scherpe_rand_clusters_opp)
} else {
  spaansevlag_bosrand2_opp <- template_Vlaanderen
}

rm(spaansevlag_bosrand1_opp)
gc()

# 10. STAP 5: FINAAL LEEFGEBIED VIA BUFFER-INTERSECTIE (100M) -------------------
message("-> Stap 5: Finaal leefgebied berekenen via 100m bosrandbuffer...")

if (!all(is.na(suppressWarnings(terra::minmax(spaansevlag_bosrand2_opp))))) {
  r_binair_rand2_opp      <- terra::ifel(!is.na(spaansevlag_bosrand2_opp) & spaansevlag_bosrand2_opp > 0, 1, NA)
  tmp_buf_rand            <- tempfile(fileext = ".tif")
  bosrand_buffer_100m_opp <- terra::buffer(r_binair_rand2_opp, width = 100, filename = tmp_buf_rand, overwrite = TRUE)
  
  spaansevlag_leefgebied1_opp <- terra::mask(spaansevlag_voortplanting2_opp, bosrand_buffer_100m_opp)
  spaansevlag_leefgebied2_opp <- terra::mask(spaansevlag_bosrand2_opp, bosrand_buffer_100m_opp)
  
  l1_clean <- terra::ifel(is.na(spaansevlag_leefgebied1_opp), 0, spaansevlag_leefgebied1_opp)
  l2_clean <- terra::ifel(is.na(spaansevlag_leefgebied2_opp), 0, spaansevlag_leefgebied2_opp)
  
  final_opp_raw <- terra::clamp(l1_clean + l2_clean, upper = 1.0)
  final_opp     <- terra::ifel(final_opp_raw > 0, final_opp_raw, NA)
  
  unlink(tmp_buf_rand)
  rm(r_binair_rand2_opp, bosrand_buffer_100m_opp, spaansevlag_leefgebied1_opp, spaansevlag_leefgebied2_opp, l1_clean, l2_clean, final_opp_raw)
} else {
  final_opp <- template_Vlaanderen
}

if (!all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  tmp_patch_cl <- tempfile(fileext = ".tif")
  cl_opp       <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE, filename = tmp_patch_cl, overwrite = TRUE)
} else {
  cl_opp <- template_Vlaanderen
}

rm(spaansevlag_voortplanting2_opp, spaansevlag_bos_opp, spaansevlag_bosrand2_opp)
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

message("🏁 LEEFGEBIEDENGRID SPAANSE VLAG SUCCESVOL GEËXPORTEERD VOOR HEEL VLAANDEREN!")
