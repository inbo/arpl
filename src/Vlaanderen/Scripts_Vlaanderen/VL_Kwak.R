# ==============================================================================
# LEEFGEBIEDENGRID VOOR KWAK - HEEL VLAANDEREN (WERKELIJKE OPPERVLAKTE)
# Ecologische analyse broedbiotoop, moeras-watercomplexen & commuterkoppeling (INBO)
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
soort <- "kwak"

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

openwater_opp     <- if("openwater" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["openwater"]] else template_Vlaanderen
voortplanting0_opp <- if("voortplanting0" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["voortplanting0"]] else template_Vlaanderen
not_eikenbos_opp  <- if("not_eikenbos" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["not_eikenbos"]] else template_Vlaanderen
moeras_opp        <- if("moeras" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["moeras"]] else template_Vlaanderen

rm(ds_vlaanderen, resultaten_gegroepeerd, df_nieuw, lijst_oppervlaktes)
gc()

# 6. VOORTPLANTINGSHABITAT FILTEREN EN CLUSTEREN (5 HA / 100M) -----------------
message("-> Broedplaatsen filteren (min. 5 ha) en eikenbos-uitsluiting toepassen...")

r_kwak_voortplanting3_opp <- template_Vlaanderen

if (!all(is.na(suppressWarnings(terra::minmax(voortplanting0_opp))))) {
  r_binair_voortp_opp <- terra::ifel(!is.na(voortplanting0_opp) & voortplanting0_opp > 0, 1, NA)
  
  cl_voortp_opp <- cluster_filter_compleet(
    masker     = r_binair_voortp_opp,
    opp_laag   = voortplanting0_opp,
    drempel_m2 = 5 * 10000, # 5 ha
    dist_m     = 100,       
    werkelijk  = TRUE
  )
  r_kwak_voortplanting2_opp <- cl_voortp_opp$raster
  
  if (!all(is.na(suppressWarnings(terra::minmax(r_kwak_voortplanting2_opp))))) {
    if (!all(is.na(suppressWarnings(terra::minmax(not_eikenbos_opp))))) {
      r_kwak_voortplanting3_opp <- terra::mask(r_kwak_voortplanting2_opp, not_eikenbos_opp > 0, inverse = TRUE)
    } else {
      r_kwak_voortplanting3_opp <- r_kwak_voortplanting2_opp
    }
  }
  rm(r_binair_voortp_opp, cl_voortp_opp, r_kwak_voortplanting2_opp)
}

gc()

# 7. FOERAGEERHABITAT FILTEREN EN PROCENTUELE MOERASTOETS (MIN. 20% WATER) -----
message("-> Foerageergebieden clusteren (min. 10 ha) en procentuele water-in-moeras toets toepassen...")

r_kwak_foerageer1_opp <- template_Vlaanderen

if (!all(is.na(suppressWarnings(terra::minmax(moeras_opp))))) {
  r_binair_moeras_opp <- terra::ifel(!is.na(moeras_opp) & moeras_opp > 0, 1, NA)
  
  cl_moeras_opp <- cluster_filter_compleet(
    masker     = r_binair_moeras_opp,
    opp_laag   = moeras_opp,
    drempel_m2 = 10 * 10000, # 10 ha
    dist_m     = 50,          
    werkelijk  = TRUE
  )
  r_moeras1_ids_opp <- cl_moeras_opp$clusters
  
  if (!all(is.na(suppressWarnings(terra::minmax(r_moeras1_ids_opp))))) {
    clean_moeras_opp <- terra::ifel(is.na(moeras_opp), 0, moeras_opp)
    clean_water_opp  <- terra::ifel(is.na(openwater_opp), 0, openwater_opp)
    
    som_moeras_opp <- terra::zonal(clean_moeras_opp, r_moeras1_ids_opp, fun = "sum", na.rm = TRUE)
    som_water_opp  <- terra::zonal(clean_water_opp, r_moeras1_ids_opp, fun = "sum", na.rm = TRUE)
    
    dt_toets_opp <- data.table(
      ClusterID = som_moeras_opp[[1]], 
      Moeras_ha = som_moeras_opp[[2]] * 0.01, 
      Water_ha  = som_water_opp[[2]] * 0.01
    )
    dt_toets_opp[, Water_Aandeel := Water_ha / Moeras_ha]
    goedgekeurde_foer_ids_opp <- dt_toets_opp[Water_Aandeel >= 0.20, ClusterID]
    
    if (length(goedgekeurde_foer_ids_opp) > 0) {
      r_kwak_foerageer1_opp <- terra::mask(moeras_opp, r_moeras1_ids_opp %in% goedgekeurde_foer_ids_opp)
    }
    rm(clean_moeras_opp, clean_water_opp, som_moeras_opp, som_water_opp, dt_toets_opp)
  }
  rm(r_binair_moeras_opp, cl_moeras_opp, r_moeras1_ids_opp)
}

gc()

# 8. KOPPELING (5 KM COMMUTER-AFSTAND) EN FINALE CLUSTERING -------------------
message("-> Wederzijdse 5 km commuter-toets toepassen en netwerk opbouwen...")

if (!all(is.na(suppressWarnings(terra::minmax(r_kwak_voortplanting3_opp)))) && !all(is.na(suppressWarnings(terra::minmax(r_kwak_foerageer1_opp))))) {
  r_binair_vp3_opp  <- terra::ifel(!is.na(r_kwak_voortplanting3_opp) & r_kwak_voortplanting3_opp > 0, 1, NA)
  r_binair_foer_opp <- terra::ifel(!is.na(r_kwak_foerageer1_opp) & r_kwak_foerageer1_opp > 0, 1, NA)
  
  tmp_dist_foer       <- tempfile(fileext = ".tif")
  dist_tot_foer_opp   <- terra::distance(r_binair_foer_opp, filename = tmp_dist_foer, overwrite = TRUE)
  
  tmp_dist_vp         <- tempfile(fileext = ".tif")
  dist_tot_voortp_opp <- terra::distance(r_binair_vp3_opp, filename = tmp_dist_vp, overwrite = TRUE)
  
  r_voortplanting1_opp <- terra::mask(r_kwak_voortplanting3_opp, dist_tot_foer_opp <= 5000)
  r_foerageer1_opp     <- terra::mask(r_kwak_foerageer1_opp, dist_tot_voortp_opp <= 5000)
  
  vpt_clean  <- terra::ifel(is.na(r_voortplanting1_opp), 0, r_voortplanting1_opp)
  foer_clean <- terra::ifel(is.na(r_foerageer1_opp), 0, r_foerageer1_opp)
  
  r_kwak_leefgebied1_opp_raw <- terra::clamp(vpt_clean + foer_clean, upper = 1.0)
  r_kwak_leefgebied1_opp     <- terra::ifel(r_kwak_leefgebied1_opp_raw > 0, r_kwak_leefgebied1_opp_raw, NA)
  r_binair_leef1_opp         <- terra::ifel(!is.na(r_kwak_leefgebied1_opp) & r_kwak_leefgebied1_opp > 0, 1, NA)
  
  cl_finaal_opp <- cluster_filter_compleet(
    masker     = r_binair_leef1_opp,
    opp_laag   = r_kwak_leefgebied1_opp,
    drempel_m2 = 10 * 10000, # 10 ha
    dist_m     = 50,          
    werkelijk  = TRUE
  )
  final_opp <- cl_finaal_opp$raster
  cl_opp    <- cl_finaal_opp$clusters
  
  unlink(c(tmp_dist_foer, tmp_dist_vp))
  rm(r_binair_vp3_opp, r_binair_foer_opp, dist_tot_foer_opp, dist_tot_voortp_opp, 
     r_voortplanting1_opp, r_foerageer1_opp, vpt_clean, foer_clean, r_kwak_leefgebied1_opp_raw, r_kwak_leefgebied1_opp, r_binair_leef1_opp, cl_finaal_opp)
} else {
  final_opp <- template_Vlaanderen
  cl_opp    <- template_Vlaanderen
}

rm(r_kwak_voortplanting3_opp, r_kwak_foerageer1_opp, openwater_opp, moeras_opp, voortplanting0_opp, not_eikenbos_opp)
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

message("🏁 LEEFGEBIEDENGRID KWAK SUCCESVOL GEËXPORTEERD VOOR HEEL VLAANDEREN!")
