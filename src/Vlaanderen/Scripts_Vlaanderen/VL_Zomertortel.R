# ==============================================================================
# LEEFGEBIEDENGRID VOOR ZOMERTORTEL - HEEL VLAANDEREN (WERKELIJKE OPPERVLAKTE)
# Ecologische analyse broedbos x LGP-landbouw x waterbuffers (300m/500m) (INBO)
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
soort <- "zomertortel"

df <- read_excel(here::here("data/input/Excel_files/Soorten_bwk_afstanden.xlsx"))
resultaat <- df %>%
  filter(tolower(trimws(Soort)) == soort) %>%
  select(Type, MinOpp_ha, AfstandBiotopen_m, Dispersiecap_m)

oppervlakte_ha       <- resultaat$MinOpp_ha[1]
afstand_m            <- resultaat$AfstandBiotopen_m[1]
buffer_m             <- resultaat$Dispersiecap_m[1]

foerageer_afstand_m  <- 300
drinkwater_afstand_m <- 500

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

df_namen_sleutel  <- read_csv(here("data/input/Excel_files/BWK_Laag_Namen_2025.csv"), show_col_types = FALSE)
gouden_namenlijst <- tolower(trimws(df_namen_sleutel$Laagnaam))

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
water_opp        <- if("water" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["water"]] else template_Vlaanderen
foerageer_opp    <- if("foerageer" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["foerageer"]] else template_Vlaanderen

rm(resultaten_gegroepeerd, df_nieuw, lijst_oppervlaktes)
gc()

# 6. EXTERNE DATASETS INLEZEN EN SYNCHRONISEREN -------------------------------
message("-> Externe kaarten inlezen (LGP, Watervlakken, Huet-zones, Groenkaart)...")

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

# 1. LANDBOUWPERCELEN (LGP 2025)
r_lgp_vlaanderen <- load_and_sync(here("data/input/Raster_Vlaanderen/vlaanderen_lgp_gwscod_2025_10m.tif"))
zomer_gewascodes <- c(201, 202, 311, 321, 601, 800) 
r_lgp_opp        <- terra::ifel(!is.na(r_lgp_vlaanderen) & r_lgp_vlaanderen %in% zomer_gewascodes, 0.01, NA)

f_bwk_clean  <- terra::ifel(is.na(foerageer_opp), 0, foerageer_opp)
f_lgp_clean  <- terra::ifel(is.na(r_lgp_opp), 0, r_lgp_opp)
foerageer_totaal_opp_raw <- terra::clamp(f_bwk_clean + f_lgp_clean, upper = 1.0)
foerageer_totaal_opp     <- terra::ifel(foerageer_totaal_opp_raw > 0, foerageer_totaal_opp_raw, NA)

rm(r_lgp_vlaanderen, r_lgp_opp, f_bwk_clean, f_lgp_clean, foerageer_totaal_opp_raw, foerageer_opp)
gc()

# 2. WATERVLAKKEN EN HUET-ZONES
r_watervlakken     <- load_and_sync(here("data/input/Raster_Vlaanderen/vlaanderen_watervlakken_2024_10m.tif"))
r_watervlakken_opp <- terra::ifel(!is.na(r_watervlakken) & r_watervlakken > 0, 0.01, NA)

r_huetzon     <- load_and_sync(here("data/input/Raster_Vlaanderen/vlaanderen_huetzon_10m.tif"))
r_huet_opp    <- terra::ifel(!is.na(r_huetzon) & r_huetzon > 0, 0.01, NA)

w_bwk_clean   <- terra::ifel(is.na(water_opp), 0, water_opp)
w_vlak_clean  <- terra::ifel(is.na(r_watervlakken_opp), 0, r_watervlakken_opp)
w_huet_clean  <- terra::ifel(is.na(r_huet_opp), 0, r_huet_opp)

water_totaal_opp_raw <- terra::clamp(w_bwk_clean + w_vlak_clean + w_huet_clean, upper = 1.0)
water_totaal_opp     <- terra::ifel(water_totaal_opp_raw > 0, water_totaal_opp_raw, NA)

rm(r_watervlakken, r_watervlakken_opp, r_huetzon, r_huet_opp, w_bwk_clean, w_vlak_clean, w_huet_clean, water_totaal_opp_raw, water_opp)
gc()

# 3. GROENKAART / BOSFILTER (HOOGGROEN > 20 ARE)
r_groenkaart    <- load_and_sync(here("data/input/ASCI Files/Groenkaart_2021.tif"))
r_hooggroen_bin <- terra::ifel(r_groenkaart == 1, 1, NA)

w_matrix        <- matrix(1, nrow = 5, ncol = 5)
tmp_focal_g     <- tempfile(fileext = ".tif")
groen_som       <- terra::focal(r_hooggroen_bin, w = w_matrix, fun = "sum", na.rm = TRUE, filename = tmp_focal_g, overwrite = TRUE)

r_bos_20are     <- terra::ifel(groen_som >= 12 & r_hooggroen_bin == 1, 1, NA)

unlink(tmp_focal_g)
rm(r_groenkaart, r_hooggroen_bin, groen_som)
gc()

# 4. GROTE DORPSKERNEN BEBOUWINGSFILTER (U% CODES > 10 HA)
bebouwings_lagen <- gouden_namenlijst[grepl("^u[a-z0-9]", gouden_namenlijst) | gouden_namenlijst == "u"]
bebouwings_lagen <- bebouwings_lagen[bebouwings_lagen != "ulm"]

tabel_bebouw <- ds_vlaanderen %>%
  filter(tolower(CODE) %in% bebouwings_lagen) %>%
  collect() %>%
  as.data.table()

r_bebouwing_schoon <- template_Vlaanderen
if (nrow(tabel_bebouw) > 0) {
  tabel_bebouw[, CODE := tolower(trimws(CODE))]
  tabel_bebouw_unique <- unique(tabel_bebouw, by = c("cel_id", "CODE"))
  tabel_bebouw_som    <- tabel_bebouw_unique[, .(Oppervlakte = pmin(sum(BWK_FRAC, na.rm = TRUE), 1.0)), by = .(cel_id)]
  r_bebouwing_schoon[tabel_bebouw_som$cel_id] <- 1
}

res_bebouw_groot <- cluster_filter_compleet(
  masker     = r_bebouwing_schoon,
  opp_laag   = r_bebouwing_schoon,
  drempel_m2 = 100000, # 10 ha
  dist_m     = 100,
  werkelijk  = FALSE
)

r_bebouw_groot <- res_bebouw_groot$raster

if (!is.null(r_bebouw_groot) && !all(is.na(suppressWarnings(terra::minmax(r_bebouw_groot))))) {
  tmp_buf_beb              <- tempfile(fileext = ".tif")
  masker_bebouw_100m_groot <- terra::buffer(r_bebouw_groot, width = 100, filename = tmp_buf_beb, overwrite = TRUE)
  masker_bebouw_100m_groot <- terra::ifel(!is.na(masker_bebouw_100m_groot) & masker_bebouw_100m_groot > 0, 1, NA)
} else {
  masker_bebouw_100m_groot <- template_Vlaanderen
}

rm(tabel_bebouw, r_bebouwing_schoon, res_bebouw_groot, r_bebouw_groot, ds_vlaanderen)
gc()

# 7. STAP 3: BIOTOOP FILTERING EN CLUSTERING -----------------------------------
message("-> Broed-, foerageer- en waterbiotopen filteren op minimale drempels...")

# Broedbiotoop zuiveren van bos >20are en bebouwing >10ha
voortplanting_schoon_opp <- terra::mask(voortplanting_opp, r_bos_20are, inverse = TRUE)

if (!all(is.na(suppressWarnings(terra::minmax(masker_bebouw_100m_groot))))) {
  voortplanting_schoon_opp <- terra::mask(voortplanting_schoon_opp, masker_bebouw_100m_groot, inverse = TRUE)
}

r_binair_broed_opp <- terra::ifel(!is.na(voortplanting_schoon_opp) & voortplanting_schoon_opp > 0, 1, NA)

res_broed_opp <- cluster_filter_compleet(
  masker     = r_binair_broed_opp,
  opp_laag   = voortplanting_schoon_opp,
  drempel_m2 = 300, # 300 m2 (0.03 ha)
  dist_m     = 50,    
  werkelijk  = TRUE
)
broed_gefilterd_opp <- res_broed_opp$raster

# Foerageerbiotoop clusteren op min. 1 ha
r_binair_foerageer_opp <- terra::ifel(!is.na(foerageer_totaal_opp) & foerageer_totaal_opp > 0, 1, NA)

res_foerageer_opp <- cluster_filter_compleet(
  masker     = r_binair_foerageer_opp,
  opp_laag   = foerageer_totaal_opp,
  drempel_m2 = 10000, # 1 ha
  dist_m     = 50,      
  werkelijk  = TRUE
)
foerageer_gefilterd_opp <- res_foerageer_opp$raster

water_gefilterd_opp     <- water_totaal_opp

rm(voortplanting_opp, voortplanting_schoon_opp, r_binair_broed_opp, res_broed_opp, 
   foerageer_totaal_opp, r_binair_foerageer_opp, res_foerageer_opp, water_totaal_opp,
   r_bos_20are, masker_bebouw_100m_groot)
gc()

# 8. STAP 4: OMGEVINGSBUFFERS BEREKENEN (300M VOEDSEL, 500M WATER) --------------
message("-> Omgevingsbuffers berekenen: 300m naar voedsel, 500m naar drinkwater...")

r_bin_foer <- terra::ifel(!is.na(foerageer_gefilterd_opp) & foerageer_gefilterd_opp > 0, 1, NA)
r_bin_wat  <- terra::ifel(!is.na(water_gefilterd_opp) & water_gefilterd_opp > 0, 1, NA)

tmp_buf_f300         <- tempfile(fileext = ".tif")
r_foerageer_buf_opp  <- terra::buffer(r_bin_foer, width = foerageer_afstand_m, filename = tmp_buf_f300, overwrite = TRUE)
r_foerageer_mask_opp <- terra::ifel(!is.na(r_foerageer_buf_opp) & r_foerageer_buf_opp > 0, 1, NA)

tmp_buf_w500         <- tempfile(fileext = ".tif")
r_water_buf_opp      <- terra::buffer(r_bin_wat, width = drinkwater_afstand_m, filename = tmp_buf_w500, overwrite = TRUE)
r_water_mask_opp     <- terra::ifel(!is.na(r_water_buf_opp) & r_water_buf_opp > 0, 1, NA)

r_omgeving_ok_opp    <- terra::mask(r_foerageer_mask_opp, r_water_mask_opp)

unlink(c(tmp_buf_f300, tmp_buf_w500))
rm(r_bin_foer, r_bin_wat, r_foerageer_buf_opp, r_water_buf_opp, r_foerageer_mask_opp, r_water_mask_opp,
   foerageer_gefilterd_opp, water_gefilterd_opp)
gc()

# 9. STAP 5: FINAAL LEEFGEBIED SAMENSTELLEN EN ID-RASTER GENEREREN -------------
message("-> Finaal leefgebied berekenen op werkelijke oppervlaktes...")

final_opp <- terra::mask(broed_gefilterd_opp, r_omgeving_ok_opp)

if (!all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  tmp_patch_cl <- tempfile(fileext = ".tif")
  cl_opp       <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE, filename = tmp_patch_cl, overwrite = TRUE)
} else {
  cl_opp <- template_Vlaanderen
}

rm(broed_gefilterd_opp, r_omgeving_ok_opp)
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

message("🏁 LEEFGEBIEDENGRID ZOMERTORTEL SUCCESVOL GEËXPORTEERD VOOR HEEL VLAANDEREN!")
