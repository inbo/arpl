# ==============================================================================
# LEEFGEBIEDENGRID VOOR BLAUWE KIEKENDIEF - HEEL VLAANDEREN (WERKELIJKE OPPERVLAKTE)
# Ecologische analyse & ruimtelijke clustering op open jaaggebied (INBO)
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
soort <- "blauwekiekendief"

resultaat <- df %>%
  filter(tolower(trimws(Soort)) == soort) %>%
  select(Type, MinOpp_ha, AfstandBiotopen_m, Dispersiecap_m)

oppervlakte_ha <- 50
buffer_m       <- 50000

print(resultaat)
rm(df, resultaat)

# 4. GEBIED VOORBEREIDEN (MASTERGRID VLAANDEREN) -------------------------------
message("-> Mastergrid Vlaanderen inlezen...")
master_grid <- rast(here("data/input/Raster_Vlaanderen/Vlaanderen_MasterGrid_10m.tif"))[[1]]

template_Vlaanderen <- master_grid
values(template_Vlaanderen) <- NA

df_namen_sleutel  <- read_csv(here("data/input/Excel_files/BWK_Laag_Namen_2025.csv"), show_col_types = FALSE)
gouden_namenlijst <- tolower(trimws(df_namen_sleutel$Laagnaam))

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

foer_raw_opp <- if("foerageer" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["foerageer"]] else template_Vlaanderen

rm(resultaten_gegroepeerd, df_nieuw, lijst_oppervlaktes)
gc()

# 6. LANDBOUWGEBRUIKSPERCELEN (LGP) EN BEBOUWINGSLAGEN -------------------------
message("-> Landbouwgebruikspercelen (LGP) inlezen...")

r_lgp_vlaanderen <- rast(here("data/input/Raster_Vlaanderen/vlaanderen_lgp_gwscod_2025_10m.tif"))
if(!terra::compareGeom(r_lgp_vlaanderen, template_Vlaanderen, stopOnError = FALSE)) {
  r_lgp_vlaanderen <- terra::resample(r_lgp_vlaanderen, template_Vlaanderen, method = "near")
}

kiekendief_winter_gewassen <- c(
  100, 101, 102, 
  201, 202,      
  311, 321,      
  601, 602, 603, 
  800, 801, 810, 
  900, 910, 931, 951 
) 

r_lgp_bin <- terra::ifel(!is.na(r_lgp_vlaanderen) & r_lgp_vlaanderen %in% kiekendief_winter_gewassen, 1, NA)

# Samenvoegen BWK-foerageerhabitat met geschikte LGP-gewassen
r_foer_bwk_clean <- terra::ifel(is.na(foer_raw_opp), 0, foer_raw_opp)
r_lgp_clean      <- terra::ifel(is.na(r_lgp_bin), 0, 1)

foer_raw_opp <- pmin(r_foer_bwk_clean + r_lgp_clean, 1.0)
foer_raw_opp <- terra::ifel(foer_raw_opp > 0, foer_raw_opp, NA)

rm(r_lgp_vlaanderen, r_lgp_bin, r_foer_bwk_clean, r_lgp_clean)
gc()

# BEBOUWINGSFILTER OPPARSEN UIT PARQUET
message("-> Bebouwingsklassen ophalen uit BWK Parquet...")
bebouwings_lagen <- gouden_namenlijst[grepl("^u[a-z0-9]", gouden_namenlijst) | gouden_namenlijst == "u" | gouden_namenlijst == "gh"]

tabel_bebouw <- ds_vlaanderen %>%
  filter(tolower(CODE) %in% bebouwings_lagen) %>%
  collect() %>%
  as.data.table()

r_bebouwing_schoon <- template_Vlaanderen

if (nrow(tabel_bebouw) > 0) {
  tabel_bebouw[, CODE := tolower(trimws(CODE))]
  tabel_bebouw_unique <- unique(tabel_bebouw, by = c("cel_id", "CODE"))
  tabel_bebouw_som    <- tabel_bebouw_unique[, .(Present = ifelse(sum(BWK_FRAC, na.rm = TRUE) >= 0.3, 1, 0)), by = .(cel_id)]
  tabel_bebouw_clean  <- tabel_bebouw_som[Present == 1]
  
  if (nrow(tabel_bebouw_clean) > 0) {
    r_bebouwing_schoon[tabel_bebouw_clean$cel_id] <- 1
  }
}

rm(ds_vlaanderen, tabel_bebouw)
gc()

# 7. VERBOSSINGSFILTER (EXTERNE GROENKAART FRACTIES) --------------------------
message("-> Verbossingsfilter berekenen: hooggroencomplexen > 20 are uitsluiten...")

r_hooggroen_frac <- rast(here("data/input/Raster_Vlaanderen/Groenkaart_2021_Fracties_5Bands_10m.tif"))[[1]]
if(!terra::compareGeom(r_hooggroen_frac, template_Vlaanderen, stopOnError = FALSE)) {
  r_hooggroen_frac <- terra::resample(r_hooggroen_frac, template_Vlaanderen, method = "near")
}

drempel_fractie     <- 0.3 
r_hooggroen_bin_10m <- terra::ifel(!is.na(r_hooggroen_frac) & r_hooggroen_frac >= drempel_fractie, 1, NA)

w_3x3     <- matrix(1, nrow = 3, ncol = 3)
r_som_3x3 <- terra::focal(r_hooggroen_bin_10m, w = w_3x3, fun = "sum", na.rm = TRUE)
r_hooggroen_schoon_10m <- terra::ifel(r_som_3x3 >= 2 & r_hooggroen_bin_10m == 1, 1, NA)

r_hooggroen_bin_40m <- terra::aggregate(r_hooggroen_schoon_10m, fact = 4, fun = "max", na.rm = TRUE)

verbossing_clusters_40m <- cluster_filter_compleet(
  masker     = r_hooggroen_bin_40m,
  opp_laag   = r_hooggroen_bin_40m,
  drempel_m2 = 2000,                
  dist_m     = 0,                    
  werkelijk  = FALSE      
)

r_bos_10m_back <- terra::disagg(verbossing_clusters_40m$raster, fact = 4)
r_bos_20are    <- terra::ifel(!is.na(r_bos_10m_back) & r_bos_10m_back > 0, 1, NA)

rm(r_hooggroen_frac, r_hooggroen_bin_10m, r_som_3x3, r_hooggroen_schoon_10m, r_hooggroen_bin_40m, verbossing_clusters_40m, r_bos_10m_back)
gc()

# 8. MODEL OP 50 HA OPEN JAAGGEBIED -------------------------------------------
message("-> Filteren van open jaaggebied (>= 50 ha) zonder dicht bos & bebouwingsbuffer (75m)...")

masker_bebouw_75m <- terra::buffer(r_bebouwing_schoon, width = 75)

r_bos_masker    <- terra::ifel(!is.na(r_bos_20are) & r_bos_20are > 0, 1, NA)
r_bebouw_masker <- terra::ifel(!is.na(masker_bebouw_75m) & masker_bebouw_75m > 0, 1, NA)

foer_schoon_opp <- terra::ifel(!is.na(r_bos_masker) | !is.na(r_bebouw_masker), NA, foer_raw_opp)

res_foer_opp <- cluster_filter_compleet(
  masker     = foer_schoon_opp,
  opp_laag   = foer_schoon_opp,
  drempel_m2 = 500000, # 50 ha
  dist_m     = 400,    
  werkelijk  = TRUE
)

final_opp <- res_foer_opp$raster
cl_opp    <- res_foer_opp$clusters

if (is.null(cl_opp) || all(is.na(terra::values(cl_opp, mat = FALSE)))) {
  if (!is.null(final_opp) && !all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
    tmp_patch <- tempfile(fileext = ".tif")
    cl_opp    <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE, 
                                filename = tmp_patch, overwrite = TRUE)
  } else {
    cl_opp <- template_Vlaanderen
  }
}

rm(r_bebouwing_schoon, masker_bebouw_75m, r_bos_20are, r_bos_masker, r_bebouw_masker, foer_raw_opp, foer_schoon_opp, res_foer_opp)
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
message(paste("    [OK] Geëxporteerd:", basename(file_path_opp)))

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

message("🏁 LEEFGEBIEDENGRID BLAUWE KIEKENDIEF SUCCESVOL GEËXPORTEERD VOOR HEEL VLAANDEREN!")
