# ==============================================================================
# LEEFGEBIEDENGRID VOOR ZWARTE SPECHT - HEEL VLAANDEREN (WERKELIJKE OPPERVLAKTE)
# Ecologische analyse boscomplexen (50-300ha / >300ha) x 4km actieradius (INBO)
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
soort <- "zwartespecht"

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

bwk_opp <- if("bwk" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["bwk"]] else template_Vlaanderen

rm(ds_vlaanderen, resultaten_gegroepeerd, df_nieuw, lijst_oppervlaktes)
gc()

# 6. STAP 1: BASISBIOTOOP EN HOOGGROEN COUPLING (100M FUZZY) --------------------
message("-> Groenkaart inlezen en bosbiotoop koppelen via 100m fuzzy netwerk...")

r_groen_raw <- rast(here("data/input/ASCI Files/Groenkaart_2021.tif"))
if(!terra::compareGeom(r_groen_raw, template_Vlaanderen, stopOnError = FALSE)) {
  r_groen <- terra::resample(r_groen_raw, template_Vlaanderen, method = "near")
} else {
  r_groen <- r_groen_raw
}

zwartespecht_groen_ok <- r_groen == 1
r_groen_opp_clean    <- terra::ifel(!is.na(zwartespecht_groen_ok) & zwartespecht_groen_ok == 1, 1.0, 0)

b_clean                <- terra::ifel(is.na(bwk_opp), 0, bwk_opp)
bos_som_raw            <- b_clean + r_groen_opp_clean
zwartespecht_basis_opp <- terra::clamp(bos_som_raw, upper = 1.0)
zwartespecht_basis_opp <- terra::ifel(zwartespecht_basis_opp > 0, zwartespecht_basis_opp, NA)

r_binair_basis_opp        <- terra::ifel(!is.na(zwartespecht_basis_opp) & zwartespecht_basis_opp > 0, 1, NA)
zwartespecht_clusters_opp <- cluster_filter_compleet(
  masker     = r_binair_basis_opp,
  opp_laag   = zwartespecht_basis_opp,
  drempel_m2 = 0,   
  dist_m     = 100,        
  werkelijk  = TRUE        
)
r_bos_patches_opp <- zwartespecht_clusters_opp$clusters

rm(r_groen_raw, r_groen, zwartespecht_groen_ok, r_groen_opp_clean, b_clean, bos_som_raw, 
   r_binair_basis_opp, zwartespecht_clusters_opp, bwk_opp)
gc()

# 7. STAP 2: CLASSIFICATIE BOSGROOTTE (>300 HA VS. 50-300 HA) ------------------
message("-> Bosclusters classificeren in grote kernen (>= 300 ha) en kleine bosjes (50-300 ha)...")

r_bos_groot_300ha_opp <- template_Vlaanderen
r_bos_klein_50ha_opp  <- template_Vlaanderen

if (!all(is.na(suppressWarnings(terra::minmax(r_bos_patches_opp))))) {
  stats_ha_opp           <- terra::zonal(zwartespecht_basis_opp, r_bos_patches_opp, fun = "sum", na.rm = TRUE)
  colnames(stats_ha_opp) <- c("ID", "Grootte_ha_sum")
  stats_ha_opp$Grootte_ha<- stats_ha_opp$Grootte_ha_sum * 0.01
  
  dt_stats_ha      <- as.data.table(stats_ha_opp)
  ids_300ha_opp    <- dt_stats_ha[Grootte_ha >= 300, ID]
  ids_50_300ha_opp <- dt_stats_ha[Grootte_ha >= 50 & Grootte_ha < 300, ID]
  
  if (length(ids_300ha_opp) > 0)    r_bos_groot_300ha_opp <- terra::mask(zwartespecht_basis_opp, r_bos_patches_opp %in% ids_300ha_opp)
  if (length(ids_50_300ha_opp) > 0) r_bos_klein_50ha_opp  <- terra::mask(zwartespecht_basis_opp, r_bos_patches_opp %in% ids_50_300ha_opp)
  
  rm(stats_ha_opp, dt_stats_ha)
}
gc()

# 8. STAP 3: RECONSTRUCTIE NETWERK (4 KM BUFFER INTERACTIE) --------------------
message("-> Netwerkrelaties controleren binnen een actieradius van 4 km...")

r_bos_klein_goedgekeurd_opp <- template_Vlaanderen

if (!all(is.na(suppressWarnings(terra::minmax(r_bos_klein_50ha_opp))))) {
  r_binair_klein_opp   <- terra::ifel(!is.na(r_bos_klein_50ha_opp) & r_bos_klein_50ha_opp > 0, 1, NA)
  
  tmp_buf_k4km         <- tempfile(fileext = ".tif")
  buffer_4km_klein_opp <- terra::buffer(r_binair_klein_opp, width = 4000, filename = tmp_buf_k4km, overwrite = TRUE)
  
  if (!all(is.na(suppressWarnings(terra::minmax(r_bos_groot_300ha_opp))))) {
    r_binair_groot_opp   <- terra::ifel(!is.na(r_bos_groot_300ha_opp) & r_bos_groot_300ha_opp > 0, 1, NA)
    
    tmp_buf_g4km         <- tempfile(fileext = ".tif")
    buffer_4km_groot_opp <- terra::buffer(r_binair_groot_opp, width = 4000, filename = tmp_buf_g4km, overwrite = TRUE)
    
    netwerk_masker_opp   <- !is.na(buffer_4km_groot_opp) | !is.na(buffer_4km_klein_opp)
    
    unlink(tmp_buf_g4km)
    rm(r_binair_groot_opp, buffer_4km_groot_opp)
  } else {
    netwerk_masker_opp   <- !is.na(buffer_4km_klein_opp)
  }
  
  r_bos_klein_goedgekeurd_opp <- terra::mask(r_bos_klein_50ha_opp, netwerk_masker_opp)
  
  unlink(tmp_buf_k4km)
  rm(r_binair_klein_opp, buffer_4km_klein_opp, netwerk_masker_opp)
}

g_clean <- terra::ifel(is.na(r_bos_groot_300ha_opp), 0, r_bos_groot_300ha_opp)
k_clean <- terra::ifel(is.na(r_bos_klein_goedgekeurd_opp), 0, r_bos_klein_goedgekeurd_opp)

som_raw         <- g_clean + k_clean
som_cl          <- terra::clamp(som_raw, upper = 1.0)
leefgebied1_opp <- terra::ifel(som_cl > 0, som_cl, NA)

rm(r_bos_groot_300ha_opp, r_bos_klein_50ha_opp, r_bos_klein_goedgekeurd_opp, g_clean, k_clean, som_raw, som_cl)
gc()

# 9. STAP 4: FINALE CLUSTERING EN OPPERVLAKTEFILTER (MIN. 50 HA) ---------------
message("-> Finale 100m clustering en 50 ha minimumoppervlakte check...")

r_binair_leef1_opp             <- terra::ifel(!is.na(leefgebied1_opp) & leefgebied1_opp > 0, 1, NA)
finaal_leefgebied_clusters_opp <- cluster_filter_compleet(
  masker     = r_binair_leef1_opp,
  opp_laag   = leefgebied1_opp,
  drempel_m2 = 500000, # 50 ha
  dist_m     = 100,     
  werkelijk  = TRUE    
)

final_opp <- finaal_leefgebied_clusters_opp$raster
cl_opp    <- finaal_leefgebied_clusters_opp$clusters

if (is.null(cl_opp) || all(is.na(terra::values(cl_opp, mat = FALSE)))) {
  if (!is.null(final_opp) && !all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
    tmp_patch <- tempfile(fileext = ".tif")
    cl_opp    <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE, filename = tmp_patch, overwrite = TRUE)
  } else {
    cl_opp <- template_Vlaanderen
  }
}

rm(zwartespecht_basis_opp, r_bos_patches_opp, leefgebied1_opp, r_binair_leef1_opp, finaal_leefgebied_clusters_opp)
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

message("🏁 LEEFGEBIEDENGRID ZWARTE SPECHT SUCCESVOL GEËXPORTEERD VOOR HEEL VLAANDEREN!")
