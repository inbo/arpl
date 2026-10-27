# ==============================================================================
# LEEFGEBIEDENGRID VOOR BOOMPIEPER - HEEL VLAANDEREN (WERKELIJKE OPPERVLAKTE)
# Ecologische analyse & bosrand/open biotoop gradiëntfiltering (INBO)
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
  memfrac = 0.6,       # Behoudend geheugenbeheer
  memmin = 1,
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
  
  tmp_buf <- tempfile(pattern = "cl_buf_", fileext = ".tif")
  on.exit(unlink(tmp_buf), add = TRUE)
  
  if (dist_m > 0) {
    r_binair <- terra::ifel(!is.na(masker) & masker > 0, 1, NA)
    
    r_buffered <- terra::buffer(
      r_binair, 
      width = dist_m, 
      filename = tmp_buf, 
      overwrite = TRUE, 
      gdal = c("COMPRESS=LZW")
    )
    
    tmp_patch  <- tempfile(pattern = "cl_patch_", fileext = ".tif")
    cl_network <- terra::patches(r_buffered, directions = 4, zeroAsNA = TRUE, filename = tmp_patch, overwrite = TRUE)
    cl_biotoop_only <- terra::mask(cl_network, masker)
    unlink(tmp_patch)
  } else {
    tmp_patch       <- tempfile(pattern = "cl_patch_", fileext = ".tif")
    cl_biotoop_only <- terra::patches(masker, directions = 8, zeroAsNA = TRUE, filename = tmp_patch, overwrite = TRUE)
    unlink(tmp_patch)
  }
  
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
  
  cl_finaal <- terra::ifel(cl_biotoop_only %in% voldoet_ids, cl_biotoop_only, NA)
  r_finaal  <- terra::mask(masker, cl_finaal)
  
  rm(df_cl, dt_calc, stats_dt, cl_biotoop_only)
  gc()
  
  return(list(raster = r_finaal, clusters = cl_finaal))
}

# 3. SOORT INFORMATIE OPHALEN --------------------------------------------------
soort <- "boompieper"

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

bwk_opp        <- if("bwk" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["bwk"]] else template_Vlaanderen
bwk_bossen_opp <- if("bwk_bossen" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["bwk_bossen"]] else template_Vlaanderen
heide_opp      <- if("heide" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["heide"]] else template_Vlaanderen
naaldbos_opp   <- if("naaldbos" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["naaldbos"]] else template_Vlaanderen

rm(ds_vlaanderen, resultaten_gegroepeerd, df_nieuw, lijst_oppervlaktes)
gc()

# 6. ONAFHANKELIJK CLUSTEREN OP 5 HA (50M AFSTAND) -----------------------------
message("-> Eerst onafhankelijk clusteren van open biotopen en bos op 5 ha...")

r_binair_open_opp <- terra::ifel(!is.na(bwk_opp) & bwk_opp > 0, 1, NA)
open_cluster_5ha_opp <- cluster_filter_compleet(
  masker     = r_binair_open_opp,
  opp_laag   = bwk_opp, 
  drempel_m2 = 5 * 10000,
  dist_m     = 50, 
  werkelijk  = TRUE
)$raster
rm(r_binair_open_opp)
gc()

r_binair_bos_opp <- terra::ifel(!is.na(bwk_bossen_opp) & bwk_bossen_opp > 0, 1, NA)
bos_cluster_5ha_opp <- cluster_filter_compleet(
  masker     = r_binair_bos_opp,
  opp_laag   = bwk_bossen_opp, 
  drempel_m2 = 5 * 10000,
  dist_m     = 50, 
  werkelijk  = TRUE
)$raster
rm(r_binair_bos_opp)
gc()

tmp_open_opp_file <- tempfile(pattern = "temp_open_opp_", fileext = ".tif")
tmp_bos_opp_file  <- tempfile(pattern = "temp_bos_opp_",  fileext = ".tif")

writeRaster(open_cluster_5ha_opp, tmp_open_opp_file, overwrite = TRUE)
writeRaster(bos_cluster_5ha_opp,  tmp_bos_opp_file,  overwrite = TRUE)

rm(open_cluster_5ha_opp, bos_cluster_5ha_opp)
gc()

open_cluster_5ha_opp <- rast(tmp_open_opp_file)
bos_cluster_5ha_opp  <- rast(tmp_bos_opp_file)

# 7. 100M BOSRAND UITSNIJDEN VAN DE 5 HA BOSCLUSTERS ---------------------------
message("-> 100m bosranden uitsnijden van de grote bosclusters via afstandskaart...")

dist_tot_rand    <- terra::distance(!is.na(bos_cluster_5ha_opp))
bos_100m_5ha_opp <- terra::mask(bos_cluster_5ha_opp, dist_tot_rand <= 100, maskvalue = FALSE)

rm(dist_tot_rand)
gc()

# 8. CONTACTZONE (20M) & EDGEN TUSSEN OPEN BIOTOOP EN BOSRAND ------------------
message("-> Edgen: Alleen bosranden behouden die grenzen aan open biotoopclusters...")

w_edge_20m <- matrix(1, nrow = 3, ncol = 3)

f_buf_open_tif <- tempfile(pattern = "buf_open_", fileext = ".tif")
buf_open_opp <- terra::focal(
  terra::ifel(!is.na(open_cluster_5ha_opp) & open_cluster_5ha_opp > 0, 1, NA), 
  w = w_edge_20m, fun = "max", na.rm = TRUE,
  filename = f_buf_open_tif, overwrite = TRUE, gdal = c("COMPRESS=LZW")
)

f_buf_bos_tif <- tempfile(pattern = "buf_bos_", fileext = ".tif")
buf_bos_opp <- terra::focal(
  terra::ifel(!is.na(bos_100m_5ha_opp) & bos_100m_5ha_opp > 0, 1, NA), 
  w = w_edge_20m, fun = "max", na.rm = TRUE,
  filename = f_buf_bos_tif, overwrite = TRUE, gdal = c("COMPRESS=LZW")
)

contactzone_groot_opp <- terra::ifel(!is.na(buf_open_opp) & !is.na(buf_bos_opp), 1, NA)

bosrand_gekoppeld_opp <- terra::mask(bos_100m_5ha_opp, contactzone_groot_opp)

bwk_cluster_opp <- terra::cover(open_cluster_5ha_opp, bosrand_gekoppeld_opp)

# 9. RANDZONES BEPALEN VOOR KLEIN LEEFGEBIED (HEIDE + NAALDBOS) ----------------
message("-> Randzones (edges) bepalen tussen heide en naaldbos voor klein leefgebied...")

r_naaldbos_bin <- terra::ifel(!is.na(naaldbos_opp) & naaldbos_opp > 0, 1, NA)
r_heide_bin    <- terra::ifel(!is.na(heide_opp) & heide_opp > 0, 1, NA)

if(!all(is.na(suppressWarnings(terra::minmax(r_naaldbos_bin)))) && !all(is.na(suppressWarnings(terra::minmax(r_heide_bin))))) {
  f_buf_naald_tif <- tempfile(pattern = "buf_naald_", fileext = ".tif")
  buf_naaldbos <- terra::focal(
    r_naaldbos_bin, 
    w = w_edge_20m, fun = "max", na.rm = TRUE,
    filename = f_buf_naald_tif, overwrite = TRUE, gdal = c("COMPRESS=LZW")
  )
  
  f_buf_heide_tif <- tempfile(pattern = "buf_heide_", fileext = ".tif")
  buf_heide <- terra::focal(
    r_heide_bin, 
    w = w_edge_20m, fun = "max", na.rm = TRUE,
    filename = f_buf_heide_tif, overwrite = TRUE, gdal = c("COMPRESS=LZW")
  )
  
  contactzone_klein  <- terra::ifel(!is.na(buf_naaldbos) & !is.na(buf_heide), 1, NA)
  heide_naaldbos_opp <- terra::mask(heide_opp, contactzone_klein)
} else {
  heide_naaldbos_opp <- template_Vlaanderen
}

gc()

# 10. CLUSTEREN EN FILTEREN OP 2 HECTARE (KLEIN LEEFGEBIED) ---------------------
message("-> Clusteren en filteren van heide-bosranden op minimale oppervlakte van 2 ha...")

r_bin_klein <- terra::ifel(!is.na(heide_naaldbos_opp) & heide_naaldbos_opp > 0, 1, NA)
bwk_cluster_klein_opp_res <- cluster_filter_compleet(r_bin_klein, heide_naaldbos_opp, dist_m = 50, drempel_m2 = 2 * 10000, werkelijk = TRUE)
bwk_cluster_klein_opp     <- bwk_cluster_klein_opp_res$raster

tijdelijke_bestanden <- c(
  tmp_open_opp_file, tmp_bos_opp_file,
  if(exists("f_buf_open_tif")) f_buf_open_tif,
  if(exists("f_buf_bos_tif")) f_buf_bos_tif,
  if(exists("f_buf_naald_tif")) f_buf_naald_tif,
  if(exists("f_buf_heide_tif")) f_buf_heide_tif
)
unlink(tijdelijke_bestanden)

rm(heide_naaldbos_opp, open_cluster_5ha_opp, bos_cluster_5ha_opp, bos_100m_5ha_opp, bosrand_gekoppeld_opp, buf_open_opp, buf_bos_opp)
gc()

# 11. COMBINEREN (GROOT OR KLEIN) -----------------------------------------------
message("-> Leefgebieden combineren (Groot OF Klein)...")

r_groot_opp_status <- terra::ifel(!is.na(bwk_cluster_opp) & bwk_cluster_opp > 0, bwk_cluster_opp, NA)
r_klein_opp_status <- terra::ifel(!is.na(bwk_cluster_klein_opp) & bwk_cluster_klein_opp > 0, bwk_cluster_klein_opp, NA)

final_opp <- terra::cover(r_groot_opp_status, r_klein_opp_status)

if (!all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  tmp_patch <- tempfile(fileext = ".tif")
  cl_opp    <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE, filename = tmp_patch, overwrite = TRUE)
} else {
  cl_opp <- template_Vlaanderen
}

rm(bwk_cluster_opp, bwk_cluster_klein_opp, r_groot_opp_status, r_klein_opp_status)
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

message("🏁 LEEFGEBIEDENGRID BOOMPIEPER SUCCESVOL GEËXPORTEERD VOOR HEEL VLAANDEREN!")
