# ==============================================================================
# LEEFGEBIEDENGRID VOOR GRUTTO - HEEL VLAANDEREN (WERKELIJKE OPPERVLAKTE)
# Ecologische analyse weidevogelbiotopen, drainage & netwerkknooppunten (INBO)
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
soort <- "grutto"

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

cultuur_0a_opp <- if("cultuur0a" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["cultuur0a"]] else template_Vlaanderen
akker_opp      <- if("akker" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["akker"]] else template_Vlaanderen
hooiland_opp   <- if("hooiland" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["hooiland"]] else template_Vlaanderen
water_bwk_opp  <- if("kleinwater_bwk" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["kleinwater_bwk"]] else template_Vlaanderen
kern1_opp      <- if("kern1" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["kern1"]] else template_Vlaanderen
kern2_opp      <- if("kern2" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["kern2"]] else template_Vlaanderen
c1b_opp        <- if("cultuurhooi1b" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["cultuurhooi1b"]] else template_Vlaanderen

rm(resultaten_gegroepeerd, df_nieuw, lijst_oppervlaktes)
gc()

# --- 6. CULTUURGRASLAND & HOOILAND ---
grutto_cultuur_opp <- terra::mask(cultuur_0a_opp, akker_opp > 0, inverse = TRUE)
grutto_hooiland_opp <- hooiland_opp

# --- 7. WATER (GRB + BWK HYBRID) ---
message("-> GRB Waterlopen/klein water verwerken...")
r_grb_raw <- rast(here("data/input/Raster_Vlaanderen/vlaanderen_grb_water_10m.tif"))
if(!terra::compareGeom(r_grb_raw, template_Vlaanderen, stopOnError = FALSE)) {
  r_grb_sync <- terra::resample(r_grb_raw, template_Vlaanderen, method = "near")
} else {
  r_grb_sync <- r_grb_raw
}

grb_water_bin <- terra::ifel(r_grb_sync > 0, 1, NA)
rm(r_grb_raw, r_grb_sync)
gc()

temp_patchfile <- tempfile(fileext = ".tif")
water_patches  <- terra::patches(grb_water_bin, directions = 8, zeroAsNA = TRUE, filename = temp_patchfile, overwrite = TRUE)

df_water_cl <- terra::as.data.frame(water_patches, cells = TRUE)
if (nrow(df_water_cl) > 0) {
  colnames(df_water_cl) <- c("cell", "ID")
  stats_water_count <- data.table(ID = df_water_cl$ID)[, .(count = .N), by = ID]
  klein_water_ids   <- stats_water_count[count <= 50, ID]
  grutto_water_grb_opp <- terra::ifel(water_patches %in% klein_water_ids, 0.01, NA)
} else {
  grutto_water_grb_opp <- template_Vlaanderen
}

w_grb_clean <- terra::ifel(is.na(grutto_water_grb_opp), 0, grutto_water_grb_opp)
w_bwk_clean <- terra::ifel(is.na(water_bwk_opp), 0, water_bwk_opp)

grutto_water_opp_raw <- terra::clamp(w_grb_clean + w_bwk_clean, upper = 1.0)
grutto_water_opp     <- terra::ifel(grutto_water_opp_raw > 0, grutto_water_opp_raw, NA)

unlink(temp_patchfile)
rm(water_patches, grutto_water_grb_opp, grb_water_bin, w_grb_clean, w_bwk_clean, grutto_water_opp_raw)
gc()

# --- 8. GROENKAART STAP ---
message("-> Groenkaart maskers berekenen...")
hooggroen_raw <- rast(here("data/input/ASCI Files/Groenkaart_2021.tif"))
if(!terra::compareGeom(hooggroen_raw, template_Vlaanderen, stopOnError = FALSE)) {
  hooggroen_sync <- terra::resample(hooggroen_raw, template_Vlaanderen, method = "near")
} else {
  hooggroen_sync <- hooggroen_raw
}

grutto_cultuur_opp_f  <- terra::mask(grutto_cultuur_opp, hooggroen_sync == 1, inverse = TRUE)
grutto_hooiland_opp_f <- terra::mask(grutto_hooiland_opp, hooggroen_sync == 1, inverse = TRUE)
grutto_water_opp_f    <- terra::mask(grutto_water_opp, hooggroen_sync == 1, inverse = TRUE)

rm(hooggroen_raw)
gc()

# --- 9. KERNEN INRICHTEN & DRAINAGEKAART ---
message("-> Kernbiotopen opbouwen x natte bodems...")
k1_clean <- terra::ifel(is.na(kern1_opp), 0, kern1_opp)
k2_clean <- terra::ifel(is.na(kern2_opp), 0, kern2_opp)

grutto_basis_kern_opp_raw <- terra::clamp(k1_clean + k2_clean, upper = 1.0)
grutto_basis_kern_opp     <- terra::ifel(grutto_basis_kern_opp_raw > 0, grutto_basis_kern_opp_raw, NA)

r_drain_raw <- rast(here("data/input/Raster_Vlaanderen/vlaanderen_drainage_n2khab_10m.tif"))
if(!terra::compareGeom(r_drain_raw, template_Vlaanderen, stopOnError = FALSE)) {
  r_drain_sync <- terra::resample(r_drain_raw, template_Vlaanderen, method = "near")
} else {
  r_drain_sync <- r_drain_raw
}

drain_cats             <- terra::cats(r_drain_sync)[[1]]
geselecteerde_letters <- c("d", "e", "f", "h", "i", "g", "e-f", "h-i", "e-i")
grutto_drain_ids      <- drain_cats$value[drain_cats$Label %in% geselecteerde_letters]

masker_drainage  <- r_drain_sync %in% grutto_drain_ids
grutto_kern3_opp <- terra::mask(grutto_basis_kern_opp, masker_drainage, maskvalues = FALSE)

rm(r_drain_raw, r_drain_sync, masker_drainage, k1_clean, k2_clean, grutto_basis_kern_opp_raw)
gc()

# --- 10. BEBOUWINGSFILTER ---
message("-> Grote bebouwingskernen (>= 10 ha) en hun 100m buffer berekenen...")
bebouwings_lagen <- gouden_namenlijst[grepl("^u[a-z0-9]", gouden_namenlijst) | gouden_namenlijst == "u"]
bebouwings_lagen <- bebouwings_lagen[bebouwings_lagen != "ulm"]

tabel_bebouw <- ds_vlaanderen %>%
  filter(tolower(CODE) %in% bebouwings_lagen) %>%
  collect() %>%
  as.data.table()

r_bebouwing_raw <- template_Vlaanderen

if (nrow(tabel_bebouw) > 0) {
  tabel_bebouw[, CODE := tolower(trimws(CODE))]
  tabel_bebouw_unique <- unique(tabel_bebouw, by = c("cel_id", "CODE"))
  tabel_bebouw_som    <- tabel_bebouw_unique[, .(Present = ifelse(sum(BWK_FRAC, na.rm = TRUE) > 0, 1, 0)), by = .(cel_id)]
  tabel_bebouw_clean  <- tabel_bebouw_som[Present == 1]
  
  if (nrow(tabel_bebouw_clean) > 0) {
    r_bebouwing_raw[tabel_bebouw_clean$cel_id] <- 1
  }
}

res_bebouw_groot <- cluster_filter_compleet(
  masker     = r_bebouwing_raw,
  opp_laag   = r_bebouwing_raw,
  drempel_m2 = 100000, # 10 ha
  dist_m     = 100,    
  werkelijk  = FALSE
)

r_bebouw_groot <- res_bebouw_groot$raster

if (!is.null(r_bebouw_groot) && !terra::global(is.na(r_bebouw_groot), "sum")[[1]] == terra::ncell(r_bebouw_groot)) {
  tmp_buf_b                <- tempfile(fileext = ".tif")
  buf_bebouw               <- terra::buffer(r_bebouw_groot, width = 100, filename = tmp_buf_b, overwrite = TRUE)
  masker_bebouw_100m_groot <- terra::ifel(!is.na(buf_bebouw) & buf_bebouw > 0, 1, NA)
  
  grutto_kern_bebouw_opp   <- terra::mask(grutto_kern3_opp, masker_bebouw_100m_groot, inverse = TRUE)
  unlink(tmp_buf_b)
} else {
  grutto_kern_bebouw_opp   <- grutto_kern3_opp
  masker_bebouw_100m_groot <- template_Vlaanderen
}

rm(tabel_bebouw, r_bebouwing_raw, res_bebouw_groot, r_bebouw_groot)
gc()

# --- 11. HOOGGROEN FILTER ---
message("-> Hooggroen-filter (20 are cluster + 50m buffer) berekenen...")
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
buf_groen                <- terra::buffer(r_groen_20are, width = 50, filename = tmp_buf_g, overwrite = TRUE)
masker_groen_50m_cluster <- terra::ifel(!is.na(buf_groen) & buf_groen > 0, 1, NA)

grutto_kern_finaal_opp   <- terra::mask(grutto_kern_bebouw_opp, masker_groen_50m_cluster, inverse = TRUE)

unlink(c(temp_groen_file, tmp_buf_g))
rm(r_groen_bin, groen_som, r_groen_20are, hooggroen_sync, buf_groen)
gc()

# --- 12. CLUSTER-FILTER KERNBIOTOOP (MIN. 10 HA) ---
drempel_kern_m2     <- 10 * 10000 
r_binair_finaal_opp <- terra::ifel(!is.na(grutto_kern_finaal_opp) & grutto_kern_finaal_opp > 0, 1, NA)

res_kern_opp <- cluster_filter_compleet(
  masker     = r_binair_finaal_opp,
  opp_laag   = grutto_kern_finaal_opp,
  drempel_m2 = drempel_kern_m2,
  dist_m     = 50,
  werkelijk  = TRUE
)
grutto_kern_10ha_opp <- res_kern_opp$raster

rm(r_binair_finaal_opp, res_kern_opp, grutto_kern_finaal_opp)
gc()

# --- 13. SAMENVOEGEN BOUWSTENEN & CULTUURHOOI1B ---
c_clean   <- terra::ifel(is.na(grutto_cultuur_opp_f), 0, grutto_cultuur_opp_f)
h_clean   <- terra::ifel(is.na(grutto_hooiland_opp_f), 0, grutto_hooiland_opp_f)
w_clean   <- terra::ifel(is.na(grutto_water_opp_f), 0, grutto_water_opp_f)
c1b_clean <- terra::ifel(is.na(c1b_opp), 0, c1b_opp)

grutto_cultuurhooi1_opp_raw <- terra::clamp(c_clean + h_clean + w_clean + c1b_clean, upper = 1.0)
grutto_cultuurhooi1_opp     <- terra::ifel(grutto_cultuurhooi1_opp_raw > 0, grutto_cultuurhooi1_opp_raw, NA)

filter_akker <- terra::ifel(!is.na(akker_opp) & akker_opp > 0, 1, NA)

grutto_c1_m1_opp     <- terra::mask(grutto_cultuurhooi1_opp, filter_akker, inverse = TRUE)
grutto_c1_m2_opp     <- terra::mask(grutto_c1_m1_opp, masker_bebouw_100m_groot, inverse = TRUE)
grutto_scenB_finaal_opp <- terra::mask(grutto_c1_m2_opp, masker_groen_50m_cluster, inverse = TRUE)

rm(c_clean, h_clean, w_clean, c1b_clean, grutto_cultuurhooi1_opp_raw, grutto_cultuurhooi1_opp, 
   grutto_c1_m1_opp, grutto_c1_m2_opp, filter_akker, masker_bebouw_100m_groot, masker_groen_50m_cluster)
gc()

# --- 14. WEDERZIJDSE KOPPELING (KERN & BOUWSTEEN ON 200M / MIN. 10 HA KERN & MIN. 25 HA TOTAAL) ---
voer_koppeling_en_drempel_uit_strikt_opp <- function(kern_rast, bouw_rast, afstand = 200, min_kern_ha = 10, min_totaal_ha = 25) {
  if (is.null(kern_rast) || terra::global(is.na(kern_rast), "sum")[[1]] == terra::ncell(kern_rast)) return(NULL)
  
  k_clean <- terra::ifel(is.na(kern_rast), 0, kern_rast)
  b_clean <- terra::ifel(is.na(bouw_rast), 0, bouw_rast)
  
  r_comb_val_raw <- terra::clamp(k_clean + b_clean, upper = 1.0)
  r_comb_val     <- terra::ifel(r_comb_val_raw > 0, r_comb_val_raw, NA)
  r_totaal_bin   <- terra::ifel(!is.na(r_comb_val) & r_comb_val > 0, 1, NA)
  
  if (terra::global(is.na(r_totaal_bin), "sum")[[1]] == terra::ncell(r_totaal_bin)) return(NULL)
  
  stral_cellen <- ceiling((afstand / 2) / 10)
  f_matrix     <- matrix(1, nrow = (2 * stral_cellen + 1), ncol = (2 * stral_cellen + 1))
  
  tmp_focal  <- tempfile(fileext = ".tif")
  r_buf      <- terra::focal(r_totaal_bin, w = f_matrix, fun = "max", na.rm = TRUE, filename = tmp_focal, overwrite = TRUE)
  r_buf      <- terra::ifel(r_buf > 0, 1, NA)
  
  tmp_patch  <- tempfile(fileext = ".tif")
  cl_netwerk <- terra::patches(r_buf, directions = 8, zeroAsNA = TRUE, filename = tmp_patch, overwrite = TRUE)
  cl_netwerk_biotoop <- terra::mask(cl_netwerk, r_totaal_bin)
  
  stats_netwerk <- terra::zonal(r_comb_val, cl_netwerk_biotoop, fun = "sum", na.rm = TRUE)
  colnames(stats_netwerk) <- c("netwerk_id", "totaal_ha_sum")
  stats_netwerk$totaal_ha  <- stats_netwerk$totaal_ha_sum * 0.01
  
  cl_kernen_only <- terra::mask(cl_netwerk_biotoop, kern_rast)
  stats_kernen   <- terra::zonal(kern_rast, cl_kernen_only, fun = "sum", na.rm = TRUE)
  colnames(stats_kernen) <- c("netwerk_id", "kern_ha_sum")
  stats_kernen$kern_ha   <- stats_kernen$kern_ha_sum * 0.01
  
  dt_net <- as.data.table(stats_netwerk)
  dt_ker <- as.data.table(stats_kernen)
  
  valide_netwerken <- dt_net[dt_ker, on = "netwerk_id"][totaal_ha >= min_totaal_ha & kern_ha >= min_kern_ha, netwerk_id]
  
  unlink(c(tmp_focal, tmp_patch))
  
  if (length(valide_netwerken) == 0) return(NULL)
  
  m_valide <- cl_netwerk_biotoop %in% valide_netwerken
  m_valide_bin <- terra::ifel(m_valide == 1, 1, NA)
  
  return(list(
    kern = terra::mask(kern_rast, m_valide_bin),
    bouw = terra::mask(bouw_rast, m_valide_bin)
  ))
}

resB_strikt <- voer_koppeling_en_drempel_uit_strikt_opp(grutto_kern_10ha_opp, grutto_scenB_finaal_opp, 200, 10, 25)

if (!is.null(resB_strikt)) {
  k_fin <- terra::ifel(is.na(resB_strikt$kern), 0, resB_strikt$kern)
  b_fin <- terra::ifel(is.na(resB_strikt$bouw), 0, resB_strikt$bouw)
  
  final_opp_raw <- terra::clamp(k_fin + b_fin, upper = 1.0)
  final_opp     <- terra::ifel(final_opp_raw > 0, final_opp_raw, NA)
} else {
  final_opp <- template_Vlaanderen
}

if (!all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  tmp_patch_cl <- tempfile(fileext = ".tif")
  cl_opp       <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE, filename = tmp_patch_cl, overwrite = TRUE)
} else {
  cl_opp <- template_Vlaanderen
}

rm(ds_vlaanderen, grutto_kern_10ha_opp, grutto_scenB_finaal_opp, resB_strikt)
gc()

# 15. EXPORT NAAR DISK (ENKEL 00_ID_Rasters EN 01_Werkelijke_Oppervlaktes) -----
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

message("🏁 LEEFGEBIEDENGRID GRUTTO SUCCESVOL GEËXPORTEERD VOOR HEEL VLAANDEREN!")
