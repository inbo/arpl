# ==============================================================================
# LEEFGEBIEDENGRID VOOR WATERSNIP - HEEL VLAANDEREN (WERKELIJKE OPPERVLAKTE)
# Ecologische analyse natte gronden, drainage, landschapskwaliteit & terrestrische kern (INBO)
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
soort <- "watersnip"

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

biotopen_opp <- if("biotopen" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["biotopen"]] else template_Vlaanderen

rm(resultaten_gegroepeerd, df_nieuw, lijst_oppervlaktes)
gc()

# 6. EXTERNE DATASETS INLEZEN EN SYNCHRONISEREN -------------------------------
message("-> Drainage, groen- en waterkaarten inlezen en synchroniseren...")

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

r_drainage <- load_and_sync(here("data/input/Raster_Vlaanderen/vlaanderen_drainage_n2khab_10m.tif"))
r_groen    <- load_and_sync(here("data/input/ASCI Files/Groenkaart_2021.tif"))
r_grb      <- load_and_sync(here("data/input/Raster_Vlaanderen/vlaanderen_grb_water_10m.tif"))

drain_cats            <- terra::cats(r_drainage)[[1]]
geselecteerde_letters <- c("d", "e", "f", "h", "i", "g", "e-f", "h-i", "e-i")
watersnip_drain_ids   <- drain_cats$value[drain_cats$Label %in% geselecteerde_letters]
watersnip_drainage_ok <- r_drainage %in% watersnip_drain_ids

# DRAINAGE-FILTER TOEPASSEN OP WERKELIJKE OPPERVLAKTE
watersnip_basis_opp <- terra::mask(biotopen_opp, terra::ifel(watersnip_drainage_ok, 1, NA))

rm(r_drainage, drain_cats, watersnip_drain_ids, watersnip_drainage_ok, ds_vlaanderen)
gc()

# 7. STAP 2: PRIMAIRE CLUSTERING (MIN. 10 HA, 100M BUFFER) ----------------------
message("-> Primaire clustering (100m fuzzy koppeling) en filteren op min. 10 ha...")

drempel_10ha_m2     <- 10 * 10000
r_binair_basis_opp  <- terra::ifel(!is.na(watersnip_basis_opp) & watersnip_basis_opp > 0, 1, NA)

watersnip_clusters_opp <- cluster_filter_compleet(
  masker     = r_binair_basis_opp,
  opp_laag   = watersnip_basis_opp,
  drempel_m2 = drempel_10ha_m2,
  dist_m     = 100, 
  werkelijk  = TRUE
)

r_cluster_ids_opp <- watersnip_clusters_opp$clusters

rm(r_binair_basis_opp)
gc()

# 8. STAP 3: LANDSCHAPSKWALITEIT CHECKS PER CLUSTER (<75% BOS & <90% WATER) -----
message("-> Landschapskwaliteit beoordelen per unieke cluster...")

ids_kwaliteit_ok_opp <- c()

if (!all(is.na(suppressWarnings(terra::minmax(r_cluster_ids_opp))))) {
  totale_opp_opp <- terra::zonal(watersnip_basis_opp, r_cluster_ids_opp, fun = "sum", na.rm = TRUE)
  colnames(totale_opp_opp) <- c("Cluster_ID", "Totaal_ha_sum")
  totale_opp_opp$Totaal_ha <- totale_opp_opp$Totaal_ha_sum * 0.01
  
  bos_in_zone_opp   <- terra::mask(terra::ifel(r_groen == 1, 1.0, NA), r_cluster_ids_opp)
  bos_opp           <- terra::zonal(bos_in_zone_opp, r_cluster_ids_opp, fun = "sum", na.rm = TRUE)
  colnames(bos_opp) <- c("Cluster_ID", "Bos_ha_sum")
  bos_opp$Bos_ha    <- bos_opp$Bos_ha_sum * 0.01
  
  water_in_zone_opp   <- terra::mask(terra::ifel(r_grb == 1, 1.0, NA), r_cluster_ids_opp)
  water_opp           <- terra::zonal(water_in_zone_opp, r_cluster_ids_opp, fun = "sum", na.rm = TRUE)
  colnames(water_opp) <- c("Cluster_ID", "Water_ha_sum")
  water_opp$Water_ha  <- water_opp$Water_ha_sum * 0.01
  
  dt_tot   <- as.data.table(totale_opp_opp)
  dt_bos   <- as.data.table(bos_opp)
  dt_water <- as.data.table(water_opp)
  
  stats_opp <- dt_tot[dt_bos, on = "Cluster_ID"][dt_water, on = "Cluster_ID"]
  stats_opp[is.na(Bos_ha), Bos_ha := 0]
  stats_opp[is.na(Water_ha), Water_ha := 0]
  
  stats_opp[, Pct_Bos := (Bos_ha / Totaal_ha) * 100]
  stats_opp[, Pct_Water := (Water_ha / Totaal_ha) * 100]
  
  ids_kwaliteit_ok_opp <- stats_opp[Pct_Bos < 75 & Pct_Water < 90, Cluster_ID]
  
  rm(totale_opp_opp, bos_in_zone_opp, bos_opp, water_in_zone_opp, water_opp, dt_tot, dt_bos, dt_water, stats_opp)
}
rm(r_groen)
gc()

# 9. STAP 4: TERRESTRISCHE KERN CHECK & FINALE SAMENVOEGING (>= 10 HA LAND) -----
message("-> Terrestrische kern check uitvoeren (min. 10 ha niet-openwater biotoop)...")

watersnip_leefgebied_opp <- template_Vlaanderen

if (length(ids_kwaliteit_ok_opp) > 0) {
  r_grb_bin            <- terra::ifel(!is.na(r_grb) & r_grb == 1, 1, NA)
  r_terr_basis_opp_raw <- terra::mask(watersnip_basis_opp, r_grb_bin, inverse = TRUE)
  r_terr_basis_opp     <- terra::ifel(!is.na(r_terr_basis_opp_raw) & r_terr_basis_opp_raw > 0, r_terr_basis_opp_raw, NA)
  
  r_binair_terr_opp    <- terra::ifel(!is.na(r_terr_basis_opp) & r_terr_basis_opp > 0, 1, NA)
  
  terr_patches_opp     <- cluster_filter_compleet(
    masker     = r_binair_terr_opp,
    opp_laag   = r_terr_basis_opp,
    drempel_m2 = 100000, # 10 ha
    dist_m     = 0,      
    werkelijk  = TRUE     
  )
  r_terr_patch_ids_opp <- terr_patches_opp$clusters
  
  if (!all(is.na(suppressWarnings(terra::minmax(r_terr_patch_ids_opp))))) {
    df_overlap_opp <- terra::as.data.frame(c(r_cluster_ids_opp, r_terr_patch_ids_opp), na.rm = TRUE)
    if (nrow(df_overlap_opp) > 0) {
      colnames(df_overlap_opp) <- c("Cluster_ID", "Terr_Patch_ID")
      dt_overlap_opp           <- as.data.table(df_overlap_opp)
      ids_finaal_opp           <- unique(dt_overlap_opp[Cluster_ID %in% ids_kwaliteit_ok_opp, Cluster_ID])
      
      if (length(ids_finaal_opp) > 0) {
        watersnip_leefgebied_opp <- terra::mask(watersnip_clusters_opp$raster, r_cluster_ids_opp %in% ids_finaal_opp)
      }
      rm(df_overlap_opp, dt_overlap_opp)
    }
  }
  rm(r_grb_bin, r_terr_basis_opp_raw, r_terr_basis_opp, r_binair_terr_opp, terr_patches_opp, r_terr_patch_ids_opp)
}

final_opp <- watersnip_leefgebied_opp

if (!all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  tmp_patch_cl <- tempfile(fileext = ".tif")
  cl_opp       <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE, filename = tmp_patch_cl, overwrite = TRUE)
} else {
  cl_opp <- template_Vlaanderen
}

rm(r_grb, watersnip_basis_opp, watersnip_clusters_opp, r_cluster_ids_opp, ids_kwaliteit_ok_opp, watersnip_leefgebied_opp)
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

message("🏁 LEEFGEBIEDENGRID WATERSNIP SUCCESVOL GEËXPORTEERD VOOR HEEL VLAANDEREN!")
