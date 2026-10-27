# ==============================================================================
# LEEFGEBIEDENGRID VOOR PORSELEINHOEN - HEEL VLAANDEREN (WERKELIJKE OPPERVLAKTE)
# Ecologische analyse natte gronden, drainage, hooggroen & meertrapsclustering (INBO)
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
soort <- "porseleinhoen"

df <- read_excel(here("data/input/Excel_files/Soorten_bwk_afstanden.xlsx"))
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

# 6. EXTERNE DATASET (GROENKAART EN TERRILS) ----------------------------------
message("-> Groenkaart inlezen en synchroniseren voor heel Vlaanderen...")

hooggroen_raw <- rast(here("data/input/ASCI Files/Groenkaart_2021.tif"))
if(!terra::compareGeom(hooggroen_raw, template_Vlaanderen, stopOnError = FALSE)) {
  hooggroen_aligned <- terra::resample(hooggroen_raw, template_Vlaanderen, method = "near")
} else {
  hooggroen_aligned <- hooggroen_raw
}

hooggroen_binair <- terra::ifel(hooggroen_aligned == 1, 1, NA)

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

stoorlaag <- terra::cover(hooggroen_binair, r_terril_opp)

rm(hooggroen_raw, hooggroen_aligned, hooggroen_binair, r_terril_opp, tabel_terril)
gc()

# 7. STAP 1: BIOTOPEN COMBINEREN EN UITSLUITINGSZONES TOEPASSEN ---------------
message("-> Stap 1: Grasland- en moeraslagen combineren en stoorlagen uitsluiten...")

g_clean <- terra::ifel(is.na(grasland_opp), 0, grasland_opp)
m_clean <- terra::ifel(is.na(moeras_opp), 0, moeras_opp)

biotoop_opp_som <- g_clean + m_clean
biotoop_opp_raw <- terra::clamp(biotoop_opp_som, upper = 1.0)
biotoop_opp     <- terra::ifel(biotoop_opp_raw > 0, biotoop_opp_raw, NA)

porseleinhoen_geschikt2_opp <- terra::mask(biotoop_opp, stoorlaag, inverse = TRUE)

rm(g_clean, m_clean, biotoop_opp_som, biotoop_opp_raw, biotoop_opp, stoorlaag)
gc()

# 8. STAP 2: NATTE GRONDEN FILTER VIA DRAINAGEKAART ---------------------------
message("-> Stap 2: Beperken tot natte gronden via Drainagekaart...")

r_drain_raw <- terra::rast(here("data/input/Raster_Vlaanderen/vlaanderen_drainage_n2khab_10m.tif"))
if(!terra::compareGeom(r_drain_raw, template_Vlaanderen, stopOnError = FALSE)) {
  r_drain_local <- terra::resample(r_drain_raw, template_Vlaanderen, method = "near")
} else {
  r_drain_local <- r_drain_raw
}

drain_cats            <- terra::cats(r_drain_local)[[1]]
geselecteerde_letters <- c("d", "e", "f", "h", "i", "g", "e-f", "h-i", "e-i")

porseleinhoen_drain_ids <- drain_cats$value[drain_cats$Label %in% geselecteerde_letters]

masker_drainage <- r_drain_local %in% porseleinhoen_drain_ids
mal_drainage    <- terra::ifel(masker_drainage == 1, 1, NA)

porseleinhoen_geschikt3_opp <- terra::mask(porseleinhoen_geschikt2_opp, mal_drainage)

rm(r_drain_raw, r_drain_local, masker_drainage, mal_drainage, drain_cats, porseleinhoen_geschikt2_opp)
gc()

# 9. STAP 3: EERSTE CLUSTERING EN FILTERING (3 HA) -----------------------------
message("-> Stap 3: Eerste clustering en filtering op minstens 3 hectare...")

drempel_3ha_m2          <- 3 * 10000
r_binair_geschikt3_opp  <- terra::ifel(!is.na(porseleinhoen_geschikt3_opp) & porseleinhoen_geschikt3_opp > 0, 1, NA)

cluster_3ha_opp <- cluster_filter_compleet(
  masker     = r_binair_geschikt3_opp,
  opp_laag   = porseleinhoen_geschikt3_opp,
  drempel_m2 = drempel_3ha_m2,
  dist_m     = 50,
  werkelijk  = TRUE
)
porseleinhoen_geschikt4_opp <- cluster_3ha_opp$raster

rm(cluster_3ha_opp, r_binair_geschikt3_opp, porseleinhoen_geschikt3_opp)
gc()

# 10. STAP 4: OPEN ZICHT FILTER (100M AFSTAND ROND BOS >= 20 ARE) --------------
message("-> Stap 4: Open zicht filter & re-evaluatie (20 are bosdrempel)...")

groenkaart_raw <- terra::rast(here("data/input/ASCI Files/Groenkaart_2021.tif"))
if(!terra::compareGeom(groenkaart_raw, template_Vlaanderen, stopOnError = FALSE)) {
  hooggroen_aligned <- terra::resample(groenkaart_raw, template_Vlaanderen, method = "near")
} else {
  hooggroen_aligned <- groenkaart_raw
}

r_groen_bin <- terra::ifel(hooggroen_aligned == 1, 1, NA)

tmp_patch_groen <- tempfile(fileext = ".tif")
groen_patches   <- terra::patches(r_groen_bin, directions = 8, zeroAsNA = TRUE, filename = tmp_patch_groen, overwrite = TRUE)

df_groen_cl <- terra::as.data.frame(groen_patches, cells = TRUE)

if (nrow(df_groen_cl) > 0) {
  colnames(df_groen_cl) <- c("cell", "ID")
  stats_groen_count     <- data.table(ID = df_groen_cl$ID)[, .(count = .N), by = ID]
  
  # Bos-patches van minstens 20 are (20 cellen van 10x10m)
  grote_bos_ids         <- stats_groen_count[count >= 20, ID]
  
  if (length(grote_bos_ids) > 0) {
    bos_groot_bin <- terra::ifel(groen_patches %in% grote_bos_ids, 1, NA)
    
    tmp_buf_g                 <- tempfile(fileext = ".tif")
    bos_buffer_100m           <- terra::buffer(bos_groot_bin, width = 100, filename = tmp_buf_g, overwrite = TRUE)
    masker_groen_100m_cluster <- terra::ifel(!is.na(bos_buffer_100m) & bos_buffer_100m > 0, 1, NA)
    
    biotoop_fase5_kaal_frac   <- terra::mask(porseleinhoen_geschikt4_opp, masker_groen_100m_cluster, inverse = TRUE)
    
    unlink(tmp_buf_g)
    
    if (!all(is.na(suppressWarnings(terra::minmax(biotoop_fase5_kaal_frac))))) {
      r_binair_kaal_frac          <- terra::ifel(!is.na(biotoop_fase5_kaal_frac) & biotoop_fase5_kaal_frac > 0, 1, NA)
      porseleinhoen_geschikt5_opp <- cluster_filter_compleet(r_binair_kaal_frac, opp_laag = biotoop_fase5_kaal_frac, drempel_m2 = 30000, dist_m = 50, werkelijk = TRUE)$raster
      rm(r_binair_kaal_frac, biotoop_fase5_kaal_frac)
    } else {
      porseleinhoen_geschikt5_opp <- template_Vlaanderen
    }
  } else {
    porseleinhoen_geschikt5_opp   <- porseleinhoen_geschikt4_opp
  }
} else {
  porseleinhoen_geschikt5_opp     <- porseleinhoen_geschikt4_opp
}

unlink(tmp_patch_groen)
rm(r_groen_bin, groen_patches, df_groen_cl, groenkaart_raw, hooggroen_aligned, porseleinhoen_geschikt4_opp)
gc()

# 11. STAP 5: FINALE CLUSTERING & EINDDREMPEL (10 HA) --------------------------
message("-> Stap 5: Finale filtering van het leefgebiedcomplex op 10 hectare...")

drempel_10ha_m2        <- 10 * 10000
r_binair_geschikt5_opp <- terra::ifel(!is.na(porseleinhoen_geschikt5_opp) & porseleinhoen_geschikt5_opp > 0, 1, NA)

finaal_10ha_opp <- cluster_filter_compleet(
  masker     = r_binair_geschikt5_opp,
  opp_laag   = porseleinhoen_geschikt5_opp,
  drempel_m2 = drempel_10ha_m2,
  dist_m     = 50,
  werkelijk  = TRUE
)

final_opp <- finaal_10ha_opp$raster
cl_opp    <- finaal_10ha_opp$clusters

if (is.null(cl_opp) || all(is.na(terra::values(cl_opp, mat = FALSE)))) {
  if (!is.null(final_opp) && !all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
    tmp_patch <- tempfile(fileext = ".tif")
    cl_opp    <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE, filename = tmp_patch, overwrite = TRUE)
  } else {
    cl_opp <- template_Vlaanderen
  }
}

rm(finaal_10ha_opp, r_binair_geschikt5_opp, porseleinhoen_geschikt5_opp, ds_vlaanderen, grasland_opp, moeras_opp)
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

message("🏁 LEEFGEBIEDENGRID PORSELEINHOEN SUCCESVOL GEËXPORTEERD VOOR HEEL VLAANDEREN!")
