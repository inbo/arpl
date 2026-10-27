# ==============================================================================
# LEEFGEBIEDENGRID VOOR GLADDESLANG - HEEL VLAANDEREN (WERKELIJKE OPPERVLAKTE)
# Ecologische analyse & dubbele-drempel netwerkfiltering met bosranden (INBO)
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

cluster_filter_twee_drempels <- function(masker, opp_laag, min_opp_ha, min_totale_opp_ha, dist_m) {
  if (terra::global(is.na(masker), "sum")[[1]] == terra::ncell(masker)) {
    return(list(raster = masker * NA, clusters = masker * NA))
  }
  
  min_opp_m2    <- min_opp_ha * 10000
  min_totale_m2 <- min_totale_opp_ha * 10000
  
  # 1. Individuele clusters identificeren en filteren op werkelijke opp (MinOpp_ha)
  tmp_patch_indiv <- tempfile(fileext = ".tif")
  cl_individueel  <- terra::patches(masker, directions = 8, zeroAsNA = TRUE, 
                                    filename = tmp_patch_indiv, overwrite = TRUE)
  
  df_indiv <- terra::as.data.frame(cl_individueel, cells = TRUE)
  if (nrow(df_indiv) == 0) {
    unlink(tmp_patch_indiv)
    return(list(raster = masker * NA, clusters = masker * NA))
  }
  colnames(df_indiv) <- c("cell", "ID")
  
  opp_val_indiv <- terra::extract(opp_laag, df_indiv$cell)[[1]]
  dt_indiv      <- data.table(ID = df_indiv$ID, Waarde = opp_val_indiv)
  
  stats_indiv       <- dt_indiv[!is.na(ID) & !is.na(Waarde), .(Area_m2 = sum(Waarde, na.rm = TRUE) * 100), by = ID]
  voldoet_indiv_ids <- stats_indiv[Area_m2 >= min_opp_m2, ID]
  
  unlink(tmp_patch_indiv)
  rm(df_indiv, dt_indiv, stats_indiv)
  gc()
  
  if (length(voldoet_indiv_ids) == 0) return(list(raster = masker * NA, clusters = masker * NA))
  
  masker_gefilterd <- terra::ifel(cl_individueel %in% voldoet_indiv_ids, masker, NA)
  
  # 2. Netwerkvorming via Dispersiecapaciteit (dist_m) met schijf-buffering
  if (dist_m > 0) {
    r_binair <- terra::ifel(!is.na(masker_gefilterd) & masker_gefilterd > 0, 1, NA)
    
    stral_cellen <- ceiling((dist_m / 2) / 10)
    f_matrix     <- matrix(1, nrow = (2 * stral_cellen + 1), ncol = (2 * stral_cellen + 1))
    
    tmp_focal  <- tempfile(fileext = ".tif")
    r_buffered <- terra::focal(r_binair, w = f_matrix, fun = "max", na.rm = TRUE, 
                               filename = tmp_focal, overwrite = TRUE)
    r_buffered <- terra::ifel(r_buffered > 0, 1, NA)
    
    tmp_patch  <- tempfile(fileext = ".tif")
    cl_network <- terra::patches(r_buffered, directions = 4, zeroAsNA = TRUE, 
                                 filename = tmp_patch, overwrite = TRUE)
    
    cl_biotoop_only <- terra::mask(cl_network, masker_gefilterd)
    
    unlink(c(tmp_focal, tmp_patch))
    rm(r_buffered, cl_network, r_binair)
  } else {
    tmp_patch       <- tempfile(fileext = ".tif")
    cl_biotoop_only <- terra::patches(masker_gefilterd, directions = 8, zeroAsNA = TRUE, 
                                      filename = tmp_patch, overwrite = TRUE)
  }
  
  # 3. Totale netwerkoppervlakte op basis van werkelijke fracties
  df_net <- terra::as.data.frame(cl_biotoop_only, cells = TRUE)
  if (nrow(df_net) == 0) return(list(raster = masker * NA, clusters = masker * NA))
  colnames(df_net) <- c("cell", "ID")
  
  opp_val_net <- terra::extract(opp_laag, df_net$cell)[[1]]
  dt_net      <- data.table(ID = df_net$ID, Waarde = opp_val_net)
  
  stats_net       <- dt_net[!is.na(ID) & !is.na(Waarde), .(Area_m2 = sum(Waarde, na.rm = TRUE) * 100), by = ID]
  voldoet_net_ids <- stats_net[Area_m2 >= min_totale_m2, ID]
  
  rm(df_net, dt_net, stats_net)
  gc()
  
  if (length(voldoet_net_ids) == 0) return(list(raster = masker * NA, clusters = masker * NA))
  
  final_network_mask <- cl_biotoop_only %in% voldoet_net_ids
  r_finaal  <- terra::mask(opp_laag, terra::ifel(final_network_mask, 1, NA))
  cl_finaal <- terra::mask(cl_biotoop_only, r_finaal)
  
  return(list(raster = r_finaal, clusters = cl_finaal))
}

# 3. SOORT INFORMATIE OPHALEN --------------------------------------------------
soort <- "gladdeslang"

df <- read_excel(here::here("data/input/Excel_files/Soorten_bwk_afstanden.xlsx"))
resultaat <- df %>%
  filter(tolower(trimws(Soort)) == soort) %>%
  select(Type, MinOpp_ha, AfstandBiotopen_m, Dispersiecap_m)

min_oppervlakte_ha <- 20
oppervlakte_ha     <- min_oppervlakte_ha
min_totale_opp_ha  <- 760
afstand_m          <- 0   
buffer_m           <- 500 

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

heide_opp <- if("biotoop_bwk" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["biotoop_bwk"]] else template_Vlaanderen
lijn_opp  <- if("lijn_bwk" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["lijn_bwk"]] else template_Vlaanderen

rm(ds_vlaanderen, resultaten_gegroepeerd, df_nieuw, lijst_oppervlaktes)
gc()

gladdeslang_biotoop1_opp <- cover(heide_opp, lijn_opp)
rm(heide_opp, lijn_opp)

# 6. EXTERNE GROENKAART INLEZEN EN SYNCHRONISEREN ------------------------------
message("-> Externe Groenkaart (hooggroen) inlezen voor heel Vlaanderen...")

hooggroen_raw <- rast(here("data/input/ASCI Files/Groenkaart_2021.tif"))
if(!terra::compareGeom(hooggroen_raw, template_Vlaanderen, stopOnError = FALSE)) {
  hooggroen_raw <- terra::resample(hooggroen_raw, template_Vlaanderen, method = "near")
}

hooggroen_binair  <- terra::classify(hooggroen_raw, matrix(c(0.5, 1.5, 1), ncol = 3, byrow = TRUE), others = 0)
hooggroen_sync    <- terra::ifel(hooggroen_binair == 1, 1, NA)

rm(hooggroen_raw, hooggroen_binair)
gc()

# 7. HABITATOPBOUW (OPEN HABITAT + 50M BOSRAND) --------------------------------
message("-> Open habitat filteren en 50m bosranden toevoegen...")

rand_breedte_m <- 50

gladdeslang_open_opp <- terra::mask(gladdeslang_biotoop1_opp, hooggroen_sync, maskvalues = 1)

r_binair_open_opp <- terra::ifel(!is.na(gladdeslang_open_opp) & gladdeslang_open_opp > 0, 1, NA)

tmp_buf_rand      <- tempfile(fileext = ".tif")
zoekzone_rand_opp <- terra::buffer(r_binair_open_opp, width = rand_breedte_m, filename = tmp_buf_rand, overwrite = TRUE)

gladdeslang_bufferbos_opp <- terra::ifel(hooggroen_sync == 1 & zoekzone_rand_opp == 1, 1, NA)

bruto_biotoop3_opp <- terra::cover(gladdeslang_open_opp, gladdeslang_bufferbos_opp)

unlink(tmp_buf_rand)
rm(gladdeslang_open_opp, r_binair_open_opp, zoekzone_rand_opp, gladdeslang_bufferbos_opp, hooggroen_sync, gladdeslang_biotoop1_opp)
gc()

# 8. DUBBELE DREMPELFILTERING (INDIVIDUEEL >= 20 HA EN NETWERK >= 760 HA) ------
message("-> Dubbele drempelfilter uitvoeren op werkelijke oppervlaktes...")

r_binair_bruto_opp <- terra::ifel(!is.na(bruto_biotoop3_opp) & bruto_biotoop3_opp > 0, 1, NA)

res_opp <- cluster_filter_twee_drempels(
  masker            = r_binair_bruto_opp,
  opp_laag          = bruto_biotoop3_opp,
  min_opp_ha        = min_oppervlakte_ha,
  min_totale_opp_ha = min_totale_opp_ha,
  dist_m            = buffer_m
)

final_opp <- res_opp$raster
cl_opp    <- res_opp$clusters

if (is.null(cl_opp) || all(is.na(terra::values(cl_opp, mat = FALSE)))) {
  if (!is.null(final_opp) && !all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
    tmp_patch <- tempfile(fileext = ".tif")
    cl_opp    <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE, filename = tmp_patch, overwrite = TRUE)
  } else {
    cl_opp <- template_Vlaanderen
  }
}

rm(bruto_biotoop3_opp, r_binair_bruto_opp, res_opp)
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

message("🏁 LEEFGEBIEDENGRID GLADDESLANG SUCCESVOL GEËXPORTEERD VOOR HEEL VLAANDEREN!")
