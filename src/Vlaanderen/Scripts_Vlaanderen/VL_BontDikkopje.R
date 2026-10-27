# ==============================================================================
# LEEFGEBIEDENGRID VOOR BONTIKKOPJE - HEEL VLAANDEREN (WERKELIJKE OPPERVLAKTE)
# Ecologische analyse & bosrand-contactvoorwaarden (INBO)
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
soort <- "bontdikkopje"

resultaat <- df %>% filter(tolower(trimws(Soort)) == soort) %>% select(Type, MinOpp_ha, AfstandBiotopen_m, Dispersiecap_m)

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
bosranden_opp         <- if("bosranden" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["bosranden"]] else template_Vlaanderen
water_bwk_opp         <- if("water_bwk" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["water_bwk"]] else template_Vlaanderen
foerageer_bwk_opp     <- if("foerageer_bwk" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["foerageer_bwk"]] else template_Vlaanderen

rm(ds_vlaanderen, resultaten_gegroepeerd, df_nieuw, lijst_oppervlaktes)
gc()

# 6. DRAINAGEFILTER INLEZEN ----------------------------------------------------
message("-> Bodemdrainage inlezen...")

r_drain_vlaanderen <- rast(here("data/input/Raster_Vlaanderen/vlaanderen_drainage_n2khab_10m.tif"))
if(!terra::compareGeom(r_drain_vlaanderen, template_Vlaanderen, stopOnError = FALSE)) {
  r_drain_vlaanderen <- terra::resample(r_drain_vlaanderen, template_Vlaanderen, method = "near")
}

drain_cats             <- terra::cats(r_drain_vlaanderen)[[1]]
geselecteerde_letters <- c("d", "e", "f", "h", "i", "g", "e-f", "h-i", "e-i")

drain_labels_char      <- as.character(drain_cats$Label)
bontdikkopje_drain_ids <- drain_cats$value[drain_labels_char %in% geselecteerde_letters]

masker_drainage         <- r_drain_vlaanderen %in% bontdikkopje_drain_ids
voortplanting_vocht_opp <- terra::mask(voortplanting_bwk_opp, masker_drainage, maskvalues = FALSE)

rm(r_drain_vlaanderen, masker_drainage)
gc()

# 7. AFSTANDSFILTER TOT OPEN WATER (200M) -------------------------------------
message("-> Waterlopen (GRB-water) inlezen en afstandsfilter 200m toepassen...")

r_grb_water_vlaanderen <- rast(here("data/input/Raster_Vlaanderen/vlaanderen_grb_water_10m.tif"))
if(!terra::compareGeom(r_grb_water_vlaanderen, template_Vlaanderen, stopOnError = FALSE)) {
  r_grb_water_vlaanderen <- terra::resample(r_grb_water_vlaanderen, template_Vlaanderen, method = "near")
}

masker_grb_water <- terra::ifel(r_grb_water_vlaanderen == 1, 1, NA)

bontdikkopjeWater <- terra::cover(masker_grb_water, water_bwk_opp)

afstand_tot_water <- terra::distance(bontdikkopjeWater)
zone_binnen_200m  <- terra::ifel(afstand_tot_water <= 200, 1, NA)

voortplanting_vocht_opp <- terra::mask(voortplanting_vocht_opp, zone_binnen_200m)

rm(r_grb_water_vlaanderen, masker_grb_water, bontdikkopjeWater, afstand_tot_water, zone_binnen_200m)
gc()

# 8. CLUSTER-ANALYSE HELPER MET BOSRANDCONTACT (WERKELIJKE OPPERVLAKTE) ------
verwerk_habitat_met_bosrand_opp <- function(r_habitat_opp, r_bos_opp, dist_m = 20, drempel_m2, afstand_m = 0) {
  
  r_habitat_bin_opp <- terra::ifel(!is.na(r_habitat_opp) & r_habitat_opp > 0, 1, NA)
  
  if (!terra::hasValues(r_habitat_bin_opp) || terra::global(r_habitat_bin_opp, "notNA")$notNA == 0) {
    return(list(opp = template_Vlaanderen, clusters_opp = template_Vlaanderen))
  }
  
  temp_mask_opp <- r_habitat_bin_opp
  if (afstand_m > 0) temp_mask_opp <- terra::buffer(r_habitat_bin_opp, width = afstand_m / 2)
  
  tmp_patch_hab       <- tempfile(fileext = ".tif")
  cl_habitat_opp      <- terra::patches(temp_mask_opp, directions = 8, zeroAsNA = TRUE, filename = tmp_patch_hab, overwrite = TRUE)
  cl_habitat_only_opp <- terra::mask(cl_habitat_opp, r_habitat_bin_opp)
  
  df_cl <- terra::as.data.frame(cl_habitat_only_opp, cells = TRUE)
  if (nrow(df_cl) == 0) {
    unlink(tmp_patch_hab)
    return(list(opp = template_Vlaanderen, clusters_opp = template_Vlaanderen))
  }
  colnames(df_cl) <- c("cell", "ID")
  
  opp_val <- terra::extract(r_habitat_opp, df_cl$cell)[[1]]
  dt_calc <- data.table(ID = df_cl$ID, Waarde = opp_val)
  
  stats_hab_opp <- dt_calc[!is.na(ID) & !is.na(Waarde), .(Area_m2 = sum(Waarde, na.rm = TRUE) * 100), by = ID]
  voldoet_opp_ids_opp <- stats_hab_opp[Area_m2 >= drempel_m2, ID]
  
  unlink(tmp_patch_hab)
  rm(df_cl, dt_calc, stats_hab_opp)
  gc()
  
  voldoet_bos_ids_opp <- c()
  if (length(voldoet_opp_ids_opp) > 0) {
    cl_habitat_groot_id_opp <- terra::ifel(cl_habitat_only_opp %in% voldoet_opp_ids_opp, cl_habitat_only_opp, NA)
    
    r_bos_bin_opp <- terra::ifel(!is.na(r_bos_opp) & r_bos_opp > 0, 1, NA)
    if (terra::hasValues(r_bos_bin_opp) && terra::global(r_bos_bin_opp, "notNA")$notNA > 0) {
      dist_naar_hab_opp       <- terra::distance(cl_habitat_groot_id_opp)
      bos_contact_pixels_opp  <- terra::ifel(dist_naar_hab_opp <= 10 & !is.na(r_bos_bin_opp), 1, NA)
      
      if (terra::hasValues(bos_contact_pixels_opp) && terra::global(bos_contact_pixels_opp, "notNA")$notNA > 0) {
        tmp_patch_con         <- tempfile(fileext = ".tif")
        cl_contact_lijnen_opp <- terra::patches(bos_contact_pixels_opp, directions = 8, zeroAsNA = TRUE, filename = tmp_patch_con, overwrite = TRUE)
        
        f_contact_opp        <- terra::freq(cl_contact_lijnen_opp)
        stats_contact_opp    <- data.frame(ID = f_contact_opp$value, Lengte_m = f_contact_opp$count * 10)   
        raaklijn_50m_ids_opp <- stats_contact_opp$ID[!is.na(stats_contact_opp$ID) & stats_contact_opp$Lengte_m >= 50]
        
        unlink(tmp_patch_con)
        
        if (length(raaklijn_50m_ids_opp) > 0) {
          r_raaklijn_50m_opp     <- terra::ifel(cl_contact_lijnen_opp %in% raaklijn_50m_ids_opp, 1, NA)
          hab_raakt_50m_lijn_opp <- terra::mask(cl_habitat_groot_id_opp, terra::buffer(r_raaklijn_50m_opp, width = 5))
          voldoet_bos_ids_opp    <- unique(na.omit(terra::values(hab_raakt_50m_lijn_opp, mat = FALSE)))
        }
      }
    }
  }
  
  if (length(voldoet_bos_ids_opp) > 0) {
    r_hab_goedgekeurd_opp    <- terra::mask(r_habitat_opp, cl_habitat_groot_id_opp %in% voldoet_bos_ids_opp)
    
    dist_naar_hab_finaal_opp <- terra::distance(r_hab_goedgekeurd_opp)
    r_bos_bin_opp            <- terra::ifel(!is.na(r_bos_opp) & r_bos_opp > 0, 1, NA)
    bos_binnen_20m_opp       <- terra::ifel(dist_naar_hab_finaal_opp <= dist_m & !is.na(r_bos_bin_opp), 1, NA)
    
    combo_potentieel_opp     <- terra::cover(r_hab_goedgekeurd_opp, terra::mask(r_bos_opp, bos_binnen_20m_opp))
    
    tmp_patch_combo          <- tempfile(fileext = ".tif")
    cl_combo_opp             <- terra::patches(terra::ifel(!is.na(combo_potentieel_opp) & combo_potentieel_opp > 0, 1, NA), 
                                               directions = 8, zeroAsNA = TRUE, filename = tmp_patch_combo, overwrite = TRUE)
    
    cl_die_habitat_raken_opp <- terra::mask(cl_combo_opp, r_hab_goedgekeurd_opp)
    aangesloten_ids_opp      <- unique(na.omit(terra::values(cl_die_habitat_raken_opp, mat = FALSE)))
    
    combo_finaal_mask_opp    <- cl_combo_opp %in% aangesloten_ids_opp
    complex_raw_opp          <- terra::mask(combo_potentieel_opp, combo_finaal_mask_opp)
    cl_finaal_opp            <- terra::mask(cl_combo_opp, complex_raw_opp)
    
    unlink(tmp_patch_combo)
  } else {
    complex_raw_opp <- template_Vlaanderen
    cl_finaal_opp   <- template_Vlaanderen
  }
  
  return(list(opp = complex_raw_opp, clusters_opp = cl_finaal_opp))
}

# 9. CLUSTER-ANALYSE VOORTPLANTING & FOERAGEREN --------------------------------
drempel_m2 <- oppervlakte_ha * 10000

message("-> VOORTPLANTING: Bosrand-voorwaarde & clusteranalyse opbouwen...")
vpt_complex <- verwerk_habitat_met_bosrand_opp(
  r_habitat_opp = voortplanting_vocht_opp,
  r_bos_opp     = bosranden_opp,
  dist_m        = 20,
  drempel_m2    = drempel_m2,
  afstand_m     = afstand_m
)

voortplanting_finaal_opp <- vpt_complex$opp
cl_voortplanting_opp     <- vpt_complex$clusters_opp

rm(vpt_complex, voortplanting_vocht_opp)
gc()

message("-> FOERAGEREN: Bosrand-voorwaarde & clusteranalyse opbouwen...")
foe_complex <- verwerk_habitat_met_bosrand_opp(
  r_habitat_opp = foerageer_bwk_opp,
  r_bos_opp     = bosranden_opp,
  dist_m        = 20,
  drempel_m2    = drempel_m2,
  afstand_m     = afstand_m
)

r_binair_foe_opp  <- terra::ifel(!is.na(foe_complex$opp) & foe_complex$opp > 0, 1, NA)
res_opp_foerageer <- cluster_filter_compleet(r_binair_foe_opp, foe_complex$opp, drempel_m2, afstand_m, werkelijk = TRUE)

foerageer_finaal_opp <- res_opp_foerageer$raster
cl_foerageer_opp     <- res_opp_foerageer$clusters

rm(res_opp_foerageer, foe_complex, r_binair_foe_opp, foerageer_bwk_opp, bosranden_opp)
gc()

# 10. TOTAAL LEEFGEBIED SAMENSTELLEN -------------------------------------------
message("-> Totale leefgebiedslagen opbouwen via terra::cover...")

final_opp <- terra::cover(voortplanting_finaal_opp, foerageer_finaal_opp)

# Unieke ID-toewijzing
if (!is.null(cl_foerageer_opp) && !all(is.na(suppressWarnings(terra::minmax(cl_foerageer_opp))))) {
  max_id_voortplanting_opp <- 0
  if (!is.null(cl_voortplanting_opp) && !all(is.na(suppressWarnings(terra::minmax(cl_voortplanting_opp))))) {
    val_opp <- terra::global(cl_voortplanting_opp, "max", na.rm = TRUE)[[1]]
    if (!is.na(val_opp)) max_id_voortplanting_opp <- val_opp
  }
  
  cl_foerageer_uniek_opp <- cl_foerageer_opp + (max_id_voortplanting_opp + 1000)
  cl_opp <- if (!is.null(cl_voortplanting_opp)) terra::cover(cl_voortplanting_opp, cl_foerageer_uniek_opp) else cl_foerageer_uniek_opp
} else {
  cl_opp <- cl_voortplanting_opp
}

# Failsafe vulling voor ID rasters
if (is.null(cl_opp) || all(is.na(terra::values(cl_opp, mat = FALSE)))) {
  if (!is.null(final_opp) && !all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
    tmp_patch <- tempfile(fileext = ".tif")
    cl_opp    <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE, filename = tmp_patch, overwrite = TRUE)
  } else {
    cl_opp <- template_Vlaanderen
  }
}

rm(voortplanting_finaal_opp, foerageer_finaal_opp, cl_voortplanting_opp, cl_foerageer_opp)
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

message("🏁 LEEFGEBIEDENGRID BONTDIKKOPJE SUCCESVOL GEËXPORTEERD VOOR HEEL VLAANDEREN!")
