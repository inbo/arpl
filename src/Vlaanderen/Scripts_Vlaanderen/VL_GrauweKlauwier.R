# ==============================================================================
# LEEFGEBIEDENGRID VOOR GRAUWE KLAUWIER - HEEL VLAANDEREN (WERKELIJKE OPPERVLAKTE)
# Ecologische analyse stoorfactoren, KLE-structuur & leefgebiedsnetwerk (INBO)
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
soort <- "grauweklauwier"

resultaat <- df %>%
  filter(tolower(trimws(Soort)) == soort) %>%
  select(Type, MinOpp_ha, AfstandBiotopen_m, Dispersiecap_m)

oppervlakte_ha <- resultaat$MinOpp_ha[1]
afstand_m      <- resultaat$AfstandBiotopen_m[1]
buffer_m       <- 10000 

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

voortplanting_bwk_opp <- if("voortplanting_bwk" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["voortplanting_bwk"]] else template_Vlaanderen
foerageerbiotoop_opp <- if("foerageerbiotoop" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["foerageerbiotoop"]] else template_Vlaanderen

rm(resultaten_gegroepeerd, df_nieuw, lijst_oppervlaktes)
gc()

# 6. BEBOUWINGS- EN WEGENFILTER (STOORFACTOREN) -------------------------------
message("-> GDX-specifieke bebouwingslagen isoleren uit Parquet...")

gdx_bebouw_codes <- c("u", "ua", "ua-", "uc", "ud", "ui", "ui-", "un", "un-", "ur", "ur-", "ur+", "uv", "uv-", "uv+")

tabel_bebouw <- ds_vlaanderen %>%
  filter(tolower(CODE) %in% gdx_bebouw_codes) %>%
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

message("-> Grote bebouwingskernen (>= 10 ha) identificeren...")
p_bebouw <- terra::as.polygons(r_bebouwing_raw, aggregate = TRUE)

if (!is.null(p_bebouw) && nrow(p_bebouw) > 0) {
  p_bebouw_buf <- terra::buffer(p_bebouw, width = 25)
  p_netwerk    <- terra::aggregate(p_bebouw_buf) %>% terra::disagg()
  
  p_netwerk$area_ha <- as.numeric(st_area(st_as_sf(p_netwerk))) / 10000
  p_groot          <- p_netwerk[p_netwerk$area_ha >= 10, ]
  
  if (nrow(p_groot) > 0) {
    tf_bebouw            <- tempfile(fileext = ".tif")
    r_bebouw_groot       <- terra::rasterize(p_groot, template_Vlaanderen, field = 1, filename = tf_bebouw, overwrite = TRUE)
    dist_bebouw          <- terra::distance(r_bebouw_groot)
    masker_bebouw_stoor  <- dist_bebouw <= 100
    unlink(tf_bebouw)
  } else {
    masker_bebouw_stoor  <- template_Vlaanderen * 0
  }
} else {
  masker_bebouw_stoor    <- template_Vlaanderen * 0
}

message("-> Wegenstoring (100m) berekenen...")
if (file.exists(here("data/input/Raster_Vlaanderen/GRB_Wegverbinding_at.tif"))) {
  r_wegen_raw <- terra::rast(here("data/input/Raster_Vlaanderen/GRB_Wegverbinding_at.tif"))
  if(!terra::compareGeom(r_wegen_raw, template_Vlaanderen, stopOnError = FALSE)) {
    r_wegen <- terra::resample(r_wegen_raw, template_Vlaanderen, method = "near")
  } else {
    r_wegen <- r_wegen_raw
  }
  
  weg_cats            <- terra::cats(r_wegen)[[1]]
  gewenste_weg_labels <- c("autosnelweg", "weg met gescheiden rijbanen die geen autosnelweg is")
  ids_wegen_ok        <- weg_cats$value[tolower(trimws(weg_cats$Label)) %in% gewenste_weg_labels]
  
  if(length(ids_wegen_ok) == 0) {
    ids_wegen_ok <- weg_cats$value[grepl("autosnelweg|gescheiden rijbanen", tolower(weg_cats$Label))]
  }
  
  r_grote_wegen    <- r_wegen %in% ids_wegen_ok
  dist_wegen       <- terra::distance(terra::ifel(r_grote_wegen == 1, 1, NA))
  masker_weg_stoor <- dist_wegen <= 100
} else {
  masker_weg_stoor   <- template_Vlaanderen * 0
}

message("-> Stoorfactoren toepassen op de bruto voortplantingsbiotoop...")
voortplanting_gezuiverd_opp <- terra::mask(voortplanting_bwk_opp, terra::ifel(masker_bebouw_stoor == 1 | masker_weg_stoor == 1, NA, 1))

rm(tabel_bebouw, r_bebouwing_raw, p_bebouw, voortplanting_bwk_opp, masker_bebouw_stoor, masker_weg_stoor)
gc()

# 7. CLUSTERING EN FILTERING KERNBIOTOOP (50M FUZZY / MINIMAAL 10 HA) ----------
message("-> Kernbiotoop clusteren op 50m en zeven op minimaal 10 ha...")

r_binair_gezuiverd_opp <- terra::ifel(!is.na(voortplanting_gezuiverd_opp) & voortplanting_gezuiverd_opp > 0, 1, NA)

cl_kern_opp <- cluster_filter_compleet(
  masker     = r_binair_gezuiverd_opp,
  opp_laag   = voortplanting_gezuiverd_opp,
  drempel_m2 = 10 * 10000,
  dist_m     = 50,
  werkelijk  = TRUE
)
r_kern_10ha_opp <- cl_kern_opp$raster

rm(voortplanting_gezuiverd_opp, r_binair_gezuiverd_opp, cl_kern_opp)
gc()

# 8. EVALUATIE VAN LINEAIRE KLE'S (MINSTENS 500 METER) -------------------------
message("-> Groenkaart inlezen en Kleine Landschapselementen (KLE) valideren...")

r_groen_raw <- terra::rast(here("data/input/ASCI Files/Groenkaart_2021.tif"))
if(!terra::compareGeom(r_groen_raw, template_Vlaanderen, stopOnError = FALSE)) {
  r_groen <- terra::resample(r_groen_raw, template_Vlaanderen, method = "near")
} else {
  r_groen <- r_groen_raw
}

r_kle_groenkaart <- (r_groen == 1) | (r_groen == 2)
r_kle_binair     <- terra::ifel(r_kle_groenkaart, 1, 0)

r_binair_kern_opp <- terra::ifel(!is.na(r_kern_10ha_opp) & r_kern_10ha_opp > 0, 1, NA)

if (!all(is.na(suppressWarnings(terra::minmax(r_binair_kern_opp))))) {
  tmp_patch_k        <- tempfile(fileext = ".tif")
  cl_patches_opp     <- terra::patches(r_binair_kern_opp, directions = 8, zeroAsNA = TRUE, filename = tmp_patch_k, overwrite = TRUE)
  kle_counts_opp     <- terra::zonal(r_kle_binair, cl_patches_opp, fun = "sum", na.rm = TRUE)
  colnames(kle_counts_opp) <- c("ID", "KLE_pixels")
  
  goedgekeurde_patch_ids_opp <- kle_counts_opp$ID[kle_counts_opp$KLE_pixels >= 50] # 50 pixels = 500m
  
  if (length(goedgekeurde_patch_ids_opp) > 0) {
    masker_kle_ok_opp          <- cl_patches_opp %in% goedgekeurde_patch_ids_opp
    r_voortplanting_finaal_opp <- terra::mask(r_kern_10ha_opp, terra::ifel(masker_kle_ok_opp, 1, NA))
  } else {
    r_voortplanting_finaal_opp <- template_Vlaanderen
  }
  unlink(tmp_patch_k)
} else {
  r_voortplanting_finaal_opp   <- template_Vlaanderen
}

rm(r_groen_raw, r_groen, r_kle_groenkaart, r_kle_binair, r_binair_kern_opp, r_kern_10ha_opp)
gc()

# 9. FOERAGEERHABITAT KOPPELEN EN GRASLANDTEST (<90% HP-GRASLAND) --------------
message("-> Foerageerhabitat koppelen (100m buffer) en graslandfilter toepassen...")

grauweklauwier_leefgebied_opp <- template_Vlaanderen

if (!all(is.na(suppressWarnings(terra::minmax(r_voortplanting_finaal_opp))))) {
  r_binair_vpt_opp     <- terra::ifel(!is.na(r_voortplanting_finaal_opp) & r_voortplanting_finaal_opp > 0, 1, NA)
  tmp_buf_100          <- tempfile(fileext = ".tif")
  buffer_100m_kern_opp <- terra::buffer(r_binair_vpt_opp, width = 100, filename = tmp_buf_100, overwrite = TRUE)
  
  r_foerageer_geclipt_opp <- terra::mask(foerageerbiotoop_opp, buffer_100m_kern_opp)
  unlink(tmp_buf_100)
  
  vpt_opp_nul  <- terra::ifel(is.na(r_voortplanting_finaal_opp), 0, r_voortplanting_finaal_opp)
  foer_opp_nul <- terra::ifel(is.na(r_foerageer_geclipt_opp), 0, r_foerageer_geclipt_opp)
  
  r_leefgebied_bruto_opp_raw <- terra::clamp(vpt_opp_nul + foer_opp_nul, upper = 1.0)
  r_leefgebied_bruto_opp     <- terra::ifel(r_leefgebied_bruto_opp_raw > 0, r_leefgebied_bruto_opp_raw, NA)
  r_binair_bruto_opp         <- terra::ifel(!is.na(r_leefgebied_bruto_opp) & r_leefgebied_bruto_opp > 0, 1, NA)
  
  cl_leefgebied_opp <- cluster_filter_compleet(
    masker     = r_binair_bruto_opp,
    opp_laag   = r_leefgebied_bruto_opp, 
    drempel_m2 = 15 * 10000, # 15 ha
    dist_m     = 50,
    werkelijk  = TRUE
  )
  r_leefgebied_netwerk_opp <- cl_leefgebied_opp$raster
  
  if (!all(is.na(suppressWarnings(terra::minmax(r_leefgebied_netwerk_opp))))) {
    tmp_patch_lg         <- tempfile(fileext = ".tif")
    cl_final_patches_opp <- terra::patches(r_leefgebied_netwerk_opp, directions = 8, zeroAsNA = TRUE, filename = tmp_patch_lg, overwrite = TRUE)
    
    opp_totaal_opp            <- terra::zonal(r_leefgebied_bruto_opp, cl_final_patches_opp, fun = "sum", na.rm = TRUE)
    colnames(opp_totaal_opp)  <- c("ID", "Totaal_ha_sum")
    opp_totaal_opp$Totaal_ha  <- opp_totaal_opp$Totaal_ha_sum * 0.01
    
    # HP-Graslanden ophalen uit Parquet
    tabel_hp <- ds_vlaanderen %>%
      filter(str_detect(tolower(CODE), "^hp")) %>%
      collect() %>%
      as.data.table()
    
    r_hp_fractie <- template_Vlaanderen
    if (nrow(tabel_hp) > 0) {
      tabel_hp[, CODE := tolower(trimws(CODE))]
      tabel_hp_unique <- unique(tabel_hp, by = c("cel_id", "CODE"))
      tabel_hp_som    <- tabel_hp_unique[, .(Opp_HP = pmin(sum(BWK_FRAC, na.rm = TRUE), 1.0)), by = .(cel_id)]
      r_hp_fractie[tabel_hp_som$cel_id] <- tabel_hp_som$Opp_HP
    }
    
    opp_hp_zonal_opp          <- terra::zonal(r_hp_fractie, cl_final_patches_opp, fun = "sum", na.rm = TRUE)
    colnames(opp_hp_zonal_opp) <- c("ID", "HP_ha_sum")
    opp_hp_zonal_opp$HP_ha    <- opp_hp_zonal_opp$HP_ha_sum * 0.01
    
    grasland_check_opp <- merge(opp_totaal_opp, opp_hp_zonal_opp, by = "ID", all.x = TRUE)
    grasland_check_opp$HP_ha[is.na(grasland_check_opp$HP_ha)] <- 0
    grasland_check_opp$GraslandAandeel <- grasland_check_opp$HP_ha / grasland_check_opp$Totaal_ha
    
    goedgekeurde_leefgebied_ids_opp <- grasland_check_opp$ID[grasland_check_opp$GraslandAandeel < 0.90]
    
    if (length(goedgekeurde_leefgebied_ids_opp) > 0) {
      masker_leefgebied_ok_opp      <- cl_final_patches_opp %in% goedgekeurde_leefgebied_ids_opp
      grauweklauwier_leefgebied_opp <- terra::mask(r_leefgebied_bruto_opp, terra::ifel(masker_leefgebied_ok_opp, 1, NA))
    }
    unlink(tmp_patch_lg)
  }
}

final_opp <- grauweklauwier_leefgebied_opp

if (!all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  tmp_patch_cl <- tempfile(fileext = ".tif")
  cl_opp       <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE, filename = tmp_patch_cl, overwrite = TRUE)
} else {
  cl_opp <- template_Vlaanderen
}

rm(ds_vlaanderen, foerageerbiotoop_opp, r_voortplanting_finaal_opp, grauweklauwier_leefgebied_opp)
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

message("🏁 LEEFGEBIEDENGRID GRAUWE KLAUWIER SUCCESVOL GEËXPORTEERD VOOR HEEL VLAANDEREN!")