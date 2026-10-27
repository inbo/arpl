# ==============================================================================
# LEEFGEBIEDENGRID VOOR BEEKROMBOUT - HEEL VLAANDEREN (WERKELIJKE OPPERVLAKTE)
# Ecologische analyse & afstandsrelaties water-land (INBO)
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
  
  # 3. Zonal statistics: Hoge snelheid via C++ terra engine
  if (werkelijk) {
    zonal_res <- terra::zonal(opp_laag, cl_biotoop_only, fun = "sum", na.rm = TRUE)
  } else {
    r_ones <- terra::ifel(!is.na(cl_biotoop_only), 1, NA)
    zonal_res <- terra::zonal(r_ones, cl_biotoop_only, fun = "sum", na.rm = TRUE)
  }
  
  if (nrow(zonal_res) == 0) return(list(raster = masker * NA, clusters = masker * NA))
  
  stats_dt <- as.data.table(zonal_res)
  setnames(stats_dt, c("ID", "Waarde"))
  stats_dt[, Area_m2 := Waarde * 100] # 10m x 10m cel = 100 m²
  
  voldoet_ids <- stats_dt[!is.na(ID) & Area_m2 >= drempel_m2, ID]
  
  if (length(voldoet_ids) == 0) return(list(raster = masker * NA, clusters = masker * NA))
  
  # 4. Resultaat filteren
  cl_finaal <- terra::ifel(cl_biotoop_only %in% voldoet_ids, cl_biotoop_only, NA)
  r_finaal  <- terra::mask(masker, cl_finaal)
  
  rm(stats_dt, cl_biotoop_only, zonal_res)
  gc()
  
  return(list(raster = r_finaal, clusters = cl_finaal))
}

# 3. SOORT INFORMATIE OPHALEN --------------------------------------------------
df <- read_excel(here("data/input/Excel_files/Soorten_bwk_afstanden.xlsx"))
soort <- "beekrombout"

resultaat <- df %>%
  filter(tolower(trimws(Soort)) == soort) %>%
  select(Type, MinOpp_ha, AfstandBiotopen_m, Dispersiecap_m)

oppervlakte_ha <- resultaat$MinOpp_ha
afstand_m      <- resultaat$AfstandBiotopen_m
buffer_m       <- resultaat$Dispersiecap_m[1]
straal_water_m <- 500  # Afsnij-afstand van land rond de beken (500m)

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

waterlopen1_bwk_opp  <- if("waterlopen1_bwk" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["waterlopen1_bwk"]] else template_Vlaanderen
landbiotoop_bos_opp  <- if("landbiotoop_bos_bwk" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["landbiotoop_bos_bwk"]] else template_Vlaanderen
landbiotoop_rest_opp <- if("landbiotoop_rest_bwk" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["landbiotoop_rest_bwk"]] else template_Vlaanderen

rm(ds_vlaanderen, resultaten_gegroepeerd, df_nieuw, lijst_oppervlaktes)
gc()

# 6. HYDROGRAFISCHE DATA (HUET-ZONES VLAANDEREN) INLADEN ----------------------
message("-> Huet-zone waterlopen inlezen voor heel Vlaanderen...")
r_huetzon_vlaanderen <- rast(here("data/input/Raster_Vlaanderen/vlaanderen_huetzon_10m.tif"))

if(!terra::compareGeom(r_huetzon_vlaanderen, template_Vlaanderen, stopOnError = FALSE)) {
  r_huetzon_vlaanderen <- terra::resample(r_huetzon_vlaanderen, template_Vlaanderen, method = "near")
}

legende_huet  <- terra::cats(r_huetzon_vlaanderen)[[1]]
target_labels <- c("Bg", "BgK", "Bk1", "Bk2", "BkK1", "BkK2")

valide_ids <- legende_huet$value[legende_huet$Label %in% target_labels]
message(paste("-> Gekoppelde Huet-zone ID's op basis van labels:", paste(valide_ids, collapse = ", ")))

r_waterlopen_bin <- r_huetzon_vlaanderen %in% valide_ids
r_waterlopen_bin <- terra::ifel(r_waterlopen_bin == 1, 1, NA)

# FUSIE WATERLOPEN (BWK + HUET-ZONES)
r_waterlopen_opp <- terra::ifel(r_waterlopen_bin == 1, 0.01, NA)
waterlopen_opp   <- terra::cover(waterlopen1_bwk_opp, r_waterlopen_opp)

rm(r_huetzon_vlaanderen, r_waterlopen_bin, r_waterlopen_opp, waterlopen1_bwk_opp)
gc()

# 7. RUIMTELIJKE CLUSTERING & OPPERVLAKTEFILTER WATERLOPEN ---------------------
message("-> Start ruimtelijke clustering waterlopen (min. 0.4 ha)...")

r_binair_water_opp <- terra::ifel(!is.na(waterlopen_opp) & waterlopen_opp > 0, 1, NA)
w_matrix <- matrix(1, nrow = 3, ncol = 3) 

tmp_focal    <- tempfile(fileext = ".tif")
r_bridge_opp <- terra::focal(r_binair_water_opp, w = w_matrix, fun = "max", na.rm = TRUE, 
                             filename = tmp_focal, overwrite = TRUE)

tmp_patch  <- tempfile(fileext = ".tif")
cl_id_opp  <- terra::patches(r_bridge_opp, directions = 8, zeroAsNA = TRUE, 
                             filename = tmp_patch, overwrite = TRUE)

stats_opp <- terra::zonal(waterlopen_opp, cl_id_opp, fun = "sum", na.rm = TRUE)
colnames(stats_opp) <- c("Cluster_ID", "Werkelijke_ha")

valide_cluster_ids_opp <- stats_opp$Cluster_ID[stats_opp$Werkelijke_ha >= 0.4]

if(length(valide_cluster_ids_opp) > 0) {
  masker_opp     <- cl_id_opp %in% valide_cluster_ids_opp
  waterlopen_opp <- terra::mask(waterlopen_opp, terra::ifel(masker_opp, 1, NA))
} else {
  waterlopen_opp <- template_Vlaanderen
}

unlink(c(tmp_focal, tmp_patch))
rm(r_binair_water_opp, r_bridge_opp, cl_id_opp, stats_opp)
gc()

# 8. LANDBIOTOOP CLUSTEREN & FILTEREN -----------------------------------------
message("-> Clusteren en filteren van landbiotoop-lagen...")

r_bin_bos   <- terra::ifel(!is.na(landbiotoop_bos_opp) & landbiotoop_bos_opp > 0, 1, NA)
bos_res_opp <- cluster_filter_compleet(r_bin_bos, opp_laag = landbiotoop_bos_opp, drempel_m2 = 50000, dist_m = 50, werkelijk = TRUE)
bos_cluster_opp <- bos_res_opp$raster

r_bin_rest   <- terra::ifel(!is.na(landbiotoop_rest_opp) & landbiotoop_rest_opp > 0, 1, NA)
rest_res_opp <- cluster_filter_compleet(r_bin_rest, opp_laag = landbiotoop_rest_opp, drempel_m2 = 500000, dist_m = 50, werkelijk = TRUE)
rest_cluster_opp <- rest_res_opp$raster

# BOS FILTEREN IN LANDMATRICE: OPTIE A (100% BINNEN MATRIX MET MAX. 10M SPELING)
filter_bos_in_landmatrice <- function(r_bos_clusters, r_rest_clusters, min_overlap_fractie = 1.0) {
  if (is.null(r_bos_clusters) || is.null(r_rest_clusters)) return(template_Vlaanderen)
  if (all(is.na(terra::values(r_bos_clusters, mat = FALSE)))) return(template_Vlaanderen)
  if (all(is.na(terra::values(r_rest_clusters, mat = FALSE)))) return(template_Vlaanderen)
  
  # 10m buffer/speling rond de goedgekeurde rest-matrix (3x3 filter op 10m resolutie)
  w_matrix  <- matrix(1, nrow = 3, ncol = 3)
  rest_zone <- terra::focal(!is.na(r_rest_clusters) & r_rest_clusters > 0, w = w_matrix, fun = "max", na.rm = TRUE)
  
  r_bos_bin <- terra::ifel(!is.na(r_bos_clusters) & r_bos_clusters > 0, 1, NA)
  cl_bos_id <- terra::patches(r_bos_bin, directions = 8, zeroAsNA = TRUE)
  
  # rest_zone is 1 waar matrix + 10m aanwezig is, 0 waar afwezig
  rest_zone_bin <- terra::ifel(!is.na(rest_zone) & rest_zone > 0, 1, 0)
  
  # Zonal sums voor overlap en clustergroottes
  overlap_som <- terra::zonal(rest_zone_bin, cl_bos_id, fun = "sum", na.rm = TRUE)
  totaal_bos  <- terra::freq(cl_bos_id)
  
  dt_overlap <- as.data.table(overlap_som)
  setnames(dt_overlap, c("ID", "Overlap_Cells"))
  
  dt_totaal <- as.data.table(totaal_bos)
  setnames(dt_totaal, c("layer", "ID", "Total_Cells"))
  
  merged_stats <- merge(dt_overlap, dt_totaal, by = "ID", all.x = TRUE)
  merged_stats[, Overlap_Frac := Overlap_Cells / Total_Cells]
  
  # OPTIE A: Het boscluster moet voor de volle 100% binnen rest_zone vallen (Overlap_Frac >= 1.0)
  # (Merk op: >= 0.9999 vangt eventuele kleine floating-point afrondingen af)
  geldige_bos_ids <- merged_stats[!is.na(ID) & Overlap_Frac >= (min_overlap_fractie - 1e-6), ID]
  
  if (length(geldige_bos_ids) == 0) return(template_Vlaanderen)
  
  return(terra::mask(r_bos_clusters, cl_bos_id %in% geldige_bos_ids, maskvalues = FALSE))
}

message("-> Filteren van boskernen op volledige (100%) inbedding in landbiotoop rest (met max. 10m speling)...")
bos_goedgekeurd_opp <- filter_bos_in_landmatrice(bos_cluster_opp, rest_cluster_opp, min_overlap_fractie = 1.0)

rm(landbiotoop_bos_opp, landbiotoop_rest_opp, bos_cluster_opp, rest_cluster_opp, bos_res_opp, rest_res_opp)
gc()

# 9. AFSTANDSFILTERS LAND EN WATER (500M) -------------------------------------
message("-> Wederzijdse afstandsfilters tussen water en boskernen (500m)...")

if (terra::hasValues(waterlopen_opp) && terra::hasValues(bos_goedgekeurd_opp) &&
    terra::global(waterlopen_opp, "notNA")$notNA > 0 && 
    terra::global(bos_goedgekeurd_opp, "notNA")$notNA > 0) {
  
  # Afstand tot goedgekeurd bos berekenen
  r_bos_bin <- terra::ifel(!is.na(bos_goedgekeurd_opp) & bos_goedgekeurd_opp > 0, 1, NA)
  dist_naar_bos <- terra::distance(r_bos_bin, maxdist = 505)
  zone_bos_500m <- terra::ifel(dist_naar_bos <= 500, 1, NA)
  
  beekrombout_waterbiotoop_opp <- terra::mask(waterlopen_opp, zone_bos_500m)
  
  if (!all(is.na(suppressWarnings(terra::minmax(beekrombout_waterbiotoop_opp))))) {
    # Snelle en stabiele rasterafstandsfilter voor land rond water
    r_water_bin <- terra::ifel(!is.na(beekrombout_waterbiotoop_opp) & beekrombout_waterbiotoop_opp > 0, 1, NA)
    dist_naar_water <- terra::distance(r_water_bin, maxdist = 505)
    zone_water_500m <- terra::ifel(dist_naar_water <= 500, 1, NA)
    
    beekrombout_landbiotoop_opp <- terra::mask(bos_goedgekeurd_opp, zone_water_500m)
    rm(r_water_bin, dist_naar_water, zone_water_500m)
  } else {
    beekrombout_landbiotoop_opp <- template_Vlaanderen
  }
  
  rm(r_bos_bin, dist_naar_bos, zone_bos_500m)
} else {
  beekrombout_waterbiotoop_opp <- template_Vlaanderen
  beekrombout_landbiotoop_opp  <- template_Vlaanderen
}

gc()

# 10. DEFINITIEVE SAMENVOEGING & NETWERK ID'S ----------------------------------
message("-> Samenvoegen water- en landleefgebied...")
final_opp <- terra::cover(beekrombout_waterbiotoop_opp, beekrombout_landbiotoop_opp)

if (!all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  tmp_patch <- tempfile(fileext = ".tif")
  cl_opp    <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE, 
                              filename = tmp_patch, overwrite = TRUE)
} else {
  cl_opp <- template_Vlaanderen
}

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

message("🏁 LEEFGEBIEDENGRID BEEKROMBOUT SUCCESVOL GEËXPORTEERD VOOR HEEL VLAANDEREN!")
