# ==============================================================================
# LEEFGEBIEDENGRID VOOR BLAUVWLEUGELSPRINKHAAN - HEEL VLAANDEREN (WERKELIJKE OPPERVLAKTE)
# Ecologische analyse droge habitat, metapopulatiestructuur & netwerkeisen (INBO)
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
soort <- "blauwvleugelsprinkhaan"

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

bwk_opp <- if(length(lijst_oppervlaktes) > 0) lijst_oppervlaktes[[1]] else template_Vlaanderen

rm(ds_vlaanderen, resultaten_gegroepeerd, df_nieuw, lijst_oppervlaktes)
gc()

# 6. STAP 1: OMGEVINGSLAAG CONFIGURATIE (DROGE DRAINAGE) ------------------------
message("-> Drainagekaart inlezen en droge gronden filteren...")

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

gewenste_drainage_letters <- c("a", "b", "c", "a-b") 

drain_cats                    <- terra::cats(r_drainage)[[1]]
ids_drainage_droog            <- drain_cats$value[drain_cats$Label %in% gewenste_drainage_letters]
blauwvleugel_drainage_ok      <- r_drainage %in% ids_drainage_droog

blauwvleugelsprinkhaan_basis_opp <- terra::mask(bwk_opp, terra::ifel(blauwvleugel_drainage_ok, 1, NA))

rm(r_drainage, drain_cats, ids_drainage_droog, blauwvleugel_drainage_ok, bwk_opp)
gc()

# 7. STAP 2: PRIMAIRE LANDSCHAPSCLUSTERING (50 METER FUZZY CLUSTER) ------------
message("-> Patches groeperen op basis van 50m onderlinge afstand (min. 1 ha)...")

r_binair_basis_opp              <- terra::ifel(!is.na(blauwvleugelsprinkhaan_basis_opp) & blauwvleugelsprinkhaan_basis_opp > 0, 1, NA)
blauwvleugelsprinkhaan_clusters <- cluster_filter_compleet(
  masker     = r_binair_basis_opp,
  opp_laag   = blauwvleugelsprinkhaan_basis_opp,
  drempel_m2 = 10000,   # 1 ha minimum
  dist_m     = 50,        
  werkelijk  = TRUE        
)
r_patches_opp <- blauwvleugelsprinkhaan_clusters$clusters

rm(r_binair_basis_opp, blauwvleugelsprinkhaan_clusters)
gc()

# 8. STAP 3: METAPOPULATIE STRUCTUUR ANALYSE -----------------------------------
message("-> Grote patches (>= 5 ha) scheiden van kleine netwerkpatches (1-5 ha) met 250m koppelafstand...")

r_leefgebied_5ha_opp     <- template_Vlaanderen
r_leefgebied_metapop_opp <- template_Vlaanderen

if (!all(is.na(suppressWarnings(terra::minmax(r_patches_opp))))) {
  stats_ha_opp            <- terra::zonal(blauwvleugelsprinkhaan_basis_opp, r_patches_opp, fun = "sum", na.rm = TRUE)
  colnames(stats_ha_opp)  <- c("ID", "Grootte_ha_sum")
  stats_ha_opp$Grootte_ha <- stats_ha_opp$Grootte_ha_sum * 0.01
  
  dt_stats_ha           <- as.data.table(stats_ha_opp)
  ids_groot_5ha_opp     <- dt_stats_ha[Grootte_ha >= 5, ID]
  ids_klein_1to5ha_opp  <- dt_stats_ha[Grootte_ha >= 1 & Grootte_ha < 5, ID]
  
  if (length(ids_groot_5ha_opp) > 0) {
    r_leefgebied_5ha_opp <- terra::mask(blauwvleugelsprinkhaan_basis_opp, r_patches_opp %in% ids_groot_5ha_opp)
  }
  
  if (length(ids_klein_1to5ha_opp) > 0) {
    r_klein_patches_opp <- terra::ifel(r_patches_opp %in% ids_klein_1to5ha_opp, r_patches_opp, NA)
    p_klein_opp         <- terra::as.polygons(r_klein_patches_opp, dissolve = FALSE)
    
    if (!is.null(p_klein_opp) && nrow(p_klein_opp) > 0) {
      poly_buffer_opp      <- terra::buffer(p_klein_opp, width = buffer_m) # 250m
      intersect_matrix_opp <- matrix(terra::is.related(poly_buffer_opp, p_klein_opp, "intersects"), 
                                     nrow = nrow(poly_buffer_opp), ncol = nrow(p_klein_opp))
      
      p_klein_opp$unieke_buren_count <- rowSums(intersect_matrix_opp)
      goedgekeurde_metapop_opp       <- p_klein_opp[p_klein_opp$unieke_buren_count >= 3, ]
      
      if (nrow(goedgekeurde_metapop_opp) > 0) {
        tmp_meta_rast            <- tempfile(fileext = ".tif")
        r_meta_mask_opp          <- terra::rasterize(goedgekeurde_metapop_opp, template_Vlaanderen, field = 1, background = NA, filename = tmp_meta_rast, overwrite = TRUE)
        r_leefgebied_metapop_opp <- terra::mask(blauwvleugelsprinkhaan_basis_opp, r_meta_mask_opp)
        unlink(tmp_meta_rast)
      }
      suppressWarnings(rm(poly_buffer_opp, intersect_matrix_opp, goedgekeurde_metapop_opp))
    }
    rm(r_klein_patches_opp, p_klein_opp)
  }
  rm(stats_ha_opp, dt_stats_ha)
}
gc()

# 9. STAP 4: FINALE SAMENVOEGING & NETWERK-ID EXPORT ---------------------------
message("-> Grote kernen en goedgekeurde netwerken samensmelten tot leefgebied...")

b_clean <- terra::ifel(is.na(r_leefgebied_5ha_opp), 0, r_leefgebied_5ha_opp)
m_clean <- terra::ifel(is.na(r_leefgebied_metapop_opp), 0, r_leefgebied_metapop_opp)

som_raw   <- b_clean + m_clean
som_cl    <- terra::clamp(som_raw, upper = 1.0)
final_opp <- terra::ifel(som_cl > 0, som_cl, NA)

if (!all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  tmp_patch_cl <- tempfile(fileext = ".tif")
  cl_opp       <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE, filename = tmp_patch_cl, overwrite = TRUE)
} else {
  cl_opp <- template_Vlaanderen
}

rm(blauwvleugelsprinkhaan_basis_opp, r_patches_opp, r_leefgebied_5ha_opp, r_leefgebied_metapop_opp, 
   b_clean, m_clean, som_raw, som_cl)
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

message("🏁 LEEFGEBIEDENGRID BLAUVWLEUGELSPRINKHAAN SUCCESVOL GEËXPORTEERD VOOR HEEL VLAANDEREN!")
