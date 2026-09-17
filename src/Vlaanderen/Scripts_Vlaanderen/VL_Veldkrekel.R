# ==============================================================================
# LEEFGEBIEDENGRID VOOR VELDKREKEL - HEEL VLAANDEREN (WERKELIJKE OPPERVLAKTE)
# Author: Bert Van Hecke
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
  memfrac = 0.8,# Dwing terra om tot max. 80% van het RAM-geheugen te gebruiken
  todisk = TRUE,
  tempdir = tempdir(),  # Geef toestemming voor automatische disk-swapping bij zware rasters
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
  
  # 2. Netwerkvorming via tijdelijke schijfbestanden (voorkomt RAM-volloop)
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
  
  # ----------------------------------------------------------------------------
  # STAP 3: SLIMME DATA-EXTRACTIE (ENKEL ACTIEVE PIXELS, GEEN 600 MILJOEN CELLEN)
  # ----------------------------------------------------------------------------
  # Haal ALLEEN cellen op waar een cluster-ID aanwezig is
  df_cl <- terra::as.data.frame(cl_biotoop_only, cells = TRUE)
  
  if (nrow(df_cl) == 0) {
    return(list(raster = masker * NA, clusters = masker * NA))
  }
  
  colnames(df_cl) <- c("cell", "ID")
  
  if (werkelijk) {
    # Haal alleen de oppervlaktes op voor die specifieke actieve cel-indexen
    opp_waarden <- terra::extract(opp_laag, df_cl$cell)[[1]]
    dt_calc <- data.table(ID = df_cl$ID, Waarde = opp_waarden)
  } else {
    dt_calc <- data.table(ID = df_cl$ID, Waarde = 1)
  }
  
  # Sommeer oppervlakte per ID (in m2)
  stats_dt <- dt_calc[!is.na(ID) & !is.na(Waarde), .(Area_m2 = sum(Waarde, na.rm = TRUE) * 100), by = ID]
  
  voldoet_ids <- stats_dt[Area_m2 >= drempel_m2, ID]
  
  if (length(voldoet_ids) == 0) {
    return(list(raster = masker * NA, clusters = masker * NA))
  }
  
  # ----------------------------------------------------------------------------
  # STAP 4: FILTEREN EN SNELE RECONSTRUCTIE
  # ----------------------------------------------------------------------------
  cl_finaal <- terra::ifel(cl_biotoop_only %in% voldoet_ids, cl_biotoop_only, NA)
  r_finaal  <- terra::mask(masker, cl_finaal)
  
  rm(df_cl, dt_calc, stats_dt, cl_biotoop_only)
  gc()
  
  return(list(raster = r_finaal, clusters = cl_finaal))
}

# 3. SOORT INFORMATIE & SCENARIO OPHALEN ----------------------------------------
df <- read_excel(here("data/input/Excel_files/Soorten_bwk_afstanden.xlsx"))
soort <- "veldkrekel"

# Bepaal scenario naam (bijv. voor mappenstructuur)
scenario_naam <- "BWK_2025"

resultaat <- df %>%
  filter(Soort == soort) %>%
  select(Type, MinOpp_ha, AfstandBiotopen_m, Dispersiecap_m)

oppervlakte_ha <- resultaat$MinOpp_ha
afstand_m      <- resultaat$AfstandBiotopen_m
buffer_m       <- resultaat$Dispersiecap_m

print(resultaat)
rm(df, resultaat)

# 4. GEBIED VOORBEREIDEN (MASTERGRID VLAANDEREN) -------------------------------
message("-> Mastergrid Vlaanderen inlezen...")
master_grid <- rast(here("data/input/Raster_Vlaanderen/Vlaanderen_MasterGrid_10m.tif"))[[1]]

template_Vlaanderen <- master_grid
values(template_Vlaanderen) <- NA

# 5. BIOTOOPFILTERING VIA PARQUET ----------------------------------------------
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

# Open Parquet (Lazy evaluation via Arrow)
ds_vlaanderen <- open_dataset(here("data/input/Raster_Vlaanderen/BWK_TidyTabel_Smal_Vlaanderen_2025.parquet"))

h_data <- resultaten_gegroepeerd %>% filter(Type == "bwk") %>% pull(Data) %>% .[[1]]

exact_codes <- tolower(trimws(h_data$Code[h_data$Match == "exact"]))
bevat_codes <- tolower(trimws(h_data$Code[h_data$Match == "bevat"]))

if (length(bevat_codes) > 0) {
  bevat_codes_escaped  <- gsub("([\\.\\^\\$\\*\\+\\?\\(\\)\\[\\{\\\\\\|])", "\\\\\\1", bevat_codes)
  bevat_codes_anchored <- paste0("^", bevat_codes_escaped)
  regex_term           <- paste0(bevat_codes_anchored, collapse = "|")
  
  tabel_vlaanderen <- ds_vlaanderen %>%
    filter(tolower(CODE) %in% exact_codes | str_detect(tolower(CODE), regex_term)) %>%
    collect() %>%
    as.data.table()
} else {
  tabel_vlaanderen <- ds_vlaanderen %>%
    filter(tolower(CODE) %in% exact_codes) %>%
    collect() %>%
    as.data.table()
}

tabel_vlaanderen[, CODE := tolower(trimws(CODE))]
tabel_unique <- unique(tabel_vlaanderen, by = c("cel_id", "CODE"))

# Som per cel_id berekenen en aftoppen op max 1.0 (100%)
tabel_cel_som <- tabel_unique[, .(Oppervlakte = pmin(sum(BWK_FRAC, na.rm = TRUE), 1.0)), by = .(cel_id)]

# Direct wegschrijven naar template via globale cel_ids
bwk_opp <- template_Vlaanderen
bwk_opp[tabel_cel_som$cel_id] <- tabel_cel_som$Oppervlakte

rm(ds_vlaanderen, tabel_vlaanderen, tabel_unique, tabel_cel_som, resultaten_gegroepeerd, df_nieuw)
gc()

# 6. DRAINAGE FILTERING --------------------------------------------------------
message("-> Bodemdrainage toepassen op Vlaams niveau...")

r_drain_raw <- rast(here("data/input/Raster_Vlaanderen/vlaanderen_drainage_10m.tif"))

# CRS expliciet synchroniseren om "CRS do not match" te voorkomen
crs(r_drain_raw) <- crs(template_Vlaanderen)

# Herlijnen op template grid
r_drain_aligned <- terra::resample(r_drain_raw, template_Vlaanderen, method = "near")

drain_cats <- terra::cats(r_drain_aligned)[[1]]
geselecteerde_letters <- c("a", "b", "c", "a-b") 
veldkrekel_drain_ids <- drain_cats$value[drain_cats$Label %in% geselecteerde_letters]

masker_drainage <- terra::ifel(r_drain_aligned %in% veldkrekel_drain_ids, 1, NA)

veldkrekel_basis_opp <- terra::mask(bwk_opp, masker_drainage)

rm(r_drain_raw, r_drain_aligned, masker_drainage, drain_cats, bwk_opp)
gc()

# 7. CLUSTERING ----------------------------------------------------------------
message("-> Patches groeperen op basis van 50m onderlinge afstand...")

r_binair_basis_opp <- terra::ifel(!is.na(veldkrekel_basis_opp) & veldkrekel_basis_opp > 0, 1, NA)

veldkrekel_clusters_opp <- cluster_filter_compleet(
  masker     = r_binair_basis_opp,
  opp_laag   = veldkrekel_basis_opp,
  drempel_m2 = 10000,   # 1 ha minimum
  dist_m     = 50,       
  werkelijk  = TRUE      
)

rm(r_binair_basis_opp)
gc()

# 8. METAPOPULATIE STRUCTUUR ANALYSE ------------------------------------------
message("-> Metapopulatiestructuur berekenen voor Vlaanderen...")

id_export_rast <- veldkrekel_clusters_opp$clusters
r_patches_opp  <- veldkrekel_clusters_opp$clusters

r_leefgebied_5ha_opp     <- template_Vlaanderen * NA
r_leefgebied_metapop_opp <- template_Vlaanderen * NA

if (!all(is.na(terra::values(r_patches_opp, mat = FALSE)))) {
  
  # SLIMME EXTRACTIE: Alleen actieve cellen ophalen
  df_patches <- terra::as.data.frame(r_patches_opp, cells = TRUE)
  colnames(df_patches) <- c("cell", "ID")
  
  opp_val <- terra::extract(veldkrekel_basis_opp, df_patches$cell)[[1]]
  
  dt_meta <- data.table(ID = df_patches$ID, Waarde = opp_val)[!is.na(ID) & !is.na(Waarde)]
  stats_ha_opp <- dt_meta[, .(Grootte_ha = sum(Waarde, na.rm = TRUE) * 0.01), by = ID]
  
  ids_groot_5ha_opp    <- stats_ha_opp[Grootte_ha >= 5, ID]
  ids_klein_1to5ha_opp <- stats_ha_opp[Grootte_ha >= 1 & Grootte_ha < 5, ID]
  
  # A. Grote clusters (>= 5 ha) direct behouden
  if (length(ids_groot_5ha_opp) > 0) {
    r_leefgebied_5ha_opp <- terra::mask(veldkrekel_basis_opp, r_patches_opp %in% ids_groot_5ha_opp)
  }
  
  # B. Kleine clusters (1 tot 5 ha): controleer op >= 6 buren binnen 500m
  if (length(ids_klein_1to5ha_opp) > 0) {
    # 1. Isoleer kleine clusters en geef elke cluster 1 representatieve cel (centrogram)
    r_klein_mask <- r_patches_opp %in% ids_klein_1to5ha_opp
    r_klein_patches <- terra::ifel(r_klein_mask, r_patches_opp, NA)
    
    # 2. Maak binaire puntenkaart van kleine patches voor snelle buren-telling
    # Bepaal een unieke cel per cluster om dubbeltellingen van pixels te voorkomen
    df_pts <- terra::as.data.frame(r_klein_patches, cells = TRUE)
    colnames(df_pts) <- c("cell", "ID")
    df_unique_pts <- df_pts[!duplicated(df_pts$ID), ]
    
    r_pts_bin <- template_Vlaanderen * NA
    r_pts_bin[df_unique_pts$cell] <- 1
    
    # 3. Tel het aantal unieke kleine buren binnen 500m (straal = 50 cellen van 10m)
    # Cirkelvormige focal matrix van 500m
    f_matrix_500m <- terra::focalMat(template_Vlaanderen, d = 500, type = "circle")
    f_matrix_500m[f_matrix_500m > 0] <- 1 # Binaire aanwezigheidsmatrix
    
    tmp_focal_count <- tempfile(fileext = ".tif")
    r_buren_count <- terra::focal(r_pts_bin, w = f_matrix_500m, fun = "sum", na.rm = TRUE,
                                  filename = tmp_focal_count, overwrite = TRUE)
    
    # 4. Filter clusters die minimaal 6 buren hebben (inclusief zichzelf)
    r_metapop_goedgekeurd <- terra::mask(r_klein_patches, r_buren_count >= 6)
    
    if (!all(is.na(terra::values(r_metapop_goedgekeurd, mat = FALSE)))) {
      valid_meta_ids <- unique(na.omit(terra::values(r_metapop_goedgekeurd, mat = FALSE)))
      r_leefgebied_metapop_opp <- terra::mask(veldkrekel_basis_opp, r_patches_opp %in% valid_meta_ids)
    }
    
    unlink(tmp_focal_count)
    rm(r_klein_mask, r_klein_patches, df_pts, df_unique_pts, r_pts_bin, r_buren_count)
    gc()
  }
  
  # Combineren van grote clusters en goedgekeurde metapopulaties
  w_clean <- terra::ifel(is.na(r_leefgebied_5ha_opp), 0, r_leefgebied_5ha_opp)
  m_clean <- terra::ifel(is.na(r_leefgebied_metapop_opp), 0, r_leefgebied_metapop_opp)
  som_raw <- w_clean + m_clean
  som_cl  <- terra::clamp(som_raw, upper = 1.0)
  
  veldkrekel_leefgebied_opp <- terra::ifel(som_cl > 0, som_cl, NA)
  rm(w_clean, m_clean, som_raw, som_cl, r_leefgebied_5ha_opp, r_leefgebied_metapop_opp)
} else {
  veldkrekel_leefgebied_opp <- template_Vlaanderen * NA
}

cat("\n----------------------------------------------------\n")
cat("Finaal Werkelijk Leefgebied Veldkrekel Vlaanderen (ha):", round(calc_ha_exact(veldkrekel_leefgebied_opp), 2), "\n")
cat("----------------------------------------------------\n\n")

# Zorg dat id_export_rast HIER NIET wordt gewist!
suppressWarnings(rm(r_patches_opp, veldkrekel_clusters_opp, veldkrekel_basis_opp))
gc()

# 9. DYNAMISCHE EXPORT MAKEN (INCL. 00_ID_RASTERS) -----------------------------
message("-> Start geformatteerde export voor Vlaanderen...")

base_dir <- here::here("data/output/Vlaanderen/Rasters_Soorten", scenario_naam)

folders <- list(
  id_raster = file.path(base_dir, "00_ID_Rasters"),
  werkelijk = file.path(base_dir, "02_Werkelijke_Oppervlaktes")
)
purrr::walk(folders, ~if (!dir.exists(.x)) dir.create(.x, showWarnings = FALSE, recursive = TRUE))

# --- A. BINAIR EXPORT WERKELIJK HABITAT ---
if (exists("veldkrekel_leefgebied_opp") && !all(is.na(terra::values(veldkrekel_leefgebied_opp, mat = FALSE)))) {
  werkelijk_export_rast <- terra::ifel(!is.na(veldkrekel_leefgebied_opp) & veldkrekel_leefgebied_opp > 0, 1, NA)
} else {
  werkelijk_export_rast <- template_Vlaanderen * NA
}

file_path_werkelijk <- file.path(folders$werkelijk, paste0("Habitat_Werkelijke_Oppervlaktes_", soort, ".tif"))

terra::writeRaster(
  werkelijk_export_rast, 
  filename = file_path_werkelijk, 
  overwrite = TRUE, 
  gdal = c("COMPRESS=LZW"), 
  datatype = "INT1U",
  NAflag = 255
)
message(paste("    [OK] Werkelijk Habitat geëxporteerd:", basename(file_path_werkelijk)))

# --- B. ANALYTISCH ID-RASTER EXPORT (VOOR SCRIPT 2 / ARPL) ---
file_path_id <- file.path(folders$id_raster, paste0("ID_Netwerken_", soort, ".tif"))

# Exporteer ID raster met INT4U datatype om unieke netwerk-IDs te behouden
terra::writeRaster(
  id_export_rast, 
  filename = file_path_id, 
  overwrite = TRUE, 
  gdal = c("COMPRESS=LZW"), 
  datatype = "INT4U",
  NAflag = 0
)
message(paste("    [OK] Analytisch ID-raster geëxporteerd:", basename(file_path_id)))

suppressWarnings(rm(werkelijk_export_rast, id_export_rast))
gc()

message(paste("🏁 SCRIPT SUCCESVOL AFGEROND VOOR:", toupper(soort)))
