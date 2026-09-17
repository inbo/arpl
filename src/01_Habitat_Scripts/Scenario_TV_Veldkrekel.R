# ==============================================================================
# LEEFGEBIEDENGRID VOOR VELDKREKEL - HEEL VLAANDEREN (INCL. ANALYTISCH ID-RASTER)
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

conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("intersect", "terra")
conflicted::conflict_prefer("any", "terra")

# Terra-geheugenbeheer instellen tegen out-of-memory crashes
terraOptions(
  memfrac = 0.6,        # Maximaal 60% RAM voor Terra
  todisk = TRUE,        # Automatische disk-swapping bij zware rasters
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

# Geheugenefficiënte clusterfilter met Focal Dilatie (Voorkomt std::bad_alloc)
cluster_filter_compleet <- function(masker, opp_laag, drempel_m2, dist_m, werkelijk = FALSE) {
  if (terra::global(is.na(masker), "sum")[[1]] == terra::ncell(masker)) {
    return(list(raster = masker * NA, clusters = masker * NA))
  }
  
  # STAP 1: NETWERKVORMING VIA FOCAL DILATIE
  if (dist_m > 0) {
    r_binair <- terra::ifel(!is.na(masker) & masker > 0, 1, NA)
    
    # 50m netwerk = 25m buffer = 2.5 cellen -> 5x5 matrix
    stral_cellen <- ceiling((dist_m / 2) / 10)
    f_matrix <- matrix(1, nrow = (2 * stral_cellen + 1), ncol = (2 * stral_cellen + 1))
    
    r_buffered <- terra::focal(r_binair, w = f_matrix, fun = "max", na.rm = TRUE)
    r_buffered <- terra::ifel(r_buffered > 0, 1, NA)
    
    cl_network <- terra::patches(r_buffered, directions = 4, zeroAsNA = TRUE)
    cl_biotoop_only <- terra::mask(cl_network, masker)
    
    rm(r_binair, r_buffered, cl_network)
    gc()
  } else {
    cl_biotoop_only <- terra::patches(masker, directions = 8, zeroAsNA = TRUE)
  }
  
  # STAP 2: OPPERVLAKTE PER NETWERK BEREKENEN
  if (werkelijk) {
    stats_df <- terra::zonal(opp_laag, cl_biotoop_only, fun = "sum", na.rm = TRUE)
    colnames(stats_df) <- c("ID", "Waarde")
    stats_df$Area_m2 <- stats_df$Waarde * 100 
  } else {
    f <- terra::freq(cl_biotoop_only)
    stats_df <- data.frame(ID = f$value, Waarde = f$count)
    stats_df$Area_m2 <- stats_df$Waarde * 100 
  }
  
  stats_df <- stats_df[!is.na(stats_df$ID), ]
  if (nrow(stats_df) == 0) return(list(raster = masker * NA, clusters = masker * NA))
  
  # STAP 3: FILTEREN OP DREMPELWAARDE
  voldoet_ids <- stats_df$ID[stats_df$Area_m2 >= drempel_m2]
  if (length(voldoet_ids) == 0) return(list(raster = masker * NA, clusters = masker * NA))
  
  masker_binair <- cl_biotoop_only %in% voldoet_ids
  final_network_mask <- terra::ifel(masker_binair == 1, 1, NA)
  
  r_finaal  <- terra::mask(masker, final_network_mask)
  cl_finaal <- terra::mask(cl_biotoop_only, r_finaal) 
  
  return(list(raster = r_finaal, clusters = cl_finaal))
}

# 3. SOORT INFORMATIE & SCENARIO OPHALEN ----------------------------------------
df <- read_excel(here::here("data/input/Excel_files/Soorten_bwk_afstanden.xlsx"))
soort <- "veldkrekel"

# Dynamic check op scenario rds
if (!exists("params") || is.null(params$scenario_rds_path)) {
  scenario_rds_path <- "data/input/Raster_Vlaanderen/BWK_TidyTabel_Smal_Vlaanderen_2025.parquet"
} else {
  scenario_rds_path <- params$scenario_rds_path
}

p_raw <- gsub("^([.][.]/)+", "", scenario_rds_path)
scenario_path <- here::here(p_raw)

scen_volledig <- basename(scenario_path)
scenario_naam <- gsub("^Vlaanderen_Scenario_|^Scenario_|.rds$|.parquet$", "", scen_volledig)

message(paste("Verwerken van soort:", soort, "binnen scenario:", scenario_naam))

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

# 5. BIOTOOPFILTERING VIA PARQUET (ARROW) --------------------------------------
message("-> BWK Parquet dataset koppelen via Arrow...")

df_nieuw <- read_csv(here("data/input/Excel_files/Resultaten_Totaal_Samengevoegd.csv"), show_col_types = FALSE)

resultaten_gegroepeerd <- df_nieuw %>%
  mutate(Soort_clean = tolower(trimws(Soort))) %>%
  filter(Soort_clean == soort) %>%
  group_by(Type) %>%
  nest(Data = c(Code, Match))

# Open Parquet via Arrow (Lazy evaluation)
ds_vlaanderen <- open_dataset(scenario_path)

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

tabel_cel_som <- tabel_unique[, .(Oppervlakte = pmin(sum(BWK_FRAC, na.rm = TRUE), 1.0)), by = .(cel_id)]

bwk_opp <- template_Vlaanderen
bwk_opp[tabel_cel_som$cel_id] <- tabel_cel_som$Oppervlakte

rm(ds_vlaanderen, tabel_vlaanderen, tabel_unique, tabel_cel_som, resultaten_gegroepeerd, df_nieuw)
gc()

# 6. DRAINAGE FILTERING --------------------------------------------------------
message("-> Bodemdrainage toepassen op Vlaams niveau...")

r_drain_raw <- rast(here("data/input/Raster_Vlaanderen/vlaanderen_drainage_10m.tif"))

# Exact herlijnen op mastergrid
r_drain_aligned <- terra::resample(r_drain_raw, template_Vlaanderen, method = "near")

drain_cats <- terra::cats(r_drain_aligned)[[1]]
geselecteerde_letters <- c("a", "b", "c", "a-b") 
veldkrekel_drain_ids <- drain_cats$value[drain_cats$Label %in% geselecteerde_letters]

masker_drainage <- terra::ifel(r_drain_aligned %in% veldkrekel_drain_ids, 1, NA)

veldkrekel_basis_opp <- terra::mask(bwk_opp, masker_drainage)

rm(r_drain_raw, r_drain_aligned, masker_drainage, drain_cats, bwk_opp)
gc()

# 7. CLUSTERING & ID-RASTER BEREKENING -----------------------------------------
message("-> Patches groeperen op basis van 50m onderlinge afstand...")

r_binair_basis_opp <- terra::ifel(!is.na(veldkrekel_basis_opp) & veldkrekel_basis_opp > 0, 1, NA)

veldkrekel_clusters_opp <- cluster_filter_compleet(
  masker     = r_binair_basis_opp,
  opp_laag   = veldkrekel_basis_opp,
  drempel_m2 = 10000,   # 1 ha minimum
  dist_m     = 50,       
  werkelijk  = TRUE      
)

# Bewaar het analytische ID-raster van de netwerken
id_export_rast <- veldkrekel_clusters_opp$clusters

rm(r_binair_basis_opp)
gc()

# 8. METAPOPULATIE STRUCTUUR ANALYSE ------------------------------------------
message("-> Metapopulatiestructuur berekenen voor Vlaanderen...")

r_patches_opp <- veldkrekel_clusters_opp$clusters
r_leefgebied_5ha_opp     <- template_Vlaanderen * NA
r_leefgebied_metapop_opp <- template_Vlaanderen * NA

if (!all(is.na(terra::values(r_patches_opp, mat = FALSE)))) {
  stats_ha_opp <- terra::zonal(veldkrekel_basis_opp, r_patches_opp, fun = "sum", na.rm = TRUE)
  colnames(stats_ha_opp) <- c("ID", "Grootte_ha")
  stats_ha_opp$Grootte_ha <- stats_ha_opp$Grootte_ha * 0.01
  
  ids_groot_5ha_opp    <- stats_ha_opp$ID[stats_ha_opp$Grootte_ha >= 5]
  ids_klein_1to5ha_opp <- stats_ha_opp$ID[stats_ha_opp$Grootte_ha >= 1 & stats_ha_opp$Grootte_ha < 5]
  
  if (length(ids_groot_5ha_opp) > 0) {
    r_leefgebied_5ha_opp <- terra::mask(veldkrekel_basis_opp, r_patches_opp %in% ids_groot_5ha_opp)
  }
  
  if (length(ids_klein_1to5ha_opp) > 0) {
    r_klein_patches_opp <- terra::ifel(r_patches_opp %in% ids_klein_1to5ha_opp, r_patches_opp, NA)
    p_klein_opp <- terra::as.polygons(r_klein_patches_opp, dissolve = FALSE)
    
    if (!is.null(p_klein_opp) && nrow(p_klein_opp) > 0) {
      poly_buffer_opp <- terra::buffer(p_klein_opp, width = 500)
      intersect_matrix_opp <- matrix(terra::is.related(poly_buffer_opp, p_klein_opp, "intersects"), 
                                     nrow = nrow(poly_buffer_opp), ncol = nrow(p_klein_opp))
      
      p_klein_opp$unieke_buren_count <- rowSums(intersect_matrix_opp)
      goedgekeurde_metapop_opp <- p_klein_opp[p_klein_opp$unieke_buren_count >= 6, ]
      
      if (nrow(goedgekeurde_metapop_opp) > 0) {
        r_meta_mask_opp <- terra::rasterize(goedgekeurde_metapop_opp, template_Vlaanderen, field = 1, background = NA)
        r_leefgebied_metapop_opp <- terra::mask(veldkrekel_basis_opp, r_meta_mask_opp)
      }
      suppressWarnings(rm(poly_buffer_opp, intersect_matrix_opp, goedgekeurde_metapop_opp, p_klein_opp))
    }
  }
  
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
