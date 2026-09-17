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

# 6. DRAINAGE FILTERING & GRADATIE OPBOUWEN -----------------------------------
message("-> Bodemdrainage toepassen en gradaties berekenen...")

r_drain_raw <- rast(here("data/input/Raster_Vlaanderen/vlaanderen_drainage_10m.tif"))
crs(r_drain_raw) <- crs(template_Vlaanderen)

r_drain_aligned <- terra::resample(r_drain_raw, template_Vlaanderen, method = "near")

drain_cats <- terra::cats(r_drain_aligned)[[1]]
geselecteerde_letters <- c("a", "b", "c", "a-b") 
veldkrekel_drain_ids <- drain_cats$value[drain_cats$Label %in% geselecteerde_letters]

# A. Masker 1: STRIKT DROOG (Enkel goedgekeurde drainageklassen)
masker_strikt <- terra::ifel(r_drain_aligned %in% veldkrekel_drain_ids, 1, NA)

# B. Masker 2: INCLUSIEF NA (Droog óf Onbekende drainage)
masker_inclusief <- terra::ifel(r_drain_aligned %in% veldkrekel_drain_ids | is.na(r_drain_aligned), 1, NA)

# C. Masker 3: GRADATIEKAART (2 = Zeker droog, 1 = Onbekend/NA, NA = Ongeschikt/Te nat)
masker_gradatie <- terra::ifel(
  r_drain_aligned %in% veldkrekel_drain_ids, 2,
  terra::ifel(is.na(r_drain_aligned), 1, NA)
)

# Biotoopoppervlaktes maskeren
veldkrekel_basis_strikt    <- terra::mask(bwk_opp, masker_strikt)
veldkrekel_basis_inclusief <- terra::mask(bwk_opp, masker_inclusief)

# Gradatiekoppeling op aanwezige BWK biotoop
veldkrekel_basis_gradatie  <- terra::mask(masker_gradatie, bwk_opp)

rm(r_drain_raw, r_drain_aligned, masker_strikt, masker_inclusief, masker_gradatie, bwk_opp)
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

# 8. METAPOPULATIE STRUCTUUR ANALYSE (VECTOR-OPTIMALISATIE) --------------------
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
  
  # B. Kleine clusters (1 tot 5 ha): supersnel testen via SF afstanden (geen raster-focal!)
  if (length(ids_klein_1to5ha_opp) > 0) {
    message(" -> Kleine clusters analyseren via snelle vector-afstand (500m)...")
    
    # 1. Isoleer kleine clusters
    r_klein_patches <- terra::ifel(r_patches_opp %in% ids_klein_1to5ha_opp, r_patches_opp, NA)
    
    # 2. Pak de unieke coördinaten (centroids) van elke kleine cluster
    df_klein_pts <- terra::as.data.frame(r_klein_patches, xy = TRUE)
    colnames(df_klein_pts) <- c("x", "y", "ID")
    
    # 1 punt per unieke cluster ID (gemiddelde X/Y voor het centrum)
    dt_centroids <- as.data.table(df_klein_pts)[, .(x = mean(x), y = mean(y)), by = ID]
    
    if (nrow(dt_centroids) > 0) {
      # Zet om naar een licht SF-puntenobject
      sf_centroids <- st_as_sf(dt_centroids, coords = c("x", "y"), crs = crs(template_Vlaanderen))
      
      # Bereken binnen 500m alle buren via st_is_within_distance (duurt < 2 seconden)
      buren_lijst <- st_is_within_distance(sf_centroids, sf_centroids, dist = 500)
      
      # Tel het aantal unieke buren (inclusief zichzelf) per cluster
      dt_centroids$buren_count <- lengths(buren_lijst)
      
      # Selecteer de cluster IDs die minimaal 6 buren hebben
      valid_meta_ids <- dt_centroids[buren_count >= 6, ID]
      
      if (length(valid_meta_ids) > 0) {
        r_leefgebied_metapop_opp <- terra::mask(veldkrekel_basis_opp, r_patches_opp %in% valid_meta_ids)
      }
      
      rm(sf_centroids, buren_lijst, dt_centroids, df_klein_pts)
    }
    rm(r_klein_patches)
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
# Isoleer de goedgekeurde netwerken op de gradatiekaart
veldkrekel_leefgebied_gradatie <- terra::mask(veldkrekel_basis_gradatie, veldkrekel_leefgebied_opp)
cat("\n----------------------------------------------------\n")
cat("Finaal Werkelijk Leefgebied Veldkrekel Vlaanderen (ha):", round(calc_ha_exact(veldkrekel_leefgebied_opp), 2), "\n")
cat("----------------------------------------------------\n\n")

suppressWarnings(rm(r_patches_opp, veldkrekel_clusters_opp, veldkrekel_basis_opp))
gc()

# 9. DYNAMISCHE EXPORT MAKEN ----------------------------------------------------
message("-> Start geformatteerde export voor Vlaanderen...")

base_dir <- here::here("data/output/Vlaanderen/Rasters_Soorten", scenario_naam)

folders <- list(
  id_raster  = file.path(base_dir, "00_ID_Rasters"),
  werkelijk  = file.path(base_dir, "02_Werkelijke_Oppervlaktes"),
  gradatie   = file.path(base_dir, "03_Geschiktheid_Gradatie")
)
purrr::walk(folders, ~if (!dir.exists(.x)) dir.create(.x, showWarnings = FALSE, recursive = TRUE))

# --- EXPORT GRADATIEKAART (0 = NA / Ongeschikt, 1 = Onbekend, 2 = Zeker geschikt) ---
file_path_gradatie <- file.path(folders$gradatie, paste0("Geschiktheid_Gradatie_", soort, ".tif"))

terra::writeRaster(
  veldkrekel_leefgebied_gradatie, 
  filename = file_path_gradatie, 
  overwrite = TRUE, 
  gdal = c("COMPRESS=LZW"), 
  datatype = "INT1U",
  NAflag = 0
)
message(paste("    [OK] Gradatiekaart (0/1/2) geëxporteerd:", basename(file_path_gradatie)))
