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

# ==============================================================================
# 6. DRAINAGE MASKERS OPBOUWEN
# ==============================================================================
message("-> Bodemdrainage maskers opbouwen (Strikt vs. Inclusief NA)...")

r_drain_raw <- rast(here("data/input/Raster_Vlaanderen/vlaanderen_drainage_10m.tif"))
crs(r_drain_raw) <- crs(template_Vlaanderen)
r_drain_aligned  <- terra::resample(r_drain_raw, template_Vlaanderen, method = "near")

drain_cats <- terra::cats(r_drain_aligned)[[1]]
veldkrekel_drain_ids <- drain_cats$value[drain_cats$Label %in% c("a", "b", "c", "a-b")]

# Masker A: STRIKT (Enkel bekende droge bodems)
masker_strikt    <- terra::ifel(r_drain_aligned %in% veldkrekel_drain_ids, 1, NA)
# Masker B: INCLUSIEF NA (Droog óf Onbekend)
masker_inclusief <- terra::ifel(r_drain_aligned %in% veldkrekel_drain_ids | is.na(r_drain_aligned), 1, NA)

veldkrekel_basis_strikt    <- terra::mask(bwk_opp, masker_strikt)
veldkrekel_basis_inclusief <- terra::mask(bwk_opp, masker_inclusief)

rm(r_drain_raw, r_drain_aligned, masker_strikt, masker_inclusief, bwk_opp)
gc()

# ==============================================================================
# HELPER FUNCTIE VOOR VOLLEDIGE PIPELINE (Clustering + Metapopulatie)
# ==============================================================================
run_volledige_pijplijn <- function(basis_laag, label_naam) {
  message(paste0("\n=== START PIPELINE VOOR SCENARIO: ", toupper(label_naam), " ==="))
  
  # 1. Clustering (50m)
  r_binair <- terra::ifel(!is.na(basis_laag) & basis_laag > 0, 1, NA)
  cl_out   <- cluster_filter_compleet(
    masker = r_binair, opp_laag = basis_laag, drempel_m2 = 10000, dist_m = 50, werkelijk = TRUE
  )
  
  r_patches <- cl_out$clusters
  if (all(is.na(terra::values(r_patches, mat = FALSE)))) {
    return(list(habitat = template_Vlaanderen * NA, ids = template_Vlaanderen * NA))
  }
  
  # 2. Metapopulatie analyse
  df_patches <- terra::as.data.frame(r_patches, cells = TRUE)
  colnames(df_patches) <- c("cell", "ID")
  
  opp_val <- terra::extract(basis_laag, df_patches$cell)[[1]]
  dt_meta <- data.table(ID = df_patches$ID, Waarde = opp_val)[!is.na(ID) & !is.na(Waarde)]
  stats_ha <- dt_meta[, .(Grootte_ha = sum(Waarde, na.rm = TRUE) * 0.01), by = ID]
  
  ids_groot_5ha    <- stats_ha[Grootte_ha >= 5, ID]
  ids_klein_1to5ha <- stats_ha[Grootte_ha >= 1 & Grootte_ha < 5, ID]
  
  r_5ha <- template_Vlaanderen * NA
  r_meta <- template_Vlaanderen * NA
  
  if (length(ids_groot_5ha) > 0) {
    r_5ha <- terra::mask(basis_laag, r_patches %in% ids_groot_5ha)
  }
  
  if (length(ids_klein_1to5ha) > 0) {
    r_klein <- terra::ifel(r_patches %in% ids_klein_1to5ha, r_patches, NA)
    df_pts  <- terra::as.data.frame(r_klein, xy = TRUE)
    colnames(df_pts) <- c("x", "y", "ID")
    dt_cent <- as.data.table(df_pts)[, .(x = mean(x), y = mean(y)), by = ID]
    
    if (nrow(dt_cent) > 0) {
      sf_cent <- st_as_sf(dt_cent, coords = c("x", "y"), crs = crs(template_Vlaanderen))
      dt_cent$buren_count <- lengths(st_is_within_distance(sf_cent, sf_cent, dist = 500))
      valid_ids <- dt_cent[buren_count >= 6, ID]
      
      if (length(valid_ids) > 0) {
        r_meta <- terra::mask(basis_laag, r_patches %in% valid_ids)
      }
    }
  }
  
  w_clean <- terra::ifel(is.na(r_5ha), 0, r_5ha)
  m_clean <- terra::ifel(is.na(r_meta), 0, r_meta)
  som_cl  <- terra::clamp(w_clean + m_clean, upper = 1.0)
  
  finaal_hab <- terra::ifel(som_cl > 0, 1, NA) # Binaire output (1/NA)
  return(list(habitat = finaal_hab, ids = r_patches))
}

# ==============================================================================
# 7 & 8. UITVOEREN VAN BEIDE SCENARIO'S
# ==============================================================================
res_strikt    <- run_volledige_pijplijn(veldkrekel_basis_strikt, "Strikt (enkel a,b,c)")
res_inclusief <- run_volledige_pijplijn(veldkrekel_basis_inclusief, "Inclusief NA")

# ==============================================================================
# BEREKENING VAN DE GECOMBINEERDE VERSCHILKAART (0/1/2)
# ==============================================================================
message("\n-> Verschilkaart berekenen (2 = Zeker leefgebied, 1 = Potentieel leefgebied door NA)...")

hab_strikt    <- terra::ifel(is.na(res_strikt$habitat), 0, 1)
hab_inclusief <- terra::ifel(is.na(res_inclusief$habitat), 0, 1)

# Logica: 
# Als Inclusief = 1 en Strikt = 1  --> 1 + 1 = 2 (Zeker geschikte bodem)
# Als Inclusief = 1 en Strikt = 0  --> 1 + 0 = 1 (Afhankelijk van NA-bodem)
# Als Inclusief = 0 en Strikt = 0  --> 0 + 0 = 0 (Geen leefgebied)
r_vergelijking <- hab_inclusief + hab_strikt
r_gradatie_finaal <- terra::ifel(r_vergelijking > 0, r_vergelijking, NA)

# ==============================================================================
# 9. DYNAMISCHE EXPORT MAKEN
# ==============================================================================
message("-> Start geformatteerde export voor Vlaanderen...")

base_dir <- here::here("data/output/Vlaanderen/Rasters_Soorten", scenario_naam)

folders <- list(
  id_raster  = file.path(base_dir, "00_ID_Rasters"),
  strikt     = file.path(base_dir, "01_Strikt_Bodem_ABC"),
  inclusief  = file.path(base_dir, "02_Inclusief_Bodem_NA"),
  gradatie   = file.path(base_dir, "03_Verschilkaart_Gradatie")
)
purrr::walk(folders, ~if (!dir.exists(.x)) dir.create(.x, showWarnings = FALSE, recursive = TRUE))

# A. Export Strikt
terra::writeRaster(res_strikt$habitat, file.path(folders$strikt, paste0("Habitat_Strikt_", soort, ".tif")),
                   overwrite = TRUE, gdal = c("COMPRESS=LZW"), datatype = "INT1U", NAflag = 255)

# B. Export Inclusief NA
terra::writeRaster(res_inclusief$habitat, file.path(folders$inclusief, paste0("Habitat_Inclusief_NA_", soort, ".tif")),
                   overwrite = TRUE, gdal = c("COMPRESS=LZW"), datatype = "INT1U", NAflag = 255)

# C. Export Gradatiekaart (2 = Kernhabitat op ABC bodem, 1 = Onzeker habitat op NA bodem)
terra::writeRaster(r_gradatie_finaal, file.path(folders$gradatie, paste0("Habitat_Gradatie_Verschil_", soort, ".tif")),
                   overwrite = TRUE, gdal = c("COMPRESS=LZW"), datatype = "INT1U", NAflag = 0)

# D. Export ID Raster (van het Inclusief scenario voor vervolganalyses)
terra::writeRaster(res_inclusief$ids, file.path(folders$id_raster, paste0("ID_Netwerken_", soort, ".tif")),
                   overwrite = TRUE, gdal = c("COMPRESS=LZW"), datatype = "INT4U", NAflag = 0)

message("🏁 BEIDE SCENARIO'S EN VERSCHILKAART SUCCESVOL GEËXPORTEERD!")
