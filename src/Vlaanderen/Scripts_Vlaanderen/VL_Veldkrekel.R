# ==============================================================================
# LEEFGEBIEDENGRID VOOR VELDKREKEL - HEEL VLAANDEREN (WERKELIJKE OPPERVLAKTE)
# Gradatie-analyse conform Ecohydrologie & Expert-advies (INBO)
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

# 3. SOORT INFORMATIE & SCENARIO OPHALEN ----------------------------------------
df <- read_excel(here("data/input/Excel_files/Soorten_bwk_afstanden.xlsx"))
soort <- "veldkrekel"
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

# 5. BIOTOOPFILTERING VIA PARQUET (INCL. DROGE SUB-CODES) -----------------------
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

h_data <- resultaten_gegroepeerd %>% filter(Type == "bwk") %>% pull(Data) %>% .[[1]]

exact_codes <- tolower(trimws(h_data$Code[h_data$Match == "exact"]))
bevat_codes <- tolower(trimws(h_data$Code[h_data$Match == "bevat"]))

# VEILIGHEIDS-INJECTIE: Voeg ontbrekende droge heide & zandcodes toe (Natte heide 'ce' NIET opnemen)
extra_droge_codes <- c("cm-", "cmb", "cg-", "cgb", "ku", "ku-", "k(ku)")
bevat_codes <- unique(c(bevat_codes, extra_droge_codes))

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

# ==============================================================================
# 6. DRAINAGE MASKERS OPBOUWEN (GEÜPDATET VOOR N2KHAB MET VEEN-CORRECTIE)
# ==============================================================================
message("-> Bodemdrainage maskers opbouwen (N2KHAB: Strikte droogte vs. Inclusief NA/Onbekend)...")

# AANGEPAST: Nieuw raster met n2khab verwerking en Veen-correctie inschakelen
r_drain_raw <- rast(here("data/input/Raster_Vlaanderen/vlaanderen_drainage_n2khab_10m.tif"))
crs(r_drain_raw) <- crs(template_Vlaanderen)
r_drain_aligned  <- terra::resample(r_drain_raw, template_Vlaanderen, method = "near")

drain_cats <- terra::cats(r_drain_aligned)[[1]]

# A. OPTIMAAL DROOG (Enkel pure a, a-b, b, c)
ids_optimaal <- drain_cats$value[drain_cats$Label %in% c("a", "a-b", "b", "c")]

# B. NATTE KLASSEN DIE EXPLICIET UITGESLOTEN WORDEN
# Omvat e t/m i, natte combinaties én gecorrigeerd Veen ('g')
ids_ongeschikt_nat <- drain_cats$value[drain_cats$Label %in% c(
  "e", "e-f", "e-i", "f", "g", "h", "h-i", "i"
)]

# --- MASKER 1: STRIKT DROOG (Enkel bewezen a, a-b, b, c) ---
masker_strikt <- terra::ifel(r_drain_aligned %in% ids_optimaal, 1, NA)

# --- MASKER 2: INCLUSIEF SUBOPTIMAAL & NA ---
# Sluit ALLEEN expliciet natte klassen (inclusief gecorrigeerd Veen 'g') uit.
# Categorieën 'd', 'a-d', 'c-d', 'Onbekend' (MILITAIR OB, KUSTDUINEN X, POLDERS) stromen door!
masker_inclusief <- terra::ifel(r_drain_aligned %in% ids_ongeschikt_nat, NA, 1)

# Biotoop-oppervlaktes maskeren voor de 2 aparte runs
veldkrekel_basis_strikt    <- terra::mask(bwk_opp, masker_strikt)
veldkrekel_basis_inclusief <- terra::mask(bwk_opp, masker_inclusief)

rm(r_drain_raw, r_drain_aligned, masker_strikt, masker_inclusief, bwk_opp)
gc()

# ==============================================================================
# GEHEUGEN-VEILIGE PIPELINE FUNCTIE (CRASH-PROOF VOOR HEEL VLAANDEREN)
# ==============================================================================
run_volledige_pijplijn <- function(basis_laag, label_naam) {
  message(paste0("\n=== START PIPELINE VOOR SCENARIO: ", toupper(label_naam), " ==="))
  
  # 1. Clustering (50m) via snelle schijf-buffering
  r_binair <- terra::ifel(!is.na(basis_laag) & basis_laag > 0, 1, NA)
  cl_out   <- cluster_filter_compleet(
    masker = r_binair, opp_laag = basis_laag, drempel_m2 = 10000, dist_m = 50, werkelijk = TRUE
  )
  
  r_patches <- cl_out$clusters
  rm(r_binair)
  gc()
  
  if (all(is.na(terra::values(r_patches, mat = FALSE)))) {
    return(list(habitat = template_Vlaanderen * NA, ids = template_Vlaanderen * NA))
  }
  
  # 2. Oppervlakte per patch berekenen
  df_patches <- terra::as.data.frame(r_patches, cells = TRUE)
  colnames(df_patches) <- c("cell", "ID")
  
  opp_val <- terra::extract(basis_laag, df_patches$cell)[[1]]
  dt_meta <- data.table(ID = df_patches$ID, Waarde = opp_val)[!is.na(ID) & !is.na(Waarde)]
  rm(opp_val, df_patches)
  gc()
  
  stats_ha <- dt_meta[, .(Grootte_ha = sum(Waarde, na.rm = TRUE) * 0.01), by = ID]
  rm(dt_meta)
  gc()
  
  ids_groot_5ha    <- stats_ha[Grootte_ha >= 5, ID]
  ids_klein_1to5ha <- stats_ha[Grootte_ha >= 1 & Grootte_ha < 5, ID]
  
  r_5ha  <- template_Vlaanderen * NA
  r_meta <- template_Vlaanderen * NA
  
  # Grote patches (>= 5 ha) direct insluiten
  if (length(ids_groot_5ha) > 0) {
    r_5ha <- terra::mask(basis_laag, r_patches %in% ids_groot_5ha)
  }
  
  # Kleine patches (1-5 ha) analyseren met Geheugen-veilige Bounding Box / Buffer
  if (length(ids_klein_1to5ha) > 0) {
    message(paste("-> Metapopulatie netwerk berekenen voor", length(ids_klein_1to5ha), "kleine patches..."))
    
    r_klein <- terra::ifel(r_patches %in% ids_klein_1to5ha, r_patches, NA)
    df_pts  <- terra::as.data.frame(r_klein, xy = TRUE)
    rm(r_klein)
    gc()
    
    colnames(df_pts) <- c("x", "y", "ID")
    dt_cent <- as.data.table(df_pts)[, .(x = mean(x), y = mean(y)), by = ID]
    rm(df_pts)
    gc()
    
    if (nrow(dt_cent) > 0) {
      sf_cent <- st_as_sf(dt_cent, coords = c("x", "y"), crs = crs(template_Vlaanderen))
      
      # SNELLE EN GEHEUGEN-ZUINIGE BUUR-TELLING (Voorkomt crash bij duizenden punten)
      buffers      <- st_buffer(sf_cent, dist = 500)
      intersecties <- st_intersects(buffers, sf_cent)
      
      dt_cent$buren_count <- lengths(intersecties)
      rm(buffers, intersecties, sf_cent)
      gc()
      
      valid_ids <- dt_cent[buren_count >= 6, ID]
      
      if (length(valid_ids) > 0) {
        r_meta <- terra::mask(basis_laag, r_patches %in% valid_ids)
      }
    }
  }
  
  w_clean <- terra::ifel(is.na(r_5ha), 0, r_5ha)
  m_clean <- terra::ifel(is.na(r_meta), 0, r_meta)
  rm(r_5ha, r_meta)
  gc()
  
  som_cl  <- terra::clamp(w_clean + m_clean, upper = 1.0)
  rm(w_clean, m_clean)
  gc()
  
  finaal_hab <- terra::ifel(som_cl > 0, 1, NA)
  return(list(habitat = finaal_hab, ids = r_patches))
}

# ==============================================================================
# 7 & 8. UITVOEREN VAN BEIDE SCENARIO'S
# ==============================================================================
res_strikt    <- run_volledige_pijplijn(veldkrekel_basis_strikt, "Strikt Optimaal (a,b,c)")
res_inclusief <- run_volledige_pijplijn(veldkrekel_basis_inclusief, "Inclusief Suboptimaal/NA (d, OB, V, NA)")

# ==============================================================================
# BEREKENING VAN DE GECOMBINEERDE VERSCHILKAART (0/1/2) - WATERDICHT
# ==============================================================================
message("\n-> Verschilkaart herberekenen met expliciete logica...")

# Haal de zuivere binaire maskers op (1 = aanwezig, NA = afwezig)
bin_strikt    <- terra::ifel(!is.na(res_strikt$habitat) & res_strikt$habitat > 0, 1, 0)
bin_inclusief <- terra::ifel(!is.na(res_inclusief$habitat) & res_inclusief$habitat > 0, 1, 0)

# Expliciete Bepaling:
# - Als aanwezig in STRIKT én INCLUSIEF -> 2 (Zeker / Optimaal droge bodem)
# - Als ALLEEN aanwezig in INCLUSIEF     -> 1 (Potentieel / NA / OB / Kustduin zonder drainage)
# - Anders                                -> NA (Geen leefgebied)

r_gradatie_finaal <- terra::ifel(
  bin_strikt == 1 & bin_inclusief == 1, 2,
  terra::ifel(bin_strikt == 0 & bin_inclusief == 1, 1, NA)
)

# ==============================================================================
# 9. DYNAMISCHE EXPORT MAKEN
# ==============================================================================
message("-> Start geformatteerde export voor Vlaanderen...")

base_dir <- here::here("data/output/Vlaanderen/Rasters_Soorten", scenario_naam)

folders <- list(
  id_raster = file.path(base_dir, "00_ID_Rasters"),
  strikt    = file.path(base_dir, "01_Strikt_Optimaal"),
  inclusief = file.path(base_dir, "02_Inclusief_Suboptimaal_NA"),
  gradatie  = file.path(base_dir, "03_Verschilkaart_Gradatie")
)
purrr::walk(folders, ~if (!dir.exists(.x)) dir.create(.x, showWarnings = FALSE, recursive = TRUE))

terra::writeRaster(res_strikt$habitat, file.path(folders$strikt, paste0("Habitat_Strikt_", soort, ".tif")),
                   overwrite = TRUE, gdal = c("COMPRESS=LZW"), datatype = "INT1U", NAflag = 255)

terra::writeRaster(res_inclusief$habitat, file.path(folders$inclusief, paste0("Habitat_Inclusief_NA_", soort, ".tif")),
                   overwrite = TRUE, gdal = c("COMPRESS=LZW"), datatype = "INT1U", NAflag = 255)

terra::writeRaster(r_gradatie_finaal, file.path(folders$gradatie, paste0("Habitat_Gradatie_Verschil_", soort, ".tif")),
                   overwrite = TRUE, gdal = c("COMPRESS=LZW"), datatype = "INT1U", NAflag = 0)

terra::writeRaster(res_inclusief$ids, file.path(folders$id_raster, paste0("ID_Netwerken_", soort, ".tif")),
                   overwrite = TRUE, gdal = c("COMPRESS=LZW"), datatype = "INT4U", NAflag = 0)

message("🏁 BEIDE SCENARIO'S EN VERSCHILKAART SUCCESVOL GEËXPORTEERD!")
