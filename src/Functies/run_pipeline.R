# ==============================================================================
# run_pipeline.R - CENTRALE PIPELINE VOOR ALLE SOORTEN (PRODUCTIEVERSIE)
# ==============================================================================

library(here)

# 0. Stabiele temp-directory voor terra opzetten
dir_terra_tmp <- file.path(here::here("data", "cache"), "terra_tmp")
dir.create(dir_terra_tmp, showWarnings = FALSE, recursive = TRUE)

terra::terraOptions(
  tempdir = dir_terra_tmp,
  memfrac = 0.7,
  todisk  = TRUE,
  verbose = FALSE
)

# ------------------------------------------------------------------------------
# 1. Modules inladen
# ------------------------------------------------------------------------------
source(here::here("src", "Functies", "00_global_setup.R"))
#source(here::here("src", "Functies", "bereken_arpl"))
source(here::here("src", "Functies", "bouw_biotoop_universeel.R"))
source(here::here("src", "Functies", "calc_ha_exact.R"))
source(here::here("src", "Functies", "cluster_filter_universeel.R"))
source(here::here("src", "Functies", "extract_bwk_laag.R"))
source(here::here("src", "Functies", "haal_categorisch_masker.R"))
source(here::here("src", "Functies", "hulp_morfologische_opening.R"))
source(here::here("src", "Functies", "hulp_wederzijdse_koppeling.R"))
source(here::here("src", "Functies", "hulpfuncties.R"))
source(here::here("src", "Functies", "verwerk_soort.R"))
source(here::here("src", "Functies", "vul_raster_vanuit_parquet.R"))


# ==============================================================================
# SUBSET TEST SCHAKELAAR VOOR SNELLE RUNS (SHAPEFILE, EXTENT OF XY-COÖRDINAAT)
# ==============================================================================
prepareer_test_grid <- function(master_grid_full, studiegebied_input, buffer_m = 2000) {
  message(sprintf("🔍 [TESTMODE] MasterGrid bijsnijden tot testgebied (buffer: %d m)...", as.integer(buffer_m)))
  
  if (inherits(studiegebied_input, "SpatExtent")) {
    e_sub <- terra::extend(studiegebied_input, buffer_m)
    
  } else if (is.numeric(studiegebied_input) && length(studiegebied_input) == 2) {
    # XY-coördinaten c(x, y): detecteer automatisch WGS84 (lon/lat) vs Lambert 72 (EPSG:31370)
    x_val <- studiegebied_input[1]
    y_val <- studiegebied_input[2]
    is_wgs84 <- (abs(x_val) <= 180 && abs(y_val) <= 90)
    bron_crs <- if (is_wgs84) "EPSG:4326" else "EPSG:31370"
    
    message(sprintf("   -> XY-punt gedetecteerd (%.5f, %.5f) in %s...", x_val, y_val, bron_crs))
    pt <- terra::vect(cbind(x_val, y_val), type = "points", crs = bron_crs)
    pt_proj <- terra::project(pt, terra::crs(master_grid_full))
    
    # Gebruik minstens 5.000 m straal als buffer_m 0 zou zijn voor een los punt
    straal <- if (buffer_m > 0) buffer_m else 5000
    e_sub <- terra::ext(terra::buffer(pt_proj, width = straal))
    
  } else if (is.character(studiegebied_input) && file.exists(studiegebied_input)) {
    v_sub <- terra::vect(studiegebied_input)
    v_sub_proj <- terra::project(v_sub, terra::crs(master_grid_full))
    e_sub <- terra::ext(terra::buffer(v_sub_proj, width = buffer_m))
    
  } else if (inherits(studiegebied_input, "SpatVector")) {
    v_sub_proj <- terra::project(studiegebied_input, terra::crs(master_grid_full))
    e_sub <- terra::ext(terra::buffer(v_sub_proj, width = buffer_m))
    
  } else {
    stop("❌ Ongeldige invoer voor testgebied. Gebruik een shapefile-pad, SpatExtent of c(x, y).")
  }
  
  # Begrens de uitsnede altijd tot de maximale randen van Vlaanderen (master_grid_full)
  e_sub <- terra::intersect(e_sub, terra::ext(master_grid_full))
  if (is.null(e_sub)) {
    stop("❌ Het opgegeven testgebied valt volledig buiten het Vlaanderen MasterGrid!")
  }
  
  master_grid_test <- terra::crop(master_grid_full, e_sub, snap = "near")
  
  message(sprintf("   -> Testgrid afmetingen: %d kolommen x %d rijen (%.2f ha)", 
                  terra::ncol(master_grid_test), terra::nrow(master_grid_test), 
                  (terra::ncell(master_grid_test) * 0.01)))
  
  return(master_grid_test)
}


# ==============================================================================
# ⚙️ CONFIGURATIE VOOR DEZE RUN
# ==============================================================================
# Kies één van de onderstaande opties voor TEST_GEBIED:
#   - NULL                                                  -> Heel Vlaanderen
#   - here::here("data", "input", "Turnhouts_Vennegebied.shp") -> Shapefile
#   - c(3.25, 51.33)                                        -> XY in GPS (lon, lat) bv. Oostkust / Zwin
#   - c(72000, 226000)                                      -> XY in Lambert 72 (EPSG:31370)
TEST_GEBIED <- NULL
#TEST_GEBIED  <- here::here("data", "input", "Turnhouts_Vennegebied.shp")
TEST_BUFFER_M <- 0   # Buffer in meters rond shapefile of XY-punt (bv. 2000, 15000, 50000)

# Overschrijven forceren (TRUE) of reeds berekende layers hergebruiken (FALSE)
FORCE_OVERWRITE <- TRUE

# Selectie van soorten die je wil uitvoeren
soorten_selectie <- c("grutto")
# ==============================================================================

# ------------------------------------------------------------------------------
# 2. Globale data inladen en MasterGrid voorbereiden (CRASH-PROOF)
# ------------------------------------------------------------------------------
message("Inladen template & databronnen...")

# Onthoud altijd het vaste pad naar het originele Vlaanderen-grid voor transformaties
PATHS$master_grid_vlaanderen <- PATHS$master_grid

if (!is.null(TEST_GEBIED)) {
  # 1. Uitsnede maken met instelbare buffer
  master_grid_full <- terra::rast(PATHS$master_grid)[[1]]
  master_grid_sub  <- prepareer_test_grid(master_grid_full, TEST_GEBIED, buffer_m = TEST_BUFFER_M)
  rm(master_grid_full)
  
  # 2. Sla het testgrid op als VAST bestand (beschermd tegen tmpFiles-cleanup)
  f_test_grid <- file.path(here::here("data", "cache_test"), "master_grid_test.tif")
  dir.create(dirname(f_test_grid), showWarnings = FALSE, recursive = TRUE)
  master_grid <- terra::writeRaster(master_grid_sub, f_test_grid, overwrite = TRUE, datatype = "INT1U")
  rm(master_grid_sub)
  
  PATHS$dir_cache  <- here::here("data", "cache_test")
  PATHS$dir_output <- here::here("data", "output_test")
} else {
  message("🌍 Volledig Vlaanderen-grid actief (koppeling naar vast bronbestand).")
  # Wijs rechtstreeks naar het permanente bronrasterbestand
  master_grid <- terra::rast(PATHS$master_grid)[[1]]
  PATHS$dir_output <- here::here("data", "output", "Vlaanderen", "Rasters_Soorten")
}

# Zorg dat de scenario-submap leeg is zodat hij direct in Rasters_Soorten schrijft
PATHS$scenario_naam <- NULL

dir.create(PATHS$dir_output, showWarnings = FALSE, recursive = TRUE)
dir.create(PATHS$dir_cache, showWarnings = FALSE, recursive = TRUE)

# Databronnen koppelen
ds_parquet   <- arrow::open_dataset(PATHS$bwk_parquet)
df_bwk_codes <- data.table::fread(PATHS$bwk_codes_csv, data.table = FALSE)
df_params    <- readxl::read_excel(PATHS$params_excel)

alle_soorten <- unique(df_params$Soort[!is.na(df_params$Soort)])
message(paste("Totaal aantal beschikbare soorten in database:", length(alle_soorten)))
message(paste("Te verwerken soorten in deze run:", paste(soorten_selectie, collapse = ", ")))


# ------------------------------------------------------------------------------
# 3. Batch loop met foutopvang en veilige geheugenopschoning
# ------------------------------------------------------------------------------
resultaten_overzicht <- list()

for (i in seq_along(soorten_selectie)) {
  huidige_soort <- tolower(trimws(soorten_selectie[i]))
  message(sprintf("\n[%d/%d] ==============================================", i, length(soorten_selectie)))
  
  status <- tryCatch({
    verwerk_soort(
      soort_naam      = huidige_soort,
      params_soort    = df_params %>% dplyr::filter(tolower(trimws(Soort)) == huidige_soort),
      master_grid     = master_grid,
      ds_parquet      = ds_parquet,
      df_bwk_codes    = df_bwk_codes,
      paths           = PATHS,
      force_overwrite = FORCE_OVERWRITE
    )
    "SUCCESS"
  }, error = function(e) {
    message(paste("❌ FOUT bij soort:", huidige_soort, "\nMelding:", e$message))
    return("FAILED")
  })
  
  resultaten_overzicht[[huidige_soort]] <- status
  
  # Veilige geheugenopschoning (wist géén actieve pointers van master_grid)
  gc()
  suppressWarnings(terra::tmpFiles(current = FALSE, orphan = TRUE, old = TRUE, remove = TRUE))
}

message("\n==============================================")
message("🏁 Batch run voltooid!")
print(table(unlist(resultaten_overzicht)))
