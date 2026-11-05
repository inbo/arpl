# ==============================================================================
# run_pipeline.R - CENTRALE PIPELINE VOOR ALLE SOORTEN
# ==============================================================================

library(here)

# 1. Modules inladen
source(here::here("src", "Functies", "00_global_setup.R"))
source(here::here("src", "Functies", "hulpfuncties.R"))
source(here::here("src", "Functies", "calc_ha_exact.R"))
source(here::here("src", "Functies", "cluster_filter_universeel.R"))
source(here::here("src", "Functies", "extract_bwk_laag.R"))
source(here::here("src", "Functies", "haal_drainage_masker.R"))
source(here::here("src", "Functies", "verwerk_soort.R"))
source(here::here("src", "Functies", "vul_raster_vanuit_parquet.R"))
#source(here::here("src", "Functies", "bereken_arpl.R"))
source(here::here("src", "Functies", "bouw_biotoop_mozaiek_met_substraat.R"))
source(here::here("src", "Functies", "bouw_biotoop_kern_met_overgang.R"))
source(here::here("src", "Functies", "bouw_biotoop_universeel.R"))
source(here::here("src", "Functies", "bouw_biotoop_weidevogel_matrix.R"))
source(here::here("src", "Functies", "bouw_biotoop_oever_gradient.R"))
source(here::here("src", "Functies", "bouw_biotoop_vallei_grasland_mozaiek.R"))
source(here::here("src", "Functies", "bouw_biotoop_open_grasland_bosrandvrij.R"))
source(here::here("src", "Functies", "bouw_biotoop_open_inundatiemoeras.R"))
source(here::here("src", "Functies", "bouw_biotoop_rietmoeras_ecotoon_complex.R"))
source(here::here("src", "Functies", "bouw_biotoop_rietmoeras_waterkern.R"))
source(here::here("src", "Functies", "bouw_biotoop_kleinschalig_cultuurland_driehoek.R"))
source(here::here("src", "Functies", "bouw_biotoop_kolonievogel_kust_polder.R"))
source(here::here("src", "Functies", "bouw_biotoop_open_zand_stuifduin.R"))
source(here::here("src", "Functies", "bouw_biotoop_kust_stern.R"))
source(here::here("src", "Functies", "bouw_biotoop_kuststeltloper_slik.R"))
source(here::here("src", "Functies", "bouw_biotoop_kleinschalig_struweel_loofbos.R"))
source(here::here("src", "Functies", "bouw_biotoop_semi_aquatisch_oeverzoogdier.R"))
source(here::here("src", "Functies", "bouw_biotoop_otter_rivierlandschap.R"))
source(here::here("src", "Functies", "bouw_biotoop_waterplant_hydrograaf.R"))
source(here::here("src", "Functies", "bouw_biotoop_natte_metapopulatie_sprinkhaan.R"))
source(here::here("src", "Functies", "bouw_biotoop_saproxyl_hellingbos.R"))
source(here::here("src", "Functies", "bouw_biotoop_saproxyl_kle_park.R"))
source(here::here("src", "Functies", "bouw_biotoop_saproxyl_alluviaal_bos.R"))
source(here::here("src", "Functies", "bouw_biotoop_stromend_water_libel.R"))
source(here::here("src", "Functies", "bouw_biotoop_verlandingsplas_libel.R"))
source(here::here("src", "Functies", "bouw_biotoop_speerwaterjuffer.R"))
source(here::here("src", "Functies", "bouw_biotoop_venglazenmaker.R"))
source(here::here("src", "Functies", "bouw_biotoop_venwitsnuitlibel.R"))
source(here::here("src", "Functies", "bouw_biotoop_zwarteheidelibel.R"))


# ==============================================================================
# SUBSET TEST SCHAKELAAR VOOR SNELLE RUNS
# ==============================================================================

prepareer_test_grid <- function(master_grid_full, studiegebied_input) {
  message("🔍 [TESTMODE] MasterGrid bijsnijden tot testgebied...")
  
  if (inherits(studiegebied_input, "SpatExtent")) {
    e_sub <- studiegebied_input
  } else if (is.character(studiegebied_input) && file.exists(studiegebied_input)) {
    v_sub <- terra::vect(studiegebied_input)
    v_sub_proj <- terra::project(v_sub, terra::crs(master_grid_full))
    # 2 km buffer rond testgebied zodat rand-clusteringen correct gevormd worden
    e_sub <- terra::ext(terra::buffer(v_sub_proj, 2000))
  } else if (inherits(studiegebied_input, "SpatVector")) {
    v_sub_proj <- terra::project(studiegebied_input, terra::crs(master_grid_full))
    e_sub <- terra::ext(terra::buffer(v_sub_proj, 2000))
  } else {
    stop("❌ Ongeldige invoer voor testgebied. Geef een shapefile-pad of SpatExtent.")
  }
  
  master_grid_test <- terra::crop(master_grid_full, e_sub, snap = "near")
  
  message(sprintf("   -> Testgrid afmetingen: %d kolommen x %d rijen (%.2f ha)", 
                  terra::ncol(master_grid_test), terra::nrow(master_grid_test), 
                  (terra::ncell(master_grid_test) * 0.01)))
  
  return(master_grid_test)
}

# ------------------------------------------------------------------------------
# 2. Globale data inladen en MasterGrid voorbereiden
# ------------------------------------------------------------------------------
message("Inladen template & databronnen...")
master_grid_full <- terra::rast(PATHS$master_grid)[[1]]

# ==============================================================================
# ⚙️ CONFIGURATIE VOOR DEZE RUN
# ==============================================================================
# Zet hier het pad naar je test-shapefile (of een SpatExtent).
# Om voor HEEL VLAANDEREN te draaien: zet TEST_GEBIED <- NULL

TEST_GEBIED <- here::here("data", "input", "Mechelse_Heide.shp")
#TEST_GEBIED <- NULL

# Overschrijven forceren (TRUE) of reeds berekende soorten overslaan (FALSE)
FORCE_OVERWRITE <- T

# Selectie van soorten die je wil testen
soorten_selectie <- c("gevlektewitsnuitlibel") 
# ==============================================================================

if (!is.null(TEST_GEBIED)) {
  # 1. Verklein het mastergrid tot het testgebied
  master_grid <- prepareer_test_grid(master_grid_full, TEST_GEBIED)
  master_grid <- master_grid * NA
  
  # 2. Isoleer testcache zodat productieve cachebestanden niet vervuild raken
  PATHS$dir_cache <- here::here("data", "cache_test")
  dir.create(PATHS$dir_cache, showWarnings = FALSE, recursive = TRUE)
} else {
  message("🌍 Volledig Vlaanderen-grid actief.")
  master_grid <- master_grid_full * NA
}
rm(master_grid_full)
gc()

ds_parquet    <- arrow::open_dataset(PATHS$bwk_parquet)
df_bwk_codes  <- data.table::fread(PATHS$bwk_codes_csv, data.table = FALSE)
df_params     <- readxl::read_excel(PATHS$params_excel)

alle_soorten <- unique(df_params$Soort[!is.na(df_params$Soort)])
message(paste("Totaal aantal beschikbare soorten in database:", length(alle_soorten)))
message(paste("Te verwerken soorten in deze run:", paste(soorten_selectie, collapse = ", ")))

# ------------------------------------------------------------------------------
# 3. Batch loop met foutopvang en geheugenbeheer
# ------------------------------------------------------------------------------
resultaten_overzicht <- list()

for (i in seq_along(soorten_selectie)) {
  huidige_soort <- soorten_selectie[i]
  message(sprintf("\n[%d/%d] ==============================================", i, length(soorten_selectie)))
  
  status <- tryCatch({
    verwerk_soort(
      soort_naam      = huidige_soort,
      params_soort    = df_params %>% dplyr::filter(Soort == huidige_soort),
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
  
  # Grondige geheugenopschoning na elke soort
  gc()
  suppressWarnings(terra::tmpFiles(remove = TRUE))
}

message("\n==============================================")
message("🏁 Batch run voltooid!")
table(unlist(resultaten_overzicht))
