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
source(here::here("src", "Functies", "bouw_biotoop_complex.R"))
#source(here::here("src", "Functies", "bereken_arpl.R"))
source(here::here("src", "Functies", "bouw_biotoop_ecotoon.R"))
source(here::here("src", "Functies", "bouw_biotoop_mozaiek_met_substraat.R"))
source(here::here("src", "Functies", "bouw_biotoop_kern_met_overgang.R"))
source(here::here("src", "Functies", "bouw_biotoop_bosrand_gradient.R"))
source(here::here("src", "Functies", "bouw_biotoop_actieradius_kern.R"))
source(here::here("src", "Functies", "bouw_biotoop_kle_mozaiek_met_storing.R"))
source(here::here("src", "Functies", "bouw_biotoop_weidevogel_matrix.R"))
source(here::here("src", "Functies", "bouw_biotoop_oever_gradient.R"))
source(here::here("src", "Functies", "bouw_biotoop_commuter_watermoeras.R"))
source(here::here("src", "Functies", "bouw_biotoop_vallei_grasland_mozaiek.R"))
source(here::here("src", "Functies", "bouw_biotoop_open_grasland_bosrandvrij.R"))
source(here::here("src", "Functies", "bouw_biotoop_drasland_steltloper.R"))
source(here::here("src", "Functies", "bouw_biotoop_bosrand_lengte.R"))
source(here::here("src", "Functies", "bouw_biotoop_open_inundatiemoeras.R"))
source(here::here("src", "Functies", "bouw_biotoop_rietmoeras_ecotoon_complex.R"))
source(here::here("src", "Functies", "bouw_biotoop_bos_landschap_actieradius.R"))
source(here::here("src", "Functies", "bouw_biotoop_rietmoeras_waterkern.R"))
source(here::here("src", "Functies", "bouw_biotoop_grootschalige_weidevogel.R"))
source(here::here("src", "Functies", "bouw_biotoop_kleinschalig_cultuurland_driehoek.R"))
source(here::here("src", "Functies", "bouw_biotoop_bos_metapopulatie_netwerk.R"))
source(here::here("src", "Functies", "bouw_biotoop_kolonievogel_kust_polder.R"))
source(here::here("src", "Functies", "bouw_biotoop_open_zand_stuifduin.R"))
source(here::here("src", "Functies", "bouw_biotoop_kust_stern.R"))
source(here::here("src", "Functies", "bouw_biotoop_grootschalig_akkercomplex.R"))
source(here::here("src", "Functies", "bouw_biotoop_kuststeltloper_slik.R"))
source(here::here("src", "Functies", "bouw_biotoop_dynamisch_kustgrasland.R"))
source(here::here("src", "Functies", "bouw_biotoop_kleinschalig_struweel_loofbos.R"))
source(here::here("src", "Functies", "bouw_biotoop_semi_aquatisch_oeverzoogdier.R"))
source(here::here("src", "Functies", "bouw_biotoop_otter_rivierlandschap.R"))

# 2. Globale data eenmalig in het geheugen houden
message("Inladen template & databronnen...")
master_grid   <- terra::rast(PATHS$master_grid)[[1]]
master_grid <- master_grid * NA

ds_parquet    <- arrow::open_dataset(PATHS$bwk_parquet)
df_bwk_codes  <- data.table::fread(PATHS$bwk_codes_csv, data.table = FALSE)
df_params     <- readxl::read_excel(PATHS$params_excel)

# Unieke soorten ophalen
alle_soorten <- unique(df_params$Soort[!is.na(df_params$Soort)])
message(paste("Totaal aantal te verwerken soorten:", length(alle_soorten)))

# 3. Batch loop met foutopvang en geheugenbeheer
# Stel in of je alles opnieuw wilt forceren of enkel ontbrekende soorten wilt berekenen
FORCE_OVERWRITE <- F

# Optioneel: selectie om eerst een subset te draaien (bijv. alle_soorten[1:5])
soorten_selectie <- c("grutto", "kwartelkoning", "watersnip") #alle_soorten 

resultaten_overzicht <- list()

for (i in seq_along(soorten_selectie)) {
  huidige_soort <- soorten_selectie[i]
  message(sprintf("\n[%d/%d] ==============================================", i, length(soorten_selectie)))
  
  status <- tryCatch({
    verwerk_soort(
      soort_naam      = huidige_soort,
      params_soort    = df_params %>% dplyr::filter(Soort == huidige_soort),
      master_grid   = master_grid,
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
