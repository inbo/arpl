# ==============================================================================
# WORK-QUEUE PARALLEL RUNNER VOOR PURE R-SCRIPTS + FINALE MAKER
# AUTEUR: Bert Van Hecke (HPC Optimized & OOM Stabilized)
# ==============================================================================

totaal_start <- Sys.time()

library(here)
library(purrr)
library(readxl)
library(dplyr)
library(readr)
library(future)
library(furrr)
library(tidyterra)
library(future.callr)

# Store absolute project root explicitly to avoid relative path breakage in workers
proj_root <- here::here()

# ------------------------------------------------------------------------------
# 1. INSTELLINGEN & SCENARIO SELECTIE PER GEBIED
# ------------------------------------------------------------------------------
# Capped at 8 workers max to prevent Linux OOM kills on heavy raster calculations.
# Gives each worker process ~8 GB RAM on a standard 64 GB node.
AANTAL_CORES <- 16

SCENARIO_SELECTIE <- list(
  De_Maten              = "DM_Scenario_BWK_2025.rds",
  Heesbossen            = "HB_Scenario_BWK_2025.rds",
  Kalmthoutse_Heide     = "KH_Scenario_BWK_2025.rds",
  Mechelse_Heide        = "MH_Scenario_BWK_2025.rds",
  Turnhouts_Vennegebied = "TV_Scenario_BWK_2025.rds",
  Voerstreek            = "VS_Scenario_BWK_2025.rds"
)

gebieden_info <- list(
  De_Maten              = list(code = "DM", col = "De_Maten",            simpel_script = "Scenario_DM_Leefgebieden_Simpel.R"),
  Heesbossen            = list(code = "HB", col = "Heesbossen",          simpel_script = "Scenario_HB_Leefgebieden_Simpel.R"),
  Kalmthoutse_Heide     = list(code = "KH", col = "Kalmthoutse_Heide",     simpel_script = "Scenario_KH_Leefgebieden_Simpel.R"),
  Mechelse_Heide        = list(code = "MH", col = "Mechelse_Heide",        simpel_script = "Scenario_MH_Leefgebieden_Simpel.R"),
  Turnhouts_Vennegebied = list(code = "TV", col = "Turnhouts_Vennegebied", simpel_script = "Scenario_TV_Leefgebieden_Simpel.R"),
  Voerstreek            = list(code = "VS", col = "Voerstreek",          simpel_script = "Scenario_VS_Leefgebieden_Simpel.R")
)

schoon_naam_op <- function(x) {
  x %>% 
    basename() %>% 
    tolower() %>% 
    gsub(" ", "", .) %>% 
    gsub("^habitat_werkelijke_oppervlaktes_|^habitat_maximale_potentie_|^id_netwerken_|^rapport_|^scenario_", "", .) %>% 
    gsub("^(tv|dm|hb|kh|mh|vs)_", "", .) %>% 
    gsub("_wv\\.tif$|_wv\\.r$|\\.tif$|\\.rds$|\\.html$|\\.r$", "", .) %>% 
    trimws()
}

# ------------------------------------------------------------------------------
# 2. VERZAMEL DYNAMISCH ALLE R-SCRIPTS PER GEBIED (INCL. SLIMME SKIP-CHECK)
# ------------------------------------------------------------------------------
taken_lijst <- list()
excel_pad <- file.path(proj_root, "data/input/Excel_files/Soortenlijst_Maatwerkgebieden_Gefilterd.xlsx")
excel_data <- if(file.exists(excel_pad)) read_excel(excel_pad) else NULL

if(!is.null(excel_data)) {
  colnames(excel_data) <- tolower(colnames(excel_data))
}

for (gb_naam in names(gebieden_info)) {
  info <- gebieden_info[[gb_naam]]
  map_scripts <- file.path(proj_root, "src", gb_naam, "Scripts_Scenario")
  
  if (!dir.exists(map_scripts)) next
  
  gekozen_rds_naam <- SCENARIO_SELECTIE[[gb_naam]]
  scenario_rds      <- file.path(proj_root, "data/input/Scenario_rds", gekozen_rds_naam)
  
  if (!file.exists(scenario_rds)) {
    warning("⚠️ Het opgegeven RDS bestand '", gekozen_rds_naam, "' bestaat niet voor ", gb_naam, "!")
    next
  }
  
  huidig_scenario <- gsub("^.*_Scenario_|^Scenario_|_wv\\.rds$|\\.rds$", "", basename(scenario_rds), ignore.case = TRUE)
  map_werkelijk   <- file.path(proj_root, "data/output", gb_naam, "Rasters_Soorten", huidig_scenario, "02_Werkelijke_Oppervlaktes")
  
  bestaande_soorten_clean <- if (dir.exists(map_werkelijk)) {
    list.files(map_werkelijk, pattern = "\\.tif$", full.names = FALSE) %>% 
      schoon_naam_op() %>% 
      unique()
  } else {
    c()
  }
  
  alle_r_scripts <- list.files(path = map_scripts, pattern = "\\.r$", full.names = TRUE, ignore.case = TRUE)
  soorten_script_pad <- alle_r_scripts[grep(paste0(info$simpel_script, "$"), alle_r_scripts, ignore.case = TRUE)]
  if(length(soorten_script_pad) > 0) soorten_script_pad <- soorten_script_pad[1] else soorten_script_pad <- file.path(map_scripts, info$simpel_script)
  
  unieke_r_scripts <- setdiff(alle_r_scripts, soorten_script_pad)
  
  for (r_script in unieke_r_scripts) {
    soort_uit_script <- schoon_naam_op(r_script)
    if (soort_uit_script %in% bestaande_soorten_clean) next
    
    taken_lijst[[length(taken_lijst) + 1]] <- list(
      Type        = "Uniek_R_Script",
      Gebied      = gb_naam,
      GebiedCode  = info$code,
      ScriptPad   = r_script,
      ScenarioRDS = scenario_rds,
      Soort       = NA
    )
  }
  
  if (!is.null(excel_data) && file.exists(soorten_script_pad)) {
    col_naam <- tolower(info$col)
    
    if (col_naam %in% colnames(excel_data)) {
      naam_kolom <- if("soort" %in% colnames(excel_data)) "soort" else "nederlandse naam"
      
      soorten <- excel_data %>% 
        mutate(
          automatisch_clean = tolower(trimws(as.character(automatisch))),
          gebied_val        = suppressWarnings(as.numeric(.data[[col_naam]]))
        ) %>% 
        filter(automatisch_clean == "ja") %>% 
        filter(gebied_val == 1) %>% 
        pull(!!sym(naam_kolom)) %>% 
        trimws() %>% 
        unique() %>% 
        na.omit()
      
      for (srt in soorten) {
        srt_clean <- schoon_naam_op(srt)
        if (srt_clean %in% bestaande_soorten_clean) next
        
        taken_lijst[[length(taken_lijst) + 1]] <- list(
          Type        = "Simpele_Soort",
          Gebied      = gb_naam,
          GebiedCode  = info$code,
          ScriptPad   = soorten_script_pad,
          ScenarioRDS = scenario_rds,
          Soort       = srt
        )
      }
    }
  }
}

message("==================================================")
message(" WORK-QUEUE KLAARGEZET")
message(" Totaal te verwerken taken na skip-check: ", length(taken_lijst))
message(" Aantal actieve cores op HPC: ", AANTAL_CORES)
message("==================================================")

# ------------------------------------------------------------------------------
# 3. VERWERKINGSFUNCTIE VOOR 1 R-SCRIPT TAAK (OOM & TEMP ISOLATION FIX)
# ------------------------------------------------------------------------------
verwerk_r_taak <- function(taak, p_root, n_cores) {
  options(here.root = p_root)
  
  # Prevent C++ library thread contention across parallel workers
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
  
  # Load packages inside isolated worker process
  suppressPackageStartupMessages({
    library(here)
    library(dplyr)
    library(purrr)
    library(readr)
    library(terra)
    library(sf)
    library(tidyterra)
    library(tidyverse)
  })
  
  # Unique isolated temporary directory on Scratch for each worker
  worker_temp <- file.path(p_root, "data/temp_workers", paste0("gis_worker_", Sys.getpid(), "_", sample(1000:9999, 1)))
  dir.create(worker_temp, showWarnings = FALSE, recursive = TRUE)
  
  # Cap terra RAM usage to 1/12th per worker to guarantee no system OOM crashes
  if (requireNamespace("terra", quietly = TRUE)) {
    terra::terraOptions(threads = 1, tempdir = worker_temp, memfrac = 0.08, verbose = FALSE)
  }
  
  tijd_start <- Sys.time()
  
  # Environment inheritance for child scripts
  run_env <- new.env(parent = globalenv())
  
  huidig_scen_naam <- gsub("^.*_Scenario_|^Scenario_|_wv\\.rds$|\\.rds$", "", basename(taak$ScenarioRDS), ignore.case = TRUE)
  
  run_env$SCENARIO_RDS_PAD  <- taak$ScenarioRDS
  run_env$scenario_rds_path <- taak$ScenarioRDS
  run_env$GEBIED_NAAM       <- taak$Gebied
  run_env$GEBIED_CODE       <- taak$GebiedCode
  run_env$huidig_scenario   <- huidig_scen_naam
  run_env$output_dir        <- file.path(p_root, "data/output", taak$Gebied, "HTML_Rapporten_Scenario", huidig_scen_naam)
  
  if (!is.na(taak$Soort)) {
    run_env$huidige_soort       <- taak$Soort
    run_env$soort_invoer        <- taak$Soort
    run_env$soort               <- taak$Soort
    run_env$soort_naam          <- taak$Soort
    run_env$SOORT               <- taak$Soort
    run_env$huidige_soort_clean <- tolower(gsub(" ", "", taak$Soort))
    item_naam                   <- taak$Soort
  } else {
    item_naam                   <- basename(taak$ScriptPad)
  }
  
  old_wd <- setwd(dirname(taak$ScriptPad))
  on.exit({
    setwd(old_wd)
    unlink(worker_temp, recursive = TRUE, force = TRUE)
  })
  
  tryCatch({
    sys.source(taak$ScriptPad, envir = run_env)
    
    if (requireNamespace("terra", quietly = TRUE)) {
      terra::tmpFiles(current = TRUE, orphan = TRUE, old = TRUE, remove = TRUE)
    }
    gc(verbose = FALSE)
    tijd_eind <- Sys.time()
    
    return(data.frame(
      Gebied = taak$Gebied, Item = item_naam, Type = taak$Type, 
      Status = "SUCCES", Start_Tijd = format(tijd_start, "%H:%M:%S"),
      Duurtijd_Min = round(as.numeric(difftime(tijd_eind, tijd_start, units="mins")), 2),
      Fout = "Geen", stringsAsFactors = FALSE
    ))
    
  }, error = function(e) {
    gc(verbose = FALSE)
    return(data.frame(
      Gebied = taak$Gebied, Item = item_naam, 
      Type = taak$Type, Status = "CRASH", Start_Tijd = format(tijd_start, "%H:%M:%S"),
      Duurtijd_Min = NA, Fout = e$message, stringsAsFactors = FALSE
    ))
  })
}

# ------------------------------------------------------------------------------
# 4. PARALLELLE WORK QUEUE UITVOEREN (FASE 1 - WITH STABILITY FIXES)
# ------------------------------------------------------------------------------
if (length(taken_lijst) > 0) {
  # 'callr' launches clean external R processes that catch crashes gracefully
  plan(callr, workers = AANTAL_CORES)
  
  eind_logboek <- furrr::future_map_dfr(
    taken_lijst, 
    ~verwerk_r_taak(.x, p_root = proj_root, n_cores = AANTAL_CORES), 
    .options = furrr_options(
      packages = c("here", "dplyr", "purrr", "readr", "terra", "sf", "tidyterra", "tidyverse"),
      seed = TRUE,
      chunk_size = 1  # Process 1 item per batch to immediately free RAM
    ),
    .progress = FALSE
  )
  
  plan(sequential)
  
  log_file_path <- file.path(proj_root, paste0("Logboek_R_WorkQueue_Totaal_", format(Sys.time(), "%Y%m%d_%H%M"), ".csv"))
  readr::write_excel_csv(eind_logboek, log_file_path)
  
  cat("\n==================================================\n")
  cat(" FASE 1: ALLE LOSSE R-SCRIPTS ZIJN AFGEROND\n")
  cat(" Totaal nieuwe onderdelen gedraaid : ", nrow(eind_logboek), "\n")
  cat(" Succesvol                           : ", sum(eind_logboek$Status == "SUCCES"), "\n")
  cat(" Gecrasht                            : ", sum(eind_logboek$Status == "CRASH"), "\n")
  cat("==================================================\n\n")
} else {
  cat("\n==================================================\n")
  cat(" FASE 1: Alle rasters voor alle gebieden bestaan al op schijf!\n")
  cat("==================================================\n\n")
}

# ------------------------------------------------------------------------------
# 5. FASE 2: RUN FINALE MAKER SCRIPT (ARPL_Actueel_Maker.R)
# ------------------------------------------------------------------------------
maker_script <- list.files(
  path = file.path(proj_root, "src"), 
  pattern = "^ARPL_Actueel_Maker\\.r$", 
  full.names = TRUE, 
  recursive = TRUE, 
  ignore.case = TRUE
)[1]

if (!is.null(maker_script) && file.exists(maker_script)) {
  cat("==================================================\n")
  cat(" STARTEN FASE 2: FINALE ARPL ACTUEEL MAKER SCRIPT\n")
  cat("==================================================\n")
  
  tijd_maker_start <- Sys.time()
  
  tryCatch({
    maker_env <- new.env(parent = globalenv())
    maker_env$SCENARIO_SELECTIE <- SCENARIO_SELECTIE
    
    sys.source(maker_script, envir = maker_env)
    
    tijd_maker_eind <- Sys.time()
    duur_maker <- round(as.numeric(difftime(tijd_maker_eind, tijd_maker_start, units = "mins")), 2)
    cat("\n✅ ARPL_Actueel_Maker.R succesvol afgerond in", duur_maker, "minuten!\n")
    
  }, error = function(e) {
    cat("\n❌ FOUT bij uitvoeren van ARPL_Actueel_Maker.R:\n", e$message, "\n")
  })
  
} else {
  warning("⚠️ Het script 'ARPL_Actueel_Maker.R' kon niet automatisch gevonden worden in src/.")
}

totaal_eind <- Sys.time()
totale_duurtijd_uur <- round(as.numeric(difftime(totaal_eind, totaal_start, units = "hours")), 2)

cat("\n==================================================\n")
cat(" 🎉 VOLLEDIGE PIPELINE GESLAAGD IN: ", totale_duurtijd_uur, " UUR\n")
cat("==================================================\n")
