# ==============================================================================
# WORK-QUEUE PARALLEL RUNNER VOOR PURE R-SCRIPTS + FINALE MAKER (6 CORES)
# AUTEUR: Bert Van Hecke
# ==============================================================================

library(here)
library(purrr)
library(readxl)
library(dplyr)
library(readr)
library(future)
library(furrr)

# ------------------------------------------------------------------------------
# 1. INSTELLINGEN & SCENARIO SELECTIE PER GEBIED
# ------------------------------------------------------------------------------
AANTAL_CORES <- 6  # Ingesteld op 6 cores voor HPC/supercomputer

SCENARIO_SELECTIE <- list(
  De_Maten              = "DM_Scenario_BWK_2025.rds",
  Heesbossen            = "HB_Scenario_BWK_2025.rds",
  Kalmthoutse_Heide     = "KH_Scenario_BWK_2025.rds",
  Mechelse_Heide        = "MH_Scenario_BWK_2025.rds",
  Turnhouts_Vennegebied = "TV_Scenario_BWK_2025.rds",
  Voerstreek            = "VS_Scenario_BWK_2025.rds"
)

gebieden_info <- list(
  De_Maten              = list(code = "DM", col = "De_Maten",              simpel_script = "Scenario_DM_Leefgebieden_Simpel.R"),
  Heesbossen            = list(code = "HB", col = "Heesbossen",            simpel_script = "Scenario_HB_Leefgebieden_Simpel.R"),
  Kalmthoutse_Heide     = list(code = "KH", col = "Kalmthoutse_Heide",     simpel_script = "Scenario_KH_Leefgebieden_Simpel.R"),
  Mechelse_Heide        = list(code = "MH", col = "Mechelse_Heide",        simpel_script = "Scenario_MH_Leefgebieden_Simpel.R"),
  Turnhouts_Vennegebied = list(code = "TV", col = "Turnhouts_Vennegebied", simpel_script = "Scenario_TV_Leefgebieden_Simpel.R"),
  Voerstreek            = list(code = "VS", col = "Voerstreek",            simpel_script = "Scenario_VS_Leefgebieden_Simpel.R")
)

# Hulpfunctie om zowel scriptnamen, soorten als rasters te herleiden tot de schone kernnaam (zonder spaties)
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
excel_pad <- here("data/input/Excel_files/Soortenlijst_Maatwerkgebieden_Gefilterd.xlsx")
excel_data <- if(file.exists(excel_pad)) read_excel(excel_pad) else NULL

if(!is.null(excel_data)) {
  colnames(excel_data) <- tolower(colnames(excel_data))
}

for (gb_naam in names(gebieden_info)) {
  info <- gebieden_info[[gb_naam]]
  map_scripts <- here("src", gb_naam, "Scripts_Scenario")
  
  if (!dir.exists(map_scripts)) next
  
  gekozen_rds_naam <- SCENARIO_SELECTIE[[gb_naam]]
  scenario_rds      <- here("data/input/Scenario_rds", gekozen_rds_naam)
  
  if (!file.exists(scenario_rds)) {
    warning("⚠️ Het opgegeven RDS bestand '", gekozen_rds_naam, "' bestaat niet voor ", gb_naam, "!")
    next
  }
  
  huidig_scenario <- gsub("^.*_Scenario_|^Scenario_|_wv\\.rds$|\\.rds$", "", basename(scenario_rds), ignore.case = TRUE)
  
  # Detecteer al verwerkte TIF-rasters voor dit specifieke gebied & scenario
  map_werkelijk <- here("data/output", gb_naam, "Rasters_Soorten", huidig_scenario, "02_Werkelijke_Oppervlaktes")
  
  bestaande_soorten_clean <- if (dir.exists(map_werkelijk)) {
    list.files(map_werkelijk, pattern = "\\.tif$", full.names = FALSE) %>% 
      schoon_naam_op() %>% 
      unique()
  } else {
    c()
  }
  
  # Zoek alle .R scripts op
  alle_r_scripts <- list.files(path = map_scripts, pattern = "\\.r$", full.names = TRUE, ignore.case = TRUE)
  soorten_script_pad <- alle_r_scripts[grep(paste0(info$simpel_script, "$"), alle_r_scripts, ignore.case = TRUE)]
  if(length(soorten_script_pad) > 0) soorten_script_pad <- soorten_script_pad[1] else soorten_script_pad <- file.path(map_scripts, info$simpel_script)
  
  unieke_r_scripts <- setdiff(alle_r_scripts, soorten_script_pad)
  
  # A. Unieke R-scripts filteren met skip-check
  for (r_script in unieke_r_scripts) {
    soort_uit_script <- schoon_naam_op(r_script)
    
    # Sla over als het raster al bestaat op schijf
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
  
  # B. Simpele Soorten uit Excel filteren met skip-check
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
        
        # Sla over als het raster al bestaat op schijf
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
# 3. VERWERKINGSFUNCTIE VOOR 1 R-SCRIPT TAAK
# ------------------------------------------------------------------------------
verwerk_r_taak <- function(taak, proj_root) {
  options(here.root = proj_root)
  library(here)
  
  # Maak een unieke tempdir per worker (voorkomt race-conditions tussen cores)
  worker_temp <- tempfile(pattern = "gis_temp_")
  dir.create(worker_temp, showWarnings = FALSE, recursive = TRUE)
  
  if (requireNamespace("terra", quietly = TRUE)) {
    # memfrac op 0.15 zodat 6 parallelle workers nooit het HPC geheugen overbelasten
    terra::terraOptions(tempdir = worker_temp, memfrac = 0.15, verbose = FALSE)
  }
  
  tijd_start <- Sys.time()
  
  # Geïsoleerde omgeving per taak
  run_env <- new.env(parent = globalenv())
  
  huidig_scen_naam <- gsub("^.*_Scenario_|^Scenario_|_wv\\.rds$|\\.rds$", "", basename(taak$ScenarioRDS), ignore.case = TRUE)
  
  run_env$SCENARIO_RDS_PAD  <- taak$ScenarioRDS
  run_env$scenario_rds_path <- taak$ScenarioRDS
  run_env$GEBIED_NAAM       <- taak$Gebied
  run_env$GEBIED_CODE       <- taak$GebiedCode
  run_env$huidig_scenario   <- huidig_scen_naam
  run_env$output_dir        <- here("data/output", taak$Gebied, "HTML_Rapporten_Scenario", huidig_scen_naam)
  
  if (!is.na(taak$Soort)) {
    # Alle mogelijke naam-varianten meegeven voor het simpele script
    run_env$huidige_soort       <- taak$Soort
    run_env$soort_invoer        <- taak$Soort
    run_env$soort              <- taak$Soort
    run_env$soort_naam         <- taak$Soort
    run_env$SOORT              <- taak$Soort
    run_env$huidige_soort_clean <- tolower(gsub(" ", "", taak$Soort))
    item_naam                   <- taak$Soort
  } else {
    item_naam                   <- basename(taak$ScriptPad)
  }
  
  old_wd <- setwd(dirname(taak$ScriptPad))
  on.exit({
    setwd(old_wd)
    unlink(worker_temp, recursive = TRUE)
  })
  
  tryCatch({
    sys.source(taak$ScriptPad, envir = run_env)
    
    terra::tmpFiles(current = TRUE, orphan = TRUE, old = TRUE, remove = TRUE)
    gc(verbose = FALSE)
    tijd_eind <- Sys.time()
    
    return(data.frame(
      Gebied = taak$Gebied, Item = item_naam, Type = taak$Type, 
      Status = "SUCCES", Duurtijd_Min = round(as.numeric(difftime(tijd_eind, tijd_start, units="mins")), 2),
      Fout = "Geen", stringsAsFactors = FALSE
    ))
    
  }, error = function(e) {
    gc(verbose = FALSE)
    return(data.frame(
      Gebied = taak$Gebied, Item = item_naam, 
      Type = taak$Type, Status = "CRASH", Duurtijd_Min = NA, Fout = e$message, stringsAsFactors = FALSE
    ))
  })
}

# ------------------------------------------------------------------------------
# 4. PARALLELLE WORK QUEUE UITVOEREN (FASE 1)
# ------------------------------------------------------------------------------
if (length(taken_lijst) > 0) {
  plan(multisession, workers = AANTAL_CORES)
  
  totaal_start <- Sys.time()
  proj_root <- here()
  
  eind_logboek <- furrr::future_map_dfr(
    taken_lijst, 
    ~verwerk_r_taak(.x, proj_root = proj_root), 
    .options = furrr_options(packages = c("here", "dplyr", "readr", "terra")),
    .progress = TRUE
  )
  
  plan(sequential)
  
  log_file_path <- here(paste0("Logboek_R_WorkQueue_Totaal_", format(Sys.time(), "%Y%m%d_%H%M"), ".csv"))
  readr::write_excel_csv(eind_logboek, log_file_path)
  
  cat("\n==================================================\n")
  cat(" FASE 1: ALLE LOSSE R-SCRIPTS ZIJN AFGEROND\n")
  cat(" Totaal nieuwe onderdelen gedraaid : ", nrow(eind_logboek), "\n")
  cat(" Succesvol                          : ", sum(eind_logboek$Status == "SUCCES"), "\n")
  cat(" Gecrasht                           : ", sum(eind_logboek$Status == "CRASH"), "\n")
  cat("==================================================\n\n")
} else {
  totaal_start <- Sys.time()
  cat("\n==================================================\n")
  cat(" FASE 1: Alle rasters voor alle gebieden bestaan al op schijf!\n")
  cat("==================================================\n\n")
}

# ------------------------------------------------------------------------------
# 5. FASE 2: RUN FINALE MAKER SCRIPT (ARPL_Actueel_Maker.R)
# ------------------------------------------------------------------------------
maker_script <- list.files(
  path = here("src"), 
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
