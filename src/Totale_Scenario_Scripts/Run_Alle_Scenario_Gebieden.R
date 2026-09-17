# ==============================================================================
# WORK-QUEUE PARALLEL RUNNER VOOR PURE R-SCRIPTS + FINALE MAKER
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
AANTAL_CORES <- 4 

# STEL HIER PER GEBIED HET GEWENSTE SCENARIO RDS-BESTAND IN:
# (Geef de exacte bestandsnaam op die in 'data/input/Scenario_rds/' staat)
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

# ------------------------------------------------------------------------------
# 2. VERZAMEL DYNAMISCH ALLE R-SCRIPTS PER GEBIED
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
  
  # A. Koppel exact het gekozen Scenario RDS-bestand uit SCENARIO_SELECTIE
  gekozen_rds_naam <- SCENARIO_SELECTIE[[gb_naam]]
  scenario_rds      <- here("data/input/Scenario_rds", gekozen_rds_naam)
  
  if (!file.exists(scenario_rds)) {
    warning("⚠️ Het opgegeven RDS bestand '", gekozen_rds_naam, "' bestaat niet in data/input/Scenario_rds/ voor ", gb_naam, "!")
    next
  }
  
  # B. Zoek alle Unieke .R scripts
  alle_r_scripts <- list.files(path = map_scripts, pattern = "\\.R$", full.names = TRUE)
  soorten_script_pad <- file.path(map_scripts, info$simpel_script)
  unieke_r_scripts <- setdiff(alle_r_scripts, soorten_script_pad)
  
  for (r_script in unieke_r_scripts) {
    taken_lijst[[length(taken_lijst) + 1]] <- list(
      Type = "Uniek_R_Script",
      Gebied = gb_naam,
      GebiedCode = info$code,
      ScriptPad = r_script,
      ScenarioRDS = scenario_rds,
      Soort = NA
    )
  }
  
  # C. Zoek de Simpele Soorten uit de Excel
  if (!is.null(excel_data) && file.exists(soorten_script_pad)) {
    col_naam <- tolower(info$col) # tolower gebruikt omdat colnames(excel_data) tolower is gemaakt
    
    if (col_naam %in% colnames(excel_data)) {
      soorten <- excel_data %>% 
        filter(automatisch == "ja") %>% 
        filter(.data[[col_naam]] == 1) %>% 
        pull(soort) %>% unique() %>% na.omit()
      
      for (srt in soorten) {
        taken_lijst[[length(taken_lijst) + 1]] <- list(
          Type = "Simpele_Soort",
          Gebied = gb_naam,
          GebiedCode = info$code,
          ScriptPad = soorten_script_pad,
          ScenarioRDS = scenario_rds,
          Soort = srt
        )
      }
    }
  }
}

# ------------------------------------------------------------------------------
# 3. VERWERKINGSFUNCTIE VOOR 1 R-SCRIPT TAAK
# ------------------------------------------------------------------------------
verwerk_r_taak <- function(taak, proj_root) {
  options(here.root = proj_root)
  library(here)
  
  tijd_start <- Sys.time()
  
  # Maak een schone, geïsoleerde R-omgeving aan per script
  run_env <- new.env(parent = globalenv())
  
  # Geef de vereiste variabelen/parameters door aan de omgeving van het script
  run_env$SCENARIO_RDS_PAD <- taak$ScenarioRDS
  run_env$GEBIED_NAAM      <- taak$Gebied
  run_env$GEBIED_CODE      <- taak$GebiedCode
  
  if (!is.na(taak$Soort)) {
    run_env$huidige_soort  <- taak$Soort
    item_naam              <- taak$Soort
  } else {
    item_naam              <- basename(taak$ScriptPad)
  }
  
  # Zet de werkmap tijdelijk om naar de map van het script
  old_wd <- setwd(dirname(taak$ScriptPad))
  on.exit(setwd(old_wd))
  
  tryCatch({
    # Snel en direct het .R script uitvoeren binnen de schone omgeving
    sys.source(taak$ScriptPad, envir = run_env)
    
    gc(verbose = FALSE) # RAM direct vrijgeven
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
plan(multisession, workers = AANTAL_CORES)

totaal_start <- Sys.time()

proj_root <- here()
eind_logboek <- future_map_dfr(
  taken_lijst, 
  ~verwerk_r_taak(.x, proj_root = proj_root), 
  .options = furrr_options(packages = c("here", "dplyr", "readr")),
  .progress = TRUE
)

# Sluit het parallelle cluster af om geheugen vrij te maken voor de finale maker
plan(sequential)

# Sla het overkoepelende logboek van de losse scripts op
log_file_path <- here(paste0("Logboek_R_WorkQueue_Totaal_", format(Sys.time(), "%Y%m%d_%H%M"), ".csv"))
readr::write_excel_csv(eind_logboek, log_file_path)

cat("\n==================================================\n")
cat(" FASE 1: ALLE LOSSE R-SCRIPTS ZIJN AFGEROND\n")
cat(" Totaal onderdelen : ", nrow(eind_logboek), "\n")
cat(" Succesvol         : ", sum(eind_logboek$Status == "SUCCES"), "\n")
cat(" Gecrasht          : ", sum(eind_logboek$Status == "CRASH"), "\n")
cat("==================================================\n\n")

# ------------------------------------------------------------------------------
# 5. FASE 2: RUN FINALE MAKER SCRIPT (ARPL_Actueel_Maker.R)
# ------------------------------------------------------------------------------
maker_script <- list.files(path = here("src/Totale_Scripts"), pattern = "^ARPL_Actueel_Maker\\.R$", full.names = TRUE, recursive = TRUE)[1]

if (!is.null(maker_script) && file.exists(maker_script)) {
  cat("==================================================\n")
  cat(" STARTEN FASE 2: FINALE ARPL ACTUEEL MAKER SCRIPT\n")
  cat("==================================================\n")
  
  tijd_maker_start <- Sys.time()
  
  tryCatch({
    # Geef de SCENARIO_SELECTIE lijst mee aan het maker script
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
  warning("⚠️ Het script 'ARPL_Actueel_Maker.R' kon niet automatisch gevonden worden.")
}

totaal_eind <- Sys.time()
totale_duurtijd_uur <- round(as.numeric(difftime(totaal_eind, totaal_start, units = "hours")), 2)

cat("\n==================================================\n")
cat(" 🎉 VOLLEDIGE PIPELINE GESLAAGD IN: ", totale_duurtijd_uur, " UUR\n")
cat(" Logboek opgeslagen in: ", log_file_path, "\n")
cat("==================================================\n")
