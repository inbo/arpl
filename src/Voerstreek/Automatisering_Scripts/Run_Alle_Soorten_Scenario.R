# ==============================================================================
# MASTER AUTOMATISERINGSSCRIPT VOOR R-SCRIPTS (.R)
# AUTEUR: Bert Van Hecke
# ==============================================================================

library(purrr)
library(readxl)
library(dplyr)
library(readr)
library(here)

# ------------------------------------------------------------------------------
# 0. GLOBALE INSTELLINGEN
# ------------------------------------------------------------------------------
SCENARIO_RDS_PAD <- here("data/input/Scenario_rds/VS_Scenario_BWK_2025.rds")

HUIDIG_SCENARIO  <- gsub("^VS_Scenario_|.rds$", "", basename(SCENARIO_RDS_PAD))
MAP_SCRIPTS      <- here("src/Voerstreek/Scripts_Scenario")

OUTPUT_DIR       <- here("data/output/Voerstreek/HTML_Rapporten_Scenario", HUIDIG_SCENARIO)

if(!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

algemeen_logboek <- list()
message("==================================================")
message(" START MASTER RUN VOOR SCENARIO: ", toupper(HUIDIG_SCENARIO))
message("==================================================")

# ------------------------------------------------------------------------------
# 1. DETECTEER EN RUN DE UNIEKE R-SCRIPTS
# ------------------------------------------------------------------------------
# Zoek naar .R bestanden (case-insensitive)
alle_scripts <- list.files(path = MAP_SCRIPTS, pattern = "\\.r$", full.names = TRUE, ignore.case = TRUE)

# Filter het generieke script eruit als dat in dezelfde map staat
soorten_script_pad <- alle_scripts[grep("Scenario_VS_Leefgebieden_Simpel\\.r$", alle_scripts, ignore.case = TRUE)]

if (length(soorten_script_pad) == 0) {
  soorten_script_pad <- file.path(MAP_SCRIPTS, "Scenario_VS_Leefgebieden_Simpel.R")
} else {
  soorten_script_pad <- soorten_script_pad[1]
}

unieke_scripts <- setdiff(alle_scripts, soorten_script_pad)

message("=> Gedetecteerde unieke R-scripts (Aantal: ", length(unieke_scripts), "):")
print(basename(unieke_scripts))
message("--------------------------------------------------")

for(script in unieke_scripts) {
  bestandsnaam <- basename(script)
  message("=> Uitvoeren van script via source(): ", bestandsnaam)
  
  tryCatch({
    # Maak een schone omgeving aan en ken variabelen/parameters toe die het script verwacht
    script_env <- new.env(parent = globalenv())
    script_env$scenario_rds_path <- SCENARIO_RDS_PAD
    script_env$output_dir        <- OUTPUT_DIR
    script_env$huidig_scenario   <- HUIDIG_SCENARIO
    
    # Voer het R-script uit
    source(script, local = script_env)
    
    algemeen_logboek[[bestandsnaam]] <- data.frame(
      Item = bestandsnaam, Type = "Soort Script", Status = "SUCCES", Fout = "Geen", stringsAsFactors = FALSE
    )
  }, error = function(e) {
    message("❌ FOUT bij: ", bestandsnaam, " - Fout: ", e$message)
    algemeen_logboek[[bestandsnaam]] <- data.frame(
      Item = bestandsnaam, Type = "Soort Script", Status = "CRASH", Fout = e$message, stringsAsFactors = FALSE
    )
  })
  
  gc(verbose = FALSE)
}

# ------------------------------------------------------------------------------
# 2. OPTIONEEL: SIMPELE SOORTEN VIA GENERIEK R-SCRIPT
# ------------------------------------------------------------------------------
if (file.exists(soorten_script_pad)) {
  excel_data <- read_excel(here("data/input/Excel_files/Soortenlijst_Maatwerkgebieden_Gefilterd.xlsx"))
  colnames(excel_data) <- tolower(colnames(excel_data))
  
  # Verwacht dat 'nederlandse naam' of 'soort' de soortnaam bevat
  naam_kolom <- if("soort" %in% colnames(excel_data)) "soort" else "nederlandse naam"
  
  soorten_lijst <- excel_data %>% 
    filter(automatisch == "ja") %>%                   
    filter(voerstreek == 1) %>%        
    pull(!!sym(naam_kolom)) %>%                                   
    unique() %>%                                     
    na.omit()                                        
  
  message("\n=> Aantal geselecteerde simpele soorten voor ", HUIDIG_SCENARIO, ": ", length(soorten_lijst))
  
  draai_leefgebied_scenario_model <- function(huidige_soort) {
    res_row <- data.frame(
      Item = huidige_soort,
      Type = "Simpele Soort",
      Status = "SUCCES",
      Fout = "Geen",
      stringsAsFactors = FALSE
    )
    
    message("   -> Starten met simulatie voor: ", toupper(HUIDIG_SCENARIO), " - ", toupper(huidige_soort))
    tryCatch({
      script_env <- new.env(parent = globalenv())
      script_env$soort_invoer      <- huidige_soort
      script_env$scenario_rds_path <- SCENARIO_RDS_PAD
      script_env$output_dir        <- OUTPUT_DIR
      script_env$huidig_scenario   <- HUIDIG_SCENARIO
      
      source(soorten_script_pad, local = script_env)
      
    }, error = function(e) {
      message("   ❌ FOUTMELDING bij ", huidige_soort, ": ", e$message)
      res_row$Status <<- "CRASH"
      res_row$Fout <<- e$message
    })
    
    gc(verbose = FALSE)
    return(res_row)
  }
  
  if (length(soorten_lijst) > 0) {
    soorten_logboek <- purrr::map_dfr(soorten_lijst, draai_leefgebied_scenario_model)
    eind_logboek    <- bind_rows(bind_rows(algemeen_logboek), soorten_logboek)
  } else {
    eind_logboek    <- bind_rows(algemeen_logboek)
  }
} else {
  eind_logboek <- bind_rows(algemeen_logboek)
}

# ------------------------------------------------------------------------------
# 3. LOGBESTAND WEGSCHRIJVEN & SAMENVATTING
# ------------------------------------------------------------------------------
log_file_path <- file.path(OUTPUT_DIR, paste0("Logboek_ScenarioRun_", HUIDIG_SCENARIO, "_", format(Sys.time(), "%Y%m%d_%H%M"), ".csv"))
readr::write_excel_csv(eind_logboek, log_file_path)
aantal_crashes <- sum(eind_logboek$Status == "CRASH")

message("\n==================================================")
message(" MASTER RUN COMPLEET VOOR ", toupper(HUIDIG_SCENARIO))
message(" Totaal onderdelen gedraaid: ", nrow(eind_logboek))
message(" Succesvol:                   ", nrow(eind_logboek) - aantal_crashes)
message(" Gecrasht:                    ", aantal_crashes)
message(" Logboek opgeslagen als:     ", log_file_path)
message("==================================================")
