# ==============================================================================
# MASTER AUTOMATISERINGSSCRIPT MET SCENARIO & DYNAMISCHE SCRIPT-DETECTIE
# AUTEUR: Bert Van Hecke
# ==============================================================================

library(rmarkdown)
library(purrr)
library(readxl)
library(dplyr)
library(readr)
library(here) # Toegevoegd om 'here' fouten te voorkomen

# ------------------------------------------------------------------------------
# 0. GLOBALE INSTELLINGEN (PAS HIER JE SCENARIO AAN)
# ------------------------------------------------------------------------------
# Geef hier het pad op naar het gewenste Scenario RDS-bestand:
SCENARIO_RDS_PAD <- here("data/input/Scenario_rds/TV_Scenario_bosbehoudss_ss31fix_tvg_vrij_2026.rds")

# Extraheer de naam van het scenario voor de logboeken en bestandsnamen
HUIDIG_SCENARIO  <- gsub("^TV_Scenario_|.rds$", "", basename(SCENARIO_RDS_PAD))

MAP_SCRIPTS      <- here("src/Turnhouts_Vennegebied/Scripts_Scenario") # Map met de 60 scenario Rmd's
OUTPUT_DIR       <- here("data/output/Turnhouts_Vennegebied/HTML_Rapporten_Scenario")

if(!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# Initialiseer een globale loglijst voor de unieke scripts
algemeen_logboek <- list()

message("==================================================")
message(" START MASTER RUN VOOR SCENARIO: ", toupper(HUIDIG_SCENARIO))
message("==================================================")


# ------------------------------------------------------------------------------
# 1. DETECTEER EN RUN DE UNIEKE SCRIPTS (Automatisch uit de map)
# ------------------------------------------------------------------------------
# 1a. Zoek alle Rmd bestanden in de map
alle_scripts <- list.files(path = MAP_SCRIPTS, pattern = "\\.Rmd$", full.names = TRUE)

# 1b. Definieer welk script GEEN uniek script is (als je nog een apart generiek script hebt)
soorten_script_pad <- file.path(MAP_SCRIPTS, "Scenario_TV_Leefgebieden_Simpel.Rmd")

# 1c. Trek het soorten-script af van de totale lijst om de unieke scripts over te houden
unieke_scripts <- setdiff(alle_scripts, soorten_script_pad)

message("=> Gedetecteerde unieke scripts om uit te voeren (Aantal: ", length(unieke_scripts), "):")
print(basename(unieke_scripts))
message("--------------------------------------------------")

# 1d. Loop door de unieke scripts en render ze
for(script in unieke_scripts) {
  bestandsnaam <- basename(script)
  output_html  <- paste0("Rapport_", HUIDIG_SCENARIO, "_", sub(".Rmd$", ".html", bestandsnaam))
  
  message("=> Renderen van script: ", bestandsnaam)
  
  tryCatch({
    rmarkdown::render(
      input = script,
      output_file = output_html,
      output_dir = OUTPUT_DIR,
      # GEFIXTE PARAMETER: Matcht met wat de 60 gegeneerde scripts verwachten!
      params = list(scenario_rds_path = SCENARIO_RDS_PAD), 
      # GEFIXTE SCHONE OMGEVING: Voorkomt vervuiling tussen zware GIS-scripts
      envir = new.env(parent = globalenv()),
      quiet = TRUE
    )
    algemeen_logboek[[bestandsnaam]] <- data.frame(
      Item = bestandsnaam, Type = "Soort Script", Status = "SUCCES", Fout = "Geen", stringsAsFactors = FALSE
    )
  }, error = function(e) {
    message("❌ FOUT bij: ", bestandsnaam, " - Fout: ", e$message)
    algemeen_logboek[[bestandsnaam]] <- data.frame(
      Item = bestandsnaam, Type = "Soort Script", Status = "CRASH", Fout = e$message, stringsAsFactors = FALSE
    )
  })
  
  # Maak het RAM-geheugen leeg na elk GIS-script om bad_alloc te voorkomen
  gc(verbose = FALSE)
}


# ------------------------------------------------------------------------------
# 2. OPTIONEEL: SIMPELE SOORTEN VIA GENERIEK SCRIPT (Indien van toepassing)
# ------------------------------------------------------------------------------
if (file.exists(soorten_script_pad)) {
  excel_data <- read_excel(here("data/input/Excel_files/Soorten_bwk_afstanden.xlsx"))
  
  soorten_lijst <- excel_data %>% 
    filter(Script == "Simpel") %>%                
    filter(Turnhouts_Vennegebied == 1) %>%        
    pull(Soort) %>%                               
    unique() %>%                                  
    na.omit()                                     
  
  message("\n=> Aantal geselecteerde simpele soorten voor ", HUIDIG_SCENARIO, ": ", length(soorten_lijst))
  
  draai_leefgebied_scenario_model <- function(huidige_soort) {
    output_file_name <- paste0(HUIDIG_SCENARIO, "_TV_", huidige_soort, ".html")
    
    res_row <- data.frame(
      Item = huidige_soort,
      Type = "Simpele Soort",
      Status = "SUCCES",
      Fout = "Geen",
      stringsAsFactors = FALSE
    )
    
    message("   -> Starten met simulatie voor: ", toupper(HUIDIG_SCENARIO), " - ", toupper(huidige_soort))
    
    tryCatch({
      rmarkdown::render(
        input = soorten_script_pad, 
        output_file = output_file_name,
        output_dir = OUTPUT_DIR,
        params = list(
          soort_invoer = huidige_soort, 
          scenario_rds_path = SCENARIO_RDS_PAD 
        ), 
        envir = new.env(parent = globalenv()),
        quiet = TRUE
      )
    }, error = function(e) {
      message("   ❌ FOUTMELDING bij ", huidige_soort)
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
message(" Succesvol:                  ", nrow(eind_logboek) - aantal_crashes)
message(" Gecrasht:                   ", aantal_crashes)
message(" Logboek opgeslagen als:     ", log_file_path)
message("==================================================")