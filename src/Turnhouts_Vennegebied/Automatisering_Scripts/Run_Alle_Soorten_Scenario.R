# ==============================================================================
# MASTER AUTOMATISERINGSSCRIPT MET SCENARIO & DYNAMISCHE SCRIPT-DETECTIE
# AUTEUR: Bert Van Hecke & AI
# ==============================================================================

library(rmarkdown)
library(purrr)
library(readxl)
library(dplyr)
library(readr)

# ------------------------------------------------------------------------------
# 0. GLOBALE INSTELLINGEN (PAS HIER JE SCENARIO AAN)
# ------------------------------------------------------------------------------
HUIDIG_SCENARIO <- "Scenario_A" # Vul hier je actieve scenario in
MAP_SCRIPTS      <- here("src/Turnhouts_Vennegebied/Scripts_Scenario's") # De map waar al je Rmd's staan
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

# 1b. Definieer welk script GEEN uniek script is (omdat het per soort draait)
soorten_script_pad <- file.path(MAP_SCRIPTS, "Leefgebieden_Simpel_TV_Scenario.Rmd")

# 1c. Trek het soorten-script af van de totale lijst om de unieke scripts over te houden
unieke_scripts <- setdiff(alle_scripts, soorten_script_pad)

message("=> Gedetecteerde unieke scripts om uit te voeren (Aantal: ", length(unieke_scripts), "):")
print(basename(unieke_scripts))
message("--------------------------------------------------")

# 1d. Loop door de unieke scripts en render ze
for(script in unieke_scripts) {
  bestandsnaam <- basename(script)
  output_html  <- paste0(HUIDIG_SCENARIO, "_", sub(".Rmd$", ".html", bestandsnaam))
  
  message("=> Renderen van uniek script: ", bestandsnaam)
  
  tryCatch({
    rmarkdown::render(
      input = script,
      output_file = output_html,
      output_dir = OUTPUT_DIR,
      params = list(scenario = HUIDIG_SCENARIO), # Het script moet deze parameter accepteren indien nodig
      quiet = TRUE
    )
    algemeen_logboek[[bestandsnaam]] <- data.frame(
      Item = bestandsnaam, Type = "Uniek Script", Status = "SUCCES", Fout = "Geen", stringsAsFactors = FALSE
    )
  }, error = function(e) {
    message("❌ FOUT bij: ", bestandsnaam)
    algemeen_logboek[[bestandsnaam]] <- data.frame(
      Item = bestandsnaam, Type = "Uniek Script", Status = "CRASH", Fout = e$message, stringsAsFactors = FALSE
    )
  })
}


# ------------------------------------------------------------------------------
# 2. TREK DE SOORTENLIJST DYNAMISCH UIT DE EXCEL
# ------------------------------------------------------------------------------
excel_data <- read_excel(here("data/input/Soorten_bwk_afstanden.xlsx"))

soorten_lijst <- excel_data %>% 
  filter(Script == "Simpel") %>%                # Alleen de simpele wasstraat
  filter(Turnhouts_Vennegebied == 1) %>%        # Alleen soorten op de TV-lijst
  pull(Soort) %>%                               # Trek de kolom 'Soort' leeg als vector
  unique() %>%                                  # Verwijder eventuele dubbele vermeldingen
  na.omit()                                     # Verwijder lege rijen

message("\n=> Aantal geselecteerde simpele soorten voor ", HUIDIG_SCENARIO, ": ", length(soorten_lijst))


# ------------------------------------------------------------------------------
# 3. RUN HET SPECIFIEKE SCRIPT PER "SIMPELE SOORT"
# ------------------------------------------------------------------------------
# De functie die per soort wordt uitgevoerd
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
      input = soorten_script_pad, # Maakt gebruik van het gedefinieerde pad naar het soorten-script
      output_file = output_file_name,
      output_dir = OUTPUT_DIR,
      params = list(
        soort_invoer = huidige_soort, 
        scenario_invoer = HUIDIG_SCENARIO 
      ), 
      quiet = TRUE
    )
  }, error = function(e) {
    message("   ❌ FOUTMELDING bij ", huidige_soort)
    res_row$Status <<- "CRASH"
    res_row$Fout <<- e$message
  })
  
  return(res_row)
}

# Uitvoeren van de soorten-loop indien er soorten zijn gevonden
if (length(soorten_lijst) > 0) {
  soorten_logboek <- purrr::map_dfr(soorten_lijst, draai_leefgebied_scenario_model)
  eind_logboek    <- bind_rows(bind_rows(algemeen_logboek), soorten_logboek)
} else {
  eind_logboek    <- bind_rows(algemeen_logboek)
  message("Waarschuwing: Geen simpele soorten gevonden die voldoen aan de filters.")
}


# ------------------------------------------------------------------------------
# 4. LOGBESTAND WEGSCHRIJVEN & SAMENVATTING
# ------------------------------------------------------------------------------
log_file_path <- paste0("Logboek_ScenarioRun_", HUIDIG_SCENARIO, "_", format(Sys.time(), "%Y%m%d_%H%M"), ".csv")
readr::write_excel_csv(eind_logboek, log_file_path)

aantal_crashes <- sum(eind_logboek$Status == "CRASH")
message("\n==================================================")
message(" MASTER RUN COMPLEET VOOR ", toupper(HUIDIG_SCENARIO))
message(" Totaal onderdelen gedraaid: ", nrow(eind_logboek))
message(" Succesvol:                  ", nrow(eind_logboek) - aantal_crashes)
message(" Gecrasht:                   ", aantal_crashes)
message(" Logboek opgeslagen als:     ", log_file_path)
message("==================================================")