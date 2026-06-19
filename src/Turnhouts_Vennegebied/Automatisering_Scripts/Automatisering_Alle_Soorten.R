# ==============================================================================
# MASTER AUTOMATISERINGSSCRIPT VOOR REGULIERE TV-SCRIPTS (KNITTEN)
# AUTEUR: Bert Van Hecke
# ==============================================================================

library(rmarkdown)
library(purrr)
library(readxl)
library(dplyr)
library(readr)
library(here)

# ------------------------------------------------------------------------------
# 0. INSTELLINGEN EN STRATEGISCHE MAPSTRUCTUUR
# ------------------------------------------------------------------------------
MAP_SCRIPTS_TV   <- here("src/Turnhouts_Vennegebied/Scripts_TV")  # Pad naar de TV-scripts (pas aan indien nodig)
OUTPUT_DIR       <- here("data/output/Turnhouts_Vennegebied/HTML_Rapporten_Soorten")

if(!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# Initialiseer de loglijst
algemeen_logboek <- list()

message("==================================================")
message(" START MASTER RUN: NORMALE TV SCRIPTS KNITTEN     ")
message("==================================================")


# ------------------------------------------------------------------------------
# 1. DETECTEER EN KNIT DE UNIEKE TV-SCRIPTS
# ------------------------------------------------------------------------------
# 1a. Zoek alle Rmd bestanden in de map Scripts_TV
alle_tv_scripts <- list.files(path = MAP_SCRIPTS_TV, pattern = "\\.Rmd$", full.names = TRUE)

# 1b. Definieer welk script per soort moet draaien (en dus niet hier 1x geknit mag worden)
soorten_script_tv_pad <- file.path(MAP_SCRIPTS_TV, "Leefgebieden_Simpel_TV.Rmd")

# 1c. Filter het soortenscript eruit om alleen de unieke scripts over te houden
unieke_tv_scripts <- setdiff(alle_tv_scripts, soorten_script_tv_pad)

message("=> Gedetecteerde unieke TV-scripts om te knitten (Aantal: ", length(unieke_tv_scripts), "):")
print(basename(unieke_tv_scripts))
message("--------------------------------------------------")

# 1d. Loop door de unieke scripts en knit ze
for(script in unieke_tv_scripts) {
  bestandsnaam <- basename(script)
  output_html  <- sub(".Rmd$", ".html", bestandsnaam) # Gewoon de originele naam, maar dan .html
  
  message("=> Knitten van uniek script: ", bestandsnaam)
  
  tryCatch({
    rmarkdown::render(
      input = script,
      output_file = output_html,
      output_dir = OUTPUT_DIR,
      quiet = TRUE
    )
    algemeen_logboek[[bestandsnaam]] <- data.frame(
      Item = bestandsnaam, Type = "Uniek TV Script", Status = "SUCCES", Fout = "Geen", stringsAsFactors = FALSE
    )
  }, error = function(e) {
    message("❌ FOUT bij: ", bestandsnaam)
    algemeen_logboek[[bestandsnaam]] <- data.frame(
      Item = bestandsnaam, Type = "Uniek TV Script", Status = "CRASH", Fout = e$message, stringsAsFactors = FALSE
    )
  })
}


# ------------------------------------------------------------------------------
# 2. TREK DE SOORTENLIJST DYNAMISCH UIT DE EXCEL
# ------------------------------------------------------------------------------
excel_data <- read_excel(here("data/input/Excel_files/Soorten_bwk_afstanden.xlsx"))

soorten_lijst <- excel_data %>% 
  filter(Script == "Simpel") %>%                # Alleen de simpele wasstraat
  filter(Turnhouts_Vennegebied == 1) %>%        # Alleen soorten op de TV-lijst
  pull(Soort) %>%                               # Trek de kolom 'Soort' leeg als vector
  unique() %>%                                  # Verwijder eventuele dubbele vermeldingen
  na.omit()                                     # Verwijder lege rijen

message("\n=> Aantal geselecteerde simpele soorten voor TV-run: ", length(soorten_lijst))


# ------------------------------------------------------------------------------
# 3. RUN HET 'SIMPEL_TV' SCRIPT PER SOORT (Knitten met parameter)
# ------------------------------------------------------------------------------
draai_leefgebied_model <- function(huidige_soort) {
  output_file_name <- paste0("TV_", huidige_soort, ".html")
  
  res_row <- data.frame(
    Item = huidige_soort,
    Type = "Simpele Soort (Normaal)",
    Status = "SUCCES",
    Fout = "Geen",
    stringsAsFactors = FALSE
  )
  
  message("   -> Knitten voor soort: ", toupper(huidige_soort))
  
  tryCatch({
    rmarkdown::render(
      input = soorten_script_tv_pad, # Knif de 'Leefgebieden_Simpel_TV.Rmd' uit de Scripts_TV map
      output_file = output_file_name,
      output_dir = OUTPUT_DIR,
      params = list(soort_invoer = huidige_soort), # Alleen de soortnaam doorgeven, geen scenario nodig
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
  soorten_logboek <- purrr::map_dfr(soorten_lijst, draai_leefgebied_model)
  eind_logboek    <- bind_rows(bind_rows(algemeen_logboek), soorten_logboek)
} else {
  eind_logboek    <- bind_rows(algemeen_logboek)
  message("Waarschuwing: Geen simpele soorten gevonden die voldoen aan de filters.")
}


# ------------------------------------------------------------------------------
# 4. LOGBESTAND WEGSCHRIJVEN & SAMENVATTING
# ------------------------------------------------------------------------------
log_file_path <- file.path(OUTPUT_DIR, paste0("Logboek_TV_NormaalRun_", format(Sys.time(), "%Y%m%d_%H%M"), ".csv"))
readr::write_excel_csv(eind_logboek, log_file_path)

aantal_crashes <- sum(eind_logboek$Status == "CRASH")
message("\n==================================================")
message(" MASTER RUN TV NORMAAL COMPLEET!                  ")
message(" Totaal onderdelen geknit:   ", nrow(eind_logboek))
message(" Succesvol:                  ", nrow(eind_logboek) - aantal_crashes)
message(" Gecrasht:                   ", aantal_crashes)
message(" Logboek opgeslagen als:     ", log_file_path)
message("==================================================")