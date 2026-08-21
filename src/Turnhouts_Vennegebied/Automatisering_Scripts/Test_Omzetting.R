library(readxl)
library(dplyr)
library(stringr)
library(fs)
library(purrr)
library(here)

# ------------------------------------------------------------------------------
# 1. PARAMETERS & GEBIEDSCONFIGURATIE INSTELLEN
# ------------------------------------------------------------------------------

input_map <- here("src/Totale_Scripts")
excel_pad <- here("data/Input/Excel_files/Soortenlijst_Maatwerkgebieden.xlsx")

# Sjabloon / Brongegevens
bron_naam <- "Turnhouts Vennegebied"
bron_snake <- "Turnhouts_Vennegebied"
bron_acroniem <- "TV"

# Configuratie-tabel met specifieke TEST-SOORTEN per gebied
gebieden_config <- tibble::tribble(
  ~naam,                   ~snake,                  ~acroniem, ~test_soorten,
  "De Maten",              "De_Maten",              "DM",      character(0),
  "Heesbossen",            "Heesbossen",            "HB",      character(0),
  "Kalmthoutse Heide",     "Kalmthoutse_Heide",     "KH",      c("Gladde slang"),
  "Mechelse Heide",        "Mechelse_Heide",        "MH",      character(0),
  "Turnhouts Vennegebied", "Turnhouts_Vennegebied", "TV",      character(0),
  "Voerstreek",            "Voerstreek",            "VS",      character(0)
)

# ------------------------------------------------------------------------------
# 2. EXCEL INLEZEN
# ------------------------------------------------------------------------------

df_excel <- read_excel(excel_pad)

# ------------------------------------------------------------------------------
# 3. VERWERKINGSFUNCTIE PER GEBIED
# ------------------------------------------------------------------------------

verwerk_gebied <- function(doel_naam, doel_snake, doel_acroniem, test_soorten) {
  
  # Sla gebieden over waarvoor geen testsoorten zijn ingesteld
  if (length(test_soorten) == 0) {
    return(NULL)
  }
  
  cat("========================================================\n")
  cat("Start testverwerking voor:", doel_naam, "(", doel_acroniem, ")\n")
  cat("Testsoorten:", paste(test_soorten, collapse = ", "), "\n")
  cat("========================================================\n")
  
  excel_kolom <- doel_snake
  
  if (!excel_kolom %in% colnames(df_excel)) {
    warning("Kolom '", excel_kolom, "' niet gevonden in de Excel! Gebied wordt overgeslagen.")
    return(NULL)
  }
  
  # Filter uitsluitend op de gekozen testsoorten uit de Excel
  maatwerk_scripts <- df_excel %>%
    filter(
      tolower(coalesce(Model, "")) == "ja",
      tolower(coalesce(Automatisch, "")) == "nee",
      coalesce(.data[[excel_kolom]], 0) == 1,
      tolower(`Nederlandse naam`) %in% tolower(test_soorten)
    ) %>%
    mutate(
      Script_Naam = paste0(
        bron_acroniem, "_", 
        gsub(" ", "", str_to_title(`Nederlandse naam`)), 
        ".Rmd"
      )
    ) %>%
    pull(Script_Naam) %>%
    na.omit() %>%
    unique()
  
  if (length(maatwerk_scripts) == 0) {
    cat("  ⚠️ Geen overeenkomende scripts gevonden in de Excel voor deze soorten.\n\n")
    return(NULL)
  }
  
  cat("Aantal weg te schrijven test-scripts:", length(maatwerk_scripts), "\n\n")
  
  # Bepaal en maak de outputmap aan (bijv. src/De_Maten/Scripts_DM)
  output_map <- here("src", doel_snake, paste0("Scripts_", doel_acroniem))
  dir_create(output_map)
  
  # Functie om afzonderlijk Rmd-bestand aan te passen en op te slaan
  verwerk_rmd <- function(script_naam) {
    
    # Check 1: Bestaat het standaard script (bijv. TV_Adder.Rmd)?
    bron_bestand <- file.path(input_map, script_naam)
    werkelijke_script_naam <- script_naam
    
    # Check 2: Indien niet gevonden, check op '_wv.Rmd'
    if (!file_exists(bron_bestand)) {
      mogelijke_wv_naam <- str_replace(script_naam, "\\.Rmd$", "_wv.Rmd")
      mogelijke_wv_pad  <- file.path(input_map, mogelijke_wv_naam)
      
      if (file_exists(mogelijke_wv_pad)) {
        bron_bestand <- mogelijke_wv_pad
        werkelijke_script_naam <- mogelijke_wv_naam
      }
    }
    
    if (!file_exists(bron_bestand)) {
      cat("  ⚠️ Overgeslagen (bronbestand niet gevonden):", script_naam, "\n")
      return(NULL)
    }
    
    # Bestandsnaam aanpassen met het nieuwe acroniem
    nieuwe_script_naam <- werkelijke_script_naam %>%
      str_replace(paste0("^", bron_acroniem, "_"), paste0(doel_acroniem, "_")) %>%
      str_replace_all(bron_snake, doel_snake)
    
    doel_bestand <- file.path(output_map, nieuwe_script_naam)
    
    # Inhoud inlezen en tekst in het Rmd-bestand vervangen
    tekst <- readLines(bron_bestand, encoding = "UTF-8", warn = FALSE)
    
    tekst <- str_replace_all(tekst, bron_naam, doel_naam)
    tekst <- str_replace_all(tekst, bron_snake, doel_snake)
    tekst <- str_replace_all(tekst, paste0("/", bron_acroniem, "_"), paste0("/", doel_acroniem, "_"))
    tekst <- str_replace_all(tekst, str_c("(?<=\\b|_)", bron_acroniem, "(?=\\b|_)"), doel_acroniem)
    
    # Opslaan
    writeLines(tekst, doel_bestand, useBytes = FALSE)
    cat("  ✓ Aangemaakt en opgeslagen:", nieuwe_script_naam, "\n")
  }
  
  walk(maatwerk_scripts, verwerk_rmd)
  cat("\nKlaar voor:", doel_naam, "\n\n")
}

# ------------------------------------------------------------------------------
# 4. RUN VOOR ALLE GEBIEDEN
# ------------------------------------------------------------------------------

pwalk(
  list(
    gebieden_config$naam,
    gebieden_config$snake,
    gebieden_config$acroniem,
    gebieden_config$test_soorten
  ),
  verwerk_gebied
)

cat("========================================================\n")
cat("TEST-BESTANDEN SUCCESVOL WEGGESCHREVEN!\n")
cat("========================================================\n")