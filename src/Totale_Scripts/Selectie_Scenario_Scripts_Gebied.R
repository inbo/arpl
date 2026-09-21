library(readxl)
library(dplyr)
library(stringr)
library(fs)
library(purrr)
library(here)

# ------------------------------------------------------------------------------
# 1. PARAMETERS & GEBIEDSCONFIGURATIE INSTELLEN
# ------------------------------------------------------------------------------

input_map <- here("src/01_Habitat_Scripts")
excel_pad <- here("data/Input/Excel_files/Soortenlijst_Maatwerkgebieden.xlsx")

# Sjabloon / Brongegevens
bron_naam     <- "Turnhouts Vennegebied"
bron_snake    <- "Turnhouts_Vennegebied"
bron_acroniem <- "TV"

# Configuratie-tabel voor alle 6 de gebieden
gebieden_config <- tibble::tribble(
  ~naam,                  ~snake,                 ~acroniem,
  "De Maten",             "De_Maten",             "DM",
  "Heesbossen",           "Heesbossen",           "HB",
  "Kalmthoutse Heide",    "Kalmthoutse_Heide",    "KH",
  "Mechelse Heide",       "Mechelse_Heide",       "MH",
  "Turnhouts Vennegebied","Turnhouts_Vennegebied","TV",
  "Voerstreek",           "Voerstreek",           "VS"
)

# ------------------------------------------------------------------------------
# HULPFUNCTIE: Check op aanwezigheid van waarnemingen-CSV
# ------------------------------------------------------------------------------
heeft_waarnemingen_bestand <- function(soort_naam, map_waarnemingen) {
  if (!dir.exists(map_waarnemingen)) return(FALSE)
  
  alle_csvs <- list.files(path = map_waarnemingen, pattern = "\\.csv$", full.names = FALSE)
  if (length(alle_csvs) == 0) return(FALSE)
  
  schoon_ras <- sub("_wv$", "", soort_naam)
  zoek_naam  <- paste0("Waarnemingen_", schoon_ras)
  
  normaliseer <- function(tekst) {
    tolower(gsub("[ _]", "", tekst))
  }
  
  zoek_naam_norm <- normaliseer(zoek_naam)
  csvs_norm      <- normaliseer(tools::file_path_sans_ext(alle_csvs))
  
  return(zoek_naam_norm %in% csvs_norm)
}

# ------------------------------------------------------------------------------
# 2. EXCEL INLEZEN
# ------------------------------------------------------------------------------

df_excel <- read_excel(excel_pad)

# ------------------------------------------------------------------------------
# 3. VERWERKINGSFUNCTIE PER GEBIED
# ------------------------------------------------------------------------------

verwerk_gebied <- function(doel_naam, doel_snake, doel_acroniem) {
  
  cat("========================================================\n")
  cat("Start verwerking voor:", doel_naam, "(", doel_acroniem, ")\n")
  cat("========================================================\n")
  
  excel_kolom <- doel_snake
  
  if (!excel_kolom %in% colnames(df_excel)) {
    warning("Kolom '", excel_kolom, "' niet gevonden in de Excel! Gebied wordt overgeslagen.")
    return(NULL)
  }
  
  map_waarnemingen_gebied <- here("data/input/Waarnemingen_Soorten", doel_snake)
  
  # A. Filter handmatige maatwerkscripts
  maatwerk_scripts <- df_excel %>%
    filter(
      tolower(coalesce(Model, "")) == "ja",
      tolower(coalesce(Automatisch, "")) == "nee",
      coalesce(.data[[excel_kolom]], 0) == 1
    ) %>%
    filter(map_lgl(`Nederlandse naam`, ~ heeft_waarnemingen_bestand(.x, map_waarnemingen_gebied))) %>%
    mutate(
      Script_Naam = paste0(
        "Scenario_", bron_acroniem, "_", 
        gsub(" ", "", str_to_title(`Nederlandse naam`)), 
        ".R"
      )
    ) %>%
    pull(Script_Naam) %>%
    na.omit() %>%
    unique()
  
  # B. Controleer op automatische soorten in dit gebied MET een CSV-bestand
  heeft_automatische_soorten <- df_excel %>%
    filter(
      tolower(coalesce(Automatisch, "")) == "ja",
      coalesce(.data[[excel_kolom]], 0) == 1
    ) %>%
    filter(map_lgl(`Nederlandse naam`, ~ heeft_waarnemingen_bestand(.x, map_waarnemingen_gebied))) %>%
    nrow() > 0
  
  # C. Algemeen automatisch script
  automatisch_script_naam <- paste0("Scenario_", bron_acroniem, "_Leefgebieden_Simpel.R")
  
  alle_te_verwerken_scripts <- maatwerk_scripts
  
  if (heeft_automatische_soorten) {
    alle_te_verwerken_scripts <- unique(c(alle_te_verwerken_scripts, automatisch_script_naam))
    cat("Inclusief automatisch script:", automatisch_script_naam, "\n")
  }
  
  cat("Aantal te verwerken scripts uit Excel (met CSV-check):", length(alle_te_verwerken_scripts), "\n\n")
  
  # D. Outputmap Schoonmaken & Aanmaken (Dwingt een 100% verse kopie af!)
  output_map <- here("src", doel_snake, "Scripts_Scenario")
  if (dir_exists(output_map)) {
    dir_delete(output_map) # Verwijder oude versies
  }
  dir_create(output_map)
  
  # E. Functie om afzonderlijk R-bestand aan te passen en op te slaan
  verwerk_r <- function(script_naam) {
    
    bron_bestand <- file.path(input_map, script_naam)
    werkelijke_script_naam <- script_naam
    
    if (!file_exists(bron_bestand)) {
      mogelijke_wv_naam <- str_replace(script_naam, "\\.R$", "_wv.R")
      mogelijke_wv_pad  <- file.path(input_map, mogelijke_wv_naam)
      
      if (file_exists(mogelijke_wv_pad)) {
        bron_bestand <- mogelijke_wv_pad
        werkelijke_script_naam <- mogelijke_wv_naam
      }
    }
    
    if (!file_exists(bron_bestand)) {
      cat("  ⚠️ Overgeslagen (nog geen R-scriptbestand op schijf):", script_naam, "\n")
      return(NULL)
    }
    
    # 1. Bestandsnaam aanpassen
    nieuwe_script_naam <- werkelijke_script_naam %>%
      str_replace(paste0("^Scenario_", bron_acroniem, "_"), paste0("Scenario_", doel_acroniem, "_")) %>%
      str_replace(paste0("^", bron_acroniem, "_"), paste0(doel_acroniem, "_")) %>%
      str_replace_all(bron_snake, doel_snake)
    
    doel_bestand <- file.path(output_map, nieuwe_script_naam)
    
    # 2. Inhoud inlezen en tekst in het R-bestand vervangen
    tekst <- readLines(bron_bestand, encoding = "UTF-8", warn = FALSE)
    
    tekst <- str_replace_all(tekst, bron_naam, doel_naam)
    tekst <- str_replace_all(tekst, bron_snake, doel_snake)
    tekst <- str_replace_all(tekst, paste0("/", bron_acroniem, "_"), paste0("/", doel_acroniem, "_"))
    tekst <- str_replace_all(tekst, paste0("Scenario_", bron_acroniem, "_"), paste0("Scenario_", doel_acroniem, "_"))
    tekst <- str_replace_all(tekst, str_c("(?<=\\b|_)", bron_acroniem, "(?=\\b|_)"), doel_acroniem)
    
    # Opslaan
    writeLines(tekst, doel_bestand, useBytes = FALSE)
    cat("  ✓ Vers gegenereerd in Scripts_Scenario:", nieuwe_script_naam, "\n")
  }
  
  walk(alle_te_verwerken_scripts, verwerk_r)
  cat("\nKlaar voor:", doel_naam, "\n\n")
}

# ------------------------------------------------------------------------------
# 4. RUN VOOR ALLE GEBIEDEN
# ------------------------------------------------------------------------------

pwalk(
  list(
    gebieden_config$naam,
    gebieden_config$snake,
    gebieden_config$acroniem
  ),
  verwerk_gebied
)

cat("========================================================\n")
cat("ALLE GEBIEDEN SUCCESVOL HERGENERERD MET DE ALLERNIEUWSTE FIXES!\n")
cat("========================================================\n")
