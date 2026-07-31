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
# 2. EXCEL INLEZEN & SCRIPTNAAM DYNAMISCH GENEREREN
# ------------------------------------------------------------------------------

df_excel <- read_excel(excel_pad)

# ------------------------------------------------------------------------------
# 3. VERWERKINGSFUNCTIE PER GEBIED
# ------------------------------------------------------------------------------

verwerk_gebied <- function(doel_naam, doel_snake, doel_acroniem) {
  
  cat("========================================================\n")
  cat("Start verwerking voor:", doel_naam, "(", doel_acroniem, ")\n")
  cat("========================================================\n")
  
  # A. Bepaal de exacte kolomnaam in de Excel (met underscore)
  excel_kolom <- doel_snake
  
  if (!excel_kolom %in% colnames(df_excel)) {
    warning("Kolom '", excel_kolom, "' niet gevonden in de Excel! Gebied wordt overgeslagen.")
    return(NULL)
  }
  
  # B. Filter handmatige maatwerkscripts & genereer de TV_...Rmd bestandsnaam
  # "Grote wolfsklauw" -> "TV_GroteWolfsklauw.Rmd"
  maatwerk_scripts <- df_excel %>%
    filter(
      tolower(coalesce(Model, "")) == "ja",
      tolower(coalesce(Automatisch, "")) == "nee",
      coalesce(.data[[excel_kolom]], 0) == 1
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
  
  # C. Controleer op automatische soorten in dit gebied
  heeft_automatische_soorten <- df_excel %>%
    filter(
      tolower(coalesce(Automatisch, "")) == "ja",
      coalesce(.data[[excel_kolom]], 0) == 1
    ) %>%
    nrow() > 0
  
  # D. Voeg algemeen automatisch script toe indien van toepassing
  automatisch_script_naam <- paste0(bron_snake, "_Leefgebieden_Simpel.Rmd")
  
  alle_te_verwerken_scripts <- maatwerk_scripts
  
  if (heeft_automatische_soorten) {
    alle_te_verwerken_scripts <- unique(c(alle_te_verwerken_scripts, automatisch_script_naam))
    cat("Inclusief automatisch script:", automatisch_script_naam, "\n")
  }
  
  cat("Aantal te verwerken scripts:", length(alle_te_verwerken_scripts), "\n")
  
  # E. Bepaal en maak de outputmap aan (bijv. src/De_Maten/Scripts_DM)
  output_map <- here("src", doel_snake, paste0("Scripts_", doel_acroniem))
  dir_create(output_map)
  
  # F. Functie om afzonderlijk Rmd-bestand aan te passen en op te slaan
  verwerk_rmd <- function(script_naam) {
    bron_bestand <- file.path(input_map, script_naam)
    
    if (!file_exists(bron_bestand)) {
      warning("❌ Bronbestand NIET gevonden in inputmap: ", bron_bestand)
      return(NULL)
    }
    
    # 1. Garandeer dat de bestandsnaam de nieuwe acroniem-prefix krijgt
    # Bijv. TV_Adder.Rmd -> MH_Adder.Rmd
    nieuwe_script_naam <- script_naam %>%
      str_replace(paste0("^", bron_acroniem, "_"), paste0(doel_acroniem, "_")) %>%
      str_replace_all(bron_snake, doel_snake)
    
    doel_bestand <- file.path(output_map, nieuwe_script_naam)
    
    # 2. Als het Turnhouts Vennegebied zelf is, kopiëren we 1-op-1
    if (doel_acroniem == bron_acroniem) {
      file_copy(bron_bestand, doel_bestand, overwrite = TRUE)
      cat("  ✓ Gekopieerd (origineel):", nieuwe_script_naam, "\n")
    } else {
      # Inhoud inlezen en tekst in het Rmd-bestand vervangen
      tekst <- readLines(bron_bestand, encoding = "UTF-8", warn = FALSE)
      
      # Vervang de gebiedsnamen en acroniemen in de R-code en Markdown tekst
      tekst <- str_replace_all(tekst, bron_naam, doel_naam)           # Turnhouts Vennegebied -> Mechelse Heide
      tekst <- str_replace_all(tekst, bron_snake, doel_snake)         # Turnhouts_Vennegebied -> Mechelse_Heide
      tekst <- str_replace_all(tekst, paste0("/", bron_acroniem, "_"), paste0("/", doel_acroniem, "_"))
      
      # Opslaan in de doelmap
      writeLines(tekst, doel_bestand, useBytes = FALSE)
      cat("  ✓ Verwerkt en opgeslagen:", nieuwe_script_naam, "\n")
    }
  }
  
  # Voer verwerking uit voor alle geselecteerde scripts
  walk(alle_te_verwerken_scripts, verwerk_rmd)
  cat("Klaar voor:", doel_naam, "\n\n")
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
cat("ALLE GEBIEDEN SUCCESVOL VERWERKT!\n")
cat("========================================================\n")