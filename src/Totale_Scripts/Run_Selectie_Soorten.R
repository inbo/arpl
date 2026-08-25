# ==============================================================================
# PIPELINE: RUN ZOOGDIEREN (MET Universele CSV-Check Voor Elk Gebied)
# ==============================================================================

library(readxl)
library(dplyr)
library(stringr)
library(fs)
library(purrr)
library(here)
library(rmarkdown)
library(readr)
library(callr)

# Directe paden naar de bron-scripts voor zoogdieren
alle_tv_scripts <- c(
  here("src/Totale_Scripts/TV_EuropeseOtter.Rmd"),
  here("src/Totale_Scripts/TV_Kwak.Rmd"),
  here("src/Totale_Scripts/TV_Roerdomp.Rmd"),
  here("src/Totale_Scripts/TV_ZwarteHeidelibel.Rmd"),
  here("src/Totale_Scripts/TV_NoordseWitsnuitlibel.Rmd"),
  here("src/Totale_Scripts/TV_Venwitsnuitlibel.Rmd"),
  here("src/Totale_Scripts/TV_Zadelsprinkhaan.Rmd"),
  here("src/Totale_Scripts/TV_Woudaap.Rmd")
)

# Sjabloongegevens
bron_naam     <- "Turnhouts Vennegebied"
bron_snake    <- "Turnhouts_Vennegebied"
bron_acroniem <- "TV"

# Alle 6 de maatwerkgebieden
gebieden_config <- tibble::tribble(
  ~naam,                  ~snake,                 ~acroniem,
  "De Maten",             "De_Maten",             "DM",
  "Heesbossen",           "Heesbossen",           "HB",
  "Kalmthoutse Heide",    "Kalmthoutse_Heide",    "KH",
  "Mechelse Heide",       "Mechelse_Heide",       "MH",
  "Turnhouts Vennegebied","Turnhouts_Vennegebied","TV",
  "Voerstreek",           "Voerstreek",           "VS"
)

# Logboek initialiseren
gecombineerd_logboek <- list()

message("==================================================")
message(" START PIPELINE FOR ZOOGDIEREN                    ")
message("==================================================")

# ------------------------------------------------------------------------------
# HULPFUNCTIE: SLIMME CHECK (IS HTML AL VAN VANDAAG?)
# ------------------------------------------------------------------------------
moet_knitten <- function(bestandsnaam, map) {
  pad_naar_bestand <- file.path(map, bestandsnaam)
  
  if (!file.exists(pad_naar_bestand)) {
    return(TRUE) # Bestaat niet -> verplicht knitten
  }
  
  info <- file.info(pad_naar_bestand)
  wijzigingsdatum <- as.Date(info$mtime)
  vandaag          <- Sys.Date()
  
  if (wijzigingsdatum == vandaag) {
    return(FALSE) # Al van vandaag -> overslaan
  } else {
    return(TRUE)  # Ouder -> opnieuw knitten
  }
}

# ------------------------------------------------------------------------------
# LUS PER GEBIED EN PER SCRIPT
# ------------------------------------------------------------------------------

for (i in 1:nrow(gebieden_config)) {
  doel_naam     <- gebieden_config$naam[i]
  doel_snake    <- gebieden_config$snake[i]
  doel_acroniem <- gebieden_config$acroniem[i]
  
  cat("\n========================================================\n")
  cat(" Start verwerking voor gebied:", doel_naam, "(", doel_acroniem, ")\n")
  cat("========================================================\n")
  
  # Mappen aanmaken
  scripts_output_map <- here("src", doel_snake, paste0("Scripts_", doel_acroniem))
  html_output_map    <- here("data/output", doel_snake, "HTML_Rapporten_Soorten")
  waarnemingen_map   <- here("data/input/Waarnemingen_Soorten", doel_snake)
  
  dir_create(scripts_output_map)
  dir_create(html_output_map)
  
  for (bron_script_pad in alle_tv_scripts) {
    
    if (!file_exists(bron_script_pad)) {
      next
    }
    
    bron_bestandsnaam <- basename(bron_script_pad)
    soort_deel        <- str_remove(str_remove(bron_bestandsnaam, "^TV_"), "\\.Rmd$")
    
    # --------------------------------------------------------------------------
    # UNIVERSELE HOOFDLETTERONGEVOELIGE CSV CHECK (VOOR ELK GEBIED, OOK TV)
    # --------------------------------------------------------------------------
    bestanden_in_map <- list.files(waarnemingen_map, pattern = "\\.csv$", full.names = FALSE)
    
    waarneming_aanwezig <- any(
      str_detect(tolower(bestanden_in_map), "waarnemingen") & 
        str_detect(tolower(bestanden_in_map), tolower(soort_deel))
    )
    
    if (!waarneming_aanwezig) {
      cat("  ⏭️ OVERSLAAN (geen waarnemingen-CSV in map):", soort_deel, "\n")
      gecombineerd_logboek[[paste0(doel_acroniem, "_", soort_deel)]] <- data.frame(
        Gebied = doel_naam, Soort = soort_deel, Type = "Knit HTML",
        Status = "OVERSLAAN_GEEN_WAARNEMINGEN", Fout = "Geen CSV waarnemingenbestand gevonden",
        stringsAsFactors = FALSE
      )
      next
    }
    
    # --------------------------------------------------------------------------
    # RMD SCRIPT OMZETTEN EN TEKST VERVANGEN
    # --------------------------------------------------------------------------
    doel_rmd_naam <- paste0(doel_acroniem, "_", soort_deel, ".Rmd")
    doel_rmd_pad  <- file.path(scripts_output_map, doel_rmd_naam)
    
    if (doel_acroniem == bron_acroniem) {
      file_copy(bron_script_pad, doel_rmd_pad, overwrite = TRUE)
    } else {
      tekst <- readLines(bron_script_pad, encoding = "UTF-8", warn = FALSE)
      
      tekst <- str_replace_all(tekst, bron_naam, doel_naam)
      tekst <- str_replace_all(tekst, bron_snake, doel_snake)
      tekst <- str_replace_all(tekst, paste0("/", bron_acroniem, "_"), paste0("/", doel_acroniem, "_"))
      tekst <- str_replace_all(tekst, str_c("(?<=\\b|_)", bron_acroniem, "(?=\\b|_)"), doel_acroniem)
      
      writeLines(tekst, doel_rmd_pad, useBytes = FALSE)
    }
    
    # --------------------------------------------------------------------------
    # HTML RENDEREN
    # --------------------------------------------------------------------------
    doel_html_naam <- paste0(doel_acroniem, "_", soort_deel, ".html")
    
    if (!moet_knitten(doel_html_naam, html_output_map)) {
      cat("  ⏭️ HTML al up-to-date (vandaag geknit):", doel_html_naam, "\n")
      gecombineerd_logboek[[paste0(doel_acroniem, "_", soort_deel)]] <- data.frame(
        Gebied = doel_naam, Soort = soort_deel, Type = "Knit HTML",
        Status = "GEKOZEN_OVERSLAAN", Fout = "Reeds vandaag geknit",
        stringsAsFactors = FALSE
      )
      next
    }
    
    cat("  -> Knitten naar HTML:", doel_html_naam, "...\n")
    
    tryCatch({
      callr::r(
        function(input, output_file, output_dir) {
          rmarkdown::render(input = input, output_file = output_file, output_dir = output_dir, quiet = FALSE)
        },
        args = list(input = doel_rmd_pad, output_file = doel_html_naam, output_dir = html_output_map),
        show = TRUE
      )
      
      cat("     🎉 SUCCESVOL GEKNIT!\n")
      gecombineerd_logboek[[paste0(doel_acroniem, "_", soort_deel)]] <- data.frame(
        Gebied = doel_naam, Soort = soort_deel, Type = "Knit HTML",
        Status = "SUCCES", Fout = "Geen",
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      cat("     ❌ CRASH BIJ RENDEREN VAN:", doel_rmd_naam, "\n")
      cat("     Foutmelding details:\n", e$message, "\n")
      gecombineerd_logboek[[paste0(doel_acroniem, "_", soort_deel)]] <- data.frame(
        Gebied = doel_naam, Soort = soort_deel, Type = "Knit HTML",
        Status = "CRASH", Fout = e$message,
        stringsAsFactors = FALSE
      )
    })
  }
}

# Logboek wegschrijven
eind_logboek <- bind_rows(gecombineerd_logboek)
log_dir <- here("data/output/Logboeken")
dir_create(log_dir)
log_pad <- file.path(log_dir, paste0("Logboek_Zoogdieren_Run_", format(Sys.time(), "%Y%m%d_%H%M"), ".csv"))
readr::write_excel_csv(eind_logboek, log_pad)

message("\n==================================================")
message(" RUN ZOOGDIEREN AFGEROND!                        ")
message("==================================================")
