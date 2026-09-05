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
  here("src/Totale_Scripts/TV_Blauwborst.Rmd"),
  here("src/Totale_Scripts/TV_Boompieper.Rmd"),
  here("src/Totale_Scripts/TV_BruineKiekendief.Rmd"),
  here("src/Totale_Scripts/TV_Grutto.Rmd"),
  here("src/Totale_Scripts/TV_Kwartelkonging.Rmd"),
  here("src/Totale_Scripts/TV_MiddelsteBonteSpecht.Rmd"),
  here("src/Totale_Scripts/TV_Nachtzwaluw.Rmd"),
  here("src/Totale_Scripts/TV_Paapje.Rmd"),
  here("src/Totale_Scripts/TV_Porseleinhoen.Rmd"),
  here("src/Totale_Scripts/TV_Roerdomp.Rmd"),
  here("src/Totale_Scripts/TV_Watersnip.Rmd"),
  here("src/Totale_Scripts/TV_Wespendief.Rmd"),
  here("src/Totale_Scripts/TV_Wulp.Rmd"),
  here("src/Totale_Scripts/TV_ZwarteSpecht.Rmd"),
  here("src/Totale_Scripts/TV_Zwartkopmeeuw.Rmd"),
  here("src/Totale_Scripts/TV_Zomertortel.Rmd"),
  here("src/Totale_Scripts/TV_Leefgebieden_Simpel.Rmd")
)

# Sjabloongegevens
bron_naam     <- "Turnhouts Vennegebied"
bron_snake    <- "Turnhouts_Vennegebied"
bron_acroniem <- "TV"

# Alle 6 de maatwerkgebieden
gebieden_config <- tibble::tribble(
  ~naam,                  ~snake,                 ~acroniem,
   #"De Maten",             "De_Maten",             "DM",
  # "Heesbossen",           "Heesbossen",           "HB",
  # "Kalmthoutse Heide",    "Kalmthoutse_Heide",    "KH",
  # "Mechelse Heide",       "Mechelse_Heide",       "MH",
  "Turnhouts Vennegebied","Turnhouts_Vennegebied","TV"
  # "Voerstreek",           "Voerstreek",           "VS"
)

# Logboek initialiseren
gecombineerd_logboek <- list()

message("==================================================")
message(" START PIPELINE FOR SOORTEN                       ")
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
    
    # 1. Definieer soort_deel_volledig DIRECT aan het begin van de lus
    soort_deel_volledig <- str_remove(str_remove(bron_bestandsnaam, "^TV_"), "\\.Rmd$")
    
    # --------------------------------------------------------------------------
    # UNIVERSELE HOOFDLETTERONGEVOELIGE CSV CHECK (MET EXPLICITE FALLBACK)
    # --------------------------------------------------------------------------
    bestanden_in_map <- list.files(waarnemingen_map, pattern = "\\.csv$", full.names = FALSE)
    bestanden_lc     <- tolower(bestanden_in_map)
    
    # Stap A: Check op exacte match met volledige naam (bijv. "wulp_wv")
    zoek_volledig <- tolower(soort_deel_volledig)
    waarneming_aanwezig <- any(
      str_detect(bestanden_lc, "waarnemingen") & str_detect(bestanden_lc, zoek_volledig)
    )
    
    # Stap B: Indien niet gevonden en er staat "_wv" in de naam, val terug op naam zonder "_wv" (bijv. "bergeend")
    if (!waarneming_aanwezig && str_detect(soort_deel_volledig, "_wv$")) {
      zoek_opgeschoond <- tolower(str_remove(soort_deel_volledig, "_wv$"))
      waarneming_aanwezig <- any(
        str_detect(bestanden_lc, "waarnemingen") & str_detect(bestanden_lc, zoek_opgeschoond)
      )
    }
    
    if (!waarneming_aanwezig) {
      cat("  ⏭️ OVERSLAAN (geen waarnemingen-CSV gevonden voor):", soort_deel_volledig, "\n")
      gecombineerd_logboek[[paste0(doel_acroniem, "_", soort_deel_volledig)]] <- data.frame(
        Gebied = doel_naam, Soort = soort_deel_volledig, Type = "Knit HTML",
        Status = "OVERSLAAN_GEEN_WAARNEMINGEN", Fout = "Geen CSV waarnemingenbestand gevonden",
        stringsAsFactors = FALSE
      )
      next
    }
    
    # --------------------------------------------------------------------------
    # RMD SCRIPT OMZETTEN EN TEKST VERVANGEN
    # --------------------------------------------------------------------------
    doel_rmd_naam <- paste0(doel_acroniem, "_", soort_deel_volledig, ".Rmd")
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
    doel_html_naam <- paste0(doel_acroniem, "_", soort_deel_volledig, ".html")
    
    if (!moet_knitten(doel_html_naam, html_output_map)) {
      cat("  ⏭️ HTML al up-to-date (vandaag geknit):", doel_html_naam, "\n")
      gecombineerd_logboek[[paste0(doel_acroniem, "_", soort_deel_volledig)]] <- data.frame(
        Gebied = doel_naam, Soort = soort_deel_volledig, Type = "Knit HTML",
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
      gecombineerd_logboek[[paste0(doel_acroniem, "_", soort_deel_volledig)]] <- data.frame(
        Gebied = doel_naam, Soort = soort_deel_volledig, Type = "Knit HTML",
        Status = "SUCCES", Fout = "Geen",
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      cat("     ❌ CRASH BIJ RENDEREN VAN:", doel_rmd_naam, "\n")
      cat("     Foutmelding details:\n", e$message, "\n")
      gecombineerd_logboek[[paste0(doel_acroniem, "_", soort_deel_volledig)]] <- data.frame(
        Gebied = doel_naam, Soort = soort_deel_volledig, Type = "Knit HTML",
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
message(" RUN SOORTEN AFGEROND!                        ")
message("==================================================")
