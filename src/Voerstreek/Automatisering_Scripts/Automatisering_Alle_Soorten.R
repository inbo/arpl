# ==============================================================================
# MASTER AUTOMATISERINGSSCRIPT VOOR REGULIERE VS-SCRIPTS (KNITTEN)
# AUTEUR: Bert Van Hecke
# UPDATE: Slimme skip-logica (bestaat al + vandaag aangemaakt)
# ==============================================================================

library(rmarkdown)
library(purrr)
library(readxl)
library(dplyr)
library(readr)
library(here)
library(callr)

# ------------------------------------------------------------------------------
# 0. INSTELLINGEN EN STRATEGISCHE MAPSTRUCTUUR
# ------------------------------------------------------------------------------
MAP_SCRIPTS_VS     <- here("src/Voerstreek/Scripts_VS/")
OUTPUT_DIR         <- here("data/output/Voerstreek/HTML_Rapporten_Soorten")
MAP_WAARNEMINGEN   <- here("data/input/Waarnemingen_Soorten/Voerstreek")

if(!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# Hulpfunctie: Controleert of er een waarnemingenbestand bestaat voor de soort
heeft_waarnemingen_bestand <- function(soort_naam, map_waarnemingen) {
  # Logica voor de bestandsnaam:
  # Als de soort eindigt op '_wv', checken we of het om de Wulp gaat
  if (grepl("_wv$", soort_naam)) {
    if (tolower(soort_naam) == "wulp_wv") {
      zoek_naam <- "Waarnemingen_Wulp_wv"
    } else {
      # Andere _wv soorten (bijv. Regenwulp_wv -> Waarnemingen_Regenwulp)
      schoon_ras <- sub("_wv$", "", soort_naam)
      zoek_naam <- paste0("Waarnemingen_", schoon_ras)
    }
  } else {
    zoek_naam <- paste0("Waarnemingen_", soort_naam)
  }
  
  # Zoek naar bestanden die beginnen met zoek_naam (gevoelig voor extensies zoals .rds, .csv, .RData)
  gevonden_bestanden <- list.files(
    path = map_waarnemingen, 
    pattern = paste0("^", zoek_naam, "\\."), 
    full.names = TRUE
  )
  
  return(length(gevonden_bestanden) > 0)
}


# ------------------------------------------------------------------------------
# HULPFUNCTIE: Moet dit bestand opnieuw geknit worden?
# ------------------------------------------------------------------------------
moet_knitten <- function(bestandsnaam, map) {
  pad_naar_bestand <- file.path(map, bestandsnaam)
  
  # Check 1: Bestaat het bestand überhaupt?
  if (!file.exists(pad_naar_bestand)) {
    return(TRUE) # Bestaat niet, dus we moeten knitten
  }
  
  # Check 2: Is het vandaag aangemaakt/gewijzigd?
  info <- file.info(pad_naar_bestand)
  wijzigingsdatum <- as.Date(info$mtime)
  vandaag         <- Sys.Date()
  
  if (wijzigingsdatum == vandaag) {
    return(FALSE) # Bestaat al en is van vandaag -> overslaan!
  } else {
    return(TRUE)  # Bestaat wel, maar is van een eerdere dag -> opnieuw knitten
  }
}


# ------------------------------------------------------------------------------
# 1. DETECTEER EN KNIT DE UNIEKE VS-SCRIPTS
# ------------------------------------------------------------------------------
alle_tv_scripts <- list.files(path = MAP_SCRIPTS_VS, pattern = "\\.Rmd$", full.names = TRUE)
soorten_script_tv_pad <- file.path(MAP_SCRIPTS_VS, "VS_Leefgebieden_Simpel.Rmd")
unieke_tv_scripts <- setdiff(alle_tv_scripts, soorten_script_tv_pad)

message("=> Gedetecteerde unieke VS-scripts om te knitten (Aantal: ", length(unieke_tv_scripts), "):")
print(basename(unieke_tv_scripts))
message("--------------------------------------------------")

for(script in unieke_tv_scripts) {
  bestandsnaam <- basename(script)
  output_html  <- sub(".Rmd$", ".html", bestandsnaam)
  
  # --- SLIMME CHECK ---
  if (!moet_knitten(output_html, OUTPUT_DIR)) {
    message("⏭️  OVERSLAAN: ", output_html, " bestaat al en is vandaag al gemaakt.")
    algemeen_logboek[[bestandsnaam]] <- data.frame(
      Item = bestandsnaam, Type = "Uniek VS Script", Status = "GEKOZEN_OVERSLAAN", Fout = "Reeds vandaag geknit", stringsAsFactors = FALSE
    )
    next # Spring direct naar de volgende in de loop
  }
  
  message("=> Knitten van uniek script: ", bestandsnaam)
  
  tryCatch({
    callr::r(
      function(input, output_file, output_dir) {
        rmarkdown::render(input = input, output_file = output_file, output_dir = output_dir, quiet = TRUE)
      },
      args = list(input = script, output_file = output_html, output_dir = OUTPUT_DIR)
    )
    
    algemeen_logboek[[bestandsnaam]] <- data.frame(
      Item = bestandsnaam, Type = "Uniek VS Script", Status = "SUCCES", Fout = "Geen", stringsAsFactors = FALSE
    )
  }, error = function(e) {
    message("❌ FOUT bij: ", bestandsnaam)
    algemeen_logboek[[bestandsnaam]] <- data.frame(
      Item = bestandsnaam, Type = "Uniek VS Script", Status = "CRASH", Fout = e$message, stringsAsFactors = FALSE
    )
  })
}


# ------------------------------------------------------------------------------
# 2. TREK DE SOORTENLIJST DYNAMISCH UIT DE EXCEL & CHECK WAARNEMINGEN
# ------------------------------------------------------------------------------
excel_data <- read_excel(here("data/input/Excel_files/Soorten_bwk_afstanden.xlsx"))

soorten_lijst <- excel_data %>% 
  filter(Script == "Simpel") %>% 
  filter(Voerstreek == 1) %>% 
  pull(Soort) %>% 
  unique() %>% 
  na.omit() %>%
  # NIEUW: Filter enkel de soorten die een waarnemingenbestand hebben
  keep(~ heeft_waarnemingen_bestand(.x, MAP_WAARNEMINGEN))

message("\n=> Aantal geselecteerde simpele soorten voor VS-run (met waarnemingenfile): ", length(soorten_lijst))


# ------------------------------------------------------------------------------
# 3. RUN HET 'SIMPEL_VS' SCRIPT PER SOORT (Knitten met parameter)
# ------------------------------------------------------------------------------
draai_leefgebied_model <- function(huidige_soort, script_pad) {
  output_file_name <- paste0("VS_", huidige_soort, ".html")
  
  # --- SLIMME CHECK ---
  if (!moet_knitten(output_file_name, OUTPUT_DIR)) {
    message("   ⏭️  OVERSLAAN: ", output_file_name, " bestaat al en is vandaag al gemaakt.")
    return(data.frame(
      Item = huidige_soort,
      Type = "Simpele Soort (Normaal)",
      Status = "GEKOZEN_OVERSLAAN",
      Fout = "Reeds vandaag geknit",
      stringsAsFactors = FALSE
    ))
  }
  
  res_row <- data.frame(
    Item = huidige_soort,
    Type = "Simpele Soort (Normaal)",
    Status = "SUCCES",
    Fout = "Geen",
    stringsAsFactors = FALSE
  )
  
  message("   -> Knitten voor soort: ", toupper(huidige_soort))
  
  tryCatch({
    callr::r(
      function(input, output_file, output_dir, soort_invoer) {
        rmarkdown::render(
          input = input,
          output_file = output_file,
          output_dir = output_dir,
          params = list(soort_invoer = soort_invoer),
          quiet = TRUE
        )
      },
      args = list(
        input = script_pad,
        output_file = output_file_name,
        output_dir = OUTPUT_DIR,
        soort_invoer = huidige_soort
      )
    )
  }, error = function(e) {
    message("   ❌ FOUTMELDING bij ", huidige_soort)
    res_row$Status <<- "CRASH"
    res_row$Fout <<- e$message
  })
  
  return(res_row)
}

if (length(soorten_lijst) > 0) {
  soorten_logboek <- purrr::map_dfr(
    soorten_lijst, 
    ~draai_leefgebied_model(huidige_soort = .x, script_pad = soorten_script_tv_pad)
  )
  eind_logboek    <- bind_rows(bind_rows(algemeen_logboek), soorten_logboek)
} else {
  eind_logboek    <- bind_rows(algemeen_logboek)
  message("Waarschuwing: Geen simpele soorten gevonden die voldoen aan de filters.")
}


# ------------------------------------------------------------------------------
# 4. LOGBESTAND WEGSCHRIJVEN & SAMENVATTING
# ------------------------------------------------------------------------------
log_file_path <- file.path(OUTPUT_DIR, paste0("Logboek_VS_NormaalRun_", format(Sys.time(), "%Y%m%d_%H%M"), ".csv"))
readr::write_excel_csv(eind_logboek, log_file_path)

aantal_crashes <- sum(eind_logboek$Status == "CRASH")
aantal_skipped <- sum(eind_logboek$Status == "GEKOZEN_OVERSLAAN")
aantal_succes  <- sum(eind_logboek$Status == "SUCCES")

message("\n==================================================")
message(" MASTER RUN VS NORMAAL COMPLEET!                  ")
message(" Totaal onderdelen geëvalueerd: ", nrow(eind_logboek))
message(" Vandaag nieuw succesvol geknit: ", aantal_succes)
message(" Vandaag overgeslagen (al up-to-date): ", aantal_skipped)
message(" Gecrasht:                     ", aantal_crashes)
message(" Logboek opgeslagen als:       ", log_file_path)
message("==================================================")
