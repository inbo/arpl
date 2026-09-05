# ==============================================================================
# MASTER AUTOMATISERINGSSCRIPT VOOR REGULIERE TV-SCRIPTS (KNITTEN)
# AUTEUR: Bert Van Hecke
# UPDATE: Slimme skip-logica + CSV-check voor waarnemingen van simpele soorten
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
MAP_SCRIPTS_TV     <- here("src/Turnhouts_Vennegebied/Scripts_TV/")
OUTPUT_DIR         <- here("data/output/Turnhouts_Vennegebied/HTML_Rapporten_Soorten")
MAP_WAARNEMINGEN   <- here("data/input/Waarnemingen_Soorten/Turnhouts_Vennegebied")

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# Initialiseer de logboek lijst (cruciaal om crashes in de unieke scripts loop te voorkomen)
algemeen_logboek <- list()

# Hulpfunctie: Controleert of er een waarnemingen-CSV bestaat voor de specifieke soort
heeft_waarnemingen_bestand <- function(soort_naam, map_waarnemingen) {
  
  # 1. Haal alle CSV-bestanden op in de waarnemingenmap
  alle_csvs <- list.files(path = map_waarnemingen, pattern = "\\.csv$", full.names = FALSE)
  if (length(alle_csvs) == 0) return(FALSE)
  
  # 2. Opschonen van de soortnaam (verwijder _wv als die aanwezig is)
  schoon_ras <- sub("_wv$", "", soort_naam)
  
  # 3. Verwachte zoeknaam opbouwen (bijv. Waarnemingen_MiddelsteBonteSpecht)
  zoek_naam <- paste0("Waarnemingen_", schoon_ras)
  
  # 4. Maak de vergelijking robuust: negeer hoofdletters, spaties en underscores
  # Dit zorgt ervoor dat "Middelste_Bonte_Specht" en "MiddelsteBonteSpecht" beiden matchen
  normaliseer <- function(tekst) {
    tolower(gsub("[ _]", "", tekst))
  }
  
  zoek_naam_norm <- normaliseer(zoek_naam)
  csvs_norm      <- normaliseer(tools::file_path_sans_ext(alle_csvs))
  
  # 5. Check of de genormaliseerde naam in de lijst met CSV's staat
  return(zoek_naam_norm %in% csvs_norm)
}


# ------------------------------------------------------------------------------
# HULPFUNCTIE: Moet dit bestand opnieuw geknit worden? (Check laatste 5 dagen)
# ------------------------------------------------------------------------------
moet_knitten <- function(bestandsnaam, map) {
  pad_naar_bestand <- file.path(map, bestandsnaam)
  
  # Check 1: Bestaat het bestand überhaupt?
  if (!file.exists(pad_naar_bestand)) {
    return(TRUE) # Bestaat niet, dus we moeten knitten
  }
  
  # Check 2: Is het in de afgelopen 5 dagen aangemaakt/gewijzigd?
  info            <- file.info(pad_naar_bestand)
  wijzigingsdatum <- as.Date(info$mtime)
  grens_datum     <- Sys.Date() - 5
  
  if (wijzigingsdatum >= grens_datum) {
    return(FALSE) # Bestaat al en is max 5 dagen oud -> overslaan!
  } else {
    return(TRUE)  # Bestaat wel, maar is ouder dan 5 dagen -> opnieuw knitten
  }
}


# ------------------------------------------------------------------------------
# 1. DETECTEER EN KNIT DE UNIEKE TV-SCRIPTS (Runt ALLE unieke scripts)
# ------------------------------------------------------------------------------
alle_tv_scripts       <- list.files(path = MAP_SCRIPTS_TV, pattern = "\\.Rmd$", full.names = TRUE)
soorten_script_tv_pad <- file.path(MAP_SCRIPTS_TV, "TV_Leefgebieden_Simpel.Rmd")
unieke_tv_scripts     <- setdiff(alle_tv_scripts, soorten_script_tv_pad)

message("=> Gedetecteerde unieke TV-scripts om te knitten (Aantal: ", length(unieke_tv_scripts), "):")
print(basename(unieke_tv_scripts))
message("--------------------------------------------------")

for(script in unieke_tv_scripts) {
  bestandsnaam <- basename(script)
  output_html  <- sub(".Rmd$", ".html", bestandsnaam)
  
  # --- SLIMME CHECK ---
  if (!moet_knitten(output_html, OUTPUT_DIR)) {
    message("⏭️  OVERSLAAN: ", output_html, " bestaat al en is de afgelopen 5 dagen al gemaakt.")
    algemeen_logboek[[bestandsnaam]] <- data.frame(
      Item = bestandsnaam, Type = "Uniek TV Script", Status = "GEKOZEN_OVERSLAAN", Fout = "Reeds vandaag geknit", stringsAsFactors = FALSE
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
# 2. TREK DE SOORTENLIJST DYNAMISCH UIT DE EXCEL & CHECK CSV WAARNEMINGEN
# ------------------------------------------------------------------------------
excel_data <- read_excel(here("data/input/Excel_files/Soortenlijst_Maatwerkgebieden.xlsx"))

soorten_lijst <- excel_data %>% 
  # Filter op automatische scripts ("ja" negeert eventuele hoofdletters via tolower)
  filter(tolower(Automatisch) == "ja") %>% 
  filter(Turnhouts_Vennegebied == 1) %>% 
  pull(Soort) %>% 
  unique() %>% 
  na.omit() %>%
  # ENKEL soorten behouden waarvoor een Waarnemingen_<Soort>.csv bestaat
  keep(~ heeft_waarnemingen_bestand(.x, MAP_WAARNEMINGEN))

message("\n=> Aantal geselecteerde simpele soorten met een .csv waarnemingenbestand: ", length(soorten_lijst))


# ------------------------------------------------------------------------------
# 3. RUN HET 'SIMPEL_TV' SCRIPT PER GEFILTERDE SOORT
# ------------------------------------------------------------------------------
draai_leefgebied_model <- function(huidige_soort, script_pad) {
  output_file_name <- paste0("TV_", huidige_soort, ".html")
  
  # --- SLIMME CHECK ---
  if (!moet_knitten(output_file_name, OUTPUT_DIR)) {
    message("   ⏭️  OVERSLAAN: ", output_file_name, " bestaat al en is de afgelopen 5 dagen al gemaakt.")
    return(data.frame(
      Item = huidige_soort,
      Type = "Simpele Soort (Normaal)",
      Status = "GEKOZEN_OVERSLAAN",
      Fout = "Reeds geknit de afgelopen 5 dagen",
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
  message("Waarschuwing: Geen simpele soorten gevonden die een waarnemingen-CSV hebben.")
}


# ------------------------------------------------------------------------------
# 4. LOGBESTAND WEGSCHRIJVEN & SAMENVATTING
# ------------------------------------------------------------------------------
log_file_path <- file.path(OUTPUT_DIR, paste0("Logboek_TV_NormaalRun_", format(Sys.time(), "%Y%m%d_%H%M"), ".csv"))
readr::write_excel_csv(eind_logboek, log_file_path)

aantal_crashes <- sum(eind_logboek$Status == "CRASH")
aantal_skipped <- sum(eind_logboek$Status == "GEKOZEN_OVERSLAAN")
aantal_succes  <- sum(eind_logboek$Status == "SUCCES")

message("\n==================================================")
message(" MASTER RUN TV NORMAAL COMPLEET!                  ")
message(" Totaal onderdelen geëvalueerd: ", nrow(eind_logboek))
message(" Vandaag nieuw succesvol geknit: ", aantal_succes)
message(" Vandaag overgeslagen (al up-to-date): ", aantal_skipped)
message(" Gecrasht:                      ", aantal_crashes)
message(" Logboek opgeslagen als:        ", log_file_path)
message("==================================================")
