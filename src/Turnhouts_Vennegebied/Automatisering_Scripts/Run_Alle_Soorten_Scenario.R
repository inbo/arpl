# ==============================================================================
# MASTER AUTOMATISERINGSSCRIPT VOOR R-SCRIPTS (.R) - PER GEBIED (INCL. TIMER & ARPL)
# AUTEUR: Bert Van Hecke
# ==============================================================================

# Start de globale timer voor het volledige script
totale_start_tijd <- Sys.time()

library(purrr)
library(readxl)
library(dplyr)
library(readr)
library(here)

# Dwing terra tot strikt geheugenbeheer
if (requireNamespace("terra", quietly = TRUE)) {
  terra::terraOptions(memfrac = 0.2, tempdir = tempdir(), verbose = FALSE)
}

# ------------------------------------------------------------------------------
# 0. GLOBALE INSTELLINGEN
# ------------------------------------------------------------------------------
GEBIED_NAAM      <- "Turnhouts_Vennegebied"
GEBIED_CODE      <- "TV"

SCENARIO_RDS_PAD <- here("data/input/Scenario_rds/TV_Scenario_BWK_2025.rds")
HUIDIG_SCENARIO  <- gsub("^TV_Scenario_|.rds$", "", basename(SCENARIO_RDS_PAD), ignore.case = TRUE)

MAP_SCRIPTS      <- here("src", GEBIED_NAAM, "Scripts_Scenario")
OUTPUT_DIR       <- here("data/output", GEBIED_NAAM, "HTML_Rapporten_Scenario", HUIDIG_SCENARIO)

if(!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# ------------------------------------------------------------------------------
# DETECTIE VAN BESTAANDE RASTERS (02_Werkelijke_Oppervlaktes)
# ------------------------------------------------------------------------------
map_werkelijk <- here("data/output", GEBIED_NAAM, "Rasters_Soorten", HUIDIG_SCENARIO, "02_Werkelijke_Oppervlaktes")

if (dir.exists(map_werkelijk)) {
  bestaande_rasters <- list.files(map_werkelijk, pattern = "\\.tif$", full.names = FALSE, ignore.case = TRUE)
} else {
  bestaande_rasters <- character(0)
}

# Hulpfunctie: stript prefixen, extensies, SPATIES, UNDERSCORES EN STREEPJES voor een 100% consistente match
schoon_naam_op <- function(x) {
  if (is.null(x) || length(x) == 0) return(character(0))
  
  x %>% 
    as.character() %>%
    basename() %>% 
    tolower() %>% 
    # Strippen van bekende extensies
    # NIEUWE REGEL (behoudt _wv als onderdeel van de unieke naam):
    gsub("\\.tif$|\\.rds$|\\.html$|\\.r$", "", .) %>%
    # Strippen van bekende prefixen
    gsub("^habitat_werkelijke_oppervlaktes_|^habitat_maximale_potentie_|^id_netwerken_|^rapport_|^scenario_", "", .) %>% 
    gsub("^(tv|dm|hb|kh|mh|vs)_", "", .) %>% 
    # Verwijder ALLE niet-alfanumerieke tekens (spaties, underscores, streepjes)
    gsub("[^a-z0-9]", "", .) %>% 
    trimws()
}

bestaande_soorten_clean <- bestaande_rasters %>% 
  schoon_naam_op() %>% 
  unique()

algemeen_logboek <- list()
message("==================================================")
message(" START MASTER RUN VOOR SCENARIO: ", toupper(HUIDIG_SCENARIO))
message(" Gebied: ", GEBIED_NAAM)
message(" Bestaande rasters (.tif) reeds gedetecteerd: ", length(bestaande_soorten_clean))
message("==================================================")

# ------------------------------------------------------------------------------
# 1. DETECTEER EN RUN DE UNIEKE R-SCRIPTS (MET SLIMME SKIP-CHECK)
# ------------------------------------------------------------------------------
alle_scripts <- list.files(path = MAP_SCRIPTS, pattern = "\\.r$", full.names = TRUE, ignore.case = TRUE)

soorten_script_pad <- alle_scripts[grep("Scenario_TV_Leefgebieden_Simpel\\.r$", alle_scripts, ignore.case = TRUE)]

if (length(soorten_script_pad) == 0) {
  soorten_script_pad <- file.path(MAP_SCRIPTS, "Scenario_TV_Leefgebieden_Simpel.R")
} else {
  soorten_script_pad <- soorten_script_pad[1]
}

alle_unieke_scripts <- setdiff(alle_scripts, soorten_script_pad)

unieke_scripts <- keep(alle_unieke_scripts, function(script_pad) {
  soort_uit_script <- schoon_naam_op(script_pad)
  is_al_gedraaid     <- soort_uit_script %in% bestaande_soorten_clean
  return(!is_al_gedraaid)
})

message("=> Gedetecteerde unieke R-scripts totaal : ", length(alle_unieke_scripts))
message("=> Reeds verwerkt (overgeslagen)         : ", length(alle_unieke_scripts) - length(unieke_scripts))
message("=> Nog uit te voeren unieke scripts      : ", length(unieke_scripts))
if (length(unieke_scripts) > 0) print(basename(unieke_scripts))
message("--------------------------------------------------")

for(script in unieke_scripts) {
  bestandsnaam <- basename(script)
  tijd_script_start <- Sys.time()
  message("=> Uitvoeren van script via source(): ", bestandsnaam)
  
  tryCatch({
    script_env <- new.env(parent = globalenv())
    script_env$SCENARIO_RDS_PAD  <- SCENARIO_RDS_PAD
    script_env$scenario_rds_path <- SCENARIO_RDS_PAD
    script_env$output_dir        <- OUTPUT_DIR
    script_env$huidig_scenario   <- HUIDIG_SCENARIO
    script_env$GEBIED_NAAM       <- GEBIED_NAAM
    script_env$GEBIED_CODE       <- GEBIED_CODE
    
    source(script, local = script_env)
    
    tijd_script_eind <- Sys.time()
    duur_min <- round(as.numeric(difftime(tijd_script_eind, tijd_script_start, units = "mins")), 2)
    
    algemeen_logboek[[bestandsnaam]] <- data.frame(
      Item = bestandsnaam, Type = "Soort Script", Status = "SUCCES", 
      Duurtijd_Min = duur_min, Fout = "Geen", stringsAsFactors = FALSE
    )
  }, error = function(e) {
    message("❌ FOUT bij: ", bestandsnaam, " - Fout: ", e$message)
    algemeen_logboek[[bestandsnaam]] <- data.frame(
      Item = bestandsnaam, Type = "Soort Script", Status = "CRASH", 
      Duurtijd_Min = NA, Fout = e$message, stringsAsFactors = FALSE
    )
  })
  
  # Agressieve opruiming van RAM en tijdelijke bestanden na elk uniek script
  if (requireNamespace("terra", quietly = TRUE)) {
    terra::tmpFiles(current = TRUE, orphan = TRUE, old = TRUE, remove = TRUE)
  }
  gc(verbose = FALSE)
}

# ------------------------------------------------------------------------------
# 2. SIMPELE SOORTEN VIA GENERIEK R-SCRIPT (MET GEFIXTE SPATIE-AFHANDELING)
# ------------------------------------------------------------------------------
if (file.exists(soorten_script_pad)) {
  excel_data <- read_excel(here("data/input/Excel_files/Soortenlijst_Maatwerkgebieden_Gefilterd.xlsx"))
  colnames(excel_data) <- tolower(trimws(colnames(excel_data)))
  
  naam_kolom <- if ("soort" %in% colnames(excel_data)) {
    "soort"
  } else if ("nederlandse naam" %in% colnames(excel_data)) {
    "nederlandse naam"
  } else {
    colnames(excel_data)[2]
  }
  
  col_gebied <- tolower(GEBIED_NAAM)
  
  if (col_gebied %in% colnames(excel_data)) {
    alle_simpele_soorten <- excel_data %>% 
      mutate(
        automatisch_clean = tolower(trimws(as.character(automatisch))),
        gebied_val        = suppressWarnings(as.numeric(.data[[col_gebied]]))
      ) %>% 
      filter(automatisch_clean == "ja") %>%                     
      filter(gebied_val == 1) %>%        
      pull(!!sym(naam_kolom)) %>%                                     
      trimws() %>% 
      unique() %>%                                     
      na.omit()                                         
  } else {
    alle_simpele_soorten <- character(0)
  }
  
  soorten_lijst <- keep(alle_simpele_soorten, function(srt) {
    srt_clean <- schoon_naam_op(srt)
    is_al_gedraaid <- srt_clean %in% bestaande_soorten_clean
    return(!is_al_gedraaid)
  })
  
  message("\n--------------------------------------------------")
  message(" Totaal geselecteerde simpele soorten in Excel : ", length(alle_simpele_soorten))
  message(" Reeds verwerkt (overgeslagen)                 : ", length(alle_simpele_soorten) - length(soorten_lijst))
  message(" Nog uit te voeren simpele soorten             : ", length(soorten_lijst))
  message("--------------------------------------------------")
  if (length(soorten_lijst) > 0) print(soorten_lijst)
  
  draai_leefgebied_scenario_model <- function(huidige_soort) {
    tijd_soort_start <- Sys.time()
    res_row <- data.frame(
      Item = huidige_soort,
      Type = "Simpele Soort",
      Status = "SUCCES",
      Duurtijd_Min = NA,
      Fout = "Geen",
      stringsAsFactors = FALSE
    )
    
    # Maak expliciet variabelen aan zonder spaties en opgeschoond
    soort_zonder_spaties <- gsub(" ", "", huidige_soort, fixed = TRUE)
    soort_clean          <- schoon_naam_op(huidige_soort)
    
    message("   -> Starten met simulatie voor: ", toupper(HUIDIG_SCENARIO), " - ", toupper(soort_zonder_spaties))
    tryCatch({
      script_env <- new.env(parent = globalenv())
      
      # Alle standaard soort-variabelen krijgen nu de waarde ZONDER spatie (e.g. "bruinesnavelbies"):
      script_env$huidige_soort        <- soort_zonder_spaties
      script_env$soort_invoer         <- soort_zonder_spaties
      script_env$soort                <- soort_zonder_spaties
      script_env$soort_naam           <- soort_zonder_spaties
      script_env$SOORT                <- soort_zonder_spaties
      script_env$huidige_soort_clean  <- soort_clean
      
      # Bewaar originele naam met spaties optioneel voor titels/labels in rapporten:
      script_env$soort_orig           <- huidige_soort
      
      script_env$SCENARIO_RDS_PAD  <- SCENARIO_RDS_PAD
      script_env$scenario_rds_path <- SCENARIO_RDS_PAD
      script_env$output_dir        <- OUTPUT_DIR
      script_env$huidig_scenario   <- HUIDIG_SCENARIO
      script_env$GEBIED_NAAM       <- GEBIED_NAAM
      script_env$GEBIED_CODE       <- GEBIED_CODE
      
      source(soorten_script_pad, local = script_env)
      
      tijd_soort_eind <- Sys.time()
      res_row$Duurtijd_Min <- round(as.numeric(difftime(tijd_soort_eind, tijd_soort_start, units = "mins")), 2)
      
    }, error = function(e) {
      message("   ❌ FOUTMELDING bij ", huidige_soort, ": ", e$message)
      res_row$Status <<- "CRASH"
      res_row$Fout <<- e$message
    })
    
    if (requireNamespace("terra", quietly = TRUE)) {
      terra::tmpFiles(current = TRUE, orphan = TRUE, old = TRUE, remove = TRUE)
    }
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
aantal_crashes <- if("Status" %in% colnames(eind_logboek)) sum(eind_logboek$Status == "CRASH") else 0

# ------------------------------------------------------------------------------
# 4. FASE 2: RUN FINALE MAKER SCRIPT VOOR DIT GEBIED (ACTIEF)
# ------------------------------------------------------------------------------
maker_script <- list.files(
  path = here("src"), 
  pattern = "^ARPL_Actueel_Maker\\.r$", 
  full.names = TRUE, 
  recursive = TRUE, 
  ignore.case = TRUE
)[1]

if (!is.null(maker_script) && file.exists(maker_script)) {
  message("\n==================================================")
  message(" STARTEN FASE 2: ARPL & ACTUEEL MAKER")
  message("==================================================")
  
  tryCatch({
    maker_env <- new.env(parent = globalenv())
    scen_lijst <- list()
    scen_lijst[[GEBIED_NAAM]] <- basename(SCENARIO_RDS_PAD)
    maker_env$SCENARIO_SELECTIE <- scen_lijst
    
    source(maker_script, local = maker_env)
    message("✅ ARPL_Actueel_Maker.R succesvol afgerond voor ", GEBIED_NAAM)
    
  }, error = function(e) {
    message("❌ FOUT bij uitvoeren van ARPL_Actueel_Maker.R: ", e$message)
  })
} else {
  warning("⚠️ 'ARPL_Actueel_Maker.R' kon niet worden gevonden in src/.")
}

# ------------------------------------------------------------------------------
# 5. EINDRAPPORTAGE EN TOTALE DUURTIJD
# ------------------------------------------------------------------------------
totale_eind_tijd <- Sys.time()
totale_duur_min  <- round(as.numeric(difftime(totale_eind_tijd, totale_start_tijd, units = "mins")), 2)

message("\n==================================================")
message(" 🎉 RUN COMPLEET VOOR ", toupper(HUIDIG_SCENARIO))
message(" Nieuwe onderdelen gedraaid : ", nrow(eind_logboek))
message(" Succesvol                    : ", nrow(eind_logboek) - aantal_crashes)
message(" Gecrasht                     : ", aantal_crashes)
message(" Totale duurtijd              : ", totale_duur_min, " minuten")
message("==================================================")
