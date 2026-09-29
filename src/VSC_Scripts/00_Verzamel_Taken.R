# ==============================================================================
# DYNAMISCHE TAKEN-VERZAMELAAR OVER ALLE 6 GEBIEDEN (MET SKIP-CHECK)
# ==============================================================================
library(here)
library(readxl)
library(dplyr)
library(purrr)

proj_root <- here::here()

SCENARIO_SELECTIE <- list(
  De_Maten              = "DM_Scenario_BWK_2025.rds",
  Heesbossen            = "HB_Scenario_BWK_2025.rds",
  Kalmthoutse_Heide     = "KH_Scenario_BWK_2025.rds",
  Mechelse_Heide        = "MH_Scenario_BWK_2025.rds",
  Turnhouts_Vennegebied = "TV_Scenario_BWK_2025.rds",
  Voerstreek            = "VS_Scenario_BWK_2025.rds"
)

gebieden_info <- list(
  De_Maten              = list(code = "DM", col = "De_Maten",            simpel_script = "Scenario_DM_Leefgebieden_Simpel.R"),
  Heesbossen            = list(code = "HB", col = "Heesbossen",          simpel_script = "Scenario_HB_Leefgebieden_Simpel.R"),
  Kalmthoutse_Heide     = list(code = "KH", col = "Kalmthoutse_Heide",   simpel_script = "Scenario_KH_Leefgebieden_Simpel.R"),
  Mechelse_Heide        = list(code = "MH", col = "Mechelse_Heide",      simpel_script = "Scenario_MH_Leefgebieden_Simpel.R"),
  Turnhouts_Vennegebied = list(code = "TV", col = "Turnhouts_Vennegebied", simpel_script = "Scenario_TV_Leefgebieden_Simpel.R"),
  Voerstreek            = list(code = "VS", col = "Voerstreek",          simpel_script = "Scenario_VS_Leefgebieden_Simpel.R")
)

schoon_naam_op <- function(x) {
  if (is.null(x) || length(x) == 0) return(character(0))
  x %>% 
    as.character() %>%
    basename() %>% 
    tolower() %>% 
    gsub("\\.tif$|\\.rds$|\\.html$|\\.r$", "", .) %>%
    gsub("^habitat_werkelijke_oppervlaktes_|^habitat_maximale_potentie_|^id_netwerken_|^rapport_|^scenario_", "", .) %>% 
    gsub("^(tv|dm|hb|kh|mh|vs)_", "", .) %>% 
    gsub("[^a-z0-9]", "", .) %>% 
    trimws()
}

excel_pad <- file.path(proj_root, "data/input/Excel_files/Soortenlijst_Maatwerkgebieden_Gefilterd.xlsx")
excel_data <- if (file.exists(excel_pad)) read_excel(excel_pad) else NULL
if (!is.null(excel_data)) colnames(excel_data) <- tolower(trimws(colnames(excel_data)))

globale_takenlijst <- list()

for (gb_naam in names(gebieden_info)) {
  info <- gebieden_info[[gb_naam]]
  map_scripts <- file.path(proj_root, "src", gb_naam, "Scripts_Scenario")
  if (!dir.exists(map_scripts)) next
  
  gekozen_rds_naam <- SCENARIO_SELECTIE[[gb_naam]]
  scenario_rds      <- file.path(proj_root, "data/input/Scenario_rds", gekozen_rds_naam)
  if (!file.exists(scenario_rds)) next
  
  huidig_scenario <- gsub("^.*_Scenario_|^Scenario_|_wv\\.rds$|\\.rds$", "", basename(scenario_rds), ignore.case = TRUE)
  map_werkelijk   <- file.path(proj_root, "data/output", gb_naam, "Rasters_Soorten", huidig_scenario, "02_Werkelijke_Oppervlaktes")
  
  bestaande_soorten_clean <- if (dir.exists(map_werkelijk)) {
    aangemaakte <- list.files(map_werkelijk, pattern = "\\.tif$", full.names = TRUE)
    geldige     <- aangemaakte[file.size(aangemaakte) > 1000]
    geldige %>% schoon_naam_op() %>% unique()
  } else {
    character(0)
  }
  
  alle_scripts <- list.files(path = map_scripts, pattern = "\\.r$", full.names = TRUE, ignore.case = TRUE)
  soorten_script_pad <- alle_scripts[grep(paste0(info$simpel_script, "$"), alle_scripts, ignore.case = TRUE)]
  soorten_script_pad <- if (length(soorten_script_pad) > 0) soorten_script_pad[1] else file.path(map_scripts, info$simpel_script)
  
  # 1. Unieke R-scripts toevoegen
  unieke_scripts <- setdiff(alle_scripts, soorten_script_pad)
  for (script in unieke_scripts) {
    if (!(schoon_naam_op(script) %in% bestaande_soorten_clean)) {
      globale_takenlijst[[length(globale_takenlijst) + 1]] <- list(
        Type = "Uniek_Script", Gebied = gb_naam, GebiedCode = info$code,
        ScriptPad = script, ScenarioRDS = scenario_rds, HuidigScenario = huidig_scenario, Soort = NA
      )
    }
  }
  
  # 2. Simpele soorten uit Excel toevoegen
  if (!is.null(excel_data) && file.exists(soorten_script_pad)) {
    col_gebied <- tolower(gb_naam)
    if (col_gebied %in% colnames(excel_data)) {
      naam_kolom <- if ("soort" %in% colnames(excel_data)) "soort" else if ("nederlandse naam" %in% colnames(excel_data)) "nederlandse naam" else colnames(excel_data)[2]
      soorten <- excel_data %>% 
        mutate(automatisch_clean = tolower(trimws(as.character(automatisch))), gebied_val = suppressWarnings(as.numeric(.data[[col_gebied]]))) %>% 
        filter(automatisch_clean == "ja", gebied_val == 1) %>% 
        pull(!!sym(naam_kolom)) %>% trimws() %>% unique() %>% na.omit()
      
      for (srt in soorten) {
        if (!(schoon_naam_op(srt) %in% bestaande_soorten_clean)) {
          globale_takenlijst[[length(globale_takenlijst) + 1]] <- list(
            Type = "Simpele_Soort", Gebied = gb_naam, GebiedCode = info$code,
            ScriptPad = soorten_script_pad, ScenarioRDS = scenario_rds, HuidigScenario = huidig_scenario, Soort = srt
          )
        }
      }
    }
  }
}

dir.create(file.path(proj_root, "data/temp"), showWarnings = FALSE, recursive = TRUE)
saveRDS(globale_takenlijst, file.path(proj_root, "data/temp/globale_takenlijst.rds"))

message("==================================================")
message(" TAKENENQUÊTE COMPLEET")
message(" Totaal te verwerken taken na skip-check: ", length(globale_takenlijst))
message("==================================================")