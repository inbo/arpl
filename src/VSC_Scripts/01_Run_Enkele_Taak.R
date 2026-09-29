# ==============================================================================
# UNIVERSELE WORKER EXECUTOR VOOR 1 ATOMAIRE TAAK
# ==============================================================================
library(here)
library(dplyr)
library(readr)

proj_root <- here::here()

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("❌ Geen taakindex meegegeven aan R-worker!")
taak_idx <- as.numeric(args[1])

takenlijst_pad <- file.path(proj_root, "data/temp/globale_takenlijst.rds")
if (!file.exists(takenlijst_pad)) stop("❌ Takenlijst niet gevonden in data/temp!")

globale_takenlijst <- readRDS(takenlijst_pad)
if (taak_idx > length(globale_takenlijst)) {
  cat("⚠️ Taakindex buiten bereik van takenlijst. Mogelijk al verwerkt.\n")
  q(status = 0)
}

taak <- globale_takenlijst[[taak_idx]]

# Geheugen- en tempdir-instellingen voor Terra
worker_temp <- file.path(proj_root, "data/temp_workers", paste0("gis_job_", Sys.getpid(), "_", sample(1000:9999, 1)))
dir.create(worker_temp, showWarnings = FALSE, recursive = TRUE)

if (requireNamespace("terra", quietly = TRUE)) {
  terra::terraOptions(memfrac = 0.40, tempdir = worker_temp, verbose = FALSE)
}

schoon_naam_op <- function(x) {
  if (is.null(x) || length(x) == 0) return(character(0))
  x %>% as.character() %>% basename() %>% tolower() %>% 
    gsub("\\.tif$|\\.rds$|\\.html$|\\.r$", "", .) %>%
    gsub("^habitat_werkelijke_oppervlaktes_|^habitat_maximale_potentie_|^id_netwerken_|^rapport_|^scenario_", "", .) %>% 
    gsub("^(tv|dm|hb|kh|mh|vs)_", "", .) %>% 
    gsub("[^a-z0-9]", "", .) %>% trimws()
}

# Omgevingsvariabelen exact opzetten volgens jouw originele script
script_env <- new.env(parent = globalenv())
script_env$SCENARIO_RDS_PAD  <- taak$ScenarioRDS
script_env$scenario_rds_path <- taak$ScenarioRDS
script_env$GEBIED_NAAM       <- taak$Gebied
script_env$GEBIED_CODE       <- taak$GebiedCode
script_env$huidig_scenario   <- taak$HuidigScenario
script_env$output_dir        <- file.path(proj_root, "data/output", taak$Gebied, "HTML_Rapporten_Scenario", taak$HuidigScenario)

if (!dir.exists(script_env$output_dir)) dir.create(script_env$output_dir, recursive = TRUE)

if (taak$Type == "Simpele_Soort") {
  huidige_soort <- taak$Soort
  soort_zonder_spaties <- gsub(" ", "", huidige_soort, fixed = TRUE)
  soort_clean          <- schoon_naam_op(huidige_soort)
  
  script_env$huidige_soort       <- soort_zonder_spaties
  script_env$soort_invoer        <- soort_zonder_spaties
  script_env$soort               <- soort_zonder_spaties
  script_env$soort_naam          <- soort_zonder_spaties
  script_env$SOORT               <- soort_zonder_spaties
  script_env$huidige_soort_clean <- soort_clean
  script_env$soort_orig          <- huidige_soort
  
  cat(sprintf("-> Running Simpele Soort: %s [%s] voor %s\n", huidige_soort, soort_zonder_spaties, taak$Gebied))
} else {
  cat(sprintf("-> Running Uniek Script: %s voor %s\n", basename(taak$ScriptPad), taak$Gebied))
}

old_wd <- setwd(dirname(taak$ScriptPad))
on.exit({
  setwd(old_wd)
  unlink(worker_temp, recursive = TRUE, force = TRUE)
})

# Voer het daadwerkelijke script uit via source
source(taak$ScriptPad, local = script_env)

if (requireNamespace("terra", quietly = TRUE)) {
  terra::tmpFiles(current = TRUE, orphan = TRUE, old = TRUE, remove = TRUE)
}
gc(verbose = FALSE)

cat(sprintf("✅ SUCCESVOL AFGEROND: Taak %d (%s - %s)\n", taak_idx, taak$Gebied, ifelse(is.na(taak$Soort), basename(taak$ScriptPad), taak$Soort)))
