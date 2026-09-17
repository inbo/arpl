# ==============================================================================
# ARPL ACTUEEL MAKER - GENERATOR PER GEBIED EN SCENARIO
# AUTEUR: Bert Van Hecke
# ==============================================================================

library(here)
library(dplyr)
library(readxl)
library(readr)
library(purrr)

# ------------------------------------------------------------------------------
# 0. CONTROLEER/STEL PARAMETERS IN (MET DEFAULT FALLBACK)
# ------------------------------------------------------------------------------
# Als het script LOS wordt gedraaid, gebruikt het deze DEFAULT scenarioselectie:
DEFAULT_SCENARIO_SELECTIE <- list(
  De_Maten              = "DM_Scenario_BWK_2025.rds",
  Heesbossen            = "HB_Scenario_BWK_2025.rds",
  Kalmthoutse_Heide     = "KH_Scenario_BWK_2025.rds",
  Mechelse_Heide        = "MH_Scenario_BWK_2025.rds",
  Turnhouts_Vennegebied = "TV_Scenario_BWK_2025.rds",
  Voerstreek            = "VS_Scenario_BWK_2025.rds"
)

# Gebruik de meegegeven SCENARIO_SELECTIE uit het WorkQueue script,
# of val terug op de DEFAULT_SCENARIO_SELECTIE als deze los wordt gerund.
if (!exists("SCENARIO_SELECTIE") || !is.list(SCENARIO_SELECTIE)) {
  message("ℹ️ Geen meegegeven SCENARIO_SELECTIE gevonden. Default scenario's worden gebruikt.")
  actieve_scenarios <- DEFAULT_SCENARIO_SELECTIE
} else {
  message("✅ Dynamische SCENARIO_SELECTIE overgenomen uit WorkQueue script.")
  actieve_scenarios <- SCENARIO_SELECTIE
}

gebieden_info <- list(
  De_Maten              = list(code = "DM", col = "De_Maten"),
  Heesbossen            = list(code = "HB", col = "Heesbossen"),
  Kalmthoutse_Heide     = list(code = "KH", col = "Kalmthoutse_Heide"),
  Mechelse_Heide        = list(code = "MH", col = "Mechelse_Heide"),
  Turnhouts_Vennegebied = list(code = "TV", col = "Turnhouts_Vennegebied"),
  Voerstreek            = list(code = "VS", col = "Voerstreek")
)

message("==================================================")
message(" STARTEN ARPL & ACTUEEL MAKER VOOR ALLE GEBIEDEN")
message("==================================================")

# ------------------------------------------------------------------------------
# 1. LUS OVER ALLE GESELECTEERDE GEBIEDEN EN SCENARIO'S
# ------------------------------------------------------------------------------
for (gb_naam in names(actieve_scenarios)) {
  rds_naam <- actieve_scenarios[[gb_naam]]
  info     <- gebieden_info[[gb_naam]]
  
  if (is.null(info)) {
    warning("⚠️ Onbekend gebied '", gb_naam, "' overgeslagen.")
    next
  }
  
  rds_pad <- here("data/input/Scenario_rds", rds_naam)
  
  if (!file.exists(rds_pad)) {
    warning("⚠️ Scenario RDS niet gevonden voor ", gb_naam, ": ", rds_pad)
    next
  }
  
  # Bepaal de opgeschoonde scenarionaam voor paden (bijv. 'bosbehoudss_ss31fix_tvg_vrij_2026')
  huidig_scenario <- gsub(paste0("^", info$code, "_Scenario_|.rds$"), "", basename(rds_pad))
  
  message("\n--------------------------------------------------")
  message("-> Verwerken gebied  : ", gb_naam, " (Code: ", info$code, ")")
  message("   Scenario RDS     : ", basename(rds_pad))
  message("   Scenario Naam    : ", huidig_scenario)
  message("--------------------------------------------------")
  
  # ----------------------------------------------------------------------------
  # 2. HIER PLAATS JE JE BESTAANDE VERWERKINGSLOGICA PER GEBIED
  # ----------------------------------------------------------------------------
  # Voorbeeld van paden die nu dynamisch opgebouwd worden per gebied/scenario:
  input_dir  <- here("data/output", gb_naam, "Scenario_output", huidig_scenario)
  output_dir <- here("data/output", gb_naam, "ARPL_Actueel_Kaarten", huidig_scenario)
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  tryCatch({
    # --- PLAATS HIER DE CORE RUN-LOGICA VAN JE ARPL_Actueel_Maker SCRIPT ---
    # Je kunt hier gebruik maken van:
    # - rds_pad          (Volledig pad naar het scenario .rds bestand)
    # - gb_naam          (bijv. "De_Maten")
    # - info$code        (bijv. "DM")
    # - huidig_scenario  (bijv. "bosbehoudss_ss31fix_tvg_vrij_2026")
    # - input_dir / output_dir
    
    message("   ✅ ARPL & Actueel kaarten gegenereerd voor ", gb_naam)
    
  }, error = function(e) {
    message("   ❌ FOUT bij verwerken van ", gb_naam, ": ", e$message)
  })
  
  gc(verbose = FALSE) # Ruim RAM op per gebied
}

message("\n==================================================")
message(" 🎉 ARPL & ACTUEEL MAKER COMPLEET VOOR ALLE GEBIEDEN")
message("==================================================")
