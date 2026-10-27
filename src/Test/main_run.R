library(here)
library(tidyverse)
library(terra)
library(readxl)
library(data.table)

# 1. HELPERFUNCTIES LADEN
source(here::here("R/01_helper_functies.R"))

conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::select)
conflicted::conflicts_prefer(terra::intersect)

# Optimale HPC/RAM instellingen
terraOptions(memfrac = 0.7, tempdir = tempdir(), verbose = FALSE)

# ------------------------------------------------------------------------------
# 2. SOORT & SCENARIO SELEKTIE
# ------------------------------------------------------------------------------
if (exists("huidige_soort") && !is.null(huidige_soort)) {
  soort <- tolower(trimws(huidige_soort))
} else {
  soort <- "aardbeivlinder"  # Standaard testsoort
}

scenario_rds_path <- "data/input/Scenario_rds/TV_Scenario_BWK_2025.rds"
scenario_path     <- here::here(gsub("^([.][.]/)+", "", scenario_rds_path))
scen_volledig     <- basename(scenario_path)
scenario_naam     <- gsub("^TV_Scenario_|^Scenario_|.rds$", "", scen_volledig)

message(paste("🚀 Verwerken van soort:", toupper(soort), "voor HEEL VLAANDEREN binnen scenario:", scenario_naam))

# ------------------------------------------------------------------------------
# 3. METADATA LEZEN UIT EXCEL
# ------------------------------------------------------------------------------
df_excel <- read_excel(here::here("data/input/Excel_files/Soorten_parameters.xlsx"))

cfg <- df_excel %>%
  mutate(Soort_clean = tolower(trimws(Soort))) %>% 
  filter(Soort_clean == soort)

if(nrow(cfg) == 0) stop(paste("❌ Soort", soort, "niet gevonden in Excel!"))

oppervlakte_ha <- cfg$MinOpp_ha[1]
afstand_m      <- cfg$AfstandBiotopen_m[1]
buffer_m       <- cfg$Dispersiecap_m[1]

# Vinkjes voor optionele stappen
gebruik_drainage   <- if ("Gebruik_Drainage" %in% colnames(cfg)) isTRUE(cfg$Gebruik_Drainage[1]) else FALSE
gebruik_groenkaart <- if ("Gebruik_Groenkaart" %in% colnames(cfg)) isTRUE(cfg$Gebruik_Groenkaart[1]) else FALSE

# ------------------------------------------------------------------------------
# 4. MASTERGRID VLAANDEREN INLADEN
# ------------------------------------------------------------------------------
master_grid <- rast(here("data/input/Raster_Vlaanderen/Vlaanderen_MasterGrid_10m.tif"))[[1]]

# ------------------------------------------------------------------------------
# 5. DATA MATCHEN OP VLAAMS NIVEAU VIA DATA.TABLE (RECHTSTREEKS & SNEL)
# ------------------------------------------------------------------------------
df_nieuw <- read_csv(here("data/input/Excel_files/Resultaten_Totaal_Samengevoegd.csv"), show_col_types = FALSE)

resultaten_gegroepeerd <- df_nieuw %>%
  mutate(Soort_clean = tolower(trimws(Soort))) %>%
  filter(Soort_clean == soort) %>%
  group_by(Type) %>%
  nest(Data = c(Code, Match))

tabel_vlaanderen <- readRDS(scenario_path)
setDT(tabel_vlaanderen)
tabel_vlaanderen[, CODE := tolower(trimws(CODE))]

r_match_vlaanderen <- master_grid * NA
r_opp_vlaanderen   <- master_grid * NA

if(nrow(resultaten_gegroepeerd) > 0) {
  # Neem de codes voor deze soort op
  h_data <- resultaten_gegroepeerd$Data[[1]]
  exact_codes <- tolower(trimws(h_data$Code[h_data$Match == "exact"]))
  bevat_codes <- tolower(trimws(h_data$Code[h_data$Match == "bevat"]))
  
  if(length(bevat_codes) > 0) {
    bevat_codes_escaped <- gsub("([\\.\\^\\$\\*\\+\\?\\(\\)\\[\\{\\\\\\|])", "\\\\\\1", bevat_codes)
    regex_term          <- paste0("^", bevat_codes_escaped, collapse = "|")
    tabel_gefilterd     <- tabel_vlaanderen[CODE %in% exact_codes | grepl(regex_term, CODE)]
  } else {
    tabel_gefilterd     <- tabel_vlaanderen[CODE %in% exact_codes]
  }
  
  # Sommeer fractie per cel_id direct over heel Vlaanderen
  tabel_cel_som <- tabel_gefilterd[, .(Oppervlakte = pmin(sum(BWK_FRAC, na.rm = TRUE), 1.0)), by = .(cel_id)]
  
  if(nrow(tabel_cel_som) > 0) {
    tabel_cel_som[, Match := ifelse(Oppervlakte >= 0.01, 1, 0)]
    
    # Direct toewijzen op het Vlaanderen raster via cel-index
    tabel_matches_clean <- tabel_cel_som[Match == 1]
    if(nrow(tabel_matches_clean) > 0) {
      r_match_vlaanderen[tabel_matches_clean$cel_id] <- 1
    }
    r_opp_vlaanderen[tabel_cel_som$cel_id] <- tabel_cel_som$Oppervlakte
  }
}

masker_bron <- r_match_vlaanderen
opp_bron    <- r_opp_vlaanderen

rm(tabel_vlaanderen, df_nieuw, resultaten_gegroepeerd); gc()

# ------------------------------------------------------------------------------
# 6. DYNAMISCHE EXECUTE OP BASIS VAN EXCEL PARAMETERS
# ------------------------------------------------------------------------------

# Check of er een 2e drempel gedefinieerd is in Excel
heeft_tweede_drempel <- !is.na(min_totale_opp_ha) && min_totale_opp_ha > oppervlakte_ha

if (heeft_tweede_drempel) {
  message("-> [Algoritme] Twee-drempels netwerkfilter toepassen (zoals bij Adder)...")
  
  # SPOOR 1: MAX (Binair op basis van aanwezigheid)
  res_max <- cluster_filter_twee_drempels_vlaanderen(
    masker            = masker_bron, 
    opp_laag          = opp_bron, 
    min_opp_ha        = oppervlakte_ha, 
    min_totale_opp_ha = min_totale_opp_ha, 
    dist_m            = buffer_m, 
    werkelijk         = FALSE
  )
  
  # SPOOR 2: OPP (Volledig parallel op basis van werkelijke fracties)
  r_binair_opp <- terra::ifel(!is.na(opp_bron) & opp_bron > 0, 1, NA)
  res_opp <- cluster_filter_twee_drempels_vlaanderen(
    masker            = r_binair_opp, 
    opp_laag          = opp_bron, 
    min_opp_ha        = oppervlakte_ha, 
    min_totale_opp_ha = min_totale_opp_ha, 
    dist_m            = buffer_m, 
    werkelijk         = TRUE
  )
  
  final_max <- res_max$raster
  final_opp <- res_opp$raster
  cl_opp    <- res_opp$clusters
  
} else {
  message("-> [Algoritme] Standaard 1-drempel clusterfilter toepassen...")
  
  drempel_m2 <- oppervlakte_ha * 10000
  
  res_max <- cluster_filter_vlaanderen(
    masker     = masker_bron, 
    opp_laag   = opp_bron, 
    drempel_m2 = drempel_m2, 
    dist_m     = afstand_m, 
    werkelijk  = FALSE
  )
  
  final_max <- res_max$raster
  cl_max    <- res_max$clusters
  
  if (!all(is.na(suppressWarnings(terra::minmax(cl_max))))) {
    zonal_werkelijk <- terra::zonal(opp_bron, cl_max, fun = "sum", na.rm = TRUE)
    colnames(zonal_werkelijk) <- c("ID", "Werk_m2")
    zonal_werkelijk$Werk_m2   <- zonal_werkelijk$Werk_m2 * 100
    
    ids_werkelijk_voldoet <- zonal_werkelijk$ID[zonal_werkelijk$Werk_m2 >= drempel_m2]
    
    if (length(ids_werkelijk_voldoet) > 0) {
      final_opp <- terra::mask(final_max, cl_max %in% ids_werkelijk_voldoet)
      cl_opp    <- terra::mask(cl_max, cl_max %in% ids_werkelijk_voldoet)
    } else {
      final_opp <- final_max * NA
      cl_opp    <- final_max * NA
    }
  } else {
    final_opp <- final_max * NA
    cl_opp    <- final_max * NA
  }
}

# ------------------------------------------------------------------------------
# 7. UNIVERSELE EXPORT VOOR VLAANDEREN
# ------------------------------------------------------------------------------
exporteer_resultaten_vlaanderen(
  final_max     = final_max,
  final_opp     = final_opp,
  cl_opp        = cl_opp,
  soort         = soort,
  scenario_naam = scenario_naam,
  master_grid   = master_grid
)

message(paste("🏁 VLAANDEREN EXPORT VOLLEDIG AFGEROND VOOR:", toupper(soort)))
