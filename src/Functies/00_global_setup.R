# ==============================================================================
# 00_global_setup.R - GLOBALE INSTELLINGEN & PADEN
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(sf)
  library(terra)
  library(readxl)
  library(data.table)
  library(arrow)
  library(conflicted)
  library(qgisprocess)
})

# Conflict resolution
conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::select)
conflicted::conflicts_prefer(dplyr::first)
conflicted::conflicts_prefer(terra::intersect)
conflicted::conflicts_prefer(terra::any)

# Terra RAM & temp instellingen
terraOptions(
  memfrac = 0.8,
  todisk  = TRUE,
  tempdir = tempdir(),
  verbose = FALSE
)

# Centrale lijst met bestandslocaties
# In 00_global_setup.R binnen PATHS:
# In 00_global_setup.R binnen PATHS:
PATHS <- list(
  # Systeemmappen & Mastergrid
  dir_cache        	  = here::here("data", "cache"),
  master_grid      	  = here::here("data", "input", "Raster_Vlaanderen", "Vlaanderen_MasterGrid_10m.tif"),
  
  # Tabellen & Parquet
  bwk_parquet         = here::here("data", "input", "Raster_Vlaanderen", "BWK_TidyTabel_Smal_Vlaanderen_2025.parquet"),
  bwk_codes_csv       = here::here("data", "input", "Excel_files", "Resultaten_Totaal_Samengevoegd.csv"),
  params_excel        = here::here("data", "input", "Excel_files", "Soorten_parameters.xlsx"),
  dir_waarnemingen    = here::here("data", "input", "Waarnemingen_Soorten", "Vlaanderen"),
  
  # Externe Omgevingsrasters
  groenkaart          = here::here("data", "input", "ASCI Files", "Groenkaart_2021.tif"),
  drainage            = here::here("data", "input", "Raster_Vlaanderen", "vlaanderen_drainage_n2khab_10m.tif"),
  profiel             = here::here("data", "input", "Raster_Vlaanderen", "Vlaanderen_profiel_10m.tif"),
  overstroming        = here::here("data", "input", "Raster_Vlaanderen", "vlaanderen_ovstrg_10m.tif"),
  grb_water           = here::here("data", "input", "Raster_Vlaanderen", "vlaanderen_grb_water_10m.tif"),
  watervlakken        = here::here("data", "input", "Raster_Vlaanderen", "vlaanderen_watervlakken_2024_10m.tif"),
  maas                = here::here("data", "input", "Raster_Vlaanderen", "maas_10m.tif"),
  huetzon             = here::here("data", "input", "Raster_Vlaanderen", "vlaanderen_huetzon_10m.tif"),
  lgp_landbouw        = here::here("data", "input", "Raster_Vlaanderen", "vlaanderen_lgp_gwscod_2025_10m.tif"),
  grb_wegen           = here::here("data", "input", "Raster_Vlaanderen", "GRB_Wegverbinding_at.tif"),
  
  # Uitvoermoppen
  out_id           = here::here("data", "output", "Vlaanderen", "Rasters_Soorten", "00_ID_Rasters"),
  out_opp          = here::here("data", "output", "Vlaanderen", "Rasters_Soorten", "01_Potentieel_Leefgebied"),
  out_arpl         = here::here("data", "output", "Vlaanderen", "Rasters_Soorten", "02_ARPL")
)

# Mappen automatisch aanmaken inclusief de submap Deellagen
dir.create(PATHS$out_id, showWarnings = FALSE, recursive = TRUE)
dir.create(PATHS$out_opp, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(PATHS$out_opp, "Deellagen"), showWarnings = FALSE, recursive = TRUE)
dir.create(PATHS$out_arpl, showWarnings = FALSE, recursive = TRUE)
dir.create(PATHS$dir_cache, showWarnings = FALSE, recursive = TRUE)
