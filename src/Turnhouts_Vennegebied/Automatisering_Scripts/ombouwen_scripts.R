library(tidyverse)
library(here)

transformeer_naar_scenario_script <- function(input_rmd_pad, output_rmd_pad) {
  
  # 1. Lees het originele bestand in als tekstregels
  lijnen <- readLines(input_rmd_pad, encoding = "UTF-8", warn = FALSE)
  tekst  <- paste(lijnen, collapse = "\n")
  
  # 2. Haal de soortnaam op (bijv. soort <- "grutto")
  m <- regexec('soort\\s*<-\\s*["\']([^"\']+)["\']', tekst)
  reg_res <- regmatches(tekst, m)[[1]]
  
  if (length(reg_res) < 2) {
    warning(paste("Geen soortnaam gevonden in:", input_rmd_pad, "- Wordt overgeslagen."))
    return(NULL)
  }
  
  soort_naam <- reg_res[2]
  soort_titel <- paste0(toupper(substr(soort_naam, 1, 1)), substr(soort_naam, 2, nchar(soort_naam)))
  
  # 3. VERVANGING 1: Nieuwe YAML Header & Setup
  nieuwe_kop <- sprintf('---
title: "Scenario Leefgebiedenscript: %s"
author: "Bert Van Hecke"
date: "`r Sys.Date()`"
output: 
  html_document:
    self_contained: true
    toc: true
    theme: flatly
params:
  scenario_rds_path: "data/input/Scenario_rds/TV_Scenario_bosbehoudss_ss31fix_tvg_vrij_2026.rds"
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  echo = FALSE,
  message = FALSE,
  warning = FALSE,
  include = FALSE
)

library(knitr)
library(here)
# Dwing R Markdown om te werken vanaf de hoofdmap (arpl/)
knitr::opts_knit$set(root.dir = here::here())

library(tidyverse)
library(sf)
library(terra)
library(readxl)
library(tidyterra)
library(leaflet)
library(kableExtra)
library(data.table)

conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::select)
conflicted::conflicts_prefer(terra::intersect)
conflicted::conflicts_prefer(terra::any)

calc_ha_exact <- function(r) {
  if(is.null(r)) return(0)
  if(all(is.na(terra::values(r, mat=FALSE)))) return(0)
  area_raster <- r * terra::cellSize(r, unit = "ha")
  val <- terra::global(area_raster, "sum", na.rm = TRUE)[[1]]
  return(as.numeric(val))
}
```', soort_titel)

# 4. VERVANGING 2: Soort- en Scenario-informatie
nieuwe_info_chunk <- sprintf('# Soort en scenario informatie ophalen
```{r Soort informatie}
soort <- "%s"

# Strip eventuele relatieve klim-paden (../) als die toch worden meegegeven
p_raw <- gsub("^(\\\\.\\\\./)+", "", params$scenario_rds_path)
scenario_path <- here::here(p_raw)

if (!file.exists(scenario_path)) {
  stop(paste("❌ FOUT: Scenario RDS bestand NIET gevonden op:", scenario_path))
}

scen_volledig <- basename(scenario_path)
scenario_naam <- gsub("^TV_Scenario_|^Scenario_|.rds$", "", scen_volledig)

message(paste("Verwerken van soort:", soort, "binnen scenario:", scenario_naam))

df <- read_excel(here::here("data/input/Excel_files/Soorten_bwk_afstanden.xlsx"))
resultaat <- df %%>%% filter(tolower(trimws(Soort)) == soort) %%>%% select(Type, MinOpp_ha, AfstandBiotopen_m, Dispersiecap_m)

oppervlakte_ha <- resultaat$MinOpp_ha[1]
afstand_m      <- resultaat$AfstandBiotopen_m[1]
buffer_m       <- resultaat$Dispersiecap_m[1]

rm(df, resultaat)
```', soort_naam)

# 5. CHUNK-GEBASEERD AFKNIPPEN VAN HET STAARTSTUK
# Vind alle regelnummers waar een codechunk start (```{r ...)
chunk_starts <- grep("^```\\{r", lijnen)
n_chunks <- length(chunk_starts)

if (n_chunks < 3) {
  stop(paste("Script heeft te weinig chunks om op te splitsen:", input_rmd_pad))
}

# Bepaal het startpunt van de romp (Chunk 3 = Gebied voorbereiden)
start_line_romp <- chunk_starts[3]

# Bepaal het knippunt (standaard 3 chunks vanaf het einde)
cut_chunk_idx <- max(3, n_chunks - 2)

# Check voor de zekerheid de laatste 5 chunks op audit/statistiek sleutelwoorden
check_indices <- max(3, n_chunks - 4):(n_chunks - 2)
for (idx in check_indices) {
  chunk_header <- lijnen[chunk_starts[idx]]
  if (grepl("audit|resultaten|statistiek|wasstraat|waarneming|export|stap 13|stap 12", chunk_header, ignore.case = TRUE)) {
    cut_chunk_idx <- idx
    break
  }
}

cut_line <- chunk_starts[cut_chunk_idx]

# Trek eventuele losse koppen (# ...) vlak boven de geknipte chunk mee
while (cut_line > start_line_romp && (grepl("^#", trimws(lijnen[cut_line - 1])) || trimws(lijnen[cut_line - 1]) == "")) {
  if (grepl("^```", trimws(lijnen[cut_line - 1]))) break
  cut_line <- cut_line - 1
}

# Isoleer de zuivere ecologische filter-romp
romp_lijnen <- lijnen[start_line_romp:(cut_line - 1)]
romp <- paste(romp_lijnen, collapse = "\n")

# 6a. WATERDICHTE VERVANGING VAN RDS-INLEZING
patroon_rds <- 'readRDS\\s*\\(\\s*(here::here\\s*\\(|here\\s*\\()?\\s*["\'][^"\']*BWK_Tabel[^"\']*["\']\\s*\\)?\\s*\\)'
romp <- gsub(patroon_rds, 'readRDS(scenario_path)', romp, perl = TRUE)

# 6b. SNELHEIDSOPTIMALISATIE: TERRA::EXTRACT VERVANGEN DOOR C++ MASK + CELLS
patroon_trage_extract <- 'extractie_ids\\s*<-\\s*terra::extract\\(id_raster_TV_globale_values,\\s*area_buffer_fix,\\s*cells\\s*=\\s*TRUE\\)\\s*\n+\\s*vertaal_df\\s*<-\\s*as\\.data\\.table\\(extractie_ids\\)\\s*\n+\\s*setnames\\(vertaal_df,\\s*c\\("cell",\\s*names\\(id_raster_TV_globale_values\\)\\),\\s*c\\("lokale_id",\\s*"globale_id"\\)\\)\\s*\n+\\s*vertaal_df\\s*<-\\s*vertaal_df\\[!is\\.na\\(globale_id\\)\\]'

vervanging_snelle_mask <- '# OPTIMALISATIE: Masker het gecropte raster met de buffer (Razendsnel in C++)
id_raster_TV_masked <- terra::mask(id_raster_TV_globale_values, area_buffer_fix)

# Haal direct de niet-NA cel-indices (lokale IDs) en waarden (globale IDs) op
lokale_ids  <- terra::cells(id_raster_TV_masked)
globale_ids <- id_raster_TV_masked[lokale_ids][[1]]

vertaal_df <- data.table(
  lokale_id  = lokale_ids,
  globale_id = globale_ids
)
vertaal_df <- vertaal_df[!is.na(globale_id)]'

romp <- gsub(patroon_trage_extract, vervanging_snelle_mask, romp, perl = TRUE)

# 7. EXPORT & TRANSITIEKAART
nieuwe_staart <- '
# Export & Transitiekaart
```{r Export}
scenario_raster_dir <- here::here("data/output/Turnhouts_Vennegebied/Rasters_Soorten", scenario_naam)
if (!dir.exists(scenario_raster_dir)) dir.create(scenario_raster_dir, recursive = TRUE)

final_opp_scen <- NULL

if (exists("resB_strikt") && !is.null(resB_strikt)) {
  final_opp_scen <- max(c(resB_strikt$kern, resB_strikt$bouw), na.rm = TRUE)
} else if (exists("final_opp")) {
  final_opp_scen <- final_opp
} else if (exists("bruinekiekendief_leefgebied1_opp")) {
  final_opp_scen <- bruinekiekendief_leefgebied1_opp
} else {
  mogelijke_vars <- ls(pattern = "(_final_opp|_leefgebied1_opp|_finaal_opp)$")
  if (length(mogelijke_vars) > 0) {
    final_opp_scen <- get(mogelijke_vars[1])
  } else {
    final_opp_scen <- terra::rast(template_TV, vals = NA)
  }
}

if (!all(is.na(terra::values(final_opp_scen, mat = FALSE)))) {
  werkelijk_scenario_rast <- terra::ifel(!is.na(final_opp_scen) & final_opp_scen > 0, 1, NA)
} else {
  werkelijk_scenario_rast <- terra::rast(template_TV, vals = NA)
}

bestandsnaam <- paste0("TV_", scenario_naam, "_", soort, ".tif")
file_path_scen <- file.path(scenario_raster_dir, bestandsnaam)

export_rast <- terra::deepcopy(werkelijk_scenario_rast)
export_rast <- terra::mask(export_rast, area_buffer_fix)

if (all(is.na(terra::values(export_rast, mat = FALSE)))) {
  export_rast <- terra::ifel(is.na(export_rast), 0, 0)
  export_rast <- terra::mask(export_rast, area_buffer_fix, updatevalue = 0)
} else {
  export_rast <- terra::mask(export_rast, area_buffer_fix, updatevalue = 0)
}

writeRaster(export_rast, 
            filename = file_path_scen, 
            overwrite = TRUE, 
            gdal = c("COMPRESS=LZW"), 
            datatype = "INT1U")

message(paste("✓ Scenario-raster succesvol weggeschreven naar:", file_path_scen))

referentie_file <- here::here("data/output/Turnhouts_Vennegebied/Rasters_Soorten", "02_Werkelijke_Oppervlakte", paste0("Habitat_Werkelijke_Oppervlaktes_", soort, ".tif"))

if (file.exists(referentie_file)) {
  r_ref <- rast(referentie_file) %>% terra::resample(template_TV, method = "near")
} else {
  r_ref <- terra::rast(template_TV, vals = NA)
}

r_scen <- werkelijk_scenario_rast

r_rgb <- terra::ifel(!is.na(r_ref) & is.na(r_scen), 1, NA)
r_rgb <- terra::cover(r_rgb, terra::ifel(is.na(r_ref) & !is.na(r_scen), 2, NA))
r_rgb <- terra::cover(r_rgb, terra::ifel(!is.na(r_ref) & !is.na(r_scen), 3, NA))

r_rgb_web <- terra::project(terra::as.int(r_rgb), "EPSG:4326", method = "near")

kleuren_rgb <- c("#E74C3C", "#2ECC71", "#3498DB") 
pal_rgb <- colorFactor(palette = kleuren_rgb, levels = 1:3, na.color = "transparent")

leaflet() %>%
  addTiles(group = "OSM (Kaart)") %>%
  addProviderTiles(providers$Esri.WorldImagery, group = "Satelliet") %>%
  addRasterImage(r_rgb_web, colors = pal_rgb, opacity = 0.7, group = "Scenario Vergelijking") %>%
  addPolygons(data = grens_web, fill = FALSE, color = "black", weight = 3, group = "Grens TV") %>%
  addLegend(
    colors = kleuren_rgb,
    labels = c("Verloren gegaan leefgebied", "Nieuw bijgekomen leefgebied", "Status Quo (Behouden)"),
    title = paste("Scenario vs Referentie:", scenario_naam),
    position = "bottomright"
  ) %>%
  addLayersControl(
    baseGroups = c("OSM (Kaart)", "Satelliet"),
    overlayGroups = c("Scenario Vergelijking", "Grens TV"),
    options = layersControlOptions(collapsed = FALSE)
  )
```'

# 8. Samenvoegen en opslaan
volledig_script <- paste(nieuwe_kop, nieuwe_info_chunk, romp, nieuwe_staart, sep = "\n\n")
volledig_script <- gsub("%%>%%", "%>%", volledig_script)

writeLines(volledig_script, output_rmd_pad, useBytes = TRUE)
message(sprintf("✓ Succesvol omgezet: %s -> %s", basename(input_rmd_pad), basename(output_rmd_pad)))
}

# ==============================================================================
# BATCH EXECUTION
# ==============================================================================

map_oud   <- here::here("src/Turnhouts_Vennegebied/Scripts_TV")
map_nieuw <- here::here("src/Turnhouts_Vennegebied/Scenario_test")

if (!dir.exists(map_nieuw)) dir.create(map_nieuw, recursive = TRUE)

alle_oude_rmds <- list.files(map_oud, pattern = "\\.[rR]md$", full.names = TRUE)

cat(sprintf("Starten met omzetten van %d scripts...\n\n", length(alle_oude_rmds)))

purrr::walk(alle_oude_rmds, function(oud_pad) {
  bestand <- basename(oud_pad)
  nieuw_pad <- file.path(map_nieuw, paste0("Scenario_", bestand))
  
  tryCatch({
    transformeer_naar_scenario_script(oud_pad, nieuw_pad)
  }, error = function(e) {
    message(sprintf("❌ Fout bij %s: %s", bestand, e$message))
  })
})
