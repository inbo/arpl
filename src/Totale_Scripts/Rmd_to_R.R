library(tidyverse)
library(here)

# ==============================================================================
# 1. MAPINSTELLINGEN
# ==============================================================================
input_dir  <- here("src/Totale_Scenario_Scripts")
output_dir <- here("src/01_Habitat_Scripts")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

rmd_files <- list.files(input_dir, pattern = "\\.Rmd$", full.names = TRUE)

message(paste("🚀 Starten met het geautomatiseerd converteren van", length(rmd_files), "Rmd-bestanden...\n"))

# ==============================================================================
# 2. BATCH CONVERSIE LUS (BEHOUDT ALLE FUNCTIES UIT SETUP-CHUNKS)
# ==============================================================================
for (file in rmd_files) {
  lines <- readLines(file, encoding = "UTF-8", warn = FALSE)
  bestand_naam <- basename(file)
  
  # A. Strip YAML header (alles tussen de eerste twee ---)
  yaml_bounds <- which(lines == "---")
  if (length(yaml_bounds) >= 2) {
    lines <- lines[(yaml_bounds[2] + 1):length(lines)]
  }
  
  # B. Vind alle Rmd chunk grenzen
  chunk_starts <- which(grepl("^```\\{r", lines))
  chunk_ends   <- which(grepl("^```$", lines))
  
  if (length(chunk_starts) != length(chunk_ends)) {
    chunk_ends <- map_int(chunk_starts, ~ chunk_ends[chunk_ends > .x][1])
  }
  
  schone_chunks <- list()
  
  # C. Verwerkt elke chunk als een geïsoleerd blok
  for (i in seq_along(chunk_starts)) {
    start_idx <- chunk_starts[i]
    end_idx   <- chunk_ends[i]
    
    header_line <- lines[start_idx]
    chunk_lines <- lines[(start_idx + 1):(end_idx - 1)]
    
    # 1. UITZONDERINGSFILTER: Knip knitr instellingen weg, maar BEHOUDT de functies en imports!
    chunk_lines <- chunk_lines[!grepl("knitr::opts_", chunk_lines)]
    chunk_lines <- chunk_lines[!grepl("opts_chunk\\$set", chunk_lines)]
    chunk_lines <- chunk_lines[!grepl("opts_knit\\$set", chunk_lines)]
    
    # Als de chunk na deze opschoning volledig leeg is, sla over
    if (length(trimws(chunk_lines[chunk_lines != ""])) == 0) {
      next
    }
    
    # 2. SLUIT ONGEWENSTE RAPPORSTAGE-CHUNKS UIT
    is_ongewenst <- grepl("statistiek|transitiekaart|leaflet|audit|wasstraat|export|waarnemingen", header_line, ignore.case = TRUE) ||
      any(grepl("leaflet\\(|kable\\(|addRasterImage\\(|addPolygons\\(", chunk_lines))
    
    if (is_ongewenst) {
      next
    }
    
    # 3. VERVANG PARAMETERS IN DE SCENARIO CHUNK
    if (any(grepl("params\\$scenario_rds_path", chunk_lines))) {
      chunk_text <- paste(chunk_lines, collapse = "\n")
      params_fix <- '
# --- DYNAMISCHE SCENARIO PARAMETER CHECK ---
if (!exists("params") || is.null(params$scenario_rds_path)) {
  scenario_rds_path <- "data/input/Scenario_rds/TV_Scenario_BWK_2025.rds"
} else {
  scenario_rds_path <- params$scenario_rds_path
}

p_raw <- gsub("^([.][.]/)+", "", scenario_rds_path)
scenario_path <- here::here(p_raw)
'
    chunk_text <- str_replace(
      chunk_text,
      "(?s)p_raw\\s*<-\\s*gsub\\(.*params\\$scenario_rds_path.*?scenario_path\\s*<-\\s*here::here\\(p_raw\\)",
      params_fix
    )
    schone_chunks[[length(schone_chunks) + 1]] <- chunk_text
    next
    }
  
  # Voeg de goedgekeurde chunk toe (inclusief alle gedefinieerde functies)
  schone_chunks[[length(schone_chunks) + 1]] <- paste(chunk_lines, collapse = "\n")
  }

# D. Voeg het universele, sluitende Exportblok toe
schoon_export_blok <- '
# ==============================================================================
# SCHONE EXPORT BIOTOOP EN ANALYTISCH ID-RASTER (VOOR SCRIPT 2 / ARPL)
# ==============================================================================
base_dir <- here::here("data/output/Turnhouts_Vennegebied/Rasters_Soorten", scenario_naam)

folders <- list(
  potentie  = file.path(base_dir, "01_Maximale_Potentie"),
  werkelijk = file.path(base_dir, "02_Werkelijke_Oppervlakte"),
  id_raster = file.path(base_dir, "00_ID_Rasters")
)
purrr::walk(folders, ~if (!dir.exists(.x)) dir.create(.x, showWarnings = FALSE, recursive = TRUE))

# 1. Bepaal Maximale Potentie Raster
potentie_export_rast <- if (exists("grutto_kaart_A") && !is.null(grutto_kaart_A) && !all(is.na(terra::values(grutto_kaart_A, mat=FALSE)))) {
  terra::ifel(grutto_kaart_A > 0, 1, NA)
} else if (exists("final_max") && !all(is.na(terra::values(final_max, mat=FALSE)))) {
  terra::ifel(!is.na(final_max) & final_max > 0, 1, NA)
} else {
  terra::rast(template_TV, vals = NA)
}

# 2. Bepaal Werkelijke Oppervlakte Raster
werkelijk_export_rast <- if (exists("resB_strikt") && !is.null(resB_strikt) && (!all(is.na(terra::values(resB_strikt$kern, mat=FALSE))) || !all(is.na(terra::values(resB_strikt$bouw, mat=FALSE))))) {
  r_net_totaal_opp <- terra::cover(resB_strikt$kern, resB_strikt$bouw)
  terra::ifel(!is.na(r_net_totaal_opp) & r_net_totaal_opp > 0, 1, NA)
} else if (exists("final_opp") && !all(is.na(terra::values(final_opp, mat=FALSE)))) {
  terra::ifel(!is.na(final_opp) & final_opp > 0, 1, NA)
} else {
  terra::rast(template_TV, vals = NA)
}

# 3. Bepaal Analytisch Metacluster ID-raster
if (exists("resB_strikt") && !is.null(resB_strikt) && (!all(is.na(terra::values(resB_strikt$kern, mat=FALSE))) || !all(is.na(terra::values(resB_strikt$bouw, mat=FALSE))))) {
  r_net_totaal <- terra::cover(resB_strikt$kern, resB_strikt$bouw)
  r_buf        <- terra::buffer(r_net_totaal, width = 100)
  id_export_rast <- terra::mask(terra::patches(r_buf, directions = 8, zeroAsNA = TRUE), r_net_totaal)
} else if (exists("cl_max") && !all(is.na(terra::values(cl_max, mat=FALSE)))) {
  id_export_rast <- cl_max
} else if (exists("cl_opp") && !all(is.na(terra::values(cl_opp, mat=FALSE)))) {
  id_export_rast <- cl_opp
} else {
  id_export_rast <- terra::rast(template_TV, vals = NA)
}

# --- EXPORT LUS VOOR VISUELE BINAIR RASTERS ---
export_config <- list(
  list(rast = potentie_export_rast,  naam = "Maximale_Potentie",        folder = folders$potentie),
  list(rast = werkelijk_export_rast, naam = "Werkelijke_Oppervlaktes", folder = folders$werkelijk)
)

for(item in export_config) {
  suffix <- if (exists("zoek_sleutel") && grepl("_wv$", zoek_sleutel)) "_wv.tif" else ".tif"
  file_path <- file.path(item$folder, paste0("Habitat_", item$naam, "_", soort, suffix))
  
  export_rast <- terra::deepcopy(item$rast)
  export_rast <- terra::extend(export_rast, master_grid, fill = NA)
  
  terra::writeRaster(
    export_rast, 
    filename = file_path, 
    overwrite = TRUE, 
    gdal = c("COMPRESS=LZW"), 
    datatype = "INT1U", 
    NAflag = 255
  )
  message(paste("    [OK] Geëxporteerd:", basename(file_path)))
}

# --- EXPORT VOOR ANALYTISCH ID-RASTER ---
suffix <- if (exists("zoek_sleutel") && grepl("_wv$", zoek_sleutel)) "_wv.tif" else ".tif"
file_path_id <- file.path(folders$id_raster, paste0("ID_Netwerken_", soort, suffix))

export_id_rast <- terra::extend(id_export_rast, master_grid, fill = NA)

terra::writeRaster(
  export_id_rast, 
  filename = file_path_id, 
  overwrite = TRUE, 
  gdal = c("COMPRESS=LZW"), 
  datatype = "INT4U", 
  NAflag = 0
)
message(paste("    [OK] Analytisch ID-raster geëxporteerd:", basename(file_path_id)))

suppressWarnings(
  rm(potentie_export_rast, werkelijk_export_rast, id_export_rast, export_config, export_rast, export_id_rast)
)
gc()

message(paste("🏁 SCENARIO EXPORT VOLLEDIG AFGEROND VOOR:", toupper(soort)))
'

schone_chunks[[length(schone_chunks) + 1]] <- schoon_export_blok

# E. Samenvoegen tot finale R-code
finaal_script <- paste(schone_chunks, collapse = "\n\n")

# F. Opslaan als schoon R-script
out_name <- paste0(tools::file_path_sans_ext(bestand_naam), ".R")
writeLines(finaal_script, file.path(output_dir, out_name), useBytes = TRUE)

message(paste("✔ Geconverteerd:", out_name))
}

message("\n🎉 Alle scenario-scripts zijn geconverteerd met behoud van alle functies uit de setup-chunk!")
