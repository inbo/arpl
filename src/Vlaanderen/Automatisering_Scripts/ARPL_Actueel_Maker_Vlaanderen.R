# ==============================================================================
# ARPL & ACTUEEL MAKER - VELDKREKEL HEEL VLAANDEREN
# Inclusief automatische detectie van dispersiecapaciteit & hoge buffercap
# Author: Bert Van Hecke / Update 2026
# ==============================================================================

library(here)
library(dplyr)
library(readxl)
library(readr)
library(sf)
library(terra)

# Geheugeninstellingen dwingen voor Vlaanderen-schaal (10m grid)
terra::terraOptions(memfrac = 0.4, tempdir = tempdir(), verbose = FALSE)

# ------------------------------------------------------------------------------
# 1. INSTELLINGEN & PADEN
# ------------------------------------------------------------------------------
soort_clean   <- "veldkrekel"
soort_format  <- "Veldkrekel"
scenario_naam <- "BWK_2025"

# Input mappen (gebaseerd op je leefgebieden-script)
base_input_dir  <- here("data/output/Vlaanderen/Rasters_Soorten", scenario_naam)
id_raster_pad   <- file.path(base_input_dir, "00_ID_Rasters", paste0("ID_Netwerken_", soort_clean, ".tif"))
werkelijk_pad  <- file.path(base_input_dir, "02_Inclusief_Suboptimaal_NA", paste0("Habitat_Inclusief_NA_", soort_clean, ".tif"))

master_grid_pad  <- here("data/input/Raster_Vlaanderen/Vlaanderen_MasterGrid_10m.tif")
waarnemingen_pad <- here("data/input/Waarnemingen_Soorten/Vlaanderen", paste0("Waarnemingen_", soort_format, ".csv"))

# Output mappen
out_dir_arpl <- file.path(base_input_dir, "03_ARPL")
out_dir_act  <- file.path(base_input_dir, "04_Actuele_Verspreiding")

if (!dir.exists(out_dir_arpl)) dir.create(out_dir_arpl, recursive = TRUE)
if (!dir.exists(out_dir_act))  dir.create(out_dir_act,  recursive = TRUE)

f_out_arpl <- file.path(out_dir_arpl, paste0("Habitat_ARPL_", soort_clean, ".tif"))
f_out_act  <- file.path(out_dir_act,  paste0("Habitat_Actuele_Verspreiding_", soort_clean, ".tif"))

message("==================================================")
message(" STARTEN ARPL & ACTUEEL VOOR: ", toupper(soort_clean), " (VLAANDEREN)")
message("==================================================")

# ------------------------------------------------------------------------------
# 2. DATA INLEZEN
# ------------------------------------------------------------------------------
master_grid <- terra::rast(master_grid_pad)[[1]]

if (!file.exists(id_raster_pad)) {
  stop("❌ ID-raster niet gevonden op: ", id_raster_pad)
}
cl_max <- terra::rast(id_raster_pad)

# Dispersiecapaciteit / Buffer ophalen uit Excel
df_afstanden <- read_excel(here("data/input/Excel_files/Soorten_bwk_afstanden.xlsx")) %>% 
  mutate(Soort_clean = tolower(gsub(" ", "", trimws(Soort))))

row_dist <- df_afstanden %>% filter(Soort_clean == soort_clean)
buffer_m <- if (nrow(row_dist) > 0) row_dist$Dispersiecap_m[1] else 500

message("-> Dispersiecapaciteit ingesteld op: ", buffer_m, " meter")

# ------------------------------------------------------------------------------
# 3. WAARNEMINGEN INLEZEN & PROCESSING
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ROBUUSTE INLEESFUNCTIE VOOR WAARNEMINGEN
# ------------------------------------------------------------------------------
laad_waarnemingen <- function(pad) {
  if (!file.exists(pad)) return(NULL)
  
  # Inlezen (probeer eerst ;, anders ,)
  pts <- tryCatch(
    readr::read_delim(pad, delim = ";", show_col_types = FALSE),
    error = function(e) readr::read_csv(pad, show_col_types = FALSE)
  )
  
  colnames(pts) <- tolower(trimws(colnames(pts)))
  
  if (all(c("x", "y") %in% colnames(pts))) {
    # 1. Zet om naar character
    x_char <- as.character(pts$x)
    y_char <- as.character(pts$y)
    
    # 2. Vervang eventuele komma's door punten
    x_clean <- gsub(",", ".", x_char)
    y_clean <- gsub(",", ".", y_char)
    
    # 3. Omzetten naar getal
    pts_x <- as.numeric(x_clean)
    pts_y <- as.numeric(y_clean)
    
    # VEILIGHEID: Als de getallen alsnog te groot zijn (> 1.000.000), herstel de komma op de juiste positie (6 cijfers voor de komma bij Lambert 72)
    if (any(na.omit(pts_x) > 1e6)) {
      # Lambert 72 X ligt tussen 20.000 en 260.000 (5-6 cijfers)
      # Lambert 72 Y ligt tussen 150.000 en 250.000 (6 cijfers)
      pts_x <- ifelse(pts_x > 1e6, pts_x / 10^(nchar(as.character(floor(pts_x))) - 6), pts_x)
      pts_y <- ifelse(pts_y > 1e6, pts_y / 10^(nchar(as.character(floor(pts_y))) - 6), pts_y)
    }
    
    clean_df <- data.frame(x = pts_x, y = pts_y) %>% filter(!is.na(x) & !is.na(y))
    
    if (nrow(clean_df) > 0) {
      message("   Geldige waarnemingen ingelezen: ", nrow(clean_df))
      return(terra::vect(as.matrix(clean_df[, c("x", "y")]), type = "points", crs = "EPSG:31370"))
    }
  }
  return(NULL)
}

v_pts <- laad_waarnemingen(waarnemingen_pad)

# ------------------------------------------------------------------------------
# 4. ARPL BEREKENEN
# ------------------------------------------------------------------------------
message("-> ARPL berekenen...")

if (!is.na(buffer_m) && buffer_m >= 250000) {
  # Veiligheidscap voor hele grote buffers (bijv. Porseleinhoen-logica)
  if (!is.null(v_pts) && length(v_pts) > 0) {
    message("   💡 Buffer >= 250km met waarnemingen gedetecteerd: Werkelijke Oppervlakte overgenomen.")
    r_werkelijk <- terra::rast(werkelijk_pad)
    arpl_export_rast <- terra::ifel(!is.na(r_werkelijk) & r_werkelijk > 0, 1, NA)
    rm(r_werkelijk)
  } else {
    arpl_export_rast <- terra::rast(master_grid, vals = NA)
  }
} else {
  # Standaard bufferverwerking
  v_blobs <- NULL
  if (!is.null(v_pts)) {
    tryCatch({
      v_blobs <- terra::buffer(v_pts, width = buffer_m)
      v_blobs <- terra::aggregate(v_blobs)
    }, error = function(e) v_blobs <<- NULL)
  }
  
  if (!is.null(v_blobs) && !all(is.na(suppressWarnings(terra::minmax(cl_max))))) {
    ext_blobs <- terra::extract(cl_max, v_blobs, ID = FALSE)
    ids_arpl  <- if (!is.null(ext_blobs) && nrow(ext_blobs) > 0) unique(na.omit(ext_blobs[[1]])) else c()
    
    if (length(ids_arpl) > 0) {
      arpl_export_rast <- cl_max %in% ids_arpl
      arpl_export_rast <- terra::ifel(arpl_export_rast == 1, 1, NA)
    } else {
      arpl_export_rast <- terra::rast(master_grid, vals = NA)
    }
    rm(v_blobs)
  } else {
    arpl_export_rast <- terra::rast(master_grid, vals = NA)
  }
}

# ------------------------------------------------------------------------------
# 5. ACTUELE VERSPREIDING BEREKENEN
# ------------------------------------------------------------------------------
message("-> Actuele verspreiding berekenen...")

if (!is.null(v_pts) && !all(is.na(suppressWarnings(terra::minmax(cl_max))))) {
  ext_points  <- terra::extract(cl_max, v_pts, ID = FALSE)
  ids_actueel <- if (!is.null(ext_points) && nrow(ext_points) > 0) unique(na.omit(ext_points[[1]])) else c()
  
  if (length(ids_actueel) > 0) {
    actueel_export_rast <- cl_max %in% ids_actueel
    actueel_export_rast <- terra::ifel(actueel_export_rast == 1, 1, NA)
  } else {
    actueel_export_rast <- terra::rast(master_grid, vals = NA)
  }
} else {
  actueel_export_rast <- terra::rast(master_grid, vals = NA)
}

if (!is.null(v_pts)) rm(v_pts)

# ------------------------------------------------------------------------------
# 6. EXPORTEN NAAR SCHIJF
# ------------------------------------------------------------------------------
message("-> Uitlijnen met MasterGrid & Exporteren met LZW compressie...")

arpl_export_rast    <- terra::extend(arpl_export_rast, master_grid, fill = NA)
actueel_export_rast <- terra::extend(actueel_export_rast, master_grid, fill = NA)

terra::writeRaster(arpl_export_rast,    f_out_arpl, overwrite = TRUE, gdal = c("COMPRESS=LZW"), datatype = "INT1U", NAflag = 255)
terra::writeRaster(actueel_export_rast, f_out_act,  overwrite = TRUE, gdal = c("COMPRESS=LZW"), datatype = "INT1U", NAflag = 255)

# Opruimen van het geheugen
rm(cl_max, arpl_export_rast, actueel_export_rast, master_grid)
terra::tmpFiles(current = TRUE, orphan = TRUE, old = TRUE, remove = TRUE)
gc(verbose = FALSE)

message("\n==================================================")
message(" 🎉 VELDKREKEL ARPL & ACTUEEL SUCCESVOL GEËXPORTEERD!")
message("==================================================")
