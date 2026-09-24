library(here)
library(tidyverse)
library(sf)
library(terra)

conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::select)

terraOptions(memfrac = 0.8, tempdir = tempdir(), verbose = FALSE)

# ==============================================================================
# 0. MAPSTRUCTUUR & PADEN INSTELLEN
# ==============================================================================
folder_heide_output  <- here::here("data/output/Heide")
folder_heide_soorten <- file.path(folder_heide_output, "Individuele_Soorten")
purrr::walk(c(folder_heide_output, folder_heide_soorten), ~if (!dir.exists(.x)) dir.create(.x, showWarnings = FALSE, recursive = TRUE))

bron_map_werkelijk <- here::here("data/output/Turnhouts_Vennegebied/Rasters_Soorten/BWK_2025/02_Werkelijke_Oppervlaktes")

area_shape  <- vect(here("data/input/Turnhouts_Vennegebied.shp"))
master_grid <- rast(here("data/input/Raster_Vlaanderen/Vlaanderen_MasterGrid_10m.tif"))[[1]]

# HERSTEL CRS MATCH: Garandeer dat vector en mastergrid exact dezelfde CRS string delen
area_shape_proj <- project(area_shape, crs(master_grid))
crs(area_shape_proj) <- crs(master_grid) # Voorkomt [mask] CRS do not match

# ==============================================================================
# 1. SELECTIE PUUR OP HABITATCODES 4010 | 4030 (+ INCLUSIES / EXCLUSIES)
# ==============================================================================
# OPLOSSING PARSING WARNING: Lees alle kolommen als character in
df_crosswalk <- read_csv(
  here::here("data/input/Excel_files/Resultaten_Totaal_Samengevoegd.csv"), 
  col_types = cols(.default = "c"),
  show_col_types = FALSE
)

# A. Automatisch heidesoorten ophalen uitsluitend via HIC-codes 4010 of 4030
soorten_4010_4030 <- df_crosswalk %>%
  filter(grepl("4010|4030", Code)) %>%
  pull(Soort) %>%
  trimws() %>%
  unique()

# B. Expliciete handmatige inclusies (nieuwe/potentiële heidesoorten)
handmatige_inclusies <- c("adder", "aardbeivlinder", "bruine eikenpage", "bruineeikenpage", "zadelsprinkhaan")

# C. Expliciete handmatige exclusies
handmatige_exclusies <- c("wespendief", "zomertortel", "watersnip")

# D. Stel de definitieve doellijst samen
kandidaat_soorten   <- unique(tolower(trimws(c(soorten_4010_4030, handmatige_inclusies))))
finale_soortenlijst <- setdiff(kandidaat_soorten, handmatige_exclusies)

# ==============================================================================
# 2. VALIDEER WELKE SOORTEN EFFECTIEF EEN .TIF HEBBEN IN DE BRONMAP
# ==============================================================================
aanwezige_tif_bestanden <- list.files(bron_map_werkelijk, pattern = "\\.tif$", full.names = TRUE)

gevalideerde_heidesoorten <- list()

for (h_soort in finale_soortenlijst) {
  h_soort_regex <- gsub(" ", "[ _]?", h_soort)
  
  matchende_file <- aanwezige_tif_bestanden[
    grepl(paste0("Habitat_Werkelijke_Oppervlaktes_", h_soort_regex, "(_wv)?\\.tif$"), aanwezige_tif_bestanden, ignore.case = TRUE)
  ]
  
  if (length(matchende_file) > 0) {
    gevalideerde_heidesoorten[[h_soort]] <- matchende_file[1]
  } else {
    message(paste("   [Nog niet aanwezig/overgeslagen]:", h_soort))
  }
}

message("\n========================================================================")
message(paste("=== TOTAAL GEVALIDEERDE HEIDESOORTEN IN HEATMAP:", length(gevalideerde_heidesoorten)))
message("========================================================================")
print(names(gevalideerde_heidesoorten))

# ==============================================================================
# 3. BESTANDEN OPHALEN, KOPIËREN EN STACKEN (GEOTIFF)
# ==============================================================================
rasters_heidesoorten_lijst <- list()

for (h_soort in names(gevalideerde_heidesoorten)) {
  tif_pad_bron  <- gevalideerde_heidesoorten[[h_soort]]
  h_soort_clean <- gsub(" ", "_", h_soort)
  
  # Laad het bestaande raster in
  r_soort <- terra::rast(tif_pad_bron)
  
  # HERSTEL CRS MATCH: Dwing de exacte CRS van het mastergrid af op het geladen raster
  crs(r_soort) <- crs(master_grid) # Voorkomt [rast] CRS do not match
  
  # Bewaar in lijst voor opstapeling (afgeknipt op studiegebied)
  r_soort_cropped <- terra::crop(r_soort, area_shape_proj, snap = "near")
  r_soort_masked  <- terra::mask(r_soort_cropped, area_shape_proj)
  rasters_heidesoorten_lijst[[h_soort]] <- r_soort_masked
  
  # Kopiëer individuele GeoTIFF naar de Heide-map
  doel_tif <- file.path(folder_heide_soorten, paste0("Werkelijk_Habitat_", h_soort_clean, ".tif"))
  terra::writeRaster(r_soort, doel_tif, overwrite = TRUE, gdal = c("COMPRESS=LZW"), datatype = "FLT4S", NAflag = -9999)
  
  message(paste("   [GeoTIFF OK]: Gekopieerd voor", toupper(h_soort)))
}

# ==============================================================================
# 4. GENEREREN FINALE HEIDE-HEATMAP GEOTIFF
# ==============================================================================
if (length(rasters_heidesoorten_lijst) > 0) {
  
  # Stapel alle lagen op en bereken de cumulatieve som per rastercel
  heide_stack <- terra::rast(rasters_heidesoorten_lijst)
  heide_heatmap_raster <- terra::app(heide_stack, fun = "sum", na.rm = TRUE)
  heide_heatmap_raster <- terra::ifel(heide_heatmap_raster > 0, heide_heatmap_raster, NA)
  
  # Export Totaal Heatmap GeoTIFF op het MasterGrid formaat
  heatmap_tif_pad <- file.path(folder_heide_output, "Heatmap_Heidesoorten_Turnhouts_Vennegebied.tif")
  export_heatmap <- terra::extend(heide_heatmap_raster, master_grid, fill = NA)
  
  terra::writeRaster(
    export_heatmap, 
    heatmap_tif_pad, 
    overwrite = TRUE, 
    gdal = c("COMPRESS=LZW"), 
    datatype = "FLT4S", 
    NAflag = -9999
  )
  
  message("\n🎉 [GeoTIFF HEATMAP OK]: Heatmap_Heidesoorten_Turnhouts_Vennegebied.tif")
  
} else {
  message("❌ Geen matchende .tif-bestanden gevonden in de bronmap.")
}

message("\n🏁 HEIDE HEATMAP GEOTIFF SCRIPT VOLLEDIG EN SCHOON AFGEROND!")
