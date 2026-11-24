library(terra)
library(sf)
library(here)

# --- CONFIGURATIE ---
base_rerun_dir   <- here("data/output/Rerun")
map_id_vlaanderen <- here("data/output/Vlaanderen/Rasters_Soorten/00_ID_Rasters") # <-- Pas aan indien anders!
map_shapes       <- here("data/input")
buffer_dist      <- 50000  # 50 km in meters
FORCE_OVERWRITE  <- FALSE  # Zet op TRUE om bestaande ID-rasters in Rerun te overschrijven

# Zoek alle gebieden die in de Rerun-map staan
gebieden <- list.dirs(base_rerun_dir, full.names = FALSE, recursive = FALSE)

message("==================================================")
message(" STARTEN BIJSNIJDEN ID_NETWERKEN (O.B.V. RERUN)")
message("==================================================")

for (gebied_naam in gebieden) {
  map_pot <- file.path(base_rerun_dir, gebied_naam, "Potentieel_Leefgebied")
  map_out <- file.path(base_rerun_dir, gebied_naam, "ID_Rasters")
  shp_pad <- file.path(map_shapes, paste0(gebied_naam, ".shp"))
  
  if (!dir.exists(map_pot)) next
  
  pot_bestanden <- list.files(map_pot, pattern = "\\.tif$", full.names = FALSE)
  if (length(pot_bestanden) == 0) next
  
  if (!file.exists(shp_pad)) {
    warning("Shapefile niet gevonden voor gebied: ", gebied_naam, " (", shp_pad, ")")
    next
  }
  
  message("\n📍 Gebied: ", gebied_naam, " (", length(pot_bestanden), " soorten gevonden)")
  dir.create(map_out, showWarnings = FALSE, recursive = TRUE)
  
  # Shapefile 1x per gebied inladen en bufferen (scheelt veel rekentijd)
  gebied_sf <- st_read(shp_pad, quiet = TRUE)
  # We transformeren de shapefile straks naar het CRS van het eerste raster
  buffer_vect <- NULL
  
  for (b_naam in pot_bestanden) {
    # Extraheer de soortnaam (inclusief eventueel _wv achtervoegsel)
    # Bv. "Potentieel_Leefgebied_grutto_Heesbossen_buffer50km.tif" -> "grutto"
    soort_suffix <- tolower(b_naam) |>
      (\(x) gsub("^potentieel_leefgebied_", "", x))() |>
      (\(x) gsub(paste0("_", tolower(gebied_naam), ".*$"), "", x))() |>
      (\(x) gsub("_buffer.*$|\\.tif$", "", x))() |>
      trimws()
    
    out_pad <- file.path(map_out, paste0("ID_Netwerken_", soort_suffix, ".tif"))
    
    if (!FORCE_OVERWRITE && file.exists(out_pad)) {
      message("   ⏩ Reeds aanwezig (overgeslagen): ", basename(out_pad))
      next
    }
    
    # Zoek het bijhorende Vlaanderen ID-raster (hoofdletterongevoelig voor de zekerheid)
    alle_vl_ids <- list.files(map_id_vlaanderen, pattern = "\\.tif$", full.names = TRUE)
    doel_naam   <- paste0("id_netwerken_", soort_suffix, ".tif")
    match_idx   <- which(tolower(basename(alle_vl_ids)) == doel_naam)
    
    # Fallback als _wv niet in de ID-netwerk bestandsnaam zit
    if (length(match_idx) == 0 && grepl("_wv$", soort_suffix)) {
      doel_naam_zonder_wv <- paste0("id_netwerken_", gsub("_wv$", "", soort_suffix), ".tif")
      match_idx <- which(tolower(basename(alle_vl_ids)) == doel_naam_zonder_wv)
    }
    
    if (length(match_idx) == 0) {
      warning("   ⚠️ Geen Vlaanderen ID-raster gevonden voor: ", soort_suffix)
      next
    }
    
    id_vl_pad <- alle_vl_ids[match_idx[1]]
    message("   ✂️ Bijsnijden: ", basename(id_vl_pad), " -> ", basename(out_pad))
    
    r_id <- rast(id_vl_pad)
    
    # Maak de 50 km buffer aan (1x per gebied)
    if (is.null(buffer_vect)) {
      if (st_crs(gebied_sf) != st_crs(r_id)) {
        gebied_sf <- st_transform(gebied_sf, st_crs(r_id))
      }
      gebied_gebufferd <- st_buffer(st_union(gebied_sf), dist = buffer_dist)
      buffer_vect <- vect(gebied_gebufferd)
    }
    
    # Snijd bij
    r_cropped <- crop(r_id, buffer_vect)
    
    # Behoud exact hetzelfde datatype (belangrijk voor grote ID-nummers in heel Vlaanderen!)
    orig_dtype <- datatype(r_id)[1]
    writeRaster(
      r_cropped, 
      out_pad, 
      overwrite = TRUE, 
      gdal = c("COMPRESS=LZW"),
      datatype = orig_dtype
    )
    
    # Geheugen opschonen
    tmpFiles(current = TRUE, orphan = TRUE, old = TRUE, remove = TRUE)
  }
}

message("\n==================================================")
message(" BIJSNIJDEN ID_NETWERKEN VOLTOOID")
message("==================================================")
