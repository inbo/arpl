library(data.table)
library(sf)
library(terra)
library(purrr) # Handig voor het proper oplijsten van bestanden
library(here)

# =========================================================================
# 1. PADEN INSTELLEN
# =========================================================================
pad_grote_rds       <- here("data/input/Raster_Vlaanderen/BWK_TidyTabel_Smal_Vlaanderen_2025.rds")
pad_master_grid     <- here("data/input/Raster_Vlaanderen/Vlaanderen_MasterGrid_10m.tif")
map_met_shapefiles  <- here("data/input/") # Pas aan naar jouw map

# =========================================================================
# 2. DATA EN MASTER GRID INLADEN
# =========================================================================
message("-> Grote Vlaanderen RDS en Master Grid inladen...")
vlaanderen_dt <- readRDS(pad_grote_rds)
setDT(vlaanderen_dt) # Hard als data.table zetten voor maximale snelheid

master_grid <- rast(pad_master_grid)

# Zoek alle shapefiles in de opgegeven map
shapefiles_lijst <- c(here("data/input/Voerstreek.shp"))

if(length(shapefiles_lijst) == 0) {
  stop(paste("GEEN shapefiles gevonden in de map:", map_met_shapefiles))
}

# =========================================================================
# 3. LOOPEN OVER ELKE SHAPEFILE (GEBIED)
# =========================================================================
message(paste("-> Start extractie voor", length(shapefiles_lijst), "gebieden..."))

for (pad_shp in shapefiles_lijst) {
  
  # Haal de pure naam van het gebied uit het bestandspad (bijv. "Turnhouts_Vennegebied")
  naam_gebied <- tools::file_path_sans_ext(basename(pad_shp))
  
  message(paste0("\n--- Verwerken: ", naam_gebied, " ---"))
  
  # 1. Lees de shapefile in van het specifieke gebied
  gebied_sf <- st_read(pad_shp, quiet = TRUE)
  
  # FIX: Als de CRS ontbreekt (NA is), dwingen we Lambert 72 (31370) af
  if (is.na(st_crs(gebied_sf))) {
    message("   -> Geen CRS gevonden in bestand. Lambert 72 (EPSG:31370) handmatig toegewezen.")
    gebied_sf <- st_set_crs(gebied_sf, 31370)
  } else if (st_crs(gebied_sf) != st_crs(master_grid)) {
    # Mocht er toch een andere projectie instaan, transformeer dan naar het grid
    gebied_sf <- st_transform(gebied_sf, st_crs(master_grid))
  }
  
  # Zorg dat de projectie exact matcht met het mastergrid (Lambert 72)
  if(st_crs(gebied_sf) != st_crs(master_grid)) {
    gebied_sf <- st_transform(gebied_sf, st_crs(master_grid))
  }
  
  # 2. Maak een buffer van 50 km (50.000 meter)
  message("   -> Buffer van 50 km aanmaken...")
  gebied_buffer <- st_buffer(gebied_sf, 50000)
  
  # 3. Vraag razendsnel de cel-ID's op die binnen deze buffer vallen
  message("   -> Cel-ID's bepalen uit master grid...")
  cellen_in_buffer <- terra::cells(master_grid, vect(gebied_buffer))[, "cell"]
  
  # 4. Filter de grote data.table op basis van deze cel_id's
  # %in% behoudt de exacte volgorde waarin de cellen in vlaanderen_dt staan
  message("   -> Tabel filteren...")
  gebied_subset_dt <- vlaanderen_dt[cel_id %in% cellen_in_buffer]
  
  # 5. Opslaan als een aparte compacte RDS file
  output_naam <- paste0("BWK_Tabel_", naam_gebied, "_Plus_50km.rds")
  saveRDS(gebied_subset_dt, output_naam)
  
  message(paste("   ✓ Succesvol opgeslagen als:", output_naam))
  
  # Geheugen netjes opschonen voor de volgende shapefile
  rm(gebied_sf, gebied_buffer, cellen_in_buffer, gebied_subset_dt)
  gc()
}

message("\n✓ KLAAR! Alle gebieden zijn succesvol uitgesneden op basis van de shapefiles + 50km buffer.")