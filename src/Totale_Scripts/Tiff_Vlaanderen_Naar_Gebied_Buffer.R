library(terra)
library(sf)
library(here)

# --- CONFIGURATIE ---
map_rasters <- here("data/output/Vlaanderen/Rasters_Soorten/01_Potentieel_Leefgebied/")
map_shapes  <- here("data/input/")
map_output  <- here("data/output/Rerun/Mechelse_Heide/Potentieel_Leefgebied/")
buffer_dist <- 50000  # 50 km in meters

dir.create(map_output, showWarnings = FALSE, recursive = TRUE)

# Definieer hier je combinaties
combinaties <- data.frame(
  soort = c("Potentieel_Leefgebied_kwartelkoning.tif",
            "Potentieel_Leefgebied_paapje.tif",
            "Potentieel_Leefgebied_wielewaal.tif",
            "Potentieel_Leefgebied_geleluzernevlinder.tif",
            "Potentieel_Leefgebied_gevlektewitsnuitlibel.tif",
            "Potentieel_Leefgebied_bruinesnavelbies.tif",
            "Potentieel_Leefgebied_grotewolfsklauw.tif",
            "Potentieel_Leefgebied_kleinblaasjeskruid.tif",
            "Potentieel_Leefgebied_kleintasjeskruid.tif",
            "Potentieel_Leefgebied_kleinwarkruid.tif",
            "Potentieel_Leefgebied_kleinewolfsklauw.tif",
            "Potentieel_Leefgebied_liggendevleugeltjesbloem.tif",
            "Potentieel_Leefgebied_rodedophei.tif",
            "Potentieel_Leefgebied_rondezonnedauw.tif",
            "Potentieel_Leefgebied_wittesnavelbies.tif",),
  shapefile = c("Mechelse_Heide.shp"),
  stringsAsFactors = FALSE
)

# combinaties <- data.frame(
#   soort = c("Potentieel_Leefgebied_kwartelkoning.tif",
#             "Potentieel_Leefgebied_paapje.tif",
#             "Potentieel_Leefgebied_matkop.tif",
#             "Potentieel_Leefgebied_watersnip.tif",
#             "Potentieel_Leefgebied_wielewaal.tif",
#             "Potentieel_Leefgebied_heivlinder.tif",
#             "Potentieel_Leefgebied_gevlektewitsnuitlibel.tif",
#             "Potentieel_Leefgebied_bruinesnavelbies.tif",
#             "Potentieel_Leefgebied_grotepimpernel.tif",
#             "Potentieel_Leefgebied_kleinwarkruid.tif",
#             "Potentieel_Leefgebied_grutto.tif",
#             "Potentieel_Leefgebied_rondezonnedauw.tif",
#             "Potentieel_Leefgebied_wittesnavelbies.tif",),
#   shapefile = c("Heesbossen.shp"),
#   stringsAsFactors = FALSE
# )

# combinaties <- data.frame(
#   soort = c("Potentieel_Leefgebied_kwartelkoning.tif",
#             "Potentieel_Leefgebied_wielewaal.tif",
#             "Potentieel_Leefgebied_kleinwarkruid.tif",
#             "Potentieel_Leefgebied_kleineratelaar.tif",
#             "Potentieel_Leefgebied_harigeratelaar.tif",
#             "Potentieel_Leefgebied_geleanemoon.tif",
#             "Potentieel_Leefgebied_grauweklauwier.tif",
#             "Potentieel_Leefgebied_geleluzernevlinder.tif"),
#   shapefile = c("Voerstreek.shp"),
#   stringsAsFactors = FALSE
# )


# --- VERWERKING ---
for (i in seq_len(nrow(combinaties))) {
  raw_soort   <- combinaties$soort[i]
  shp_naam    <- combinaties$shapefile[i]
  
  # Verwijder eventuele extensies voor schone bestandsnamen
  soort_naam  <- tools::file_path_sans_ext(raw_soort)
  gebied_naam <- tools::file_path_sans_ext(shp_naam)
  
  tiff_pad <- file.path(map_rasters, paste0(soort_naam, ".tif"))
  shp_pad  <- file.path(map_shapes, paste0(gebied_naam, ".shp"))
  out_pad  <- file.path(map_output, paste0(soort_naam, "_", gebied_naam, "_buffer50km.tif"))
  
  if (!file.exists(tiff_pad) || !file.exists(shp_pad)) {
    warning(paste("Bestand ontbreekt voor:", tiff_pad, "of", shp_pad))
    next
  }
  
  message(paste("Verwerken:", soort_naam, "voor", gebied_naam, "..."))
  
  # 1. Lees het raster en de shapefile in
  r <- rast(tiff_pad)
  gebied <- st_read(shp_pad, quiet = TRUE)
  
  # 2. Zorg voor gelijke projectie (CRS)
  if (st_crs(gebied) != st_crs(r)) {
    gebied <- st_transform(gebied, st_crs(r))
  }
  
  # 3. Voeg de 50 km buffer toe
  gebied_gebufferd <- st_buffer(st_union(gebied), dist = buffer_dist)
  
  # 4. Converteer naar SpatVector voor terra
  buffer_vect <- vect(gebied_gebufferd)
  
  # 5. Snijd bij naar de bounding box van de buffer
  r_cropped <- crop(r, buffer_vect)
  
  # 6. Schrijf weg naar schijf
  writeRaster(r_cropped, out_pad, overwrite = TRUE, gdal = c("COMPRESS=LZW"))
  message(paste("Opgeslagen als:", out_pad))
}
