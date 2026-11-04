library(here)
library(sf)
library(qgisprocess)
library(terra)
library(data.table)

rasterize_layer <- function(f_in,
                            f_out,
                            f_master,
                            field = NULL,            # NULL voor binaire aanwezigheid (bijv. Maas)
                            binary_value = 1,        # Waarde voor aanwezige geometrieën als field = NULL
                            nodata = 0,              # Achtergrond/NoData-waarde (vaak 0 voor binair, 255 voor categorieën)
                            write_rat = TRUE) {
  
  # 1. Extent en resolutie ophalen uit mastergrid
  message("📐 Mastergrid-eigenschappen ophalen...")
  r_master <- terra::rast(f_master)[[1]]
  e <- terra::ext(r_master)
  res_x <- terra::xres(r_master)
  res_y <- terra::yres(r_master)
  ext_str <- sprintf("%f,%f,%f,%f", e$xmin, e$xmax, e$ymin, e$ymax)
  rm(r_master)
  
  # 2. Vector inlezen
  message("📖 Vectorbestand inlezen: ", basename(f_in))
  v <- sf::st_read(f_in, quiet = TRUE)
  
  rat <- NULL
  target_field <- "val_id"
  
  if (is.null(field)) {
    # SCENARIO A: Binair masker (bijv. Maas -> 1 = water, nodata = achtergrond)
    message("ℹ️ Geen kolom opgegeven: binaire laag wordt gegenereerd (waarde = ", binary_value, ")")
    v[[target_field]] <- as.integer(binary_value)
    
  } else {
    # SCENARIO B: Kolom opgegeven (tekst of numeriek)
    if (!field %in% names(v)) {
      stop(sprintf("Kolom '%s' niet gevonden in %s", field, basename(f_in)))
    }
    
    col_data <- v[[field]]
    
    if (is.character(col_data) || is.factor(col_data)) {
      # Automatische lookup maken voor unieke categorieën
      unique_vals <- sort(unique(na.omit(as.character(col_data))))
      lookup <- setNames(seq_along(unique_vals), unique_vals)
      
      v[[target_field]] <- unname(lookup[as.character(col_data)])
      v[[target_field]][is.na(v[[target_field]])] <- nodata
      
      rat <- data.table::data.table(
        Value = as.integer(lookup),
        Label = names(lookup)
      )
    } else {
      # Al numeriek
      v[[target_field]] <- as.integer(col_data)
      v[[target_field]][is.na(v[[target_field]])] <- nodata
    }
  }
  
  # Bepaal datatype (Byte tot 255 waarden, anders Int32)
  max_val <- max(v[[target_field]], na.rm = TRUE)
  data_type <- if (max_val <= 255 && nodata <= 255) 0 else 2 # 0 = Byte, 2 = Int32
  
  # 3. Temp GPKG schrijven
  f_temp_gpkg <- tempfile(fileext = ".gpkg")
  sf::st_write(v[, target_field], f_temp_gpkg, quiet = TRUE)
  rm(v)
  gc()
  
  # 4. GDAL Rasterize via qgisprocess
  message("⚙️ GDAL Rasterize uitvoeren...")
  qgisprocess::qgis_run_algorithm(
    "gdal:rasterize",
    INPUT        = f_temp_gpkg,
    FIELD        = target_field,
    UNITS        = 1,             # Georeferenced meters
    WIDTH        = res_x,
    HEIGHT       = res_y,
    EXTENT       = ext_str,
    NODATA       = nodata,
    DATA_TYPE    = data_type,
    INIT         = nodata,
    OPTIONS      = "COMPRESS=LZW",
    OUTPUT       = f_out,
    .quiet       = TRUE
  )
  unlink(f_temp_gpkg)
  
  # 5. Labels insluiten in GeoTIFF & CSV wegschrijven
  if (!is.null(rat)) {
    r_out <- terra::rast(f_out)
    levels(r_out) <- as.data.frame(rat)
    terra::writeRaster(r_out, f_out, overwrite = TRUE, gdal = "COMPRESS=LZW")
    
    if (write_rat) {
      rat_path <- sub("\\.tif$", "_RAT.csv", f_out)
      data.table::fwrite(rat, rat_path)
      message("📋 RAT weggeschreven naar: ", basename(rat_path))
    }
  }
  
  message("✅ Gereed: ", f_out)
}
