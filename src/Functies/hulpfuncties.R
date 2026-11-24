# ==============================================================================
# hulpfuncties.R - GDAL HELPERS VOOR STREAMING & PROXIMITY
# ==============================================================================

hulp_align_gdal <- function(f_in, master_grid) {
  e <- terra::ext(master_grid)
  f_out <- tempfile(fileext = ".tif")
  
  qgisprocess::qgis_run_algorithm(
    "gdal:warpreproject",
    INPUT             = f_in,
    TARGET_CRS        = "EPSG:31370",
    TARGET_EXTENT     = sprintf("%f,%f,%f,%f", e$xmin, e$xmax, e$ymin, e$ymax),
    TARGET_RESOLUTION = terra::res(master_grid)[1],
    RESAMPLING        = 0, # Nearest Neighbour
    CREATION_OPTIONS  = "COMPRESS=LZW",
    OUTPUT            = f_out,
    .quiet            = TRUE
  )
  return(f_out)
}

hulp_proximity_gdal <- function(f_in, max_dist_m) {
  f_out <- tempfile(fileext = ".tif")
  
  qgisprocess::qgis_run_algorithm(
    "gdal:proximity",
    INPUT        = f_in,
    BAND         = 1,
    VALUES       = "1",
    UNITS        = 0,
    MAX_DISTANCE = max_dist_m,
    NODATA       = 65535,
    DATA_TYPE    = 2, # UInt16 (tot 65.535 m)
    OUTPUT       = f_out,
    .quiet       = TRUE
  )
  return(f_out)
}
