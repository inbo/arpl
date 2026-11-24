clump_fast <- function(r_in, master_grid = NULL, diag = FALSE, titel = "cl_tmp") {
  ref_grid <- if (!is.null(master_grid)) master_grid else r_in
  
  f_in  <- tempfile(fileext = ".tif")
  f_out <- tempfile(fileext = ".tif")
  terra::writeRaster(terra::ifel(!is.na(r_in) & r_in > 0, 1L, NA_integer_), f_in, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
  
  e <- terra::ext(ref_grid)
  qgisprocess::qgis_run_algorithm(
    "grass:r.clump",
    input                           = f_in,
    output                          = f_out,
    title                           = titel,
    "-d"                            = diag,
    GRASS_REGION_PARAMETER          = sprintf("%f,%f,%f,%f", e$xmin, e$xmax, e$ymin, e$ymax),
    GRASS_REGION_CELLSIZE_PARAMETER = terra::res(ref_grid)[1],
    GRASS_RASTER_FORMAT_OPT         = "COMPRESS=LZW",
    .quiet                          = TRUE
  )
  unlink(f_in)
  
  cl <- terra::rast(f_out)
  try(terra::set.cats(cl, layer = 1, value = NULL), silent = TRUE)
  try(terra::coltab(cl) <- NULL, silent = TRUE)
  attr(cl, "tmp_file") <- f_out
  return(cl)
}
