voeg_toe_geclamped <- function(r_basis, r_nieuw, master_grid = NULL) {
  has_b <- !is_leeg(r_basis)
  has_n <- !is_leeg(r_nieuw)
  
  if (!has_b && !has_n) {
    ref <- if (!is.null(master_grid)) master_grid else if (!is.null(r_basis)) r_basis else r_nieuw
    return(if (!is.null(ref)) ref * NA_real_ else NULL)
  }
  if (!has_b) return(r_nieuw)
  if (!has_n) return(r_basis)
  
  if (!terra::compareGeom(r_basis, r_nieuw, stopOnError = FALSE)) {
    r_nieuw <- terra::resample(r_nieuw, r_basis, method = "near")
  }
  
  f_tmp <- tempfile(fileext = ".tif")
  r_som <- terra::ifel(is.na(r_basis), 0, r_basis) + terra::ifel(is.na(r_nieuw), 0, r_nieuw)
  r_clamped <- terra::clamp(r_som, upper = 1.0)
  r_out <- terra::ifel(r_clamped > 0, r_clamped, NA_real_)
  
  r_vast <- terra::writeRaster(r_out, f_tmp, overwrite = TRUE, datatype = "FLT4S", gdal = "COMPRESS=NONE")
  rm(r_som, r_clamped, r_out)
  return(r_vast)
}
