# ==============================================================================
# hulp_wederzijdse_koppeling.R
# Wederzijdse interactiefilter via GDAL Proximity (symmetrisch én asymmetrisch)
# ==============================================================================

hulp_wederzijdse_koppeling <- function(laag_a, laag_b, dist_a_m, dist_b_m = dist_a_m, master_grid) {
  
  heeft_a <- !all(is.na(suppressWarnings(terra::minmax(laag_a))))
  heeft_b <- !all(is.na(suppressWarnings(terra::minmax(laag_b))))
  
  if (!heeft_a || !heeft_b || dist_a_m <= 0 || dist_b_m <= 0) {
    return(list(
      laag_a_gekoppeld = master_grid * NA_real_,
      laag_b_gekoppeld = master_grid * NA_real_
    ))
  }
  
  # ----------------------------------------------------------------------------
  # 1. Proximity vanaf B naar A (met zoekafstand dist_a_m)
  # ----------------------------------------------------------------------------
  f_b_bin  <- tempfile(fileext = ".tif")
  terra::writeRaster(terra::ifel(!is.na(laag_b) & laag_b > 0, 1, NA), f_b_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
  f_prox_b <- hulp_proximity_gdal(f_b_bin, dist_a_m)
  unlink(f_b_bin)
  
  p_b <- terra::rast(f_prox_b)
  mask_a <- terra::ifel(!is.na(p_b) & p_b >= 0 & p_b <= dist_a_m, 1, NA)
  unlink(f_prox_b)
  rm(p_b)
  
  laag_a_gekoppeld <- terra::mask(laag_a, mask_a)
  rm(mask_a)
  
  # Guardrail: controleer of Laag A cellen overhoudt
  if (all(is.na(suppressWarnings(terra::minmax(laag_a_gekoppeld))))) {
    return(list(
      laag_a_gekoppeld = master_grid * NA_real_,
      laag_b_gekoppeld = master_grid * NA_real_
    ))
  }
  
  # ----------------------------------------------------------------------------
  # 2. Proximity vanaf overgebleven A naar B (met zoekafstand dist_b_m)
  # ----------------------------------------------------------------------------
  f_a_bin  <- tempfile(fileext = ".tif")
  terra::writeRaster(terra::ifel(!is.na(laag_a_gekoppeld) & laag_a_gekoppeld > 0, 1, NA), f_a_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
  f_prox_a <- hulp_proximity_gdal(f_a_bin, dist_b_m)
  unlink(f_a_bin)
  
  p_a <- terra::rast(f_prox_a)
  mask_b <- terra::ifel(!is.na(p_a) & p_a >= 0 & p_a <= dist_b_m, 1, NA)
  unlink(f_prox_a)
  rm(p_a)
  
  laag_b_gekoppeld <- terra::mask(laag_b, mask_b)
  rm(mask_b)
  gc()
  
  return(list(
    laag_a_gekoppeld = laag_a_gekoppeld,
    laag_b_gekoppeld = laag_b_gekoppeld
  ))
}
