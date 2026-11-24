# ==============================================================================
# morfologische_opening.R
# Morfologische opening (erosie -> dilatatie) om smalle linten/randen te elimineren
# diameter_m: breedte waaronder structuren worden geëlimineerd (bv. 40m = 20m straal)
# ==============================================================================

hulp_morfologische_opening <- function(r_input, diameter_m = 40, master_grid) {
  
  r_bin <- terra::ifel(!is.na(r_input) & r_input > 0, 1L, NA_integer_)
  if (is_leeg(r_bin)) {
    return(master_grid * NA_real_)
  }
  
  f_eroded  <- tempfile(fileext = ".tif")
  f_dilated <- tempfile(fileext = ".tif")
  on.exit(unlink(c(f_eroded, f_dilated)), add = TRUE)
  
  # Venstergrootte bepalen op basis van werkelijke resolutie
  res_m       <- terra::res(master_grid)[1]
  cel_straal  <- max(1L, as.integer(round((diameter_m / 2) / res_m)))
  venster_dim <- (cel_straal * 2L) + 1L
  w_mat       <- matrix(1, nrow = venster_dim, ncol = venster_dim)
  
  # Stap A: Erosie (alle cellen in het venster moeten 1 zijn)
  r_eroded <- terra::focal(
    r_bin,
    w         = w_mat,
    fun       = "min",
    na.rm     = FALSE,
    filename  = f_eroded,
    overwrite = TRUE,
    gdal      = c("COMPRESS=LZW")
  )
  r_eroded_bin <- terra::ifel(!is.na(r_eroded) & r_eroded == 1, 1L, NA_integer_)
  rm(r_bin, r_eroded)
  
  # Early exit als er geen enkele brede kern overblijft
  if (is_leeg(r_eroded_bin)) {
    return(master_grid * NA_real_)
  }
  
  # Stap B: Dilatatie (terug uitbreiden vanaf behouden robuuste kernen)
  r_dilated <- terra::focal(
    r_eroded_bin,
    w         = w_mat,
    fun       = "max",
    na.rm     = TRUE,
    filename  = f_dilated,
    overwrite = TRUE,
    gdal      = c("COMPRESS=LZW")
  )
  
  mask_kern   <- terra::ifel(!is.na(r_dilated) & r_dilated == 1, 1L, NA_integer_)
  f_out       <- tempfile(fileext = ".tif")
  r_gezuiverd <- terra::mask(r_input, mask_kern, filename = f_out, overwrite = TRUE, datatype = "FLT4S", gdal = "COMPRESS=NONE")
  
  rm(r_eroded_bin, r_dilated, mask_kern)
  gc()
  
  return(r_gezuiverd)
}
