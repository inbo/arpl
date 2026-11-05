# ==============================================================================
# hulp_morfologische_opening.R
# Morfologische opening (erosie -> dilatatie) om smalle linten/randen te elimineren
# diameter_m: breedte waaronder structuren worden geëlimineerd (bv. 40m = 20m straal)
# ==============================================================================

hulp_morfologische_opening <- function(r_input, diameter_m = 40, master_grid) {
  
  r_bin <- terra::ifel(!is.na(r_input) & r_input > 0, 1, NA)
  if (all(is.na(suppressWarnings(terra::minmax(r_bin))))) {
    return(master_grid * NA_real_)
  }

# Aantal cellen bepalen voor focal window (10m resolutie)
# Voor diameter 40m (straal 20m) = venster van 5x5 cellen (2 cellen links, 2 rechts, 1 centrum)
cel_straal <- round((diameter_m / 2) / 10)
venster_dim <- (cel_straal * 2) + 1
w_mat <- matrix(1, nrow = venster_dim, ncol = venster_dim)

# Stap A: Erosie (alle cellen in het venster moeten 1 zijn)
f_eroded <- tempfile(fileext = ".tif")
r_eroded <- terra::focal(
  r_bin, 
  w        = w_mat, 
  fun      = "min", 
  na.rm    = FALSE, 
  filename = f_eroded, 
  overwrite = TRUE, 
  gdal     = c("COMPRESS=LZW")
)
r_eroded_bin <- terra::ifel(!is.na(r_eroded) & r_eroded == 1, 1, NA)

# Stap B: Dilatatie (terug uitbreiden vanaf behouden robuuste kernen)
f_dilated <- tempfile(fileext = ".tif")
r_dilated <- terra::focal(
  r_eroded_bin, 
  w        = w_mat, 
  fun      = "max", 
  na.rm    = TRUE, 
  filename = f_dilated, 
  overwrite = TRUE, 
  gdal     = c("COMPRESS=LZW")
)

mask_kern <- terra::ifel(!is.na(r_dilated) & r_dilated == 1, 1, NA)
r_gezuiverd <- terra::mask(r_input, mask_kern)

unlink(c(f_eroded, f_dilated))
rm(r_bin, r_eroded, r_eroded_bin, r_dilated, mask_kern)
gc()

return(r_gezuiverd)
}