haal_drainage_masker <- function(pad_drainage, master_grid, dir_cache, 
                                 drain_letters = c("d", "e", "f", "h", "i", "g", "e-f", "h-i", "e-i")) {
  
  f_cache <- file.path(dir_cache, "masker_drainage_nat.tif")
  if (file.exists(f_cache)) {
    return(terra::rast(f_cache))
  }
  
  if (!file.exists(pad_drainage)) return(master_grid * NA)
  
  # 1. Haal categorieën uit origineel bestand (vóór GDAL warp)
  r_orig   <- terra::rast(pad_drainage)
  leg      <- terra::cats(r_orig)[[1]]
  val_kol  <- if (!is.null(leg) && "Value" %in% colnames(leg)) "Value" else 1
  
  if (!is.null(leg) && "Label" %in% colnames(leg)) {
    ok_ids <- leg[[val_kol]][tolower(trimws(as.character(leg$Label))) %in% tolower(drain_letters)]
  } else {
    ok_ids <- integer(0)
  }
  rm(r_orig)
  
  # 2. Warp met GDAL naar mastergrid
  f_aligned <- hulp_align_gdal(pad_drainage, master_grid)
  r_drain   <- terra::rast(f_aligned)
  
  if (length(ok_ids) > 0) {
    mask_drain <- terra::ifel(!is.na(r_drain) & r_drain %in% ok_ids, 1, NA)
  } else {
    mask_drain <- master_grid * NA
  }
  
  unlink(f_aligned)
  rm(r_drain)
  
  # 3. Cache bewaren voor de volgende soorten
  terra::writeRaster(mask_drain, f_cache, overwrite = TRUE, datatype = "INT1U", gdal = "COMPRESS=LZW")
  return(mask_drain)
}