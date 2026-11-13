# ==============================================================================
# laad_of_extraheer_parquet_lagen.R
# Robuuste Parquet extractie met fysieke GeoTIFF caching per soort en laagtype
# ==============================================================================

laad_of_extraheer_parquet_lagen <- function(soort_naam, types_te_laden, ds_parquet, df_bwk_codes, master_grid, soort_cache) {
  lagen <- list()
  
  for (tp in types_te_laden) {
    f_tif <- file.path(soort_cache, paste0("laag_", tp, ".tif"))
    
    if (file.exists(f_tif)) {
      message(paste0("    📦 [CACHE] Laag '", tp, "' ingeladen van schijf..."))
      lagen[[tp]] <- terra::rast(f_tif)
    } else {
      message(paste0("    ⏳ Parquet-laag '", tp, "' extraheren voor ", soort_naam, "..."))
      r_tp <- extract_bwk_laag(soort_naam, tp, ds_parquet, df_bwk_codes, master_grid)
      if (is.null(r_tp)) r_tp <- master_grid * NA_real_
      
      terra::writeRaster(
        r_tp, 
        filename  = f_tif, 
        overwrite = TRUE, 
        datatype  = "FLT4S", 
        gdal      = c("COMPRESS=LZW")
      )
      lagen[[tp]] <- terra::rast(f_tif)
    }
  }
  
  return(lagen)
}