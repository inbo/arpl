maak_bebouwingsmasker <- function(min_opp_ha = 10, buf_m = 100, afstand_cl_m = 50, ds_parquet, master_grid, paths, force_overwrite = FALSE) {
  f_cache <- file.path(
    paths$dir_cache,
    sprintf("masker_bebouwing_%.1fha_%dm_%dm.tif", min_opp_ha, as.integer(afstand_cl_m), as.integer(buf_m))
  )
  
  if (file.exists(f_cache) && !force_overwrite) {
    test_r <- suppressWarnings(terra::rast(f_cache))
    if (terra::compareGeom(test_r, master_grid, stopOnError = FALSE)) {
      message(sprintf("       -> [CACHE] Bebouwingsbuffer (>= %.1f ha, %dm) ingeladen...", min_opp_ha, as.integer(buf_m)))
      return(test_r)
    }
  }
  
  message(sprintf("       -> Bebouwingsbuffer genereren (GDX >= %.1f ha op %dm, %dm buffer)...", 
                  min_opp_ha, as.integer(afstand_cl_m), as.integer(buf_m)))
  gdx_codes <- c("u", "ua", "ua-", "uc", "ud", "ui", "ui-", "un", "un-", "ur", "ur-", "ur+", "uv", "uv-", "uv+")
  
  r_bebouw_bin <- vul_raster_vanuit_parquet(
    ds_parquet  = ds_parquet,
    codes_exact = gdx_codes,
    master_grid = master_grid,
    waarde_vast = 1
  )
  
  res_b <- cluster_filter_universeel(
    masker            = r_bebouw_bin,
    opp_laag          = r_bebouw_bin,
    min_opp_ha        = min_opp_ha,
    min_totale_opp_ha = NA,
    afstand_m         = afstand_cl_m,
    dispersie_m       = 0,
    werkelijk         = FALSE
  )
  
  if (!is_leeg(res_b$raster)) {
    f_b_bin <- tempfile(fileext = ".tif")
    terra::writeRaster(res_b$raster, f_b_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    f_prox_b <- hulp_proximity_gdal(f_b_bin, buf_m)
    unlink(f_b_bin)
    
    p_b <- terra::rast(f_prox_b)
    m_bebouw_stoor <- terra::ifel(!is.na(p_b) & p_b >= 0 & p_b <= buf_m, 1L, NA_integer_)
    unlink(f_prox_b)
    rm(p_b)
  } else {
    m_bebouw_stoor <- master_grid * NA_real_
  }
  rm(res_b, r_bebouw_bin)
  
  terra::writeRaster(m_bebouw_stoor, f_cache, overwrite = TRUE, datatype = "INT1U", NAflag = 255, gdal = "COMPRESS=LZW")
  gc()
  return(m_bebouw_stoor)
}
