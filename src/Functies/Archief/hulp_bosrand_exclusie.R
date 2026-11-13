# ==============================================================================
# hulp_bosrand_exclusie.R
# Bouwt gecombineerde boslaag (Groenkaart + BWK-bos), filtert op minimumoppervlakte
# en snijdt de proximity-buffer weg uit het open doelbiotoop.
# ==============================================================================

hulp_bosrand_exclusie <- function(r_doel, 
                                  r_bos_bwk, 
                                  pad_groenkaart, 
                                  bos_drempel_m2 = 2000, 
                                  buffer_dist_m  = 100, 
                                  master_grid) {
  
  bos_min_cellen <- max(1, round(bos_drempel_m2 / 100))
  
  # A. Groenkaart hooggroen (waarde 1)
  if (!is.null(pad_groenkaart) && file.exists(pad_groenkaart)) {
    f_groen_al     <- hulp_align_gdal(pad_groenkaart, master_grid)
    r_groen        <- terra::rast(f_groen_al)
    groenkaart_bos <- terra::ifel(!is.na(r_groen) & r_groen == 1, 1, NA)
    unlink(f_groen_al)
    rm(r_groen)
  } else {
    groenkaart_bos <- master_grid * NA_real_
  }
  
  # B. BWK-bos (overlap met doelbiotoop uitsluiten)
  if (!is.null(r_bos_bwk) && !all(is.na(suppressWarnings(terra::minmax(r_bos_bwk))))) {
    bwkbos_corr <- terra::ifel(!is.na(r_doel), NA, r_bos_bwk)
    bwkbos_bin  <- terra::ifel(!is.na(bwkbos_corr) & bwkbos_corr > 0, 1, NA)
    rm(bwkbos_corr)
  } else {
    bwkbos_bin  <- master_grid * NA_real_
  }
  
  # C. Samenvoegen
  bos_totaal_ruw <- terra::cover(bwkbos_bin, groenkaart_bos)
  rm(bwkbos_bin, groenkaart_bos)
  gc()
  
  if (all(is.na(suppressWarnings(terra::minmax(bos_totaal_ruw))))) {
    return(r_doel)
  }
  
  # D. Clumpen via GRASS om boskernen >= drempel te behouden
  f_bos_bin   <- tempfile(fileext = ".tif")
  f_clump_bos <- tempfile(fileext = ".tif")
  terra::writeRaster(bos_totaal_ruw, f_bos_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
  
  e <- terra::ext(master_grid)
  qgisprocess::qgis_run_algorithm(
    "grass:r.clump", 
    input                           = f_bos_bin, 
    output                          = f_clump_bos, 
    title                           = "cl_bos_exclusie",
    "-d"                            = FALSE, 
    GRASS_REGION_PARAMETER          = sprintf("%f,%f,%f,%f", e$xmin, e$xmax, e$ymin, e$ymax),
    GRASS_REGION_CELLSIZE_PARAMETER = 10, 
    .quiet                          = TRUE
  )
  unlink(f_bos_bin)
  
  cl_bos <- terra::rast(f_clump_bos)
  try(terra::set.cats(cl_bos, layer = 1, value = NULL), silent = TRUE)
  try(terra::coltab(cl_bos) <- NULL, silent = TRUE)
  cl_bos_only <- terra::mask(cl_bos, bos_totaal_ruw)
  unlink(f_clump_bos)
  rm(cl_bos, bos_totaal_ruw)
  gc()
  
  dt_bos <- data.table::as.data.table(terra::freq(cl_bos_only))
  data.table::setnames(dt_bos, c("layer", "ID", "count"))
  
  grote_bos_ids <- dt_bos[!is.na(ID) & ID > 0 & count >= bos_min_cellen, ID]
  rm(dt_bos)
  
  if (length(grote_bos_ids) == 0) {
    rm(cl_bos_only)
    return(r_doel)
  }
  
  bos_gefilterd <- terra::ifel(cl_bos_only %in% as.integer(grote_bos_ids), 1, NA)
  rm(cl_bos_only)
  gc()
  
  # E. GDAL Proximity buffer (100m)
  f_b_bin <- tempfile(fileext = ".tif")
  terra::writeRaster(bos_gefilterd, f_b_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
  rm(bos_gefilterd)
  
  f_prox_bos <- hulp_proximity_gdal(f_b_bin, buffer_dist_m)
  unlink(f_b_bin)
  
  p_bos <- terra::rast(f_prox_bos)
  mask_storing <- terra::ifel(!is.na(p_bos) & p_bos >= 0 & p_bos <= buffer_dist_m, 1, NA)
  unlink(f_prox_bos)
  rm(p_bos)
  
  # F. Wegsnijden uit doelbiotoop
  r_doel_open <- terra::mask(r_doel, mask_storing, inverse = TRUE)
  rm(mask_storing)
  gc()
  
  return(r_doel_open)
}