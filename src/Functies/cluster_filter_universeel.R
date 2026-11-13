# ==============================================================================
# cluster_filter_universeel.R - GEOPTIMALISEERDE PRODUCTIEVERSIE (QGIS/GRASS)
# ==============================================================================

cluster_filter_universeel <- function(masker, opp_laag, min_opp_ha, min_totale_opp_ha = NA, 
                                      afstand_m = 0, dispersie_m = 0, werkelijk = TRUE) {
  
  if (is.null(masker)) return(list(raster = NULL, clusters = NULL))
  if (!terra::hasMinMax(masker)) suppressWarnings(terra::setMinMax(masker))
  if (all(is.na(terra::minmax(masker)))) {
    return(list(raster = masker * NA_real_, clusters = masker * NA_real_))
  }
  
  min_opp_ha        <- as.numeric(min_opp_ha)
  min_totale_opp_ha <- as.numeric(min_totale_opp_ha)
  afstand_m         <- as.numeric(ifelse(is.na(afstand_m), 0, afstand_m))
  dispersie_m       <- as.numeric(ifelse(is.na(dispersie_m), 0, dispersie_m))
  
  # Houd tijdelijke GRASS-bestanden bij en ruim ze automatisch op bij afsluiten functie
  tmp_bestanden <- character(0)
  on.exit(unlink(tmp_bestanden[file.exists(tmp_bestanden)]), add = TRUE)
  
  e <- terra::ext(masker)
  extent_str <- sprintf("%f,%f,%f,%f", e$xmin, e$xmax, e$ymin, e$ymax)
  
  # ----------------------------------------------------------------------------
  # Hulpopdracht: GRASS r.clump uitvoeren (gecomprimeerd + -c vlag)
  # ----------------------------------------------------------------------------
  run_grass_clump <- function(r_input_bin) {
    f_bin   <- tempfile(pattern = "clump_in_", fileext = ".tif")
    f_clump <- tempfile(pattern = "clump_out_", fileext = ".tif")
    tmp_bestanden <<- c(tmp_bestanden, f_bin, f_clump)
    
    terra::writeRaster(r_input_bin, f_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    
    qgisprocess::qgis_run_algorithm(
      "grass:r.clump",
      input                           = f_bin,
      output                          = f_clump,
      title                           = "clumps",
      "-d"                            = FALSE,
      "-c"                            = TRUE,
      GRASS_REGION_PARAMETER          = extent_str,
      GRASS_REGION_CELLSIZE_PARAMETER = 10,
      GRASS_RASTER_FORMAT_OPT         = "COMPRESS=LZW",
      .quiet                          = TRUE
    )
    
    unlink(f_bin)
    
    cl_res <- terra::rast(f_clump)
    terra::set.cats(cl_res, layer = 1, value = NULL)
    terra::coltab(cl_res) <- NULL
    return(cl_res)
  }
  
  # ----------------------------------------------------------------------------
  # Hulpopdracht: Clusters vormen via GDAL proximity (Byte/PACKBITS indien <= 250m)
  # ----------------------------------------------------------------------------
  vorm_clusters <- function(r_in, dist) {
    dist <- as.numeric(dist)
    if (!is.na(dist) && dist > 0) {
      straal_m <- dist / 2
      f_in     <- tempfile(pattern = "prox_in_", fileext = ".tif")
      f_prox   <- tempfile(pattern = "prox_out_", fileext = ".tif")
      tmp_bestanden <<- c(tmp_bestanden, f_in, f_prox)
      
      terra::writeRaster(r_in, f_in, overwrite = TRUE, datatype = "INT1U", NAflag = 0)
      
      use_byte <- straal_m <= 250
      qgisprocess::qgis_run_algorithm(
        "gdal:proximity",
        INPUT        = f_in,
        BAND         = 1,
        VALUES       = "1",
        UNITS        = 0,
        MAX_DISTANCE = straal_m,
        NODATA       = if (use_byte) 255 else 9999,
        DATA_TYPE    = if (use_byte) 0 else 2,
        OPTIONS      = "COMPRESS=PACKBITS",
        OUTPUT       = f_prox,
        .quiet       = TRUE
      )
      
      r_prox     <- terra::rast(f_prox)
      r_buffered <- terra::ifel(!is.na(r_prox) & r_prox <= straal_m, 1, NA)
      
      cl_raw <- run_grass_clump(r_buffered)
      unlink(c(f_in, f_prox))
      return(cl_raw)
    } else {
      return(run_grass_clump(r_in))
    }
  }
  
  # ----------------------------------------------------------------------------
  # Hulpopdracht: Oppervlaktefiltering + directe sparse raster-opbouw
  # ----------------------------------------------------------------------------
  filter_en_bouw_uitvoer <- function(cl_raster, r_filter_in, min_ha, return_cluster_ids = TRUE) {
    min_m2 <- as.numeric(min_ha) * 10000
    
    cellen_idx <- terra::cells(r_filter_in)
    if (length(cellen_idx) == 0) {
      return(list(raster = masker * NA_real_, clusters = masker * NA_real_, n_ok = 0L))
    }
    
    c_ids <- as.integer(cl_raster[cellen_idx][[1]])
    
    if (werkelijk) {
      c_vals <- as.numeric(opp_laag[cellen_idx][[1]])
      dt <- data.table::data.table(ID = c_ids, Waarde = c_vals)
      stats_dt <- dt[!is.na(ID) & ID > 0 & !is.na(Waarde), .(Area_m2 = sum(Waarde, na.rm = TRUE) * 100), by = ID]
    } else {
      dt <- data.table::data.table(ID = c_ids)
      stats_dt <- dt[!is.na(ID) & ID > 0, .(Area_m2 = .N * 100), by = ID]
    }
    
    voldoet_ids <- stats_dt[Area_m2 >= min_m2, ID]
    rm(dt, stats_dt)
    
    if (length(voldoet_ids) == 0) {
      return(list(raster = masker * NA_real_, clusters = masker * NA_real_, n_ok = 0L))
    }
    
    keep <- !is.na(c_ids) & (c_ids %in% voldoet_ids)
    valid_cells <- cellen_idx[keep]
    
    cl_out <- masker * NA_integer_
    cl_out[valid_cells] <- if (return_cluster_ids) c_ids[keep] else 1L
    
    r_out <- masker * NA_real_
    r_out[valid_cells] <- if (werkelijk) c_vals[keep] else as.numeric(opp_laag[valid_cells][[1]])
    
    return(list(raster = r_out, clusters = cl_out, n_ok = length(voldoet_ids)))
  }
  
  # 1. Binaire laag opzetten
  r_binair <- terra::ifel(!is.na(masker) & masker > 0, 1, NA)
  
  # ----------------------------------------------------------------------------
  # SCENARIO 1: Dubbele drempel (MinTotaleOpp_ha aanwezig)
  # ----------------------------------------------------------------------------
  if (!is.na(min_totale_opp_ha) && min_totale_opp_ha > 0) {
    cl_indiv <- vorm_clusters(r_binair, dist = afstand_m)
    res_s1   <- filter_en_bouw_uitvoer(cl_indiv, r_binair, min_opp_ha, return_cluster_ids = FALSE)
    if (res_s1$n_ok == 0L) return(list(raster = masker * NA_real_, clusters = masker * NA_real_))
    
    masker_stap1 <- res_s1$clusters
    cl_netwerk   <- vorm_clusters(masker_stap1, dist = dispersie_m)
    res_s2       <- filter_en_bouw_uitvoer(cl_netwerk, masker_stap1, min_totale_opp_ha, return_cluster_ids = TRUE)
    
    return(list(raster = res_s2$raster, clusters = res_s2$clusters))
  } 
  # ----------------------------------------------------------------------------
  # SCENARIO 2: Enkele drempel (Standaard soorten)
  # ----------------------------------------------------------------------------
  else {
    cl_clusters <- vorm_clusters(r_binair, dist = afstand_m)
    res <- filter_en_bouw_uitvoer(cl_clusters, r_binair, min_opp_ha, return_cluster_ids = TRUE)
    return(list(raster = res$raster, clusters = res$clusters))
  }
}
