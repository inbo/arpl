# ==============================================================================
# cluster_filter_universeel.R - DEFINITIEVE PRODUCTIEVERSIE
# ==============================================================================

cluster_filter_universeel <- function(masker, opp_laag, min_opp_ha, min_totale_opp_ha = NA, 
                                      afstand_m = 0, dispersie_m = 0, werkelijk = TRUE) {
  
  if (all(is.na(suppressWarnings(terra::minmax(masker))))) {
    return(list(raster = masker * NA, clusters = masker * NA))
  }
  
  min_opp_ha        <- as.numeric(min_opp_ha)
  min_totale_opp_ha <- as.numeric(min_totale_opp_ha)
  afstand_m         <- as.numeric(ifelse(is.na(afstand_m), 0, afstand_m))
  dispersie_m       <- as.numeric(ifelse(is.na(dispersie_m), 0, dispersie_m))
  
  # Bounds en extent string voor GRASS
  e <- terra::ext(masker)
  extent_str <- sprintf("%f,%f,%f,%f", e$xmin, e$xmax, e$ymin, e$ymax)
  
  # Hulpopdracht: GRASS r.clump uitvoeren en direct opschonen van categories/RAT
  run_grass_clump <- function(r_input_bin) {
    f_bin   <- file.path(tempdir(), paste0("clump_in_", as.numeric(Sys.time()), ".tif"))
    f_clump <- file.path(tempdir(), paste0("clump_out_", as.numeric(Sys.time()), ".tif"))
    
    terra::writeRaster(r_input_bin, f_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    
    qgisprocess::qgis_run_algorithm(
      "grass:r.clump",
      input                           = f_bin,
      output                          = f_clump,
      title                           = "clumps",
      "-d"                            = FALSE,
      GRASS_REGION_PARAMETER          = extent_str,
      GRASS_REGION_CELLSIZE_PARAMETER = 10,
      .quiet                          = TRUE
    )
    
    cl_res <- terra::rast(f_clump)
    # Verwijder categories/colortab zodat terra de laag als zuivere Int32 behandelt
    terra::set.cats(cl_res, layer = 1, value = NULL)
    terra::coltab(cl_res) <- NULL
    
    unlink(f_bin)
    return(cl_res)
  }
  
  # Hulpopdracht: Clusters vormen via veilige GDAL proximity (Float32, NODATA=9999)
  vorm_clusters <- function(r_in, dist) {
    dist <- as.numeric(dist)
    if (!is.na(dist) && dist > 0) {
      straal_m <- dist / 2
      f_in     <- file.path(tempdir(), paste0("prox_in_", as.numeric(Sys.time()), ".tif"))
      f_prox   <- file.path(tempdir(), paste0("prox_out_", as.numeric(Sys.time()), ".tif"))
      
      terra::writeRaster(r_in, f_in, overwrite = TRUE, datatype = "INT1U", NAflag = 0)
      
      qgisprocess::qgis_run_algorithm(
        "gdal:proximity",
        INPUT        = f_in,
        BAND         = 1,
        VALUES       = "1",
        UNITS        = 0,
        MAX_DISTANCE = straal_m,
        NODATA       = 9999,
        DATA_TYPE    = 2,            # Float32 voorkomt byte wrapping
        OPTIONS      = "COMPRESS=LZW",
        OUTPUT       = f_prox,
        .quiet       = TRUE
      )
      
      r_prox     <- terra::rast(f_prox)
      r_buffered <- terra::ifel(!is.na(r_prox) & r_prox >= 0 & r_prox <= straal_m, 1, NA)
      
      cl_raw <- run_grass_clump(r_buffered)
      
      unlink(c(f_in, f_prox))
      return(cl_raw)
    } else {
      cl_raw <- run_grass_clump(r_in)
      return(cl_raw)
    }
  }
  
  # Hulpopdracht: Oppervlaktefiltering via directe cel-extractie
  filter_op_oppervlakte <- function(cl_raster, r_filter_in, min_ha) {
    min_m2 <- as.numeric(min_ha) * 10000
    
    # Haal direct de actieve biotoopcellen op
    cellen_idx <- terra::cells(r_filter_in)
    if (length(cellen_idx) == 0) return(NULL)
    
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
    gc()
    return(voldoet_ids)
  }
  
  # 1. Binaire laag opzetten
  r_binair <- terra::ifel(!is.na(masker) & masker > 0, 1, NA)
  
  # ----------------------------------------------------------------------------
  # SCENARIO 1: Dubbele drempel (MinTotaleOpp_ha aanwezig)
  # ----------------------------------------------------------------------------
  if (!is.na(min_totale_opp_ha) && min_totale_opp_ha > 0) {
    cl_indiv  <- vorm_clusters(r_binair, dist = afstand_m)
    ids_indiv <- filter_op_oppervlakte(cl_indiv, r_binair, min_opp_ha)
    if (length(ids_indiv) == 0) return(list(raster = masker * NA, clusters = masker * NA))
    
    masker_stap1 <- terra::ifel(cl_indiv %in% ids_indiv, 1, NA)
    
    cl_netwerk <- vorm_clusters(masker_stap1, dist = dispersie_m)
    ids_net    <- filter_op_oppervlakte(cl_netwerk, masker_stap1, min_totale_opp_ha)
    if (length(ids_net) == 0) return(list(raster = masker * NA, clusters = masker * NA))
    
    cl_finaal <- terra::ifel(cl_netwerk %in% ids_net, cl_netwerk, NA)
    r_finaal  <- terra::mask(opp_laag, cl_finaal)
    r_finaal  <- terra::mask(r_finaal, masker_stap1)
    
    return(list(raster = r_finaal, clusters = cl_finaal))
  } 
  # ----------------------------------------------------------------------------
  # SCENARIO 2: Enkele drempel (Standaard soorten)
  # ----------------------------------------------------------------------------
  else {
    cl_clusters <- vorm_clusters(r_binair, dist = afstand_m)
    ids_ok      <- filter_op_oppervlakte(cl_clusters, r_binair, min_opp_ha)
    if (length(ids_ok) == 0) return(list(raster = masker * NA, clusters = masker * NA))
    
    cl_finaal <- terra::ifel(cl_clusters %in% ids_ok, cl_clusters, NA)
    # Maskeer terug naar de werkelijke biotoopcellen
    cl_finaal <- terra::mask(cl_finaal, r_binair)
    r_finaal  <- terra::mask(opp_laag, cl_finaal)
    
    return(list(raster = r_finaal, clusters = cl_finaal))
  }
}
