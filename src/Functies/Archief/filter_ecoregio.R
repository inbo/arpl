# ==============================================================================
# filter_ecoregio.R - GENERIEKE FILTER VOOR ECOREGIO'S / ECODISTRICTEN
# ==============================================================================

filter_ecoregio <- function(raster_input, 
                            pad_ecoregio, 
                            master_grid, 
                            toegelaten_regios = NULL, 
                            uitgesloten_regios = NULL,
                            cluster_raster = NULL) {
  
  if (is.null(raster_input) || all(is.na(suppressWarnings(terra::minmax(raster_input))))) {
    return(raster_input)
  }
  
  if (!file.exists(pad_ecoregio)) {
    message("    ⚠️ [ECOREGIO] Bestand niet gevonden op pad: ", pad_ecoregio, ". Filter overgeslagen.")
    return(raster_input)
  }
  
  message("    🗺️ Ecoregio-filter toepassen...")
  
  # 1. Ecoregio-laag aligneren via GDAL warp
  f_eco_al <- hulp_align_gdal(pad_ecoregio, master_grid)
  r_eco    <- terra::rast(f_eco_al)
  
  # 2. Inspecteer categorieën / labels
  leg_eco  <- terra::cats(r_eco)[[1]]
  
  # Identificeer numerieke ID's op basis van tekstuele labels of numerieke invoer
  bepaal_ids <- function(selectie, leg) {
    if (is.numeric(selectie)) return(selectie)
    if (is.null(leg)) return(as.numeric(selectie))
    
    kol_val   <- colnames(leg)[tolower(colnames(leg)) %in% c("value", "id")][1]
    kol_label <- colnames(leg)[tolower(colnames(leg)) %in% c("label", "naam", "ecoregio", "ecodistrict")][1]
    
    if (is.na(kol_val)) kol_val <- 1
    if (is.na(kol_label)) return(as.numeric(selectie))
    
    match_ids <- leg[[kol_val]][tolower(trimws(as.character(leg[[kol_label]]))) %in% tolower(trimws(selectie))]
    return(match_ids)
  }
  
  # 3. Bouw binair ecoregio-masker
  if (!is.null(toegelaten_regios)) {
    ok_ids <- bepaal_ids(toegelaten_regios, leg_eco)
    mask_eco <- terra::ifel(!is.na(r_eco) & r_eco %in% ok_ids, 1, NA)
  } else if (!is.null(uitgesloten_regios)) {
    nok_ids <- bepaal_ids(uitgesloten_regios, leg_eco)
    mask_eco <- terra::ifel(!is.na(r_eco) & !(r_eco %in% nok_ids), 1, NA)
  } else {
    unlink(f_eco_al)
    rm(r_eco)
    return(raster_input)
  }
  
  unlink(f_eco_al)
  rm(r_eco)
  
  # 4. Toepassing: direct op cellen óf op hele clusters
  if (is.null(cluster_raster)) {
    # Direct cell-by-cell maskeren
    res <- terra::mask(raster_input, mask_eco)
  } else {
    # Cluster-niveau: als een cluster minstens 1 cel in de toegelaten ecoregio heeft, behoud de hele cluster
    df_over <- terra::zonal(mask_eco, cluster_raster, fun = "max", na.rm = TRUE)
    colnames(df_over) <- c("Cluster_ID", "In_Regio")
    ok_clusters <- df_over$Cluster_ID[!is.na(df_over$In_Regio) & df_over$In_Regio == 1]
    
    if (length(ok_clusters) > 0) {
      cl_mask <- terra::ifel(cluster_raster %in% ok_clusters, 1, NA)
      res     <- terra::mask(raster_input, cl_mask)
    } else {
      res     <- master_grid * NA
    }
  }
  
  ha_na_eco <- calc_ha_exact(res)
  message(sprintf("       -> Biotoop na ecoregio-filter: %.2f ha", ha_na_eco))
  
  return(res)
}
