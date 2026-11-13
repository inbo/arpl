# ==============================================================================
# bouw_biotoop_bosrand_lengte.R
# Archetype: Halfopen habitat x boskerndrempel x minimale bosrandlengte (>= 1.5 km)
# Geschikt voor: Nachtzwaluw, etc.
# ==============================================================================

bouw_biotoop_bosrand_lengte <- function(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
  bestandsnaam <- tolower(gsub(" ", "", soort_naam))
  soort_cache  <- file.path(paths$dir_cache, bestandsnaam)
  dir.create(soort_cache, showWarnings = FALSE, recursive = TRUE)
  
  f_cache_rds <- file.path(soort_cache, "01_raw_parquet_layers.rds")
  
  # ----------------------------------------------------------------------------
  # 1. PARQUET EXTRACTIE (met persistente cache)
  # ----------------------------------------------------------------------------
  types_in_tabel <- unique(params_soort$Type)
  
  if (file.exists(f_cache_rds)) {
    message("    📦 [CACHE] Ruwe Parquet-lagen ingeladen van schijf...")
    lagen_packed <- readRDS(f_cache_rds)
    lagen <- lapply(lagen_packed, terra::unwrap)
  } else {
    message(paste0("    ⏳ Parquet-lagen extraheren voor ", soort_naam, "..."))
    lagen <- list()
    for (tp in types_in_tabel) {
      r_tp <- extract_bwk_laag(soort_naam, tp, ds_parquet, df_bwk_codes, master_grid)
      lagen[[tp]] <- if (!is.null(r_tp)) r_tp else (master_grid * NA)
    }
    # Veilig inpakken vóór wegschrijven naar RDS
    lagen_packed <- lapply(lagen, terra::wrap)
    saveRDS(lagen_packed, f_cache_rds)
    message("    💾 [CACHE] Parquet-lagen veilig bewaard.")
  }
  
  r_vp_raw  <- if ("voortplanting" %in% names(lagen)) lagen[["voortplanting"]] else (master_grid * NA)
  r_bos_raw <- if ("bos"           %in% names(lagen)) lagen[["bos"]]           else (master_grid * NA)
  
  # ----------------------------------------------------------------------------
  # 2. BOSCLUSTERING & 5 HA KERNDREMPEL (50M ZOEFAFSTAND)
  # ----------------------------------------------------------------------------
  message("    🌲 Boscomplexen clusteren en filteren op >= 5 ha (50m netwerk)...")
  res_bos <- cluster_filter_universeel(
    masker            = r_bos_raw,
    opp_laag          = r_bos_raw,
    min_opp_ha        = 5,
    min_totale_opp_ha = NA,
    afstand_m         = 50,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_bos_5ha <- res_bos$raster
  rm(res_bos, r_bos_raw)
  gc()
  
  ha_bos_5ha <- calc_ha_exact(r_bos_5ha)
  message(sprintf("       -> Goedgekeurde boscomplexen (>= 5 ha): %.2f ha", ha_bos_5ha))
  
  # ----------------------------------------------------------------------------
  # 3. BOSRAND DETECTIE & CONTACTLENGTE-TOETS (MINSTENS 1.500 M / 150 PIXELS)
  # ----------------------------------------------------------------------------
  message("    📏 Bosranden detecteren en open biotoop toetsen op >= 1.500m bosrand...")
  
  heeft_bos <- !all(is.na(suppressWarnings(terra::minmax(r_bos_5ha))))
  heeft_vp  <- !all(is.na(suppressWarnings(terra::minmax(r_vp_raw))))
  
  if (heeft_bos && heeft_vp) {
    # Directe bosrand identificeren via 3x3 window (10m overlap)
    r_bos_bin <- terra::ifel(!is.na(r_bos_5ha) & r_bos_5ha > 0, 1, 0)
    w_3x3     <- matrix(1, nrow = 3, ncol = 3)
    
    f_focal_b <- tempfile(fileext = ".tif")
    r_bos_buf <- terra::focal(r_bos_bin, w = w_3x3, fun = "max", na.rm = TRUE, filename = f_focal_b, overwrite = TRUE)
    unlink(f_focal_b)
    rm(r_bos_bin)
    
    r_vp_bin   <- terra::ifel(!is.na(r_vp_raw) & r_vp_raw > 0, 1, 0)
    r_rand_bin <- terra::ifel(r_vp_bin == 1 & r_bos_buf == 1, 1, 0)
    rm(r_bos_buf)
    
    # Open biotoop groeperen met 25m buffer om functionele heide/kapvlaktecomplexen te vormen
    r_vp_alleen_bin <- terra::ifel(r_vp_bin == 1, 1, NA)
    rm(r_vp_bin)
    
    f_vp_bin   <- tempfile(fileext = ".tif")
    terra::writeRaster(r_vp_alleen_bin, f_vp_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 0)
    f_prox_vp  <- hulp_proximity_gdal(f_vp_bin, 25)
    unlink(f_vp_bin)
    
    p_vp       <- terra::rast(f_prox_vp)
    mask_vp25  <- terra::ifel(!is.na(p_vp) & p_vp >= 0 & p_vp <= 25, 1, NA)
    unlink(f_prox_vp)
    rm(p_vp)
    gc()
    
    # Clump van de open complexen via GRASS
    f_buf_bin <- tempfile(fileext = ".tif")
    f_cl_vp   <- tempfile(fileext = ".tif")
    terra::writeRaster(mask_vp25, f_buf_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    rm(mask_vp25)
    gc()
    
    e <- terra::ext(master_grid)
    qgisprocess::qgis_run_algorithm(
      "grass:r.clump", input = f_buf_bin, output = f_cl_vp, title = "cl_vp25m",
      "-d" = FALSE, GRASS_REGION_PARAMETER = sprintf("%f,%f,%f,%f", e$xmin, e$xmax, e$ymin, e$ymax),
      GRASS_REGION_CELLSIZE_PARAMETER = 10, .quiet = TRUE
    )
    unlink(f_buf_bin)
    
    cl_complex <- terra::rast(f_cl_vp)
    terra::set.cats(cl_complex, layer = 1, value = NULL)
    terra::coltab(cl_complex) <- NULL
    
    # Beperk complex-ID's tot de werkelijke biotoopcellen
    cl_complex_vp <- terra::mask(cl_complex, r_vp_alleen_bin)
    unlink(f_cl_vp)
    rm(r_vp_alleen_bin)
    gc()
    
    # Tel bosrandpixels per open complex via data.table
    cells_c <- terra::cells(!is.na(cl_complex_vp) & cl_complex_vp > 0)
    dt_rand <- data.table::data.table(
      ID         = as.integer(cl_complex_vp[cells_c][[1]]),
      Rand_Pixel = as.integer(r_rand_bin[cells_c][[1]])
    )
    rm(cells_c, r_rand_bin)
    gc()
    
    dt_stats <- dt_rand[, .(Totaal_Rand = sum(Rand_Pixel, na.rm = TRUE)), by = ID]
    rm(dt_rand)
    
    # Minstens 150 pixels = 1.500m bosrandlengte
    ok_complex_ids <- dt_stats[Totaal_Rand >= 150, ID]
    rm(dt_stats)
    
    message(sprintf("       -> Open biotoopcomplexen met >= 1.500m bosrand: %d", length(ok_complex_ids)))
    
    if (length(ok_complex_ids) > 0) {
      mask_ok <- terra::ifel(cl_complex_vp %in% as.integer(ok_complex_ids), 1, NA)
      r_vp_gefilterd <- terra::mask(r_vp_raw, mask_ok)
      rm(mask_ok)
    } else {
      r_vp_gefilterd <- master_grid * NA
    }
    rm(cl_complex_vp)
    gc()
  } else {
    r_vp_gefilterd <- master_grid * NA
  }
  rm(r_bos_5ha, r_vp_raw)
  gc()
  
  ha_vp_geflt <- calc_ha_exact(r_vp_gefilterd)
  message(sprintf("       -> Voortplantingsbiotoop met voldoende bosrand: %.2f ha", ha_vp_geflt))
  
  # ----------------------------------------------------------------------------
  # 4. FINALE CLUSTERING OP 25 HA (50M ZOEFAFSTAND)
  # ----------------------------------------------------------------------------
  message("    🧩 Finale leefgebiedsclustering (>= 25 ha, 50m netwerk)...")
  res_finaal <- cluster_filter_universeel(
    masker            = r_vp_gefilterd,
    opp_laag          = r_vp_gefilterd,
    min_opp_ha        = 25,
    min_totale_opp_ha = NA,
    afstand_m         = 50,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  final_opp <- res_finaal$raster
  cl_opp    <- res_finaal$clusters
  rm(res_finaal, r_vp_gefilterd)
  gc()
  
  return(list(
    final_opp = final_opp,
    cl_opp    = cl_opp
  ))
}