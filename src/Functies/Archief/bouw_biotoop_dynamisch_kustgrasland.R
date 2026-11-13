# ==============================================================================
# bouw_biotoop_dynamisch_kustgrasland.R
# Archetype voor kust-/weidevogels met dynamische oppervlaktedrempel (5-10 ha)
# en uitsluiting van dominante vennen (>95% water)
# Geschikt voor: Tureluur, etc.
# ==============================================================================

bouw_biotoop_dynamisch_kustgrasland <- function(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
  bestandsnaam <- tolower(gsub(" ", "", soort_naam))
  soort_cache  <- file.path(paths$dir_cache, bestandsnaam)
  dir.create(soort_cache, showWarnings = FALSE, recursive = TRUE)
  
  f_cache_rds <- file.path(soort_cache, "01_raw_parquet_layers.rds")
  
  # Parameters ophalen
  rij_bwk       <- params_soort %>% dplyr::filter(tolower(Type) == "bwk") %>% dplyr::slice(1)
  zoekafstand_m <- as.numeric(ifelse(nrow(rij_bwk) > 0 && !is.na(rij_bwk$AfstandBiotopen_m), rij_bwk$AfstandBiotopen_m, 100))
  
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
    lagen_packed <- lapply(lagen, terra::wrap)
    saveRDS(lagen_packed, f_cache_rds)
    message("    💾 [CACHE] Parquet-lagen veilig bewaard.")
  }
  
  r_bwk_raw    <- if ("bwk"          %in% names(lagen)) lagen[["bwk"]]          else (master_grid * NA)
  r_slik_raw   <- if ("slik_schor"   %in% names(lagen)) lagen[["slik_schor"]]   else (master_grid * NA)
  r_vennen_raw <- if ("vennen_water" %in% names(lagen)) lagen[["vennen_water"]] else (master_grid * NA)
  
  ha_totaal_ruw <- calc_ha_exact(r_bwk_raw)
  message(sprintf("       -> Doelbiotoop vóór clustering: %.2f ha", ha_totaal_ruw))
  
  # ----------------------------------------------------------------------------
  # 2. CLUSTEREN OP 100M NETWERK (GRASS R.CLUMP)
  # ----------------------------------------------------------------------------
  message(sprintf("    🏝️ Biotoop clusteren op %dm zoekafstand...", as.integer(zoekafstand_m)))
  
  r_bin <- terra::ifel(!is.na(r_bwk_raw) & r_bwk_raw > 0, 1, NA)
  
  if (all(is.na(suppressWarnings(terra::minmax(r_bin))))) {
    return(list(final_opp = master_grid * NA, cl_opp = master_grid * NA))
  }
  
  f_bin   <- tempfile(fileext = ".tif")
  f_clump <- tempfile(fileext = ".tif")
  
  # Proximity buffer voor zoekafstand
  if (zoekafstand_m > 0) {
    terra::writeRaster(r_bin, f_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    f_prox <- hulp_proximity_gdal(f_bin, zoekafstand_m / 2)
    unlink(f_bin)
    
    p_grid <- terra::rast(f_prox)
    r_fuzzy <- terra::ifel(!is.na(p_grid) & p_grid >= 0 & p_grid <= (zoekafstand_m / 2), 1, NA)
    unlink(f_prox)
    rm(p_grid)
    
    terra::writeRaster(r_fuzzy, f_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    rm(r_fuzzy)
  } else {
    terra::writeRaster(r_bin, f_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
  }
  
  e <- terra::ext(master_grid)
  qgisprocess::qgis_run_algorithm(
    "grass:r.clump", input = f_bin, output = f_clump, title = paste0("cl_", bestandsnaam),
    "-d" = FALSE, GRASS_REGION_PARAMETER = sprintf("%f,%f,%f,%f", e$xmin, e$xmax, e$ymin, e$ymax),
    GRASS_REGION_CELLSIZE_PARAMETER = 10, .quiet = TRUE
  )
  unlink(f_bin)
  
  cl_raw <- terra::rast(f_clump)
  terra::set.cats(cl_raw, layer = 1, value = NULL)
  if (terra::has.colors(cl_raw)) {
    try(terra::coltab(cl_raw) <- NULL, silent = TRUE)
  }
  cl_biotoop <- terra::mask(cl_raw, r_bin)
  unlink(f_clump)
  rm(cl_raw, r_bin)
  gc()
  
  # ----------------------------------------------------------------------------
  # 3. DYNAMISCHE OPPERVLAKTEFILTER (5 HA SLIK <-> 10 HA GRASLAND) & VENNENTOETS (<95%)
  # ----------------------------------------------------------------------------
  message("    ⚖️ Dynamische oppervlaktetoets (5-10 ha) en vennenfilter (<95%) berekenen per cluster...")
  
  if (!all(is.na(suppressWarnings(terra::minmax(cl_biotoop))))) {
    cells_idx <- terra::cells(!is.na(cl_biotoop) & cl_biotoop > 0)
    
    dt_eval <- data.table::data.table(
      ID         = as.integer(cl_biotoop[cells_idx][[1]]),
      Bwk_Val    = as.numeric(r_bwk_raw[cells_idx][[1]]),
      Slik_Val   = as.numeric(terra::ifel(is.na(r_slik_raw[cells_idx][[1]]), 0, r_slik_raw[cells_idx][[1]])),
      Vennen_Val = as.numeric(terra::ifel(is.na(r_vennen_raw[cells_idx][[1]]), 0, r_vennen_raw[cells_idx][[1]]))
    )
    rm(cells_idx)
    gc()
    
    dt_stats <- dt_eval[, .(
      Totaal_ha = sum(Bwk_Val, na.rm = TRUE) * 0.01,
      Slik_ha   = sum(Slik_Val, na.rm = TRUE) * 0.01,
      Vennen_ha = sum(Vennen_Val, na.rm = TRUE) * 0.01
    ), by = ID]
    
    # 1. Dynamische drempel: 10 - 5 * (Slik_ha / Totaal_ha)
    dt_stats[, Slik_Fractie := Slik_ha / Totaal_ha]
    dt_stats[, MinOpp_Req   := 10.0 - (5.0 * pmin(1.0, Slik_Fractie))]
    
    # 2. Vennenfractie
    dt_stats[, Vennen_Pct   := (Vennen_ha / Totaal_ha) * 100]
    
    # Selectiecriterium: Voldoen aan dynamische drempel ÉN vennen < 95%
    ok_cluster_ids <- dt_stats[!is.na(ID) & (Totaal_ha >= MinOpp_Req) & (Vennen_Pct < 95), ID]
    rm(dt_eval, dt_stats)
    gc()
    
    if (length(ok_cluster_ids) > 0) {
      cl_opp    <- terra::ifel(cl_biotoop %in% as.integer(ok_cluster_ids), cl_biotoop, NA)
      final_opp <- terra::mask(r_bwk_raw, cl_opp)
    } else {
      cl_opp    <- master_grid * NA
      final_opp <- master_grid * NA
    }
  } else {
    cl_opp    <- master_grid * NA
    final_opp <- master_grid * NA
  }
  rm(cl_biotoop, r_bwk_raw, r_slik_raw, r_vennen_raw)
  gc()
  
  ha_finaal <- calc_ha_exact(final_opp)
  message(sprintf("       -> Goedgekeurd leefgebied na dynamische filtering: %.2f ha", ha_finaal))
  
  return(list(
    final_opp = final_opp,
    cl_opp    = cl_opp
  ))
}
