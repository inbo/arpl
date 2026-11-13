# ==============================================================================
# bouw_biotoop_open_grasland_bosrandvrij.R
# Archetype: Open graslanden & venen met strikte bosrandexclusie (>= 100m buffer)
# Exacte vertaling van werkend script Turnhouts Vennegebied naar heel Vlaanderen
# Geschikt voor: Paapje, Graspieper, etc.
# ==============================================================================

bouw_biotoop_open_grasland_bosrandvrij <- function(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
  bestandsnaam <- tolower(gsub(" ", "", soort_naam))
  soort_cache  <- file.path(paths$dir_cache, bestandsnaam)
  dir.create(soort_cache, showWarnings = FALSE, recursive = TRUE)
  
  f_cache_rds <- file.path(soort_cache, "01_raw_parquet_layers.rds")
  
  # Parameters ophalen
  rij_params     <- params_soort %>% dplyr::slice(1)
  min_opp_ha     <- as.numeric(ifelse(is.na(rij_params$MinOpp_ha), 5, rij_params$MinOpp_ha))
  zoekafstand_m  <- as.numeric(ifelse(is.na(rij_params$AfstandBiotopen_m), 50, rij_params$AfstandBiotopen_m))
  bos_drempel_m2 <- as.numeric(ifelse(is.na(rij_params$Groenkaart_Drempel_m2), 2000, rij_params$Groenkaart_Drempel_m2))
  bos_min_cellen <- max(1, round(bos_drempel_m2 / 100)) # 2000 m2 = 20 cellen van 10x10m
  
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
  
  r_doel_raw <- if ("bwk"    %in% names(lagen)) lagen[["bwk"]]    else (master_grid * NA)
  r_bos_bwk  <- if ("bwkbos" %in% names(lagen)) lagen[["bwkbos"]] else (master_grid * NA)
  
  ha_doel_ruw <- calc_ha_exact(r_doel_raw)
  message(sprintf("       -> Doelbiotoop (grasland/veen) vóór filtering: %.2f ha", ha_doel_ruw))
  
  # ----------------------------------------------------------------------------
  # 2. BOS OPBOUWEN EN FILTEREN OP >= 20 ARE (20 CELLEN)
  # ----------------------------------------------------------------------------
  message(sprintf("    🌳 Boslagen samenvoegen en filteren op boomopstanden >= %d are...", as.integer(bos_min_cellen)))
  
  # A. Groenkaart hooggroen (waarde 1)
  f_groen_al     <- hulp_align_gdal(paths$groenkaart, master_grid)
  r_groen        <- terra::rast(f_groen_al)
  groenkaart_bos <- terra::ifel(!is.na(r_groen) & r_groen == 1, 1, NA)
  unlink(f_groen_al)
  rm(r_groen)
  
  # B. BWK-bos: doelbiotoop maskeren (zoals in het originele script)
  if (!all(is.na(suppressWarnings(terra::minmax(r_bos_bwk))))) {
    bwkbos_gecorrigeerd <- terra::ifel(!is.na(r_doel_raw), NA, r_bos_bwk)
    bwkbos_bin          <- terra::ifel(!is.na(bwkbos_gecorrigeerd) & bwkbos_gecorrigeerd > 0, 1, NA)
    rm(bwkbos_gecorrigeerd)
  } else {
    bwkbos_bin <- master_grid * NA
  }
  rm(r_bos_bwk)
  
  # C. Samenvoegen via cover (beide lagen bevatten enkel 1 of NA)
  paapje_bos_ruw <- terra::cover(bwkbos_bin, groenkaart_bos)
  rm(bwkbos_bin, groenkaart_bos)
  gc()
  
  # D. Clumpen via GRASS om patches van >= 20 cellen te selecteren
  f_bos_bin   <- tempfile(fileext = ".tif")
  f_clump_bos <- tempfile(fileext = ".tif")
  terra::writeRaster(paapje_bos_ruw, f_bos_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
  
  e <- terra::ext(master_grid)
  qgisprocess::qgis_run_algorithm(
    "grass:r.clump", 
    input = f_bos_bin, 
    output = f_clump_bos, 
    title = "cl_bos_paapje",
    "-d" = FALSE, 
    GRASS_REGION_PARAMETER = sprintf("%f,%f,%f,%f", e$xmin, e$xmax, e$ymin, e$ymax),
    GRASS_REGION_CELLSIZE_PARAMETER = 10, 
    .quiet = TRUE
  )
  unlink(f_bos_bin)
  
  cl_bos <- terra::rast(f_clump_bos)
  terra::set.cats(cl_bos, layer = 1, value = NULL)
  if (terra::has.colors(cl_bos)) try(terra::coltab(cl_bos) <- NULL, silent = TRUE)
  
  # Alleen patches binnen het bos behouden
  cl_bos_only <- terra::mask(cl_bos, paapje_bos_ruw)
  unlink(f_clump_bos)
  rm(cl_bos, paapje_bos_ruw)
  gc()
  
  dt_bos <- data.table::as.data.table(terra::freq(cl_bos_only))
  setnames(dt_bos, c("layer", "ID", "count"))
  
  # Filter strictly op geldige ID's met >= 20 cellen
  grote_bos_ids <- dt_bos[!is.na(ID) & ID > 0 & count >= bos_min_cellen, ID]
  rm(dt_bos)
  
  if (length(grote_bos_ids) > 0) {
    paapje_bos <- terra::ifel(cl_bos_only %in% as.integer(grote_bos_ids), 1, NA)
  } else {
    paapje_bos <- master_grid * NA
  }
  rm(cl_bos_only)
  gc()
  
  ha_bos_storend <- calc_ha_exact(paapje_bos)
  message(sprintf("       -> Storende boomopstanden (>= %d cellen): %.2f ha", as.integer(bos_min_cellen), ha_bos_storend))
  
  # ----------------------------------------------------------------------------
  # 3. 100 METER BOSRAND-BUFFER TOEPASSEN (GDAL PROXIMITY)
  # ----------------------------------------------------------------------------
  message("    🚫 Bosrand-exclusiezone (100m buffer) toepassen op doelbiotoop...")
  
  if (ha_bos_storend > 0) {
    f_b_bin <- tempfile(fileext = ".tif")
    terra::writeRaster(paapje_bos, f_b_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    rm(paapje_bos)
    gc()
    
    # GDAL Proximity: berekent afstand tot cellen met waarde 1
    f_prox_bos <- hulp_proximity_gdal(f_b_bin, 100)
    unlink(f_b_bin)
    
    p_bos <- terra::rast(f_prox_bos)
    
    # Storingszone: cellen binnen 0 t.e.m. 100 meter van bos
    mask_storing_100m <- terra::ifel(!is.na(p_bos) & p_bos >= 0 & p_bos <= 100, 1, NA)
    unlink(f_prox_bos)
    rm(p_bos)
    gc()
    
    # Exclusie: snijd doelbiotoop binnen de storingszone weg
    r_doel_open <- terra::mask(r_doel_raw, mask_storing_100m, inverse = TRUE)
    rm(mask_storing_100m)
    gc()
  } else {
    r_doel_open <- r_doel_raw
  }
  rm(r_doel_raw)
  gc()
  
  ha_doel_open <- calc_ha_exact(r_doel_open)
  message(sprintf("       -> Doelbiotoop buiten 100m bosrandzone: %.2f ha", ha_doel_open))
  
  # ----------------------------------------------------------------------------
  # 4. CLUSTERING EN OPPERVLAKTEFILTER (>= 5 HA OP 50M ZOEFAFSTAND)
  # ----------------------------------------------------------------------------
  message(sprintf("    🏝️ Clusteren op minimumoppervlakte (>= %.1f ha, zoekafstand %dm)...", min_opp_ha, as.integer(zoekafstand_m)))
  
  res_cluster <- cluster_filter_universeel(
    masker            = r_doel_open,
    opp_laag          = r_doel_open,
    min_opp_ha        = min_opp_ha,
    min_totale_opp_ha = NA,
    afstand_m         = zoekafstand_m,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  
  final_opp <- res_cluster$raster
  cl_opp    <- res_cluster$clusters
  rm(res_cluster, r_doel_open)
  gc()
  
  ha_finaal <- calc_ha_exact(final_opp)
  message(sprintf("       -> Definitief leefgebied Paapje: %.2f ha", ha_finaal))
  
  return(list(
    final_opp = final_opp,
    cl_opp    = cl_opp
  ))
}
