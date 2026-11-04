# ==============================================================================
# bouw_biotoop_open_inundatiemoeras.R
# Archetype: Inundatiegraslanden & zeggenmoerassen x drainage x hooggroen 100m x meertrapscluster
# Geschikt voor: Porseleinhoen, Kleinst waterhoen, etc.
# ==============================================================================

bouw_biotoop_open_inundatiemoeras <- function(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
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
  
  r_gras_raw   <- if ("grasland" %in% names(lagen)) lagen[["grasland"]] else (master_grid * NA)
  r_moeras_raw <- if ("moeras"   %in% names(lagen)) lagen[["moeras"]]   else (master_grid * NA)
  
  # Terrils uitsluiten uit Parquet
  t_terril <- ds_parquet %>%
    dplyr::filter(stringr::str_detect(tolower(CODE), "^kg")) %>%
    dplyr::collect() %>%
    data.table::as.data.table()
  
  r_terril_bin <- master_grid * NA
  if (nrow(t_terril) > 0) {
    t_terril[, CODE := tolower(trimws(CODE))]
    t_terril_cel <- unique(t_terril, by = c("cel_id", "CODE"))[, .(Val = 1), by = cel_id]
    r_terril_bin[t_terril_cel$cel_id] <- 1
    rm(t_terril_cel)
  }
  rm(t_terril)
  gc()
  
  # ----------------------------------------------------------------------------
  # 2. SAMENVOEGING & DIRECTE UITSLUITING (HOOGGROEN & TERRILS)
  # ----------------------------------------------------------------------------
  message("    🌾 Grasland en moeras combineren, directe hooggroen- en terril-uitsluiting...")
  g_c <- terra::ifel(is.na(r_gras_raw),   0, r_gras_raw)
  m_c <- terra::ifel(is.na(r_moeras_raw), 0, r_moeras_raw)
  r_bio_bruto <- terra::ifel((g_c + m_c) > 0, terra::clamp(g_c + m_c, upper = 1.0), NA)
  rm(g_c, m_c)
  gc()
  
  f_groen_al <- hulp_align_gdal(paths$groenkaart, master_grid)
  r_groen    <- terra::rast(f_groen_al)
  mask_direct_fout <- terra::ifel((!is.na(r_groen) & (r_groen == 1 | r_groen > 0)) | (!is.na(r_terril_bin) & r_terril_bin == 1), 1, NA)
  unlink(f_groen_al)
  rm(r_groen, r_terril_bin)
  
  r_bio_geschoond <- terra::mask(r_bio_bruto, mask_direct_fout, inverse = TRUE)
  rm(r_bio_bruto, mask_direct_fout)
  gc()
  
  # ----------------------------------------------------------------------------
  # 3. NATTE GRONDEN FILTER (BODEMDRAINAGEKAART)
  # ----------------------------------------------------------------------------
  message("    💧 Beperken tot natte drainageklassen (d, e, f, h, i, g, e-f, h-i, e-i)...")
  if (file.exists(paths$drainage)) {
    f_drain_al   <- hulp_align_gdal(paths$drainage, master_grid)
    r_drain      <- terra::rast(f_drain_al)
    leg_drain    <- terra::cats(r_drain)[[1]]
    drain_let    <- c("d", "e", "f", "h", "i", "g", "e-f", "h-i", "e-i")
    ok_drain_ids <- leg_drain$value[as.character(leg_drain$Label) %in% drain_let]
    
    mask_drain <- terra::ifel(!is.na(r_drain) & r_drain %in% ok_drain_ids, 1, NA)
    unlink(f_drain_al)
    rm(r_drain, leg_drain)
    
    r_bio_nat <- terra::mask(r_bio_geschoond, mask_drain)
    rm(mask_drain)
  } else {
    r_bio_nat <- r_bio_geschoond
  }
  rm(r_bio_geschoond)
  gc()
  
  ha_nat <- calc_ha_exact(r_bio_nat)
  message(sprintf("       -> Biotoop na drainagefilter: %.2f ha", ha_nat))
  
  # ----------------------------------------------------------------------------
  # 4. EERSTE CLUSTERING & FILTERING (>= 3 HA OP 50M NETWERK)
  # ----------------------------------------------------------------------------
  message("    🏝️ Eerste clustertoets uitvoeren (>= 3 ha, 50m netwerk)...")
  res_3ha <- cluster_filter_universeel(
    masker            = r_bio_nat,
    opp_laag          = r_bio_nat,
    min_opp_ha        = 3,
    min_totale_opp_ha = NA,
    afstand_m         = 50,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_bio_3ha <- res_3ha$raster
  rm(res_3ha, r_bio_nat)
  gc()
  
  # ----------------------------------------------------------------------------
  # 5. OPEN ZICHT FILTER (100M EXCLUSIE ROND HOOGGROEN >= 20 ARE & HERKEURING 3 HA)
  # ----------------------------------------------------------------------------
  message("    🌲 Open zicht filter: 100m buffer rond hooggroenclusters (>= 20 are)...")
  
  f_cache_groen100 <- file.path(paths$dir_cache, "masker_groen_100m_cluster.tif")
  if (file.exists(f_cache_groen100)) {
    message("       -> [CACHE] Hooggroen-100m storingsbuffer ingeladen van schijf...")
    masker_groen_100m <- terra::rast(f_cache_groen100)
  } else {
    f_groen_al   <- hulp_align_gdal(paths$groenkaart, master_grid)
    r_groen_raw  <- terra::rast(f_groen_al)
    r_groen_bin  <- terra::ifel(!is.na(r_groen_raw) & (r_groen_raw == 1 | r_groen_raw == 2 | r_groen_raw > 0), 1, NA)
    unlink(f_groen_al)
    rm(r_groen_raw)
    
    w_5x5 <- matrix(1, nrow = 5, ncol = 5)
    f_focal_g <- tempfile(fileext = ".tif")
    r_groen_som <- terra::focal(r_groen_bin, w = w_5x5, fun = "sum", na.rm = TRUE, filename = f_focal_g, overwrite = TRUE)
    r_groen_20are <- terra::ifel(!is.na(r_groen_bin) & r_groen_som >= 12, 1, NA)
    unlink(f_focal_g)
    rm(r_groen_bin, r_groen_som)
    
    if (!all(is.na(suppressWarnings(terra::minmax(r_groen_20are))))) {
      f_g20_bin <- tempfile(fileext = ".tif")
      terra::writeRaster(r_groen_20are, f_g20_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 0)
      f_prox_g  <- hulp_proximity_gdal(f_g20_bin, 100)
      unlink(f_g20_bin)
      
      p_g <- terra::rast(f_prox_g)
      masker_groen_100m <- terra::ifel(!is.na(p_g) & p_g >= 0 & p_g <= 100, 1, NA)
      unlink(f_prox_g)
      rm(p_g)
    } else {
      masker_groen_100m <- master_grid * NA
    }
    terra::writeRaster(masker_groen_100m, f_cache_groen100, overwrite = TRUE, datatype = "INT1U", gdal = "COMPRESS=LZW")
    rm(r_groen_20are)
    gc()
  }
  
  # Aftrek van de 100m buffer
  r_bio_open <- terra::mask(r_bio_3ha, masker_groen_100m, inverse = TRUE)
  rm(r_bio_3ha, masker_groen_100m)
  gc()
  
  # Directe herkeuring op de 3 ha drempel na aftrek van de bosrandstoring
  res_3ha_post <- cluster_filter_universeel(
    masker            = r_bio_open,
    opp_laag          = r_bio_open,
    min_opp_ha        = 3,
    min_totale_opp_ha = NA,
    afstand_m         = 50,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_bio_open_3ha <- res_3ha_post$raster
  rm(res_3ha_post, r_bio_open)
  gc()
  
  # ----------------------------------------------------------------------------
  # 6. FINALE CLUSTERING & EINDDREMPEL (>= 10 HA OP 50M NETWERK)
  # ----------------------------------------------------------------------------
  message("    🧩 Finale clustering op functioneel complex (>= 10 ha, 50m netwerk)...")
  res_10ha <- cluster_filter_universeel(
    masker            = r_bio_open_3ha,
    opp_laag          = r_bio_open_3ha,
    min_opp_ha        = 10,
    min_totale_opp_ha = NA,
    afstand_m         = 50,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  final_opp <- res_10ha$raster
  cl_opp    <- res_10ha$clusters
  rm(res_10ha, r_bio_open_3ha)
  gc()
  
  # Deellagen isoleren binnen de definitieve clusters
  gras_uitvoer  <- terra::mask(r_gras_raw, cl_opp)
  moeras_uitvoer<- terra::mask(r_moeras_raw, cl_opp)
  rm(r_gras_raw, r_moeras_raw)
  gc()
  
  return(list(
    final_opp = final_opp,
    cl_opp    = cl_opp,
    grasland  = gras_uitvoer,
    moeras    = moeras_uitvoer
  ))
}
