# ==============================================================================
# bouw_biotoop_saproxyl_hellingbos.R
# Archetype voor saproxyle bosrandkevers (loofbos <= 40% naald x oud bos 1km x helling 4° 1km x rand 15m/40m x 5ha)
# Geschikt voor: Vliegend hert (Lucanus cervus), etc.
# ==============================================================================

bouw_biotoop_saproxyl_hellingbos <- function(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
  bestandsnaam <- tolower(gsub(" ", "", soort_naam))
  soort_cache  <- file.path(paths$dir_cache, bestandsnaam)
  dir.create(soort_cache, showWarnings = FALSE, recursive = TRUE)
  
  f_cache_rds <- file.path(soort_cache, "01_raw_parquet_layers.rds")
  
  # Parameters ophalen uit Excel
  rij_params    <- params_soort %>% dplyr::slice(1)
  min_opp_ha    <- as.numeric(ifelse(is.na(rij_params$MinOpp_ha), 5.0, rij_params$MinOpp_ha))
  bos_opp_ha    <- as.numeric(ifelse(is.na(rij_params$MinTotaleOpp_ha), 10.0, rij_params$MinTotaleOpp_ha))
  zoekafstand_m <- as.numeric(ifelse(is.na(rij_params$AfstandBiotopen_m), 50, rij_params$AfstandBiotopen_m))
  
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
  
  r_naald_raw <- if ("naaldbos" %in% names(lagen)) lagen[["naaldbos"]] else (master_grid * NA)
  
  # ----------------------------------------------------------------------------
  # 2. BASISLAGEN INLADEN (GROENKAART, HISTORISCH BOS & HELLINGEN)
  # ----------------------------------------------------------------------------
  message("    🌲 Groenkaart, oud bos en hellingen inlezen en synchroniseren...")
  
  # A. Groenkaart: Hooggroen (1) vs Open groen (2, 3, 4)
  f_groen_al <- hulp_align_gdal(paths$groenkaart, master_grid)
  r_groen    <- terra::rast(f_groen_al)
  
  r_groen_hoog <- terra::ifel(!is.na(r_groen) & r_groen == 1, 1, NA)
  r_groen_open <- terra::ifel(!is.na(r_groen) & (r_groen %in% c(2, 3, 4)), 1, NA)
  unlink(f_groen_al)
  rm(r_groen)
  
  # B. Historisch oud bos (klassen 1 en 2) binnen 1.000 m
  f_cache_oud_bos <- file.path(paths$dir_cache, "masker_oud_bos_1000m.tif")
  
  if (file.exists(f_cache_oud_bos)) {
    mask_oud_bos_1km <- terra::rast(f_cache_oud_bos)
  } else if (!is.null(paths$bosleeftijd) && file.exists(paths$bosleeftijd)) {
    f_oud_al  <- hulp_align_gdal(paths$bosleeftijd, master_grid)
    r_oud_raw <- terra::rast(f_oud_al)
    r_oud_bin <- terra::ifel(!is.na(r_oud_raw) & (r_oud_raw %in% c(1, 2)), 1, NA)
    unlink(f_oud_al)
    rm(r_oud_raw)
    
    if (!all(is.na(suppressWarnings(terra::minmax(r_oud_bin))))) {
      f_ob_bin <- tempfile(fileext = ".tif")
      terra::writeRaster(r_oud_bin, f_ob_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
      f_prox_ob <- hulp_proximity_gdal(f_ob_bin, 1000)
      unlink(f_ob_bin)
      
      p_ob <- terra::rast(f_prox_ob)
      mask_oud_bos_1km <- terra::ifel(!is.na(p_ob) & p_ob >= 0 & p_ob <= 1000, 1, NA)
      unlink(f_prox_ob)
      rm(p_ob)
    } else {
      mask_oud_bos_1km <- master_grid * 1
    }
    rm(r_oud_bin)
    terra::writeRaster(mask_oud_bos_1km, f_cache_oud_bos, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
  } else {
    mask_oud_bos_1km <- master_grid * 1
  }
  
  # C. Hellingen >= 4 graden binnen 1.000 m
  f_cache_helling <- file.path(paths$dir_cache, "masker_helling_4gr_1000m.tif")
  
  if (file.exists(f_cache_helling)) {
    mask_helling_1km <- terra::rast(f_cache_helling)
  } else if (!is.null(paths$helling_4gr) && file.exists(paths$helling_4gr)) {
    f_hel_al  <- hulp_align_gdal(paths$helling_4gr, master_grid)
    r_hel_raw <- terra::rast(f_hel_al)
    r_hel_bin <- terra::ifel(!is.na(r_hel_raw) & r_hel_raw == 1, 1, NA)
    unlink(f_hel_al)
    rm(r_hel_raw)
    
    if (!all(is.na(suppressWarnings(terra::minmax(r_hel_bin))))) {
      f_h_bin <- tempfile(fileext = ".tif")
      terra::writeRaster(r_hel_bin, f_h_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
      f_prox_h <- hulp_proximity_gdal(f_h_bin, 1000)
      unlink(f_h_bin)
      
      p_h <- terra::rast(f_prox_h)
      mask_helling_1km <- terra::ifel(!is.na(p_h) & p_h >= 0 & p_h <= 1000, 1, NA)
      unlink(f_prox_h)
      rm(p_h)
    } else {
      mask_helling_1km <- master_grid * 1
    }
    rm(r_hel_bin)
    terra::writeRaster(mask_helling_1km, f_cache_helling, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
  } else {
    mask_helling_1km <- master_grid * 1
  }
  
  # ----------------------------------------------------------------------------
  # 3. LOOFBOSSELECTIE MET MAX. 40% NAALDHOUT EN OMGEVINGSFILTERS
  # ----------------------------------------------------------------------------
  message("    🌳 Loofbos afbakenen (<= 40% naaldhout) binnen 1 km van oud bos en hellingen...")
  
  r_naald_clean <- terra::ifel(is.na(r_naald_raw), 0, r_naald_raw)
  r_loof_frac   <- terra::clamp(1.0 - r_naald_clean, lower = 0.0, upper = 1.0)
  
  # Exclusie van opstanden met > 40% naaldhout
  r_loof_dominant <- terra::ifel(!is.na(r_groen_hoog) & r_naald_clean <= 0.40, r_loof_frac, NA)
  rm(r_naald_clean, r_groen_hoog, r_naald_raw)
  gc()
  
  # Filteren op 1 km van oud bos en 1 km van hellingen
  r_bos_gefilterd <- terra::mask(r_loof_dominant, mask_oud_bos_1km)
  r_bos_gefilterd <- terra::mask(r_bos_gefilterd, mask_helling_1km)
  rm(r_loof_dominant, mask_oud_bos_1km, mask_helling_1km)
  gc()
  
  # Clusteren van het geschikte loofbosareaal (>= 10 ha op 50 m)
  message(sprintf("    🏝️ Loofboscomplexen clusteren (>= %.1f ha, %dm netwerk)...", bos_opp_ha, as.integer(zoekafstand_m)))
  res_bos10 <- cluster_filter_universeel(
    masker            = r_bos_gefilterd,
    opp_laag          = r_bos_gefilterd,
    min_opp_ha        = bos_opp_ha,
    min_totale_opp_ha = NA,
    afstand_m         = zoekafstand_m,
    dispersie_m       = 0,
    werkelijk         = FALSE
  )
  vh_bossen_10ha <- res_bos10$raster
  rm(res_bos10, r_bos_gefilterd)
  gc()
  
  ha_bos_10ha <- calc_ha_exact(vh_bossen_10ha)
  message(sprintf("       -> Gekwalificeerde loofboscomplexen (>= 10 ha): %.2f ha", ha_bos_10ha))
  
  # ----------------------------------------------------------------------------
  # 4. AANPALEND OPEN GROEN & SCHERPE BOSRANDEN (<= 15 M)
  # ----------------------------------------------------------------------------
  message("    🌾 Aanpalend open groen clusteren (>= 5 ha op 50 m)...")
  res_open5 <- cluster_filter_universeel(
    masker            = r_groen_open,
    opp_laag          = r_groen_open,
    min_opp_ha        = 5.0,
    min_totale_opp_ha = NA,
    afstand_m         = 50,
    dispersie_m       = 0,
    werkelijk         = FALSE
  )
  vh_open_5ha <- res_open5$raster
  rm(res_open5, r_groen_open)
  gc()
  
  message("    📏 Scherpe bosranden detecteren (loofbos binnen 15 m van open groen)...")
  
  heeft_bos  <- !all(is.na(suppressWarnings(terra::minmax(vh_bossen_10ha))))
  heeft_open <- !all(is.na(suppressWarnings(terra::minmax(vh_open_5ha))))
  
  if (heeft_bos && heeft_open) {
    # 15m buffer vanaf open groen via GDAL proximity
    f_op_bin  <- tempfile(fileext = ".tif")
    terra::writeRaster(terra::ifel(!is.na(vh_open_5ha) & vh_open_5ha > 0, 1, NA), f_op_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    f_prox_op <- hulp_proximity_gdal(f_op_bin, 15)
    unlink(f_op_bin)
    
    p_op <- terra::rast(f_prox_op)
    mask_op15 <- terra::ifel(!is.na(p_op) & p_op >= 0 & p_op <= 15, 1, NA)
    unlink(f_prox_op)
    rm(p_op)
    
    r_bosrand_ruw <- terra::mask(vh_bossen_10ha, mask_op15)
    rm(mask_op15)
    
    # Filteren op aaneengesloten randlengte (>= 0.25 ha op 10 m koppeling)
    res_rand025 <- cluster_filter_universeel(
      masker            = r_bosrand_ruw,
      opp_laag          = r_bosrand_ruw,
      min_opp_ha        = 0.25,
      min_totale_opp_ha = NA,
      afstand_m         = 10,
      dispersie_m       = 0,
      werkelijk         = FALSE
    )
    vh_bosrand_valide <- res_rand025$raster
    rm(res_rand025, r_bosrand_ruw)
    gc()
  } else {
    vh_bosrand_valide <- master_grid * NA
  }
  
  ha_rand_val <- calc_ha_exact(vh_bosrand_valide)
  message(sprintf("       -> Kwalitatieve bosranden (>= 0.25 ha continue lengte): %.2f ha", ha_rand_val))
  
  # ----------------------------------------------------------------------------
  # 5. FUNCTIONELE BOSDIEPTE (40 M DIEPTEBUFFER BOSINWAARTS)
  # ----------------------------------------------------------------------------
  message("    🎯 40 m functionele bosdiepte berekenen vanaf gevalideerde bosranden...")
  
  if (!all(is.na(suppressWarnings(terra::minmax(vh_bosrand_valide))))) {
    f_rd_bin  <- tempfile(fileext = ".tif")
    terra::writeRaster(terra::ifel(!is.na(vh_bosrand_valide) & vh_bosrand_valide > 0, 1, NA), f_rd_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    f_prox_rd <- hulp_proximity_gdal(f_rd_bin, 40)
    unlink(f_rd_bin)
    
    p_rd <- terra::rast(f_prox_rd)
    mask_rd40 <- terra::ifel(!is.na(p_rd) & p_rd >= 0 & p_rd <= 40, 1, NA)
    unlink(f_prox_rd)
    rm(p_rd)
    
    # 40m diepte binnen het gekwalificeerde loofbos
    vh_bos_40m <- terra::mask(vh_bossen_10ha, mask_rd40)
    rm(mask_rd40)
    
    # Bruto leefgebied = rand + 40m diepte, gewogen met de werkelijke looffractie
    vh_leefgebied_bruto <- terra::cover(vh_bosrand_valide, vh_bos_40m)
    vh_leefgebied_bruto <- terra::mask(r_loof_frac, vh_leefgebied_bruto)
    rm(vh_bos_40m)
  } else {
    vh_leefgebied_bruto <- master_grid * NA
  }
  rm(r_loof_frac)
  gc()
  
  # ----------------------------------------------------------------------------
  # 6. FINALE CLUSTERING OP WERKELIJK HABITAT (MIN. 5 HA OP 50 M)
  # ----------------------------------------------------------------------------
  message(sprintf("    🧩 Finale clustering op leefgebiedscomplex (>= %.1f ha werkelijk, %dm netwerk)...", min_opp_ha, as.integer(zoekafstand_m)))
  
  res_finaal <- cluster_filter_universeel(
    masker            = vh_leefgebied_bruto,
    opp_laag          = vh_leefgebied_bruto,
    min_opp_ha        = min_opp_ha,
    min_totale_opp_ha = NA,
    afstand_m         = zoekafstand_m,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  final_opp <- res_finaal$raster
  cl_opp    <- res_finaal$clusters
  rm(res_finaal, vh_leefgebied_bruto)
  gc()
  
  ha_finaal <- calc_ha_exact(final_opp)
  message(sprintf("       -> Definitief leefgebied Vliegend hert: %.2f ha", ha_finaal))
  
  # Deellagen isoleren binnen goedgekeurde netwerken
  if (!all(is.na(suppressWarnings(terra::minmax(cl_opp))))) {
    bos_uitvoer  <- terra::mask(vh_bossen_10ha, cl_opp)
    rand_uitvoer <- terra::mask(vh_bosrand_valide, cl_opp)
    open_uitvoer <- terra::mask(vh_open_5ha, cl_opp)
  } else {
    bos_uitvoer  <- master_grid * NA
    rand_uitvoer <- master_grid * NA
    open_uitvoer <- master_grid * NA
  }
  rm(vh_bossen_10ha, vh_bosrand_valide, vh_open_5ha)
  gc()
  
  return(list(
    final_opp        = final_opp,
    cl_opp           = cl_opp,
    loofbos_complex  = bos_uitvoer,
    bosrand_ecotoon  = rand_uitvoer,
    open_groen       = open_uitvoer
  ))
}
