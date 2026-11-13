# ==============================================================================
# bouw_biotoop_kust_stern.R
# Archetype voor kust-/vissterns (broedlocatie >= 20ha x viswater >= 1ha binnen 5km)
# Geschikt voor: Dwergstern, Visdief, etc.
# ==============================================================================

bouw_biotoop_kust_stern <- function(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
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
    lagen_packed <- lapply(lagen, terra::wrap)
    saveRDS(lagen_packed, f_cache_rds)
    message("    💾 [CACHE] Parquet-lagen veilig bewaard.")
  }
  
  r_vp_raw   <- if ("voortplanting" %in% names(lagen)) lagen[["voortplanting"]] else (master_grid * NA)
  r_foer_bwk <- if ("foerageer"     %in% names(lagen)) lagen[["foerageer"]]     else (master_grid * NA)
  
  # ----------------------------------------------------------------------------
  # 2. SAMENSTELLEN EN FILTEREN BASISBIOTOPEN
  # ----------------------------------------------------------------------------
  message("    🌊 Waterlopen (Huet-zones) integreren en broedgebied zuiveren van open water...")
  
  # A. Broedgebied zuiveren van direct open water / moeras (not openwater)
  r_vp_gezuiverd <- terra::mask(r_vp_raw, r_foer_bwk > 0, inverse = TRUE)
  rm(r_vp_raw)
  gc()
  
  # B. Viswater samenstellen (BWK open water + Huet-zones 17 & 19)
  r_huet_fractie <- master_grid * NA
  if (file.exists(paths$huetzon)) {
    f_hz_al <- hulp_align_gdal(paths$huetzon, master_grid)
    r_hz    <- terra::rast(f_hz_al)
    r_huet_fractie <- terra::ifel(!is.na(r_hz) & r_hz %in% c(17, 19), 0.01, NA)
    unlink(f_hz_al)
    rm(r_hz)
  }
  
  w_bwk_c  <- terra::ifel(is.na(r_foer_bwk), 0, r_foer_bwk)
  w_huet_c <- terra::ifel(is.na(r_huet_fractie), 0, r_huet_fractie)
  r_water_totaal <- terra::ifel((w_bwk_c + w_huet_c) > 0, terra::clamp(w_bwk_c + w_huet_c, upper = 1.0), NA)
  rm(w_bwk_c, w_huet_c, r_foer_bwk, r_huet_fractie)
  gc()
  
  # ----------------------------------------------------------------------------
  # 3. BASISCLUSTERING: BROEDGEBIED (>= 20 HA) EN WATER (>= 1 HA)
  # ----------------------------------------------------------------------------
  message("    🏝️ Broedgebieden clusteren op >= 20 ha en viswater op >= 1 ha (50m netwerk)...")
  
  # Broedgebied >= 20 ha
  res_vp <- cluster_filter_universeel(
    masker            = r_vp_gezuiverd,
    opp_laag          = r_vp_gezuiverd,
    min_opp_ha        = 20,
    min_totale_opp_ha = NA,
    afstand_m         = 50,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_vp_20ha <- res_vp$raster
  rm(res_vp, r_vp_gezuiverd)
  gc()
  
  # Viswater >= 1 ha
  res_water <- cluster_filter_universeel(
    masker            = r_water_totaal,
    opp_laag          = r_water_totaal,
    min_opp_ha        = 1,
    min_totale_opp_ha = NA,
    afstand_m         = 50,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_water_1ha <- res_water$raster
  rm(res_water, r_water_totaal)
  gc()
  
  ha_vp_20ha <- calc_ha_exact(r_vp_20ha)
  ha_wat_1ha <- calc_ha_exact(r_water_1ha)
  message(sprintf("       -> Broedgebieden (>= 20 ha): %.2f ha", ha_vp_20ha))
  message(sprintf("       -> Viswateren (>= 1 ha): %.2f ha", ha_wat_1ha))
  
  # ----------------------------------------------------------------------------
  # 4. WEDERZIJDSE KOPPELING OP 5 KM (GDAL PROXIMITY)
  # ----------------------------------------------------------------------------
  message("    🎯 5 km wederzijdse actieradiuskoppeling tussen broedgebied en viswater...")
  
  heeft_vp    <- !all(is.na(suppressWarnings(terra::minmax(r_vp_20ha))))
  heeft_water <- !all(is.na(suppressWarnings(terra::minmax(r_water_1ha))))
  
  if (heeft_vp && heeft_water) {
    # 5 km proximity vanaf viswater
    f_w_bin <- tempfile(fileext = ".tif")
    terra::writeRaster(terra::ifel(!is.na(r_water_1ha) & r_water_1ha > 0, 1, NA), f_w_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    f_prox_w5k <- hulp_proximity_gdal(f_w_bin, 5000)
    unlink(f_w_bin)
    
    p_w5k <- terra::rast(f_prox_w5k)
    mask_water5k <- terra::ifel(!is.na(p_w5k) & p_w5k >= 0 & p_w5k <= 5000, 1, NA)
    unlink(f_prox_w5k)
    rm(p_w5k)
    gc()
    
    vp_finaal <- terra::mask(r_vp_20ha, mask_water5k)
    rm(mask_water5k)
    
    # 5 km proximity vanaf broedgebied
    if (!all(is.na(suppressWarnings(terra::minmax(vp_finaal))))) {
      f_vp_bin <- tempfile(fileext = ".tif")
      terra::writeRaster(terra::ifel(!is.na(vp_finaal) & vp_finaal > 0, 1, NA), f_vp_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
      f_prox_vp5k <- hulp_proximity_gdal(f_vp_bin, 5000)
      unlink(f_vp_bin)
      
      p_vp5k <- terra::rast(f_prox_vp5k)
      mask_vp5k <- terra::ifel(!is.na(p_vp5k) & p_vp5k >= 0 & p_vp5k <= 5000, 1, NA)
      unlink(f_prox_vp5k)
      rm(p_vp5k)
      gc()
      
      water_finaal <- terra::mask(r_water_1ha, mask_vp5k)
      rm(mask_vp5k)
    } else {
      water_finaal <- master_grid * NA
    }
  } else {
    vp_finaal    <- master_grid * NA
    water_finaal <- master_grid * NA
  }
  rm(r_vp_20ha, r_water_1ha)
  gc()
  
  # ----------------------------------------------------------------------------
  # 5. SAMENVOEGEN COMPOSIET LEEFGEBIED & FINALE 20 HA CLUSTERING
  # ----------------------------------------------------------------------------
  message("    🧩 Composiet leefgebied samenstellen en filteren op >= 20 ha complex...")
  
  r_leef_bruto <- terra::cover(vp_finaal, water_finaal)
  
  res_finaal <- cluster_filter_universeel(
    masker            = r_leef_bruto,
    opp_laag          = r_leef_bruto,
    min_opp_ha        = 20,
    min_totale_opp_ha = NA,
    afstand_m         = 50,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  final_opp <- res_finaal$raster
  cl_opp    <- res_finaal$clusters
  rm(res_finaal, r_leef_bruto)
  gc()
  
  # Deellagen netjes isoleren binnen goedgekeurde clusters
  vp_uitvoer  <- terra::mask(vp_finaal, cl_opp)
  wat_uitvoer <- terra::mask(water_finaal, cl_opp)
  rm(vp_finaal, water_finaal)
  gc()
  
  return(list(
    final_opp     = final_opp,
    cl_opp        = cl_opp,
    voortplanting = vp_uitvoer,
    foerageren    = wat_uitvoer
  ))
}