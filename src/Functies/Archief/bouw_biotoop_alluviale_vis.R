# ==============================================================================
# bouw_biotoop_alluviale_vis.R
# Archetype voor alluviale laaglandvissen (BWK + watervlakken - ap% + Huet P1/P2 x ovstrg 1-2)
# Geschikt voor: Grote modderkruiper (Misgurnus fossilis), etc.
# ==============================================================================

bouw_biotoop_alluviale_vis <- function(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
  bestandsnaam <- tolower(gsub(" ", "", soort_naam))
  soort_cache  <- file.path(paths$dir_cache, bestandsnaam)
  dir.create(soort_cache, showWarnings = FALSE, recursive = TRUE)
  
  f_cache_rds <- file.path(soort_cache, "01_raw_parquet_layers.rds")
  
  # Parameters ophalen uit Excel
  rij_params    <- params_soort %>% dplyr::slice(1)
  min_opp_ha    <- as.numeric(ifelse(is.na(rij_params$MinOpp_ha), 0.04, rij_params$MinOpp_ha))
  zoekafstand_m <- as.numeric(ifelse(is.na(rij_params$AfstandBiotopen_m), 10, rij_params$AfstandBiotopen_m))
  
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
  
  r_bwk_raw     <- if ("bwk"     %in% names(lagen)) lagen[["bwk"]]     else (master_grid * NA)
  r_not_bwk_raw <- if ("not_bwk" %in% names(lagen)) lagen[["not_bwk"]] else (master_grid * NA)
  
  # ----------------------------------------------------------------------------
  # 2. EXTERNE WATERVLAKKEN INLADEN EN SCHONEN VAN AP% (VISVIJVERS)
  # ----------------------------------------------------------------------------
  message("    🌊 Externe watervlakken schonen van intensieve visvijvers (ap%)...")
  
  if (!is.null(paths$watervlakken_raster) && file.exists(paths$watervlakken_raster)) {
    f_wv_al <- hulp_align_gdal(paths$watervlakken_raster, master_grid)
    r_wv    <- terra::rast(f_wv_al)
    r_wv_bin <- terra::ifel(!is.na(r_wv) & r_wv > 0, 1.0, NA)
    unlink(f_wv_al)
    rm(r_wv)
    
    if (!all(is.na(suppressWarnings(terra::minmax(r_not_bwk_raw))))) {
      r_not_zero  <- terra::ifel(is.na(r_not_bwk_raw), 0, r_not_bwk_raw)
      r_rest_opp  <- terra::clamp(1.0 - r_not_zero, lower = 0, upper = 1.0)
      r_plassen_schoon <- terra::ifel(!is.na(r_wv_bin) & r_rest_opp > 0, r_rest_opp, NA)
      rm(r_not_zero, r_rest_opp)
    } else {
      r_plassen_schoon <- r_wv_bin
    }
    rm(r_wv_bin)
  } else {
    r_plassen_schoon <- master_grid * NA
  }
  rm(r_not_bwk_raw)
  gc()
  
  # ----------------------------------------------------------------------------
  # 3. HYDROGRAFISCHE HUET-ZONES (P1 & P2: BRASEM- EN KARPERZONES)
  # ----------------------------------------------------------------------------
  message("    🐟 Huet-waterlopen (P1 & P2) inladen en filteren...")
  
  f_cache_huet_modder <- file.path(paths$dir_cache, "masker_huetzon_modderkruiper.tif")
  
  if (file.exists(f_cache_huet_modder)) {
    message("       -> [CACHE] Huetzone P1/P2 laag ingeladen van schijf...")
    r_huet_p1p2 <- terra::rast(f_cache_huet_modder)
  } else if (!is.null(paths$huetzon) && file.exists(paths$huetzon)) {
    f_huet_al <- hulp_align_gdal(paths$huetzon, master_grid)
    r_huet    <- terra::rast(f_huet_al)
    leg_huet  <- terra::cats(r_huet)[[1]]
    
    target_p <- c("p1", "p2")
    ok_p_ids <- integer(0)
    
    if (!is.null(leg_huet)) {
      kol_val   <- colnames(leg_huet)[tolower(colnames(leg_huet)) %in% c("value", "id")][1]
      kol_label <- colnames(leg_huet)[tolower(colnames(leg_huet)) %in% c("label", "omschrijving", "huet")][1]
      if (is.na(kol_val)) kol_val <- 1
      if (!is.na(kol_label)) {
        huet_labels_clean <- tolower(trimws(as.character(leg_huet[[kol_label]])))
        ok_p_ids <- leg_huet[[kol_val]][huet_labels_clean %in% target_p]
      }
    }
    
    if (length(ok_p_ids) > 0) {
      r_huet_p1p2 <- terra::ifel(!is.na(r_huet) & (r_huet %in% ok_p_ids), 1.0, NA)
    } else {
      r_huet_p1p2 <- master_grid * NA
    }
    unlink(f_huet_al)
    rm(r_huet)
    terra::writeRaster(r_huet_p1p2, f_cache_huet_modder, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
  } else {
    r_huet_p1p2 <- master_grid * NA
  }
  
  # ----------------------------------------------------------------------------
  # 4. SAMENVOEGING WATERBASIS
  # ----------------------------------------------------------------------------
  message("    💧 Waterbasis samenstellen (BWK + Plassen + Huet P1/P2)...")
  
  r_water_tot <- terra::cover(r_bwk_raw, r_plassen_schoon)
  r_water_tot <- terra::cover(r_water_tot, r_huet_p1p2)
  r_water_tot <- terra::clamp(r_water_tot, upper = 1.0)
  
  ha_water_tot <- calc_ha_exact(r_water_tot)
  message(sprintf("       -> Bruto waterbasis vóór overstromingsfilter: %.2f ha", ha_water_tot))
  
  # ----------------------------------------------------------------------------
  # 5. ALLUVIALE OVERSTROMINGSFILTER (OVSTRG KLASSE 1 EN 2)
  # ----------------------------------------------------------------------------
  message("    🌊 Beperken tot overstromingsgevoelige gronden (ovstrg klasse 1 en 2)...")
  
  f_cache_ovstr <- file.path(paths$dir_cache, "masker_overstroming_1_2.tif")
  
  if (file.exists(f_cache_ovstr)) {
    mask_ovstr <- terra::rast(f_cache_ovstr)
  } else if (!is.null(paths$overstroming) && file.exists(paths$overstroming)) {
    f_ovstr_al <- hulp_align_gdal(paths$overstroming, master_grid)
    r_ovstr    <- terra::rast(f_ovstr_al)
    mask_ovstr <- terra::ifel(!is.na(r_ovstr) & (r_ovstr %in% c(1, 2)), 1, NA)
    unlink(f_ovstr_al)
    rm(r_ovstr)
    terra::writeRaster(mask_ovstr, f_cache_ovstr, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
  } else {
    mask_ovstr <- master_grid * NA
  }
  
  if (!all(is.na(suppressWarnings(terra::minmax(mask_ovstr))))) {
    r_water_overstroombaar <- terra::mask(r_water_tot, mask_ovstr)
  } else {
    message("       ⚠️ Waarschuwing: Geen overstromingsgevoelige gronden gevonden!")
    r_water_overstroombaar <- master_grid * NA
  }
  rm(r_water_tot, mask_ovstr)
  gc()
  
  ha_overstroombaar <- calc_ha_exact(r_water_overstroombaar)
  message(sprintf("       -> Waterbiotoop binnen alluviaal overstromingsgebied: %.2f ha", ha_overstroombaar))
  
  # ----------------------------------------------------------------------------
  # 6. RUIMTELIJKE CLUSTERING & OPPERVLAKTEDREMPEL (>= 0.04 HA OP 10M NETWERK)
  # ----------------------------------------------------------------------------
  message(sprintf("    🏝️ Duikers overbruggen (%dm netwerk) en filteren op >= %.2f ha (400 m²)...", as.integer(zoekafstand_m), min_opp_ha))
  
  res_cluster <- cluster_filter_universeel(
    masker            = r_water_overstroombaar,
    opp_laag          = r_water_overstroombaar,
    min_opp_ha        = min_opp_ha,
    min_totale_opp_ha = NA,
    afstand_m         = zoekafstand_m,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  
  final_opp <- res_cluster$raster
  cl_opp    <- res_cluster$clusters
  rm(res_cluster, r_water_overstroombaar)
  gc()
  
  ha_finaal <- calc_ha_exact(final_opp)
  message(sprintf("       -> Definitief leefgebied Grote modderkruiper: %.2f ha", ha_finaal))
  
  # Deellagen isoleren binnen goedgekeurde netwerken
  if (!all(is.na(suppressWarnings(terra::minmax(cl_opp))))) {
    bwk_uitvoer     <- terra::mask(r_bwk_raw, cl_opp)
    plassen_uitvoer <- terra::mask(r_plassen_schoon, cl_opp)
    huet_uitvoer    <- terra::mask(r_huet_p1p2, cl_opp)
  } else {
    bwk_uitvoer     <- master_grid * NA
    plassen_uitvoer <- master_grid * NA
    huet_uitvoer    <- master_grid * NA
  }
  rm(r_bwk_raw, r_plassen_schoon, r_huet_p1p2)
  gc()
  
  return(list(
    final_opp             = final_opp,
    cl_opp                = cl_opp,
    bwk_water_moeras      = bwk_uitvoer,
    watervlakken_geschoond = plassen_uitvoer,
    huet_waterlopen       = huet_uitvoer
  ))
}
