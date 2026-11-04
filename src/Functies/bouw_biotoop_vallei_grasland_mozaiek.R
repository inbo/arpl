# ==============================================================================
# bouw_biotoop_vallei_grasland_mozaiek.R
# Archetype: Alluviale valleibodems x morfologische compactheid x graslanddominantie >= 75%
# Geschikt voor: Kwartelkoning, Watersnip, etc.
# ==============================================================================

bouw_biotoop_vallei_grasland_mozaiek <- function(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
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
  
  r_grasland_raw <- if ("grasland" %in% names(lagen)) lagen[["grasland"]] else (master_grid * NA)
  r_moeras_raw   <- if ("moeras"   %in% names(lagen)) lagen[["moeras"]]   else (master_grid * NA)
  
  # Terrils uitsluiten uit Parquet (kg codes)
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
  # 2. INITIËLE SAMENVOEGING & DIRECTE UITSLUITINGEN (HOOGGROEN & TERRILS)
  # ----------------------------------------------------------------------------
  message("    🌾 Grasland en moeras samenvoegen, hooggroen en terrils direct uitsluiten...")
  g_c <- terra::ifel(is.na(r_grasland_raw), 0, r_grasland_raw)
  m_c <- terra::ifel(is.na(r_moeras_raw),   0, r_moeras_raw)
  r_bio_bruto <- terra::ifel((g_c + m_c) > 0, terra::clamp(g_c + m_c, upper = 1.0), NA)
  rm(g_c, m_c)
  gc()
  
  # Directe uitsluiting van bomen/hooggroen en terrils
  f_groen_aligned <- hulp_align_gdal(paths$groenkaart, master_grid)
  r_groen         <- terra::rast(f_groen_aligned)
  mask_direct_fout <- terra::ifel((!is.na(r_groen) & (r_groen == 1 | r_groen > 0)) | (!is.na(r_terril_bin) & r_terril_bin == 1), 1, NA)
  unlink(f_groen_aligned)
  rm(r_groen, r_terril_bin)
  
  r_bio_fase1 <- terra::mask(r_bio_bruto, mask_direct_fout, inverse = TRUE)
  rm(r_bio_bruto, mask_direct_fout)
  gc()
  
  ha_fase1 <- calc_ha_exact(r_bio_fase1)
  message(sprintf("       -> Biotoop na uitsluiting van hooggroen en terrils: %.2f ha", ha_fase1))
  
  # ----------------------------------------------------------------------------
  # 3. VALLEIGRONDEN STUURKAART (PROFIEL p/p+x OF DRAINAGE d-i OF OVERSTROMING 1-2)
  # ----------------------------------------------------------------------------
  message("    🌊 Alluviale valleisturing toepassen (profiel, drainage en overstroming)...")
  
  # A. Bodemprofiel p, p+x
  f_cache_prof <- file.path(paths$dir_cache, "masker_profiel_p_px.tif")
  
  if (file.exists(f_cache_prof)) {
    mask_profiel <- terra::rast(f_cache_prof)
  } else if (file.exists(paths$profiel)) {
    r_prof_orig <- terra::rast(paths$profiel)
    leg_prof    <- terra::cats(r_prof_orig)[[1]]
    
    ok_prof_ids <- integer(0)
    if (!is.null(leg_prof)) {
      kol_label <- colnames(leg_prof)[tolower(colnames(leg_prof)) %in% c("label", "omschrijving", "naam", "profiel")][1]
      kol_val   <- colnames(leg_prof)[tolower(colnames(leg_prof)) %in% c("value", "id")][1]
      if (is.na(kol_val)) kol_val <- 1
      
      if (!is.na(kol_label)) {
        ok_prof_ids <- leg_prof[[kol_val]][tolower(trimws(as.character(leg_prof[[kol_label]]))) %in% c("p", "p+x")]
      }
    }
    rm(r_prof_orig)
    
    f_prof_al <- hulp_align_gdal(paths$profiel, master_grid)
    r_prof    <- terra::rast(f_prof_al)
    
    if (length(ok_prof_ids) > 0) {
      mask_profiel <- terra::ifel(!is.na(r_prof) & r_prof %in% ok_prof_ids, 1, NA)
    } else {
      # Fallback: indien profiel p numeriek gecodeerd is of geen cats heeft
      message("       ℹ️ Geen tekstuele profielcodes 'p'/'p+x' gevonden in cats; controleer profielcodering.")
      mask_profiel <- master_grid * NA
    }
    unlink(f_prof_al)
    rm(r_prof)
    
    terra::writeRaster(mask_profiel, f_cache_prof, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
  } else {
    mask_profiel <- master_grid * NA
  }
  
  # B. Drainage d, e, f, g, h, i, e-f, h-i, e-i
  mask_drainage <- master_grid * NA
  if (exists("haal_drainage_masker", mode = "function")) {
    mask_drainage <- haal_drainage_masker(paths$drainage, master_grid, paths$dir_cache)
  } else if (file.exists(paths$drainage)) {
    f_drain_al <- hulp_align_gdal(paths$drainage, master_grid)
    r_drain    <- terra::rast(f_drain_al)
    leg_drain  <- terra::cats(r_drain)[[1]]
    
    drain_codes_ok <- c("d", "e", "f", "g", "h", "i", "e-f", "h-i", "e-i")
    ok_d_ids <- integer(0)
    if (!is.null(leg_drain)) {
      kol_label_d <- colnames(leg_drain)[tolower(colnames(leg_drain)) %in% c("label", "omschrijving", "drainage")][1]
      kol_val_d   <- colnames(leg_drain)[tolower(colnames(leg_drain)) %in% c("value", "id")][1]
      if (is.na(kol_val_d)) kol_val_d <- 1
      if (!is.na(kol_label_d)) {
        ok_d_ids <- leg_drain[[kol_val_d]][tolower(trimws(as.character(leg_drain[[kol_label_d]]))) %in% drain_codes_ok]
      }
    }
    if (length(ok_d_ids) > 0) {
      mask_drainage <- terra::ifel(!is.na(r_drain) & r_drain %in% ok_d_ids, 1, NA)
    }
    unlink(f_drain_al)
    rm(r_drain)
  }
  
  # C. Overstromingsgevoeligheid 1 en 2
  f_cache_ovstr <- file.path(paths$dir_cache, "masker_overstroming_1_2.tif")
  
  if (file.exists(f_cache_ovstr)) {
    mask_overstroming <- terra::rast(f_cache_ovstr)
  } else if (file.exists(paths$overstroming)) {
    f_ovstr_al <- hulp_align_gdal(paths$overstroming, master_grid)
    r_ovstr    <- terra::rast(f_ovstr_al)
    mask_overstroming <- terra::ifel(!is.na(r_ovstr) & (r_ovstr == 1 | r_ovstr == 2), 1, NA)
    unlink(f_ovstr_al)
    rm(r_ovstr)
    terra::writeRaster(mask_overstroming, f_cache_ovstr, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
  } else {
    mask_overstroming <- master_grid * NA
  }
  
  # OR-koppeling over de drie valleistuurkaarten (profiel OF drainage OF overstroming)
  w_p <- terra::ifel(!is.na(mask_profiel), 1, 0)
  w_d <- terra::ifel(!is.na(mask_drainage), 1, 0)
  w_o <- terra::ifel(!is.na(mask_overstroming), 1, 0)
  
  vallei_geschikt <- terra::ifel((w_p + w_d + w_o) > 0, 1, NA)
  rm(w_p, w_d, w_o, mask_profiel, mask_drainage, mask_overstroming)
  gc()
  
  ha_vallei <- calc_ha_exact(vallei_geschikt)
  message(sprintf("       -> Totale geschikte valleigronden (profiel|drainage|overstr): %.2f ha", ha_vallei))
  
  if (ha_vallei > 0) {
    r_bio_fase2 <- terra::mask(r_bio_fase1, vallei_geschikt)
  } else {
    message("       ⚠️ Waarschuwing: Geen alluviale valleicellen gedetecteerd. Controleer bronrasters.")
    r_bio_fase2 <- master_grid * NA
  }
  rm(r_bio_fase1, vallei_geschikt)
  gc()
  
  ha_fase2 <- calc_ha_exact(r_bio_fase2)
  message(sprintf("       -> Biotoop binnen alluviale valleigronden: %.2f ha", ha_fase2))
  
  # ----------------------------------------------------------------------------
  # 4. MORFOLOGISCHE COMPACTHEID (MAJORITY FILTER CIRKEL 20M)
  # ----------------------------------------------------------------------------
  message("    🔬 Morfologische compactheid filteren (majority filter cirkel 20m)...")
  r_bin_f2 <- terra::ifel(!is.na(r_bio_fase2) & r_bio_fase2 > 0, 1, 0)
  
  f_mat_cirkel <- terra::focalMat(master_grid, d = 20, type = "circle")
  f_mat_cirkel[f_mat_cirkel > 0] <- 1
  drempel_cirkel <- ceiling(sum(f_mat_cirkel) / 2)
  
  f_focal_maj  <- tempfile(fileext = ".tif")
  r_som_cirkel <- terra::focal(r_bin_f2, w = f_mat_cirkel, fun = "sum", na.rm = TRUE, filename = f_focal_maj, overwrite = TRUE)
  mask_compact <- terra::ifel(!is.na(r_som_cirkel) & r_som_cirkel >= drempel_cirkel, 1, NA)
  unlink(f_focal_maj)
  rm(r_bin_f2, r_som_cirkel)
  
  r_bio_fase3 <- terra::mask(r_bio_fase2, mask_compact)
  rm(r_bio_fase2, mask_compact)
  gc()
  
  ha_fase3 <- calc_ha_exact(r_bio_fase3)
  message(sprintf("       -> Biotoop na morfologische compactheid: %.2f ha", ha_fase3))
  
  # ----------------------------------------------------------------------------
  # 5. TUSSENTIJDSE DREMPEL (>= 3 HA OP 50M ZOEFAFSTAND)
  # ----------------------------------------------------------------------------
  message("    🏝️ Tussentijdse drempel toepassen (>= 3 ha, 50m netwerk)...")
  res_3ha <- cluster_filter_universeel(
    masker            = r_bio_fase3,
    opp_laag          = r_bio_fase3,
    min_opp_ha        = 3,
    min_totale_opp_ha = NA,
    afstand_m         = 50,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_bio_fase4 <- res_3ha$raster
  rm(res_3ha, r_bio_fase3)
  gc()
  
  ha_fase4 <- calc_ha_exact(r_bio_fase4)
  message(sprintf("       -> Biotoop na 3 ha clusterdrempel: %.2f ha", ha_fase4))
  
  # ----------------------------------------------------------------------------
  # 6. NEGATIEVE HOOGGROEN BUFFER (100M BUFFER ROND CLUSTERS >= 20 ARE)
  # ----------------------------------------------------------------------------
  message("    🌲 Predatiestoring uitsluiten (100m buffer rond hooggroenclusters >= 20 are)...")
  
  f_cache_groen100 <- file.path(paths$dir_cache, "masker_groen_100m_cluster.tif")
  
  if (file.exists(f_cache_groen100)) {
    message("       -> [CACHE] Hooggroen-100m storingsbuffer ingeladen van schijf...")
    masker_groen_100m <- terra::rast(f_cache_groen100)
  } else {
    f_groen_al   <- hulp_align_gdal(paths$groenkaart, master_grid)
    r_groen_raw  <- terra::rast(f_groen_al)
    r_groen_bin  <- terra::ifel(!is.na(r_groen_raw) & (r_groen_raw == 1 | r_groen_raw > 0), 1, NA)
    unlink(f_groen_al)
    rm(r_groen_raw)
    
    # Clump via GRASS om clusters >= 20 are (20 cellen) te detecteren
    f_g_bin   <- tempfile(fileext = ".tif")
    f_cl_g    <- tempfile(fileext = ".tif")
    terra::writeRaster(r_groen_bin, f_g_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    rm(r_groen_bin)
    gc()
    
    e <- terra::ext(master_grid)
    qgisprocess::qgis_run_algorithm(
      "grass:r.clump", input = f_g_bin, output = f_cl_g, title = "cl_groen_20are",
      "-d" = FALSE, GRASS_REGION_PARAMETER = sprintf("%f,%f,%f,%f", e$xmin, e$xmax, e$ymin, e$ymax),
      GRASS_REGION_CELLSIZE_PARAMETER = 10, .quiet = TRUE
    )
    unlink(f_g_bin)
    
    cl_g <- terra::rast(f_cl_g)
    terra::set.cats(cl_g, layer = 1, value = NULL)
    if (terra::has.colors(cl_g)) try(terra::coltab(cl_g) <- NULL, silent = TRUE)
    
    dt_g <- data.table::as.data.table(terra::freq(cl_g))
    setnames(dt_g, c("layer", "ID", "count"))
    g20_ids <- dt_g[!is.na(ID) & ID > 0 & count >= 20, ID]
    rm(dt_g)
    
    if (length(g20_ids) > 0) {
      r_g20_bin <- terra::ifel(cl_g %in% as.integer(g20_ids), 1, NA)
      f_g20_bin <- tempfile(fileext = ".tif")
      terra::writeRaster(r_g20_bin, f_g20_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
      rm(r_g20_bin)
      
      f_prox_g <- hulp_proximity_gdal(f_g20_bin, 100)
      unlink(f_g20_bin)
      
      p_g <- terra::rast(f_prox_g)
      masker_groen_100m <- terra::ifel(!is.na(p_g) & p_g >= 0 & p_g <= 100, 1, NA)
      unlink(f_prox_g)
      rm(p_g)
    } else {
      masker_groen_100m <- master_grid * NA
    }
    unlink(f_cl_g)
    rm(cl_g)
    gc()
    
    terra::writeRaster(masker_groen_100m, f_cache_groen100, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
  }
  
  r_bio_fase5 <- terra::mask(r_bio_fase4, masker_groen_100m, inverse = TRUE)
  rm(r_bio_fase4, masker_groen_100m)
  gc()
  
  ha_fase5 <- calc_ha_exact(r_bio_fase5)
  message(sprintf("       -> Biotoop buiten 100m hooggroenbuffer: %.2f ha", ha_fase5))
  
  # ----------------------------------------------------------------------------
  # 7. FUNCTIONEEL LEEFGEBIED (>= 10 HA) & GRASLANDAANDEEL (>= 75%)
  # ----------------------------------------------------------------------------
  message("    📊 Functioneel complex clusteren (>= 10 ha) en aandeel grasland (>= 75%) toetsen...")
  
  res_10ha <- cluster_filter_universeel(
    masker            = r_bio_fase5,
    opp_laag          = r_bio_fase5,
    min_opp_ha        = 10,
    min_totale_opp_ha = NA,
    afstand_m         = 50,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_bio_10ha <- res_10ha$raster
  cl_10ha    <- res_10ha$clusters
  rm(res_10ha, r_bio_fase5)
  gc()
  
  if (!all(is.na(suppressWarnings(terra::minmax(cl_10ha))))) {
    cells_10 <- terra::cells(!is.na(cl_10ha) & cl_10ha > 0)
    dt_gras_check <- data.table::data.table(
      ID       = as.integer(cl_10ha[cells_10][[1]]),
      Tot_Val  = as.numeric(r_bio_10ha[cells_10][[1]]),
      Gras_Val = as.numeric(terra::ifel(is.na(r_grasland_raw[cells_10][[1]]), 0, r_grasland_raw[cells_10][[1]]))
    )
    rm(cells_10)
    gc()
    
    dt_agg <- dt_gras_check[, .(
      Tot_ha  = sum(Tot_Val, na.rm = TRUE) * 0.01,
      Gras_ha = sum(Gras_Val, na.rm = TRUE) * 0.01
    ), by = ID]
    dt_agg[, Grasland_Pct := (Gras_ha / Tot_ha) * 100]
    
    ok_complex_ids <- dt_agg[!is.na(ID) & Grasland_Pct >= 75, ID]
    rm(dt_gras_check, dt_agg)
    gc()
    
    if (length(ok_complex_ids) > 0) {
      cl_finaal      <- terra::ifel(cl_10ha %in% as.integer(ok_complex_ids), cl_10ha, NA)
      final_opp      <- terra::mask(r_bio_10ha, cl_finaal)
      gras_uitvoer   <- terra::mask(r_grasland_raw, cl_finaal)
      moeras_uitvoer <- terra::mask(r_moeras_raw, cl_finaal)
    } else {
      final_opp      <- master_grid * NA
      cl_finaal      <- master_grid * NA
      gras_uitvoer   <- master_grid * NA
      moeras_uitvoer <- master_grid * NA
    }
  } else {
    final_opp      <- master_grid * NA
    cl_finaal      <- master_grid * NA
    gras_uitvoer   <- master_grid * NA
    moeras_uitvoer <- master_grid * NA
  }
  rm(cl_10ha, r_bio_10ha, r_grasland_raw, r_moeras_raw)
  gc()
  
  return(list(
    final_opp = final_opp,
    cl_opp    = cl_finaal,
    grasland  = gras_uitvoer,
    moeras    = moeras_uitvoer
  ))
}
