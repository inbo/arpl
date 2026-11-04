# ==============================================================================
# bouw_biotoop_weidevogel_matrix.R
# Archetype voor kritische weidevogels (open polder/grasland x drainage x storingsvrij)
# Geschikt voor: Grutto, Wulp, etc.
# ==============================================================================

bouw_biotoop_weidevogel_matrix <- function(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
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
  
  r_kern1     <- if ("kern1"          %in% names(lagen)) lagen[["kern1"]]          else (master_grid * NA)
  r_kern2     <- if ("kern2"          %in% names(lagen)) lagen[["kern2"]]          else (master_grid * NA)
  r_cultuur   <- if ("cultuur0a"      %in% names(lagen)) lagen[["cultuur0a"]]      else (master_grid * NA)
  r_hooiland  <- if ("hooiland"       %in% names(lagen)) lagen[["hooiland"]]       else (master_grid * NA)
  r_c1b       <- if ("cultuurhooi1b"  %in% names(lagen)) lagen[["cultuurhooi1b"]]  else (master_grid * NA)
  r_akker     <- if ("akker"          %in% names(lagen)) lagen[["akker"]]          else (master_grid * NA)
  r_water_bwk <- if ("kleinwater_bwk" %in% names(lagen)) lagen[["kleinwater_bwk"]] else (master_grid * NA)
  
  # ----------------------------------------------------------------------------
  # 2. KLEIN WATER OPBOUWEN (GRB-WATER <= 50 CELLEN + BWK)
  # ----------------------------------------------------------------------------
  message("    💧 Klein water & plas-dras opbouwen (GRB sloten/poelen + BWK)...")
  if (!is.null(paths$grb_water) && file.exists(paths$grb_water)) {
    f_grb_aligned <- hulp_align_gdal(paths$grb_water, master_grid)
    r_grb         <- terra::rast(f_grb_aligned)
    r_grb_bin     <- terra::ifel(!is.na(r_grb) & r_grb > 0, 1, NA)
    unlink(f_grb_aligned)
    rm(r_grb)
    
    f_bin_w   <- tempfile(fileext = ".tif")
    f_clump_w <- tempfile(fileext = ".tif")
    terra::writeRaster(r_grb_bin, f_bin_w, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    
    e <- terra::ext(master_grid)
    qgisprocess::qgis_run_algorithm(
      "grass:r.clump", input = f_bin_w, output = f_clump_w, title = "cl_water",
      "-d" = FALSE, GRASS_REGION_PARAMETER = sprintf("%f,%f,%f,%f", e$xmin, e$xmax, e$ymin, e$ymax),
      GRASS_REGION_CELLSIZE_PARAMETER = 10, .quiet = TRUE
    )
    unlink(f_bin_w)
    
    cl_water <- terra::rast(f_clump_w)
    try(terra::set.cats(cl_water, layer = 1, value = NULL), silent = TRUE)
    
    dt_w <- data.table::as.data.table(terra::freq(cl_water))
    setnames(dt_w, c("layer", "ID", "count"))
    klein_ids <- dt_w[!is.na(ID) & ID > 0 & count <= 50, ID]
    
    if (length(klein_ids) > 0) {
      r_grb_klein <- terra::ifel(cl_water %in% as.integer(klein_ids), 0.01, NA)
    } else {
      r_grb_klein <- master_grid * NA
    }
    unlink(f_clump_w)
    rm(r_grb_bin, cl_water, dt_w)
    gc()
  } else {
    r_grb_klein <- master_grid * NA
  }
  
  w_grb <- terra::ifel(is.na(r_grb_klein), 0, r_grb_klein)
  w_bwk <- terra::ifel(is.na(r_water_bwk), 0, r_water_bwk)
  r_water_totaal <- terra::ifel((w_grb + w_bwk) > 0, terra::clamp(w_grb + w_bwk, upper = 1.0), NA)
  rm(w_grb, w_bwk, r_grb_klein, r_water_bwk)
  gc()
  
  # ----------------------------------------------------------------------------
  # 3. DIRECTE BOOMUISTLUITING EN STORINGSBUFFERS CONFORM BERT
  # ----------------------------------------------------------------------------
  message("    🚫 Storingsbuffers berekenen (bebouwing 100m & hooggroen 50m conform Bert)...")
  
  # A. Bebouwingsbuffer 100m (Grote kernen >= 10 ha)
  f_cache_bebouw <- file.path(paths$dir_cache, "masker_bebouwing_groot_100m.tif")
  if (file.exists(f_cache_bebouw)) {
    message("       -> [CACHE] Bebouwingsbuffer ingeladen van schijf...")
    masker_bebouw_stoor <- terra::rast(f_cache_bebouw)
  } else {
    message("       -> Bebouwingscodes ophalen en clusteren op >= 10 ha...")
    gdx_codes <- c("u", "ua", "ua-", "uc", "ud", "ui", "ui-", "un", "un-", "ur", "ur-", "ur+", "uv", "uv-", "uv+")
    r_bebouw_bin <- vul_raster_vanuit_parquet(
      ds_parquet  = ds_parquet,
      codes_exact = gdx_codes,
      master_grid = master_grid,
      waarde_vast = 1
    )
    
    if (!all(is.na(suppressWarnings(terra::minmax(r_bebouw_bin))))) {
      res_bebouw <- cluster_filter_universeel(
        masker            = r_bebouw_bin,
        opp_laag          = r_bebouw_bin,
        min_opp_ha        = 10,
        min_totale_opp_ha = NA,
        afstand_m         = 25,
        dispersie_m       = 0,
        werkelijk         = FALSE
      )
      r_bebouw_groot <- res_bebouw$raster
      
      if (!all(is.na(suppressWarnings(terra::minmax(r_bebouw_groot))))) {
        f_b_bin  <- tempfile(fileext = ".tif")
        terra::writeRaster(r_bebouw_groot, f_b_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
        f_prox_b <- hulp_proximity_gdal(f_b_bin, 100)
        unlink(f_b_bin)
        
        p_b <- terra::rast(f_prox_b)
        masker_bebouw_stoor <- terra::ifel(!is.na(p_b) & p_b >= 0 & p_b <= 100, 1, NA)
        unlink(f_prox_b)
        rm(p_b)
      } else {
        masker_bebouw_stoor <- master_grid * NA
      }
      rm(res_bebouw, r_bebouw_groot)
    } else {
      masker_bebouw_stoor <- master_grid * NA
    }
    rm(r_bebouw_bin)
    terra::writeRaster(masker_bebouw_stoor, f_cache_bebouw, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    gc()
  }
  
  # B. Hooggroen: Directe uitsluiting + 5x5 focal filter (>= 12 cellen) en 50m buffer
  f_cache_groen <- file.path(paths$dir_cache, "masker_groen_50m_cluster.tif")
  
  f_groen_aligned <- hulp_align_gdal(paths$groenkaart, master_grid)
  r_groen_raw     <- terra::rast(f_groen_aligned)
  # Zuiver hooggroen (bomen en bos = klasse 1 in Groenkaart_2021)
  r_bomen_direct  <- terra::ifel(!is.na(r_groen_raw) & r_groen_raw == 1, 1, NA)
  unlink(f_groen_aligned)
  rm(r_groen_raw)
  
  if (file.exists(f_cache_groen)) {
    message("       -> [CACHE] Hooggroen-storingsbuffer ingeladen van schijf...")
    masker_groen_stoor <- terra::rast(f_cache_groen)
  } else {
    message("       -> Houtkanten & boomclusters detecteren via 5x5 focal (>= 12 cellen) conform Bert...")
    
    r_groen_bin <- terra::ifel(!is.na(r_bomen_direct), 1, 0)
    w_mat <- matrix(1, nrow = 5, ncol = 5)
    
    f_focal_groen <- tempfile(fileext = ".tif")
    groen_som <- terra::focal(r_groen_bin, w = w_mat, fun = "sum", na.rm = TRUE, filename = f_focal_groen, overwrite = TRUE)
    rm(r_groen_bin)
    
    # Exact de drempel van Bert: >= 12 cellen in 5x5 venster en zelf boom zijn
    r_groen_cluster <- terra::ifel(!is.na(groen_som) & groen_som >= 12 & !is.na(r_bomen_direct), 1, NA)
    unlink(f_focal_groen)
    rm(groen_som)
    
    if (!all(is.na(suppressWarnings(terra::minmax(r_groen_cluster))))) {
      f_g_bin  <- tempfile(fileext = ".tif")
      terra::writeRaster(r_groen_cluster, f_g_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
      rm(r_groen_cluster)
      
      f_prox_g <- hulp_proximity_gdal(f_g_bin, 50)
      unlink(f_g_bin)
      
      p_g <- terra::rast(f_prox_g)
      masker_groen_stoor <- terra::ifel(!is.na(p_g) & p_g >= 0 & p_g <= 50, 1, NA)
      unlink(f_prox_g)
      rm(p_g)
    } else {
      masker_groen_stoor <- master_grid * NA
    }
    terra::writeRaster(masker_groen_stoor, f_cache_groen, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    gc()
  }
  
  # Gecombineerd storingsmasker
  stoor_totaal <- terra::cover(masker_bebouw_stoor, masker_groen_stoor)
  rm(masker_bebouw_stoor, masker_groen_stoor)
  gc()
  
  # ----------------------------------------------------------------------------
  # 4. SPOOR A: KERNBIOTOOP (DRAINAGEFILTER + STORINGEN + >= 10 HA CLUSTERING)
  # ----------------------------------------------------------------------------
  message("    🌱 Spoor A: Kernbiotoop filteren op natte drainageklassen & storingsvrij maken...")
  
  # Kernbiotoop = max(kern1, kern2)
  r_kern_bruto <- terra::cover(r_kern1, r_kern2)
  
  # DRAINAGEFILTER
  drain_letters <- c("d", "e", "f", "h", "i", "g", "e-f", "h-i", "e-i")
  if (!is.null(paths$drainage) && file.exists(paths$drainage)) {
    r_drain_orig <- terra::rast(paths$drainage)
    leg_drain    <- terra::cats(r_drain_orig)[[1]]
    
    ok_drain_ids <- integer(0)
    if (!is.null(leg_drain)) {
      kol_val   <- colnames(leg_drain)[tolower(colnames(leg_drain)) %in% c("value", "id")][1]
      kol_label <- colnames(leg_drain)[tolower(colnames(leg_drain)) %in% c("label", "omschrijving", "drainage")][1]
      if (!is.na(kol_val) && !is.na(kol_label)) {
        drain_labels_clean <- tolower(trimws(as.character(leg_drain[[kol_label]])))
        ok_drain_ids <- leg_drain[[kol_val]][drain_labels_clean %in% tolower(drain_letters)]
      }
    }
    if (length(ok_drain_ids) == 0) ok_drain_ids <- c(4, 5, 6, 7, 8, 9, 45, 56, 89)
    rm(r_drain_orig)
    
    f_drain_aligned <- hulp_align_gdal(paths$drainage, master_grid)
    r_drain         <- terra::rast(f_drain_aligned)
    mask_drain      <- terra::ifel(!is.na(r_drain) & (r_drain %in% ok_drain_ids), 1, NA)
    r_kern_nat      <- terra::mask(r_kern_bruto, mask_drain)
    unlink(f_drain_aligned)
    rm(r_drain, mask_drain)
  } else {
    r_kern_nat <- r_kern_bruto
  }
  rm(r_kern_bruto)
  gc()
  
  # Kern zuiveren van alle directe bomen én de storingszone
  r_kern_gezuiverd <- terra::mask(r_kern_nat, r_bomen_direct, inverse = TRUE)
  r_kern_gezuiverd <- terra::mask(r_kern_gezuiverd, stoor_totaal, inverse = TRUE)
  rm(r_kern_nat)
  gc()
  
  ha_kern_gezuiverd <- calc_ha_exact(r_kern_gezuiverd)
  message(sprintf("       -> Nat kernbiotoop storingsvrij: %.2f ha", ha_kern_gezuiverd))
  
  # Kern clusteren op >= 10 ha (50m netwerk)
  message("    🏝️ Kernbiotoop clusteren op >= 10 ha (50m netwerk)...")
  res_kern <- cluster_filter_universeel(
    masker            = r_kern_gezuiverd,
    opp_laag          = r_kern_gezuiverd,
    min_opp_ha        = 10,
    min_totale_opp_ha = NA,
    afstand_m         = 50,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_kern_10ha  <- res_kern$raster
  ha_kern_10ha <- calc_ha_exact(r_kern_10ha)
  message(sprintf("       -> Goedgekeurde kernen (>= 10 ha): %.2f ha", ha_kern_10ha))
  rm(res_kern, r_kern_gezuiverd)
  gc()
  
  # ----------------------------------------------------------------------------
  # 5. SPOOR B: MATRIX / BOUWSTENEN (CONFORM BERT: DIRECTE BOOM- & AKKERUISTLUITING)
  # ----------------------------------------------------------------------------
  message("    🌾 Spoor B: Matrix & bouwstenen opbouwen (conform Bert)...")
  
  # 1. Strikte akkerfilter: elke cel met akkerfractie uitsluiten van cultuurgrasland
  filter_akker <- terra::ifel(!is.na(r_akker) & r_akker > 0, 1, NA)
  r_cultuur_zonder_akker <- terra::mask(r_cultuur, filter_akker, inverse = TRUE)
  rm(r_cultuur, r_akker)
  
  # 2. Direct bomen uitsluiten van ELKE bouwsteen vóór samenvoeging
  r_cultuur_schoon  <- terra::mask(r_cultuur_zonder_akker, r_bomen_direct, inverse = TRUE)
  r_hooiland_schoon <- terra::mask(r_hooiland, r_bomen_direct, inverse = TRUE)
  r_water_schoon    <- terra::mask(r_water_totaal, r_bomen_direct, inverse = TRUE)
  r_c1b_schoon      <- terra::mask(r_c1b, r_bomen_direct, inverse = TRUE)
  rm(r_cultuur_zonder_akker, r_hooiland, r_water_totaal, r_c1b, r_bomen_direct)
  
  # Samenvoegen via cover conform Bert
  r_bouwsteen_tot <- terra::cover(r_cultuur_schoon, r_hooiland_schoon)
  r_bouwsteen_tot <- terra::cover(r_bouwsteen_tot, r_water_schoon)
  r_matrix_bruto  <- terra::cover(r_bouwsteen_tot, r_c1b_schoon)
  rm(r_cultuur_schoon, r_hooiland_schoon, r_water_schoon, r_c1b_schoon, r_bouwsteen_tot)
  gc()
  
  # Storingsbuffers en eventuele resterende akkers uitsluiten van de matrix
  r_matrix_gezuiverd <- terra::mask(r_matrix_bruto, filter_akker, inverse = TRUE)
  r_matrix_gezuiverd <- terra::mask(r_matrix_gezuiverd, stoor_totaal, inverse = TRUE)
  rm(r_matrix_bruto, filter_akker, stoor_totaal)
  gc()
  
  ha_matrix_gezuiverd <- calc_ha_exact(r_matrix_gezuiverd)
  message(sprintf("       -> Matrix & bouwstenen storingsvrij: %.2f ha", ha_matrix_gezuiverd))
  
  # ----------------------------------------------------------------------------
  # 6. WEDERZIJDSE KOPPELING KERN & MATRIX (200M BUFFER, NETWERK >= 25 HA, KERN >= 10 HA)
  # ----------------------------------------------------------------------------
  message("    🔗 Wederzijdse koppelings- en netwerkeis (200m buffer, >= 25 ha totaal & >= 10 ha kern)...")
  
  heeft_kern <- (ha_kern_10ha > 0)
  
  if (heeft_kern) {
    r_comb_val <- terra::cover(r_kern_10ha, r_matrix_gezuiverd)
    r_comb_bin <- terra::ifel(!is.na(r_comb_val) & r_comb_val > 0, 1, NA)
    
    # 200m netwerkvorming via GDAL proximity (100m straal)
    f_comb_bin <- tempfile(fileext = ".tif")
    terra::writeRaster(r_comb_bin, f_comb_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    f_prox_200 <- hulp_proximity_gdal(f_comb_bin, 100)
    unlink(f_comb_bin)
    
    p_200     <- terra::rast(f_prox_200)
    r_buf_200 <- terra::ifel(!is.na(p_200) & p_200 >= 0 & p_200 <= 100, 1, NA)
    unlink(f_prox_200)
    rm(p_200)
    gc()
    
    f_buf_bin <- tempfile(fileext = ".tif")
    f_cl_net  <- tempfile(fileext = ".tif")
    terra::writeRaster(r_buf_200, f_buf_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    rm(r_buf_200)
    gc()
    
    e <- terra::ext(master_grid)
    qgisprocess::qgis_run_algorithm(
      "grass:r.clump", input = f_buf_bin, output = f_cl_net, title = "cl_grutto",
      "-d" = FALSE, GRASS_REGION_PARAMETER = sprintf("%f,%f,%f,%f", e$xmin, e$xmax, e$ymin, e$ymax),
      GRASS_REGION_CELLSIZE_PARAMETER = 10, .quiet = TRUE
    )
    unlink(f_buf_bin)
    
    cl_netwerk <- terra::rast(f_cl_net)
    try(terra::set.cats(cl_netwerk, layer = 1, value = NULL), silent = TRUE)
    
    cl_net_biotoop <- terra::mask(cl_netwerk, r_comb_bin)
    unlink(f_cl_net)
    rm(cl_netwerk, r_comb_bin)
    gc()
    
    # OUT-OF-CORE ZONALE STATISTIEKEN
    df_net <- terra::zonal(r_comb_val, cl_net_biotoop, fun = "sum", na.rm = TRUE)
    colnames(df_net) <- c("netwerk_id", "totaal_ha_sum")
    dt_net <- data.table::as.data.table(df_net)
    dt_net[, totaal_ha := totaal_ha_sum * 0.01]
    rm(df_net)
    
    cl_kernen_only <- terra::mask(cl_net_biotoop, r_kern_10ha)
    df_ker <- terra::zonal(r_kern_10ha, cl_kernen_only, fun = "sum", na.rm = TRUE)
    colnames(df_ker) <- c("netwerk_id", "kern_ha_sum")
    dt_ker <- data.table::as.data.table(df_ker)
    dt_ker[, kern_ha := kern_ha_sum * 0.01]
    rm(cl_kernen_only, df_ker)
    
    valide_netwerken <- dt_net[dt_ker, on = "netwerk_id"][totaal_ha >= 25 & kern_ha >= 10, netwerk_id]
    rm(dt_net, dt_ker)
    gc()
    
    if (length(valide_netwerken) > 0) {
      cl_valide_mask <- terra::ifel(cl_net_biotoop %in% as.integer(valide_netwerken), 1, NA)
      final_opp      <- terra::mask(r_comb_val, cl_valide_mask)
      cl_opp         <- terra::mask(cl_net_biotoop, cl_valide_mask)
      kern_uitvoer   <- terra::mask(r_kern_10ha, cl_valide_mask)
      matrix_uitvoer <- terra::mask(r_matrix_gezuiverd, cl_valide_mask)
      rm(cl_valide_mask)
    } else {
      final_opp      <- master_grid * NA
      cl_opp         <- master_grid * NA
      kern_uitvoer   <- master_grid * NA
      matrix_uitvoer <- master_grid * NA
    }
    rm(cl_net_biotoop, r_comb_val)
  } else {
    final_opp      <- master_grid * NA
    cl_opp         <- master_grid * NA
    kern_uitvoer   <- master_grid * NA
    matrix_uitvoer <- master_grid * NA
  }
  rm(r_kern_10ha, r_matrix_gezuiverd)
  gc()
  
  ha_finaal <- calc_ha_exact(final_opp)
  message(sprintf("       -> Definitief leefgebied Grutto: %.2f ha", ha_finaal))
  
  return(list(
    final_opp = final_opp,
    cl_opp    = cl_opp,
    kern      = kern_uitvoer,
    matrix    = matrix_uitvoer
  ))
}
