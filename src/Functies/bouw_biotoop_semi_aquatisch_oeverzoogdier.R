# ==============================================================================
# bouw_biotoop_semi_aquatisch_oeverzoogdier.R
# Archetype voor semi-aquatische oeverzoogdieren (water x oevercontact 30% x bebouwingsvrij)
# Geschikt voor: Waterspitsmuis, etc.
# ==============================================================================

bouw_biotoop_semi_aquatisch_oeverzoogdier <- function(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
  bestandsnaam <- tolower(gsub(" ", "", soort_naam))
  soort_cache  <- file.path(paths$dir_cache, bestandsnaam)
  dir.create(soort_cache, showWarnings = FALSE, recursive = TRUE)
  
  f_cache_rds <- file.path(soort_cache, "01_raw_parquet_layers.rds")
  
  # Parameters ophalen
  rij_params    <- params_soort %>% dplyr::slice(1)
  min_opp_ha    <- as.numeric(ifelse(is.na(rij_params$MinOpp_ha), 25, rij_params$MinOpp_ha))
  zoekafstand_m <- as.numeric(ifelse(is.na(rij_params$AfstandBiotopen_m), 250, rij_params$AfstandBiotopen_m))
  
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
  
  r_wat1 <- if ("water_bwk"       %in% names(lagen)) lagen[["water_bwk"]]       else (master_grid * NA)
  r_wat2 <- if ("waterlopen1_bwk" %in% names(lagen)) lagen[["waterlopen1_bwk"]] else (master_grid * NA)
  r_land <- if ("bwk_land"        %in% names(lagen)) lagen[["bwk_land"]]        else (master_grid * NA)
  
  # ----------------------------------------------------------------------------
  # 2. HUETZONES WATERLOPEN INTEGREREN IN SPOOR A (WATER)
  # ----------------------------------------------------------------------------
  message("    🌊 Hydrografische Huetzones (Bg, Bk, Rk) koppelen met BWK-water...")
  
  f_huet_al <- hulp_align_gdal(paths$huetzon, master_grid)
  r_huet    <- terra::rast(f_huet_al)
  leg_huet  <- terra::cats(r_huet)[[1]]
  
  gewenste_letters <- c("Bg", "BgK", "Bk1", "Bk2", "BkK1", "BkK2", "Rk")
  ok_huet_ids <- integer(0)
  if (!is.null(leg_huet)) {
    kol_val   <- colnames(leg_huet)[tolower(colnames(leg_huet)) %in% c("value", "id")][1]
    kol_label <- colnames(leg_huet)[tolower(colnames(leg_huet)) %in% c("label", "omschrijving", "huet")][1]
    if (is.na(kol_val)) kol_val <- 1
    if (!is.na(kol_label)) {
      ok_huet_ids <- leg_huet[[kol_val]][tolower(trimws(as.character(leg_huet[[kol_label]]))) %in% tolower(gewenste_letters)]
    }
  }
  
  if (length(ok_huet_ids) > 0) {
    r_huet_ok <- terra::ifel(!is.na(r_huet) & r_huet %in% ok_huet_ids, 1.0, NA)
  } else {
    r_huet_ok <- master_grid * NA
  }
  unlink(f_huet_al)
  rm(r_huet)
  
  # Waterbasis samenvoegen
  w1_c <- terra::ifel(is.na(r_wat1), 0, r_wat1)
  w2_c <- terra::ifel(is.na(r_wat2), 0, r_wat2)
  wh_c <- terra::ifel(is.na(r_huet_ok), 0, r_huet_ok)
  
  r_water_tot <- terra::ifel((w1_c + w2_c + wh_c) > 0, terra::clamp(w1_c + w2_c + wh_c, upper = 1.0), NA)
  rm(w1_c, w2_c, wh_c, r_wat1, r_wat2, r_huet_ok)
  gc()
  
  # ----------------------------------------------------------------------------
  # 3. BEBOUWINGSFILTER (100M BUFFER ROND WOONCLUSTERS >= 1.5 HA)
  # ----------------------------------------------------------------------------
  message("    🚫 Bebouwingsfilter berekenen (100m buffer rond woonclusters >= 1.5 ha)...")
  
  f_cache_bebouw <- file.path(paths$dir_cache, "masker_bebouwing_1_5ha_100m.tif")
  
  if (file.exists(f_cache_bebouw)) {
    mask_bebouw_100m <- terra::rast(f_cache_bebouw)
  } else {
    harde_codes <- c("u", "ua", "ub", "uc", "ud", "ui", "un", "ur")
    t_bebouw <- ds_parquet %>%
      dplyr::filter(tolower(CODE) %in% harde_codes) %>%
      dplyr::collect() %>%
      data.table::as.data.table()
    
    if (nrow(t_bebouw) > 0) {
      t_bebouw[, CODE := tolower(trimws(CODE))]
      t_cel <- unique(t_bebouw, by = c("cel_id", "CODE"))[, .(Val = 1), by = cel_id]
      r_bebouw_bin <- master_grid * NA
      r_bebouw_bin[t_cel$cel_id] <- 1
      rm(t_cel)
      
      # Clump via GRASS op 50m koppeling om clusters >= 1.5 ha (150 cellen) te detecteren
      f_b_bin <- tempfile(fileext = ".tif")
      f_cl_b  <- tempfile(fileext = ".tif")
      terra::writeRaster(r_bebouw_bin, f_b_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
      rm(r_bebouw_bin)
      
      e <- terra::ext(master_grid)
      qgisprocess::qgis_run_algorithm(
        "grass:r.clump", input = f_b_bin, output = f_cl_b, title = "cl_bebouw_15ha",
        "-d" = FALSE, GRASS_REGION_PARAMETER = sprintf("%f,%f,%f,%f", e$xmin, e$xmax, e$ymin, e$ymax),
        GRASS_REGION_CELLSIZE_PARAMETER = 10, .quiet = TRUE
      )
      unlink(f_b_bin)
      
      cl_b <- terra::rast(f_cl_b)
      terra::set.cats(cl_b, layer = 1, value = NULL)
      if (terra::has.colors(cl_b)) try(terra::coltab(cl_b) <- NULL, silent = TRUE)
      
      dt_b <- data.table::as.data.table(terra::freq(cl_b))
      setnames(dt_b, c("layer", "ID", "count"))
      b15_ids <- dt_b[!is.na(ID) & ID > 0 & count >= 150, ID] # 150 cellen = 1.5 ha
      rm(dt_b)
      
      if (length(b15_ids) > 0) {
        r_b15_bin <- terra::ifel(cl_b %in% as.integer(b15_ids), 1, NA)
        f_b15_bin <- tempfile(fileext = ".tif")
        terra::writeRaster(r_b15_bin, f_b15_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
        rm(r_b15_bin)
        
        f_prox_b <- hulp_proximity_gdal(f_b15_bin, 100)
        unlink(f_b15_bin)
        
        p_b <- terra::rast(f_prox_b)
        mask_bebouw_100m <- terra::ifel(!is.na(p_b) & p_b >= 0 & p_b <= 100, 1, NA)
        unlink(f_prox_b)
        rm(p_b)
      } else {
        mask_bebouw_100m <- master_grid * NA
      }
      unlink(f_cl_b)
      rm(cl_b)
      gc()
    } else {
      mask_bebouw_100m <- master_grid * NA
    }
    rm(t_bebouw)
    terra::writeRaster(mask_bebouw_100m, f_cache_bebouw, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
  }
  
  # ----------------------------------------------------------------------------
  # 4. SYMMETRISCHE 30% OEVERCHECK & AANGESLOTEN 50M OEVERSTROOK
  # ----------------------------------------------------------------------------
  message("    🐸 Symmetrische oevercheck (30% contact) & 50m oeverstrook berekenen...")
  
  r_bin_w <- terra::ifel(!is.na(r_water_tot) & r_water_tot > 0, 1, 0)
  r_bin_l <- terra::ifel(!is.na(r_land)      & r_land > 0,      1, 0)
  
  m_buren <- matrix(c(1, 1, 1,
                      1, 0, 1,
                      1, 1, 1), nrow = 3, ncol = 3)
  
  # A. Waterclusters toetsen op 30% landcontact langs oever
  land_buren_w <- terra::focal(r_bin_l, w = m_buren, fun = "sum", na.rm = TRUE)
  land_buren_w <- terra::mask(land_buren_w, r_bin_w == 1, maskvalues = FALSE)
  
  niet_water    <- terra::ifel(r_bin_w == 0, 1, 0)
  oever_tot_w   <- terra::focal(niet_water, w = m_buren, fun = "sum", na.rm = TRUE)
  oever_tot_w   <- terra::mask(oever_tot_w, r_bin_w == 1, maskvalues = FALSE)
  rm(niet_water)
  
  cl_w <- terra::patches(r_water_tot, directions = 8, zeroAsNA = TRUE)
  
  df_tot_w  <- terra::zonal(oever_tot_w, cl_w, fun = "sum", na.rm = TRUE)
  df_land_w <- terra::zonal(land_buren_w, cl_w, fun = "sum", na.rm = TRUE)
  colnames(df_tot_w)  <- c("ID", "Tot_Oever")
  colnames(df_land_w) <- c("ID", "Land_Contact")
  
  dt_w <- merge(data.table::as.data.table(df_tot_w), data.table::as.data.table(df_land_w), by = "ID", all.x = TRUE)
  dt_w[is.na(Land_Contact), Land_Contact := 0]
  dt_w[, Ratio := ifelse(Tot_Oever > 0, Land_Contact / Tot_Oever, 0)]
  ok_w_ids <- dt_w[!is.na(ID) & Ratio >= 0.30, ID]
  rm(land_buren_w, oever_tot_w, df_tot_w, df_land_w, dt_w)
  gc()
  
  if (length(ok_w_ids) > 0) {
    water_ok <- terra::mask(r_water_tot, cl_w %in% as.integer(ok_w_ids))
  } else {
    water_ok <- master_grid * NA
  }
  rm(cl_w)
  
  # B. Landclusters toetsen op 30% watercontact direct
  water_buren_l <- terra::focal(r_bin_w, w = m_buren, fun = "sum", na.rm = TRUE)
  water_buren_l <- terra::mask(water_buren_l, r_bin_l == 1, maskvalues = FALSE)
  
  niet_land   <- terra::ifel(r_bin_l == 0, 1, 0)
  rand_tot_l  <- terra::focal(niet_land, w = m_buren, fun = "sum", na.rm = TRUE)
  rand_tot_l  <- terra::mask(rand_tot_l, r_bin_l == 1, maskvalues = FALSE)
  rm(niet_land, r_bin_w, r_bin_l)
  
  cl_l <- terra::patches(r_land, directions = 8, zeroAsNA = TRUE)
  
  df_tot_l  <- terra::zonal(rand_tot_l, cl_l, fun = "sum", na.rm = TRUE)
  df_wat_l  <- terra::zonal(water_buren_l, cl_l, fun = "sum", na.rm = TRUE)
  colnames(df_tot_l) <- c("ID", "Tot_Rand")
  colnames(df_wat_l) <- c("ID", "Wat_Contact")
  
  dt_l <- merge(data.table::as.data.table(df_tot_l), data.table::as.data.table(df_wat_l), by = "ID", all.x = TRUE)
  dt_l[is.na(Wat_Contact), Wat_Contact := 0]
  dt_l[, Ratio := ifelse(Tot_Rand > 0, Wat_Contact / Tot_Rand, 0)]
  ok_l_ids <- dt_l[!is.na(ID) & Ratio >= 0.30, ID]
  rm(water_buren_l, rand_tot_l, df_tot_l, df_wat_l, dt_l)
  gc()
  
  if (length(ok_l_ids) > 0) {
    land_ok_direct <- terra::mask(r_land, cl_l %in% as.integer(ok_l_ids))
  } else {
    land_ok_direct <- master_grid * NA
  }
  rm(cl_l)
  
  # C. 50 m oeverbuffer rond goedgekeurd water met fysiek contact (<= 10 m)
  if (!all(is.na(suppressWarnings(terra::minmax(water_ok))))) {
    f_w_ok <- tempfile(fileext = ".tif")
    terra::writeRaster(terra::ifel(!is.na(water_ok), 1, NA), f_w_ok, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    f_prox_w50 <- hulp_proximity_gdal(f_w_ok, 50)
    unlink(f_w_ok)
    
    p_w50 <- terra::rast(f_prox_w50)
    mask_w50 <- terra::ifel(!is.na(p_w50) & p_w50 >= 0 & p_w50 <= 50, 1, NA)
    unlink(f_prox_w50)
    rm(p_w50)
    
    land_oever_50m_raw <- terra::mask(r_land, mask_w50)
    rm(mask_w50)
    
    # Raak-check: enkel land dat binnen 10m van water ligt
    cl_land_oever <- terra::patches(land_oever_50m_raw, directions = 8, zeroAsNA = TRUE)
    
    f_w_touch <- tempfile(fileext = ".tif")
    terra::writeRaster(terra::ifel(!is.na(water_ok), 1, NA), f_w_touch, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    f_prox_w10 <- hulp_proximity_gdal(f_w_touch, 10)
    unlink(f_w_touch)
    
    p_w10 <- terra::rast(f_prox_w10)
    mask_w10 <- terra::ifel(!is.na(p_w10) & p_w10 >= 0 & p_w10 <= 10, 1, NA)
    unlink(f_prox_w10)
    rm(p_w10)
    
    land_touching <- terra::mask(cl_land_oever, mask_w10)
    rm(mask_w10)
    ids_touch <- unique(na.omit(as.integer(land_touching[])))
    rm(land_touching)
    
    if (length(ids_touch) > 0) {
      land_oever_50m <- terra::mask(land_oever_50m_raw, cl_land_oever %in% ids_touch)
    } else {
      land_oever_50m <- master_grid * NA
    }
    rm(cl_land_oever, land_oever_50m_raw)
    gc()
  } else {
    land_oever_50m <- master_grid * NA
  }
  
  # Minicluster samenstellen: Water + Direct Land + 50m Oeverstrook
  w_c  <- terra::ifel(is.na(water_ok), 0, water_ok)
  l1_c <- terra::ifel(is.na(land_ok_direct), 0, land_ok_direct)
  l2_c <- terra::ifel(is.na(land_oever_50m), 0, land_oever_50m)
  
  r_mini_raw <- terra::ifel((w_c + l1_c + l2_c) > 0, terra::clamp(w_c + l1_c + l2_c, upper = 1.0), NA)
  rm(w_c, l1_c, l2_c, water_ok, land_ok_direct, land_oever_50m, r_water_tot, r_land)
  gc()
  
  # Bebouwingszone (100m) uitsluiten van minicluster
  r_mini_schoon <- terra::mask(r_mini_raw, mask_bebouw_100m, inverse = TRUE)
  rm(r_mini_raw, mask_bebouw_100m)
  gc()
  
  ha_mini <- calc_ha_exact(r_mini_schoon)
  message(sprintf("       -> Miniclusters (oever + water) na bebouwingsuitsluiting: %.2f ha", ha_mini))
  
  # ----------------------------------------------------------------------------
  # 5. LANDSCHAPPELIJKE CLUSTERING (>= 25 HA OP 250M ZOEFAFSTAND)
  # ----------------------------------------------------------------------------
  message(sprintf("    🏝️ Miniclusters koppelen op %dm en filteren op >= %.1f ha...", as.integer(zoekafstand_m), min_opp_ha))
  
  res_cluster <- cluster_filter_universeel(
    masker            = r_mini_schoon,
    opp_laag          = r_mini_schoon,
    min_opp_ha        = min_opp_ha,
    min_totale_opp_ha = NA,
    afstand_m         = zoekafstand_m,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  
  final_opp <- res_cluster$raster
  cl_opp    <- res_cluster$clusters
  rm(res_cluster, r_mini_schoon)
  gc()
  
  ha_finaal <- calc_ha_exact(final_opp)
  message(sprintf("       -> Definitief leefgebied Waterspitsmuis: %.2f ha", ha_finaal))
  
  return(list(
    final_opp = final_opp,
    cl_opp    = cl_opp
  ))
}
