# ==============================================================================
# bouw_biotoop_verlandingsplas_libel.R
# Archetype voor verlandingslibellen (water <= 0.4ha met >= 4 plassen op 500m x land >= 50ha met 70% struweel op 800m x complex >= 35ha)
# Geschikt voor: Sierlijke witsnuitlibel (Leucorrhinia pectoralis), etc.
# ==============================================================================

bouw_biotoop_verlandingsplas_libel <- function(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
  bestandsnaam <- tolower(gsub(" ", "", soort_naam))
  soort_cache  <- file.path(paths$dir_cache, bestandsnaam)
  dir.create(soort_cache, showWarnings = FALSE, recursive = TRUE)
  
  f_cache_rds <- file.path(soort_cache, "01_raw_parquet_layers.rds")
  
  # Parameters ophalen uit Excel
  rij_wat <- params_soort %>% dplyr::filter(grepl("voortplanting", tolower(Type))) %>% dplyr::slice(1)
  min_opp_water <- as.numeric(ifelse(nrow(rij_wat) > 0 && !is.na(rij_wat$MinOpp_ha), rij_wat$MinOpp_ha, 0.01))
  
  rij_foer <- params_soort %>% dplyr::filter(grepl("foerageer", tolower(Type))) %>% dplyr::slice(1)
  min_opp_foer <- as.numeric(ifelse(nrow(rij_foer) > 0 && !is.na(rij_foer$MinOpp_ha), rij_foer$MinOpp_ha, 50.0))
  zoekafstand_foer_m <- as.numeric(ifelse(nrow(rij_foer) > 0 && !is.na(rij_foer$AfstandBiotopen_m), rij_foer$AfstandBiotopen_m, 25))
  
  koppel_m       <- 800 # 800 m interactie- en foerageerstraal
  min_complex_ha <- 35.0 # Minimaal 35 ha gecombineerd leefgebied
  
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
  
  r_wat_raw   <- if ("voortplanting_bwk" %in% names(lagen)) lagen[["voortplanting_bwk"]] else (master_grid * NA)
  r_foer_raw  <- if ("foerageer0"        %in% names(lagen)) lagen[["foerageer0"]]        else (master_grid * NA)
  r_struw_raw <- if ("struweel1"         %in% names(lagen)) lagen[["struweel1"]]         else (master_grid * NA)
  
  # ----------------------------------------------------------------------------
  # 2. KLEINE WATERVLAKKEN (<= 4.000 M²) INTEGREREN EN PLASSENCLUSTEREN
  # ----------------------------------------------------------------------------
  message("    💧 Watervlakken koppelen en filteren op netwerkconnectiviteit (>= 4 plassen binnen 500m)...")
  
  if (!is.null(paths$watervlakken_raster) && file.exists(paths$watervlakken_raster)) {
    f_wv_al <- hulp_align_gdal(paths$watervlakken_raster, master_grid)
    r_wv    <- terra::rast(f_wv_al)
    r_wv_bin <- terra::ifel(!is.na(r_wv) & r_wv > 0, 1.0, NA)
    unlink(f_wv_al)
    rm(r_wv)
    r_wat_tot <- terra::cover(r_wat_raw, r_wv_bin)
    rm(r_wv_bin)
  } else {
    r_wat_tot <- r_wat_raw
  }
  r_wat_tot <- terra::clamp(r_wat_tot, upper = 1.0)
  rm(r_wat_raw)
  gc()
  
  # Afzonderlijke plassen identificeren (8-richtingen buurcellen)
  r_wat_bin <- terra::ifel(!is.na(r_wat_tot) & r_wat_tot > 0, 1, NA)
  cl_water_ids <- terra::patches(r_wat_bin, directions = 8, zeroAsNA = TRUE)
  rm(r_wat_bin)
  
  if (!all(is.na(suppressWarnings(terra::minmax(cl_water_ids))))) {
    p_wat_poly <- terra::as.polygons(cl_water_ids, dissolve = TRUE)
    
    if (!is.null(p_wat_poly) && nrow(p_wat_poly) > 0) {
      pts_wat  <- terra::centroids(p_wat_poly, inside = TRUE)
      buf_500m <- terra::buffer(pts_wat, width = 500)
      rel_mat  <- terra::is.related(buf_500m, pts_wat, "intersects")
      
      pts_wat$buren_count <- rowSums(rel_mat)
      ok_pts <- pts_wat[pts_wat$buren_count >= 4, ] # Minstens 4 plassen binnen 500 m
      
      if (nrow(ok_pts) > 0) {
        id_col_name <- names(p_wat_poly)[1]
        valide_poel_ids <- ok_pts[[id_col_name]]
        r_water_netwerk <- terra::mask(r_wat_tot, cl_water_ids %in% as.integer(valide_poel_ids))
      } else {
        r_water_netwerk <- master_grid * NA
      }
      rm(pts_wat, buf_500m, rel_mat, p_wat_poly)
    } else {
      r_water_netwerk <- master_grid * NA
    }
  } else {
    r_water_netwerk <- master_grid * NA
  }
  rm(cl_water_ids, r_wat_tot)
  gc()
  
  ha_wat_net <- calc_ha_exact(r_water_netwerk)
  message(sprintf("       -> Goedgekeurde voortplantingsplassen in netwerk (>= 4 plassen / 500m): %.2f ha", ha_wat_net))
  
  # ----------------------------------------------------------------------------
  # 3. FOERAGEERHABITAT FILTEREN (>= 50 HA OP 25 M, >= 70% STRUWEEL)
  # ----------------------------------------------------------------------------
  message(sprintf("    🌳 Foerageerhabitat clusteren (>= %.1f ha op %dm) en filteren op >= 70%% struweel...", min_opp_foer, as.integer(zoekafstand_foer_m)))
  
  res_foer <- cluster_filter_universeel(
    masker            = r_foer_raw,
    opp_laag          = r_foer_raw,
    min_opp_ha        = min_opp_foer,
    min_totale_opp_ha = NA,
    afstand_m         = zoekafstand_foer_m,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_foer_50ha  <- res_foer$raster
  cl_foer_50ha <- res_foer$clusters
  rm(res_foer, r_foer_raw)
  gc()
  
  # Toetsen op minstens 70% struweeldekking in en rond het foerageercluster (50m perimeter)
  if (!all(is.na(suppressWarnings(terra::minmax(cl_foer_50ha))))) {
    r_struw_clean <- terra::ifel(!is.na(r_struw_raw) & r_struw_raw > 0, r_struw_raw, 0)
    
    # Circulaire buffer van 50m rondom de foerageerclusters voor de perimetercheck
    w_50m <- terra::focalMat(cl_foer_50ha, 50, type = 'circle')
    w_50m[w_50m > 0] <- 1
    cl_foer_buf <- terra::focal(cl_foer_50ha, w = w_50m, fun = "max", na.rm = TRUE)
    rm(w_50m)
    
    df_foer_tot <- terra::zonal(r_foer_50ha, cl_foer_50ha, fun = "sum", na.rm = TRUE)
    colnames(df_foer_tot) <- c("ID", "Totaal_Sum")
    
    struweel_in_zone <- terra::mask(r_struw_clean, cl_foer_buf)
    df_struw_tot <- terra::zonal(struweel_in_zone, cl_foer_buf, fun = "sum", na.rm = TRUE)
    colnames(df_struw_tot) <- c("ID", "Struweel_Sum")
    rm(cl_foer_buf, struweel_in_zone, r_struw_clean)
    
    dt_foer <- merge(data.table::as.data.table(df_foer_tot), data.table::as.data.table(df_struw_tot), by = "ID", all.x = TRUE)
    dt_foer[is.na(Struweel_Sum), Struweel_Sum := 0]
    dt_foer[, Pct_Struw := (Struweel_Sum / Totaal_Sum) * 100]
    
    valide_foer_ids <- dt_foer[!is.na(ID) & Pct_Struw >= 70, ID]
    rm(df_foer_tot, df_struw_tot, dt_foer)
    
    if (length(valide_foer_ids) > 0) {
      r_foer_gevalideerd <- terra::mask(r_struw_raw, cl_foer_50ha %in% as.integer(valide_foer_ids))
    } else {
      r_foer_gevalideerd <- master_grid * NA
    }
  } else {
    r_foer_gevalideerd <- master_grid * NA
  }
  rm(r_foer_50ha, cl_foer_50ha, r_struw_raw)
  gc()
  
  ha_foer_val <- calc_ha_exact(r_foer_gevalideerd)
  message(sprintf("       -> Goedgekeurd jachtstruweel (>= 70%% struweel/bos): %.2f ha", ha_foer_val))
  
  # ----------------------------------------------------------------------------
  # 4. WEDERZIJDSE KOPPELING BINNEN 800 M
  # ----------------------------------------------------------------------------
  message(sprintf("    🎯 Wederzijdse interactie binnen %dm berekenen via GDAL proximity...", as.integer(koppel_m)))
  
  heeft_water <- !all(is.na(suppressWarnings(terra::minmax(r_water_netwerk))))
  heeft_foer  <- !all(is.na(suppressWarnings(terra::minmax(r_foer_gevalideerd))))
  
  if (heeft_water && heeft_foer) {
    # Proximity vanaf water naar foerageergebied
    f_w_bin  <- tempfile(fileext = ".tif")
    terra::writeRaster(terra::ifel(!is.na(r_water_netwerk) & r_water_netwerk > 0, 1, NA), f_w_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    f_prox_w <- hulp_proximity_gdal(f_w_bin, koppel_m)
    unlink(f_w_bin)
    
    p_w <- terra::rast(f_prox_w)
    mask_w800 <- terra::ifel(!is.na(p_w) & p_w >= 0 & p_w <= koppel_m, 1, NA)
    unlink(f_prox_w)
    rm(p_w)
    
    # Proximity vanaf foerageergebied naar water
    f_f_bin  <- tempfile(fileext = ".tif")
    terra::writeRaster(terra::ifel(!is.na(r_foer_gevalideerd) & r_foer_gevalideerd > 0, 1, NA), f_f_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    f_prox_f <- hulp_proximity_gdal(f_f_bin, koppel_m)
    unlink(f_f_bin)
    
    p_f <- terra::rast(f_prox_f)
    mask_f800 <- terra::ifel(!is.na(p_f) & p_f >= 0 & p_f <= koppel_m, 1, NA)
    unlink(f_prox_f)
    rm(p_f)
    
    r_wat_gekoppeld  <- terra::mask(r_water_netwerk, mask_f800)
    r_foer_gekoppeld <- terra::mask(r_foer_gevalideerd, mask_w800)
    rm(mask_w800, mask_f800)
  } else {
    r_wat_gekoppeld  <- master_grid * NA
    r_foer_gekoppeld <- master_grid * NA
  }
  rm(r_water_netwerk, r_foer_gevalideerd)
  gc()
  
  # ----------------------------------------------------------------------------
  # 5. COMPLEX-OMVANG (>= 35 HA OP 25 M) & WATER-FIRST 800 M AFSNIJDING
  # ----------------------------------------------------------------------------
  message(sprintf("    🧩 Gecombineerd complex clusteren (>= %.1f ha op 25m) & land afkappen op %dm rond water...", min_complex_ha, as.integer(koppel_m)))
  
  r_complex_ruw <- terra::cover(r_wat_gekoppeld, r_foer_gekoppeld)
  r_complex_ruw <- terra::clamp(r_complex_ruw, upper = 1.0)
  
  res_complex <- cluster_filter_universeel(
    masker            = r_complex_ruw,
    opp_laag          = r_complex_ruw,
    min_opp_ha        = min_complex_ha,
    min_totale_opp_ha = NA,
    afstand_m         = 25,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_complex_35ha  <- res_complex$raster
  cl_complex_35ha <- res_complex$clusters
  rm(res_complex, r_complex_ruw)
  gc()
  
  if (!all(is.na(suppressWarnings(terra::minmax(r_complex_35ha))))) {
    water_in_complex <- terra::mask(r_wat_gekoppeld, cl_complex_35ha)
    land_in_complex  <- terra::mask(r_foer_gekoppeld, cl_complex_35ha)
    
    # Water-first afsnijding: land strikt binnen 800m van het goedgekeurde water
    if (!all(is.na(suppressWarnings(terra::minmax(water_in_complex))))) {
      f_wc_bin  <- tempfile(fileext = ".tif")
      terra::writeRaster(terra::ifel(!is.na(water_in_complex) & water_in_complex > 0, 1, NA), f_wc_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
      f_prox_wc <- hulp_proximity_gdal(f_wc_bin, koppel_m)
      unlink(f_wc_bin)
      
      p_wc <- terra::rast(f_prox_wc)
      mask_wc800 <- terra::ifel(!is.na(p_wc) & p_wc >= 0 & p_wc <= koppel_m, 1, NA)
      unlink(f_prox_wc)
      rm(p_wc)
      
      water_finaal <- water_in_complex
      land_finaal  <- terra::mask(land_in_complex, mask_wc800)
      rm(mask_wc800)
    } else {
      water_finaal <- master_grid * NA
      land_finaal  <- master_grid * NA
    }
    rm(water_in_complex, land_in_complex)
  } else {
    water_finaal <- master_grid * NA
    land_finaal  <- master_grid * NA
  }
  rm(r_wat_gekoppeld, r_foer_gekoppeld, r_complex_35ha, cl_complex_35ha)
  gc()
  
  ha_wat_fin  <- calc_ha_exact(water_finaal)
  ha_land_fin <- calc_ha_exact(land_finaal)
  message(sprintf("       -> Voortplantingswater binnen goedgekeurd complex: %.2f ha", ha_wat_fin))
  message(sprintf("       -> Jachtstruweel binnen 800m van water: %.2f ha", ha_land_fin))
  
  # ----------------------------------------------------------------------------
  # 6. DEFINITIEF LEEFGEBIED & NETWERKEN
  # ----------------------------------------------------------------------------
  final_opp <- terra::cover(water_finaal, land_finaal)
  final_opp <- terra::clamp(final_opp, upper = 1.0)
  
  ha_tot <- calc_ha_exact(final_opp)
  message(sprintf("       -> Definitief leefgebied Sierlijke witsnuitlibel: %.2f ha", ha_tot))
  
  r_fin_bin <- terra::ifel(!is.na(final_opp) & final_opp > 0, 1, NA)
  
  if (!all(is.na(suppressWarnings(terra::minmax(r_fin_bin))))) {
    f_bin   <- tempfile(fileext = ".tif")
    f_clump <- tempfile(fileext = ".tif")
    terra::writeRaster(r_fin_bin, f_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    
    e <- terra::ext(master_grid)
    qgisprocess::qgis_run_algorithm(
      "grass:r.clump", input = f_bin, output = f_clump, title = paste0("cl_", bestandsnaam),
      "-d" = FALSE, GRASS_REGION_PARAMETER = sprintf("%f,%f,%f,%f", e$xmin, e$xmax, e$ymin, e$ymax),
      GRASS_REGION_CELLSIZE_PARAMETER = 10, .quiet = TRUE
    )
    unlink(f_bin)
    
    cl_opp <- terra::rast(f_clump)
    try(terra::set.cats(cl_opp, layer = 1, value = NULL), silent = TRUE)
    try(terra::coltab(cl_opp) <- NULL, silent = TRUE)
    cl_opp <- terra::mask(cl_opp, r_fin_bin)
    unlink(f_clump)
  } else {
    cl_opp <- master_grid * NA
  }
  rm(r_fin_bin)
  gc()
  
  return(list(
    final_opp        = final_opp,
    cl_opp           = cl_opp,
    water            = water_finaal,
    struweel_800m    = land_finaal
  ))
}
