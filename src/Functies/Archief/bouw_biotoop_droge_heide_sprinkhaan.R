# ==============================================================================
# bouw_biotoop_droge_heide_sprinkhaan.R
# Archetype voor vleugellamme heidesprinkhanen (droge heide x drainage a-c x >= 5ha direct OF 1-5ha met >= 6 buren op 500m)
# Geschikt voor: Zadelsprinkhaan (Ephippiger ephippiger), etc.
# ==============================================================================

bouw_biotoop_droge_heide_sprinkhaan <- function(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
  bestandsnaam <- tolower(gsub(" ", "", soort_naam))
  soort_cache  <- file.path(paths$dir_cache, bestandsnaam)
  dir.create(soort_cache, showWarnings = FALSE, recursive = TRUE)
  
  f_cache_rds <- file.path(soort_cache, "01_raw_parquet_layers.rds")
  
  # Parameters ophalen uit Excel
  rij_params    <- params_soort %>% dplyr::slice(1)
  min_opp_ha    <- as.numeric(ifelse(is.na(rij_params$MinOpp_ha), 1.0, rij_params$MinOpp_ha))
  groot_opp_ha  <- as.numeric(ifelse(is.na(rij_params$MinTotaleOpp_ha), 5.0, rij_params$MinTotaleOpp_ha))
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
  
  r_bwk_raw <- if ("bwk" %in% names(lagen)) lagen[["bwk"]] else (master_grid * NA)
  
  ha_bwk_ruw <- calc_ha_exact(r_bwk_raw)
  message(sprintf("       -> Ruw droog open biotoop: %.2f ha", ha_bwk_ruw))
  
  # ----------------------------------------------------------------------------
  # 2. BODEMDRAINAGE-CONFIGURATIE: STRIKTE DROOGTE (KLASSEN A, B, C, A-B)
  # ----------------------------------------------------------------------------
  message("    💧 Bodemdrainage filteren op strikt droge bodems (a, b, c, a-b)...")
  
  f_cache_drain_strikt <- file.path(paths$dir_cache, "masker_drainage_strikt_droog.tif")
  
  if (file.exists(f_cache_drain_strikt)) {
    message("       -> [CACHE] Strikt droge drainagebuffer ingeladen van schijf...")
    mask_drain_strikt <- terra::rast(f_cache_drain_strikt)
  } else if (!is.null(paths$drainage) && file.exists(paths$drainage)) {
    f_drain_al <- hulp_align_gdal(paths$drainage, master_grid)
    r_drain    <- terra::rast(f_drain_al)
    leg_drain  <- terra::cats(r_drain)[[1]]
    
    strikt_letters <- c("a", "b", "c", "a-b")
    ok_strikt_ids  <- integer(0)
    
    if (!is.null(leg_drain)) {
      kol_val   <- colnames(leg_drain)[tolower(colnames(leg_drain)) %in% c("value", "id")][1]
      kol_label <- colnames(leg_drain)[tolower(colnames(leg_drain)) %in% c("label", "omschrijving", "drainage")][1]
      if (is.na(kol_val)) kol_val <- 1
      if (!is.na(kol_label)) {
        drain_labels_clean <- tolower(trimws(as.character(leg_drain[[kol_label]])))
        ok_strikt_ids      <- leg_drain[[kol_val]][drain_labels_clean %in% strikt_letters]
      }
    }
    
    if (length(ok_strikt_ids) == 0) ok_strikt_ids <- c(1, 2, 3, 12)
    
    mask_drain_strikt <- terra::ifel(!is.na(r_drain) & (r_drain %in% ok_strikt_ids), 1, NA)
    unlink(f_drain_al)
    rm(r_drain)
    terra::writeRaster(mask_drain_strikt, f_cache_drain_strikt, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
  } else {
    mask_drain_strikt <- master_grid * NA
  }
  
  r_bio_strikt <- terra::mask(r_bwk_raw, mask_drain_strikt)
  rm(mask_drain_strikt, r_bwk_raw)
  gc()
  
  ha_strikt <- calc_ha_exact(r_bio_strikt)
  message(sprintf("       -> Doelbiotoop op strikt droge bodems: %.2f ha", ha_strikt))
  
  # ----------------------------------------------------------------------------
  # 3. PRIMAIRE CLUSTERING OP 50 M (MIN. 1 HA)
  # ----------------------------------------------------------------------------
  message(sprintf("    🏝️ Primaire clustering op landschapsniveau (>= %.1f ha, %dm netwerk)...", min_opp_ha, as.integer(zoekafstand_m)))
  
  res_1ha <- cluster_filter_universeel(
    masker            = r_bio_strikt,
    opp_laag          = r_bio_strikt,
    min_opp_ha        = min_opp_ha,
    min_totale_opp_ha = NA,
    afstand_m         = zoekafstand_m,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_1ha  <- res_1ha$raster
  cl_1ha <- res_1ha$clusters
  rm(res_1ha, r_bio_strikt)
  gc()
  
  # ----------------------------------------------------------------------------
  # 4. METAPOPULATIE STRUCTUUR ANALYSE
  # ----------------------------------------------------------------------------
  message(sprintf("    📊 Patchgrootte differentiëren: Kernen >= %.1f ha vs. Netwerk (1-5 ha met >= 6 buren op 500m)...", groot_opp_ha))
  
  if (!all(is.na(suppressWarnings(terra::minmax(cl_1ha))))) {
    # Oppervlakte per unieke patch via terra::zonal
    df_opp <- terra::zonal(r_1ha, cl_1ha, fun = "sum", na.rm = TRUE)
    colnames(df_opp) <- c("ID", "Opp_Sum")
    dt_opp <- data.table::as.data.table(df_opp)
    dt_opp[, Opp_ha := Opp_Sum * 0.01]
    rm(df_opp)
    
    ids_groot <- dt_opp[!is.na(ID) & Opp_ha >= groot_opp_ha, ID]
    ids_klein <- dt_opp[!is.na(ID) & Opp_ha >= min_opp_ha & Opp_ha < groot_opp_ha, ID]
    rm(dt_opp)
    
    # Grote kernen (>= 5 ha): direct goedgekeurd
    if (length(ids_groot) > 0) {
      r_groot <- terra::mask(r_1ha, cl_1ha %in% as.integer(ids_groot))
    } else {
      r_groot <- master_grid * NA
    }
    
    # Kleine patches (1-5 ha): moeten >= 6 buren hebben binnen 500m
    if (length(ids_klein) > 0) {
      r_klein_mask <- cl_1ha %in% as.integer(ids_klein)
      cl_klein     <- terra::mask(cl_1ha, r_klein_mask)
      rm(r_klein_mask)
      
      p_poly <- terra::as.polygons(cl_klein, dissolve = TRUE)
      
      if (!is.null(p_poly) && nrow(p_poly) > 0) {
        pts_cent <- terra::centroids(p_poly, inside = TRUE)
        buf_500m <- terra::buffer(pts_cent, width = 500)
        rel_mat  <- terra::is.related(buf_500m, pts_cent, "intersects")
        
        pts_cent$buren_count <- rowSums(rel_mat)
        ok_pts <- pts_cent[pts_cent$buren_count >= 6, ]
        
        if (nrow(ok_pts) > 0) {
          id_col_name <- names(p_poly)[1]
          valide_klein_ids <- ok_pts[[id_col_name]]
          r_klein_ok <- terra::mask(r_1ha, cl_1ha %in% as.integer(valide_klein_ids))
        } else {
          r_klein_ok <- master_grid * NA
        }
        rm(pts_cent, buf_500m, rel_mat, p_poly)
      } else {
        r_klein_ok <- master_grid * NA
      }
      rm(cl_klein)
    } else {
      r_klein_ok <- master_grid * NA
    }
    rm(cl_1ha, r_1ha)
    gc()
    
    # Samenvoegen
    final_opp <- terra::cover(r_groot, r_klein_ok)
    final_opp <- terra::clamp(final_opp, upper = 1.0)
  } else {
    r_groot    <- master_grid * NA
    r_klein_ok <- master_grid * NA
    final_opp  <- master_grid * NA
  }
  
  ha_groot <- calc_ha_exact(r_groot)
  ha_klein <- calc_ha_exact(r_klein_ok)
  ha_fin   <- calc_ha_exact(final_opp)
  
  message(sprintf("       -> Grote autonome kernen (>= 5 ha): %.2f ha", ha_groot))
  message(sprintf("       -> Goedgekeurde metapopulatiepatches (1-5 ha): %.2f ha", ha_klein))
  message(sprintf("       -> Definitief leefgebied Zadelsprinkhaan: %.2f ha", ha_fin))
  
  # ----------------------------------------------------------------------------
  # 5. NETWERK CLUSTERING (CL_OPP)
  # ----------------------------------------------------------------------------
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
    final_opp              = final_opp,
    cl_opp                 = cl_opp,
    grote_kernen           = r_groot,
    metapopulatie_patches  = r_klein_ok
  ))
}
