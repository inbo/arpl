# ==============================================================================
# bouw_biotoop_droge_metapopulatie_sprinkhaan.R
# Archetype: Droogteminnende krekels/sprinkhanen (droge heide/duin x drainage a-c x patches >= 5ha OF 1-5ha met >= 6 buren op 500m)
# Geschikt voor: Veldkrekel (Gryllus campestris), etc.
# ==============================================================================

bouw_biotoop_droge_metapopulatie_sprinkhaan <- function(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
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
  
  # Ondersteun zowel generieke 'bwk' als scenario-specifieke 'bwk_strikt'
  r_bwk_raw <- if ("bwk" %in% names(lagen)) {
    lagen[["bwk"]]
  } else if ("bwk_strikt" %in% names(lagen)) {
    lagen[["bwk_strikt"]]
  } else if (length(lagen) > 0) {
    lagen[[1]]
  } else {
    master_grid * NA
  }
  
  ha_bwk_ruw <- calc_ha_exact(r_bwk_raw)
  message(sprintf("       -> Ruw droog open biotoop: %.2f ha", ha_bwk_ruw))
  
  # ----------------------------------------------------------------------------
  # 2. DRAINAGE-DIFFERENTIATIE (STRIKT VS. INCLUSIEF)
  # ----------------------------------------------------------------------------
  message("    💧 Bodemdrainage filteren: Strikt droog (a-c) vs Inclusief (!e-i)...")
  
  f_cache_drain_strikt <- file.path(paths$dir_cache, "masker_drainage_strikt_droog.tif")
  f_cache_drain_incl   <- file.path(paths$dir_cache, "masker_drainage_inclusief.tif")
  
  if (file.exists(f_cache_drain_strikt) && file.exists(f_cache_drain_incl)) {
    mask_drain_strikt <- terra::rast(f_cache_drain_strikt)
    mask_drain_incl   <- terra::rast(f_cache_drain_incl)
  } else if (!is.null(paths$drainage) && file.exists(paths$drainage)) {
    f_drain_al <- hulp_align_gdal(paths$drainage, master_grid)
    r_drain    <- terra::rast(f_drain_al)
    leg_drain  <- terra::cats(r_drain)[[1]]
    
    strikt_letters <- c("a", "b", "c", "a-b")
    nat_letters    <- c("e", "e-f", "e-i", "f", "g", "h", "h-i", "i")
    
    ok_strikt_ids <- integer(0)
    nat_ids       <- integer(0)
    
    if (!is.null(leg_drain)) {
      kol_val   <- colnames(leg_drain)[tolower(colnames(leg_drain)) %in% c("value", "id")][1]
      kol_label <- colnames(leg_drain)[tolower(colnames(leg_drain)) %in% c("label", "omschrijving", "drainage")][1]
      if (is.na(kol_val)) kol_val <- 1
      if (!is.na(kol_label)) {
        drain_labels_clean <- tolower(trimws(as.character(leg_drain[[kol_label]])))
        ok_strikt_ids      <- leg_drain[[kol_val]][drain_labels_clean %in% strikt_letters]
        nat_ids            <- leg_drain[[kol_val]][drain_labels_clean %in% nat_letters]
      }
    }
    
    # Fallback-ID's conform INBO bodemraster
    if (length(ok_strikt_ids) == 0) ok_strikt_ids <- c(1, 2, 3, 12)
    if (length(nat_ids) == 0)       nat_ids       <- c(5, 6, 7, 8, 9, 56, 89)
    
    # Strikt: enkel a-c
    mask_drain_strikt <- terra::ifel(!is.na(r_drain) & (r_drain %in% ok_strikt_ids), 1, NA)
    
    # Inclusief: alles behalve expliciet nat (ook NA, duinen en OB toegelaten)
    mask_drain_incl   <- terra::ifel(!is.na(r_drain) & (r_drain %in% nat_ids), NA, 1)
    
    unlink(f_drain_al)
    rm(r_drain)
    terra::writeRaster(mask_drain_strikt, f_cache_drain_strikt, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    terra::writeRaster(mask_drain_incl,   f_cache_drain_incl,   overwrite = TRUE, datatype = "INT1U", NAflag = 255)
  } else {
    mask_drain_strikt <- master_grid * NA
    mask_drain_incl   <- master_grid * 1
  }
  
  # Basislagen per scenario
  r_bio_strikt <- terra::mask(r_bwk_raw, mask_drain_strikt)
  r_bio_incl   <- terra::mask(r_bwk_raw, mask_drain_incl)
  rm(mask_drain_strikt, mask_drain_incl, r_bwk_raw)
  gc()
  
  ha_strikt_ruw <- calc_ha_exact(r_bio_strikt)
  ha_incl_ruw   <- calc_ha_exact(r_bio_incl)
  message(sprintf("       -> Biotoop na drainagefilter: %.2f ha (Strikt) | %.2f ha (Inclusief)", ha_strikt_ruw, ha_incl_ruw))
  
  # ----------------------------------------------------------------------------
  # 3. INTERNE HELPER: METAPOPULATIE-PIPELINE (>= 5 HA DIRECT, 1-5 HA MET >= 6 BUREN OP 500M)
  # ----------------------------------------------------------------------------
  filter_metapopulatie_veldkrekel <- function(r_bio_scenario, naam_scenario) {
    if (all(is.na(suppressWarnings(terra::minmax(r_bio_scenario))))) {
      return(list(leefgebied = master_grid * NA, grote_kernen = master_grid * NA))
    }
    
    # Stap 1: Clusteren op 50m, minimum 1 ha
    res_1ha <- cluster_filter_universeel(
      masker            = r_bio_scenario,
      opp_laag          = r_bio_scenario,
      min_opp_ha        = min_opp_ha,
      min_totale_opp_ha = NA,
      afstand_m         = zoekafstand_m,
      dispersie_m       = 0,
      werkelijk         = TRUE
    )
    r_1ha  <- res_1ha$raster
    cl_1ha <- res_1ha$clusters
    rm(res_1ha)
    gc()
    
    if (all(is.na(suppressWarnings(terra::minmax(cl_1ha))))) {
      return(list(leefgebied = master_grid * NA, grote_kernen = master_grid * NA))
    }
    
    # Stap 2: Oppervlakte per patch berekenen via terra::zonal
    df_opp <- terra::zonal(r_bio_scenario, cl_1ha, fun = "sum", na.rm = TRUE)
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
      
      # Centroïden bepalen per unieke patch via terra::patches
      p_poly <- terra::as.polygons(cl_klein, dissolve = TRUE)
      
      if (!is.null(p_poly) && nrow(p_poly) > 0) {
        pts_cent <- terra::centroids(p_poly, inside = TRUE)
        
        # Buffer 500m rondom elk punt en tel aantal snijdende centroïden
        buf_500m <- terra::buffer(pts_cent, width = 500)
        rel_mat  <- terra::is.related(buf_500m, pts_cent, "intersects")
        
        # Aantal buren per patch (inclusief zichzelf, dus drempel >= 6 buren = minstens 6 punten)
        pts_cent$buren_count <- rowSums(rel_mat)
        ok_pts <- pts_cent[pts_cent$buren_count >= 6, ]
        
        if (nrow(ok_pts) > 0) {
          # Haal de winnende ID's op
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
    
    # Samenvoegen grote en goedgekeurde kleine patches
    r_totaal_scenario <- terra::cover(r_groot, r_klein_ok)
    r_totaal_scenario <- terra::clamp(r_totaal_scenario, upper = 1.0)
    rm(r_klein_ok)
    gc()
    
    return(list(
      leefgebied  = r_totaal_scenario,
      grote_kernen = r_groot
    ))
  }
  
  # ----------------------------------------------------------------------------
  # 4. UITVOEREN PER SCENARIO
  # ----------------------------------------------------------------------------
  message("    🔬 Scenario 1: Strikt droog analyseren...")
  res_strikt <- filter_metapopulatie_veldkrekel(r_bio_strikt, "Strikt")
  r_lg_strikt <- res_strikt$leefgebied
  ha_lg_strikt <- calc_ha_exact(r_lg_strikt)
  message(sprintf("       -> Goedgekeurd leefgebied Strikt Droog: %.2f ha", ha_lg_strikt))
  
  message("    🔬 Scenario 2: Inclusief (suboptimaal/duin) analyseren...")
  res_incl   <- filter_metapopulatie_veldkrekel(r_bio_incl, "Inclusief")
  r_lg_incl  <- res_incl$leefgebied
  ha_lg_incl <- calc_ha_exact(r_lg_incl)
  message(sprintf("       -> Goedgekeurd leefgebied Inclusief: %.2f ha", ha_lg_incl))
  
  rm(r_bio_strikt, r_bio_incl, res_strikt, res_incl)
  gc()
  
  # ----------------------------------------------------------------------------
  # 5. GRADATIE SAMENSTELLEN: WAARDE 2 (STRIKT DROOG) & WAARDE 1 (SUBOPTIMAAL)
  # ----------------------------------------------------------------------------
  message("    🎨 Gradatieraster samenstellen (Klasse 2 = Zeker droog, Klasse 1 = Potentieel/Inclusief)...")
  
  r_bin_strikt <- terra::ifel(!is.na(r_lg_strikt) & r_lg_strikt > 0, 1, 0)
  r_bin_incl   <- terra::ifel(!is.na(r_lg_incl)   & r_lg_incl > 0,   1, 0)
  
  # Gradatie: als cel in strikt zit -> 2, als cel enkel in inclusief zit -> 1
  r_gradatie <- terra::ifel(r_bin_strikt == 1, 2, terra::ifel(r_bin_incl == 1, 1, NA))
  rm(r_bin_strikt, r_bin_incl)
  
  # Definitief leefgebied: unie van beide scenario's
  final_opp <- terra::cover(r_lg_strikt, r_lg_incl)
  final_opp <- terra::clamp(final_opp, upper = 1.0)
  
  ha_finaal <- calc_ha_exact(final_opp)
  message(sprintf("       -> Definitief totaal leefgebied Veldkrekel: %.2f ha", ha_finaal))
  
  # ----------------------------------------------------------------------------
  # 6. NETWERK CLUSTERING (CL_OPP)
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
    final_opp    = final_opp,
    cl_opp       = cl_opp,
    gradatie     = r_gradatie,
    strikt_droog = r_lg_strikt,
    inclusief    = r_lg_incl
  ))
}
