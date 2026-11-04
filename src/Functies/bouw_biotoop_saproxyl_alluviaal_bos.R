# ==============================================================================
# bouw_biotoop_saproxyl_alluviaal_bos.R
# Archetype voor alluviale zachthoutkevers (populieren/ooibos + Groenkaart klasse 1 x wetlands/waterloop 100m x >= 1 ha op 50m)
# Geschikt voor: Vermiljoenkever (Cucujus cinnaberinus), etc.
# ==============================================================================

bouw_biotoop_saproxyl_alluviaal_bos <- function(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
  bestandsnaam <- tolower(gsub(" ", "", soort_naam))
  soort_cache  <- file.path(paths$dir_cache, bestandsnaam)
  dir.create(soort_cache, showWarnings = FALSE, recursive = TRUE)
  
  f_cache_rds <- file.path(soort_cache, "01_raw_parquet_layers.rds")
  
  # Parameters ophalen uit Excel
  rij_params    <- params_soort %>% dplyr::slice(1)
  min_opp_ha    <- as.numeric(ifelse(is.na(rij_params$MinOpp_ha), 1.0, rij_params$MinOpp_ha))
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
  
  r_bwk_populier <- if ("bwk" %in% names(lagen)) lagen[["bwk"]] else (master_grid * NA)
  
  # ----------------------------------------------------------------------------
  # 2. HOUTIGE OPSTANDEN (BWK POPULIEREN/OOIBOS + GROENKAART HOOGGROEN)
  # ----------------------------------------------------------------------------
  message("    🌳 Houtige opstanden samenstellen (BWK zachthout + Groenkaart klasse 1)...")
  
  f_groen_al <- hulp_align_gdal(paths$groenkaart, master_grid)
  r_groen    <- terra::rast(f_groen_al)
  # Strikt klasse 1 = bomen/opgaand hout
  r_groen_bomen <- terra::ifel(!is.na(r_groen) & r_groen == 1, 1.0, NA)
  unlink(f_groen_al)
  rm(r_groen)
  
  # Combinatie van populierenbos en bomenrijen
  r_hout_tot <- terra::cover(r_bwk_populier, r_groen_bomen)
  r_hout_tot <- terra::clamp(r_hout_tot, upper = 1.0)
  rm(r_groen_bomen)
  gc()
  
  ha_hout <- calc_ha_exact(r_hout_tot)
  message(sprintf("       -> Totaal bomenareaal vóór wetlandfiltering: %.2f ha", ha_hout))
  
  # ----------------------------------------------------------------------------
  # 3. MILIEUMASKER: ALLUVIALE WETLANDS / BEEKBEGELEIDENDE NATTE GRONDEN
  # ----------------------------------------------------------------------------
  message("    🌊 Alluviale wetlandconditie toepassen (valleigronden & natte beekoevers)...")
  
  f_cache_wetland <- file.path(paths$dir_cache, "masker_alluviaal_wetland.tif")
  
  if (file.exists(f_cache_wetland)) {
    message("       -> [CACHE] Wetlandmasker ingeladen van schijf...")
    mask_wetland <- terra::rast(f_cache_wetland)
  } else {
    # Strategie A: Bestaat er een expliciete wetlandkaart (actwet)?
    heeft_wetland_file <- !is.null(paths$wetland) && file.exists(paths$wetland)
    
    if (heeft_wetland_file) {
      f_wet_al     <- hulp_align_gdal(paths$wetland, master_grid)
      r_wet_raw    <- terra::rast(f_wet_al)
      # Typische wetlandcodes (50 t/m 81, 120 of > 0)
      mask_wetland <- terra::ifel(!is.na(r_wet_raw) & (r_wet_raw %in% c(50:81, 120) | r_wet_raw > 0), 1, NA)
      unlink(f_wet_al)
      rm(r_wet_raw)
    } else {
      # Strategie B: Robuuste fallback conform Rust-script (Waterloop 100m buffer OF natte drainage e-i)
      message("       -> Geen 'actwet.tif' gevonden; fallback via 100m waterloopbuffer en natte drainage...")
      
      # 1. 100m buffer rond waterlopen (Huetzones of watervlakken)
      mask_water_100m <- master_grid * NA
      if (!is.null(paths$huetzon) && file.exists(paths$huetzon)) {
        f_hz_al  <- hulp_align_gdal(paths$huetzon, master_grid)
        r_hz_raw <- terra::rast(f_hz_al)
        r_hz_bin <- terra::ifel(!is.na(r_hz_raw) & r_hz_raw > 0, 1, NA)
        unlink(f_hz_al)
        rm(r_hz_raw)
        
        if (!all(is.na(suppressWarnings(terra::minmax(r_hz_bin))))) {
          f_hz_tmp <- tempfile(fileext = ".tif")
          terra::writeRaster(r_hz_bin, f_hz_tmp, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
          f_prox_hz <- hulp_proximity_gdal(f_hz_tmp, 100)
          unlink(f_hz_tmp)
          
          p_hz <- terra::rast(f_prox_hz)
          mask_water_100m <- terra::ifel(!is.na(p_hz) & p_hz >= 0 & p_hz <= 100, 1, NA)
          unlink(f_prox_hz)
          rm(p_hz)
        }
        rm(r_hz_bin)
      }
      
      # 2. Natte drainageklassen (e en natter: e, f, g, h, i, e-f, h-i, e-i)
      mask_drain_nat <- master_grid * NA
      if (!is.null(paths$drainage) && file.exists(paths$drainage)) {
        f_dr_al  <- hulp_align_gdal(paths$drainage, master_grid)
        r_dr_raw <- terra::rast(f_dr_al)
        leg_dr   <- terra::cats(r_dr_raw)[[1]]
        
        nat_codes <- c("e", "f", "g", "h", "i", "e-f", "h-i", "e-i")
        ok_d_ids  <- integer(0)
        if (!is.null(leg_dr)) {
          kol_l <- colnames(leg_dr)[tolower(colnames(leg_dr)) %in% c("label", "omschrijving", "drainage")][1]
          kol_v <- colnames(leg_dr)[tolower(colnames(leg_dr)) %in% c("value", "id")][1]
          if (is.na(kol_v)) kol_v <- 1
          if (!is.na(kol_l)) {
            ok_d_ids <- leg_dr[[kol_v]][tolower(trimws(as.character(leg_dr[[kol_l]]))) %in% nat_codes]
          }
        }
        if (length(ok_d_ids) == 0) ok_d_ids <- c(5, 6, 7, 8, 9, 56, 89)
        mask_drain_nat <- terra::ifel(!is.na(r_dr_raw) & (r_dr_raw %in% ok_d_ids), 1, NA)
        unlink(f_dr_al)
        rm(r_dr_raw)
      }
      
      mask_wetland <- terra::cover(mask_water_100m, mask_drain_nat)
      rm(mask_water_100m, mask_drain_nat)
    }
    
    terra::writeRaster(mask_wetland, f_cache_wetland, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
  }
  
  # Bomen maskeren met de alluviale wetlandzone
  r_bos_alluviaal <- terra::mask(r_hout_tot, mask_wetland)
  rm(r_hout_tot, mask_wetland)
  gc()
  
  ha_alluviaal <- calc_ha_exact(r_bos_alluviaal)
  message(sprintf("       -> Houtige opstanden binnen alluviaal wetlandgebied: %.2f ha", ha_alluviaal))
  
  # ----------------------------------------------------------------------------
  # 4. RUIMTELIJKE CLUSTERING & OPPERVLAKTEDREMPEL (>= 1.0 HA OP 50 M)
  # ----------------------------------------------------------------------------
  message(sprintf("    🏝️ Valleibossen clusteren (%dm netwerk) en filteren op >= %.2f ha...", as.integer(zoekafstand_m), min_opp_ha))
  
  res_cluster <- cluster_filter_universeel(
    masker            = r_bos_alluviaal,
    opp_laag          = r_bos_alluviaal,
    min_opp_ha        = min_opp_ha,
    min_totale_opp_ha = NA,
    afstand_m         = zoekafstand_m,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  
  final_opp <- res_cluster$raster
  cl_opp    <- res_cluster$clusters
  rm(res_cluster, r_bos_alluviaal)
  gc()
  
  ha_finaal <- calc_ha_exact(final_opp)
  message(sprintf("       -> Definitief leefgebied Vermiljoenkever: %.2f ha", ha_finaal))
  
  # Deellagen isoleren binnen goedgekeurde netwerken
  if (!all(is.na(suppressWarnings(terra::minmax(cl_opp))))) {
    populier_uitvoer <- terra::mask(r_bwk_populier, cl_opp)
  } else {
    populier_uitvoer <- master_grid * NA
  }
  rm(r_bwk_populier)
  gc()
  
  return(list(
    final_opp               = final_opp,
    cl_opp                  = cl_opp,
    bwk_populier_ooibos     = populier_uitvoer
  ))
}
