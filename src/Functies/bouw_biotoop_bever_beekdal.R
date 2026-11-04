# ==============================================================================
# bouw_biotoop_bever_beekdal.R
# Archetype voor Europese bever: 50 ha beekdal/water x 5 ha houtige opslag (50m koppeling)
# Conform officiële INBO-standaard ANB2019_EuropeseBever2018.asc
# ==============================================================================

bouw_biotoop_bever_beekdal <- function(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
  bestandsnaam <- tolower(gsub(" ", "", soort_naam))
  soort_cache  <- file.path(paths$dir_cache, bestandsnaam)
  dir.create(soort_cache, showWarnings = FALSE, recursive = TRUE)
  
  f_cache_rds <- file.path(soort_cache, "01_raw_parquet_layers.rds")
  
  # Parameters ophalen
  rij_alg       <- params_soort %>% dplyr::filter(tolower(Type) == "leefgebied_algemeen") %>% dplyr::slice(1)
  rij_hout      <- params_soort %>% dplyr::filter(tolower(Type) == "houtig") %>% dplyr::slice(1)
  
  min_opp_alg   <- as.numeric(ifelse(nrow(rij_alg) > 0 && !is.na(rij_alg$MinOpp_ha), rij_alg$MinOpp_ha, 50))
  min_opp_hout  <- as.numeric(ifelse(nrow(rij_hout) > 0 && !is.na(rij_hout$MinOpp_ha), rij_hout$MinOpp_ha, 5))
  zoekafstand_m <- as.numeric(ifelse(nrow(rij_alg) > 0 && !is.na(rij_alg$AfstandBiotopen_m), rij_alg$AfstandBiotopen_m, 50))
  koppel_dist_m <- as.numeric(ifelse(nrow(rij_alg) > 0 && !is.na(rij_alg$Wederzijdse_Afstand_m), rij_alg$Wederzijdse_Afstand_m, 50))
  
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
  
  r_alg_bwk  <- if ("leefgebied_algemeen" %in% names(lagen)) lagen[["leefgebied_algemeen"]] else (master_grid * NA)
  r_hout_bwk <- if ("houtig"              %in% names(lagen)) lagen[["houtig"]]              else (master_grid * NA)
  
  # ----------------------------------------------------------------------------
  # 2. HYDROGRAFISCHE WATERLOPEN EN MAAS KOPPELEN MET SPOOR A
  # ----------------------------------------------------------------------------
  message("    🌊 Hydrografische waterlopen (Huetzon 1, 17, 18, 19) en Maas toevoegen aan Spoor A...")
  
  # A. Huet-zones klassen 1, 17, 18, 19
  if (file.exists(paths$huetzon)) {
    f_huet_al <- hulp_align_gdal(paths$huetzon, master_grid)
    r_huet    <- terra::rast(f_huet_al)
    r_huet_ok <- terra::ifel(!is.na(r_huet) & (r_huet %in% c(1, 17, 18, 19)), 1.0, NA)
    unlink(f_huet_al)
    rm(r_huet)
  } else {
    r_huet_ok <- master_grid * NA
  }
  
  # B. Optionele extra Maas-laag
  if (!is.null(paths$maas) && file.exists(paths$maas)) {
    f_maas_al <- hulp_align_gdal(paths$maas, master_grid)
    r_maas    <- terra::rast(f_maas_al)
    r_maas_ok <- terra::ifel(!is.na(r_maas) & r_maas >= 1, 1.0, NA)
    unlink(f_maas_al)
    rm(r_maas)
  } else {
    r_maas_ok <- master_grid * NA
  }
  
  # Spoor A samenstellen: BWK + Huetzones + Maas
  a_bwk <- terra::ifel(is.na(r_alg_bwk), 0, r_alg_bwk)
  a_hue <- terra::ifel(is.na(r_huet_ok), 0, r_huet_ok)
  a_mas <- terra::ifel(is.na(r_maas_ok), 0, r_maas_ok)
  
  r_alg_totaal <- terra::ifel((a_bwk + a_hue + a_mas) > 0, terra::clamp(a_bwk + a_hue + a_mas, upper = 1.0), NA)
  rm(a_bwk, a_hue, a_mas, r_alg_bwk, r_huet_ok, r_maas_ok)
  gc()
  
  ha_alg_ruw  <- calc_ha_exact(r_alg_totaal)
  ha_hout_ruw <- calc_ha_exact(r_hout_bwk)
  message(sprintf("       -> Algemeen beekdal-/watercomplex vóór clustering: %.2f ha", ha_alg_ruw))
  message(sprintf("       -> Houtige opslag vóór clustering: %.2f ha", ha_hout_ruw))
  
  # ----------------------------------------------------------------------------
  # 3. ZELFSTANDIGE CLUSTERING VAN BEIDE SPOREN
  # ----------------------------------------------------------------------------
  message(sprintf("    🏝️ Spoor A clusteren op >= %.1f ha (%dm zoekafstand)...", min_opp_alg, as.integer(zoekafstand_m)))
  res_alg <- cluster_filter_universeel(
    masker            = r_alg_totaal,
    opp_laag          = r_alg_totaal,
    min_opp_ha        = min_opp_alg,
    min_totale_opp_ha = NA,
    afstand_m         = zoekafstand_m,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_alg_50ha <- res_alg$raster
  rm(res_alg, r_alg_totaal)
  gc()
  
  message(sprintf("    🪵 Spoor B clusteren op >= %.1f ha (%dm zoekafstand)...", min_opp_hout, as.integer(zoekafstand_m)))
  res_hout <- cluster_filter_universeel(
    masker            = r_hout_bwk,
    opp_laag          = r_hout_bwk,
    min_opp_ha        = min_opp_hout,
    min_totale_opp_ha = NA,
    afstand_m         = zoekafstand_m,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_hout_5ha <- res_hout$raster
  rm(res_hout, r_hout_bwk)
  gc()
  
  ha_alg_cl  <- calc_ha_exact(r_alg_50ha)
  ha_hout_cl <- calc_ha_exact(r_hout_5ha)
  message(sprintf("       -> Geclusterd water-/dalcomplex (>= %.0f ha): %.2f ha", min_opp_alg, ha_alg_cl))
  message(sprintf("       -> Geclusterde houtige opslag (>= %.0f ha): %.2f ha", min_opp_hout, ha_hout_cl))
  
  # ----------------------------------------------------------------------------
  # 4. WEDERZIJDSE KOPPELING WATERCOMPLEX <-> HOUTIGE OPSLAG (<= 50 METER)
  # ----------------------------------------------------------------------------
  message(sprintf("    🔗 Wederzijdse oeverkoppeling berekenen (<= %d m via GDAL proximity)...", as.integer(koppel_dist_m)))
  
  heeft_alg  <- (ha_alg_cl > 0)
  heeft_hout <- (ha_hout_cl > 0)
  
  if (heeft_alg && heeft_hout) {
    # A. Proximity tot Hout berekenen (max. 50 m)
    r_h_bin <- terra::ifel(!is.na(r_hout_5ha) & r_hout_5ha > 0, 1, NA)
    f_h_bin <- tempfile(fileext = ".tif")
    terra::writeRaster(r_h_bin, f_h_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    rm(r_h_bin)
    
    f_prox_h <- hulp_proximity_gdal(f_h_bin, koppel_dist_m)
    unlink(f_h_bin)
    
    p_h <- terra::rast(f_prox_h)
    mask_alg_ok <- terra::ifel(!is.na(p_h) & p_h >= 0 & p_h <= koppel_dist_m, 1, NA)
    unlink(f_prox_h)
    rm(p_h)
    
    # Filter Spoor A: watercomplex binnen 50m van hout
    r_alg_gekoppeld <- terra::mask(r_alg_50ha, mask_alg_ok)
    rm(mask_alg_ok, r_alg_50ha)
    gc()
    
    # B. Proximity tot Watercomplex berekenen (max. 50 m)
    r_a_bin <- terra::ifel(!is.na(r_alg_gekoppeld) & r_alg_gekoppeld > 0, 1, NA)
    f_a_bin <- tempfile(fileext = ".tif")
    terra::writeRaster(r_a_bin, f_a_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    rm(r_a_bin)
    
    f_prox_a <- hulp_proximity_gdal(f_a_bin, koppel_dist_m)
    unlink(f_a_bin)
    
    p_a <- terra::rast(f_prox_a)
    mask_hout_ok <- terra::ifel(!is.na(p_a) & p_a >= 0 & p_a <= koppel_dist_m, 1, NA)
    unlink(f_prox_a)
    rm(p_a)
    
    # Filter Spoor B: houtige opslag binnen 50m van goedgekeurd watercomplex
    r_hout_gekoppeld <- terra::mask(r_hout_5ha, mask_hout_ok)
    rm(mask_hout_ok, r_hout_5ha)
    gc()
    
    # C. Samensmelten tot definitief beverleefgebied
    val_a <- terra::ifel(is.na(r_alg_gekoppeld), 0, r_alg_gekoppeld)
    val_h <- terra::ifel(is.na(r_hout_gekoppeld), 0, r_hout_gekoppeld)
    
    final_opp <- terra::ifel((val_a + val_h) > 0, terra::clamp(val_a + val_h, upper = 1.0), NA)
    rm(val_a, val_h)
    gc()
    
    # D. Unieke netwerk ID's via GRASS r.clump
    f_fin_bin <- tempfile(fileext = ".tif")
    f_cl_fin  <- tempfile(fileext = ".tif")
    terra::writeRaster(terra::ifel(!is.na(final_opp), 1, NA), f_fin_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    
    e <- terra::ext(master_grid)
    qgisprocess::qgis_run_algorithm(
      "grass:r.clump", input = f_fin_bin, output = f_cl_fin, title = "cl_bever",
      "-d" = TRUE, # directions = 8
      GRASS_REGION_PARAMETER = sprintf("%f,%f,%f,%f", e$xmin, e$xmax, e$ymin, e$ymax),
      GRASS_REGION_CELLSIZE_PARAMETER = 10, .quiet = TRUE
    )
    unlink(f_fin_bin)
    
    cl_opp <- terra::rast(f_cl_fin)
    terra::set.cats(cl_opp, layer = 1, value = NULL)
    if (terra::has.colors(cl_opp)) try(terra::coltab(cl_opp) <- NULL, silent = TRUE)
    cl_opp <- terra::mask(cl_opp, final_opp)
    unlink(f_cl_fin)
    
  } else {
    message("    ⚠️ Geen overlappend beekdal- of houtcomplex na koppeling.")
    final_opp <- master_grid * NA
    cl_opp    <- master_grid * NA
  }
  
  ha_finaal <- calc_ha_exact(final_opp)
  message(sprintf("       -> Definitief leefgebied Europese bever: %.2f ha", ha_finaal))
  
  return(list(
    final_opp = final_opp,
    cl_opp    = cl_opp
  ))
}
