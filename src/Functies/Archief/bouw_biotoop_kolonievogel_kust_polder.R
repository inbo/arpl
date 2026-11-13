# ==============================================================================
# bouw_biotoop_kolonievogel_kust_polder.R
# Archetype voor kust-/kolonievogels (broedlocatie >= 30ha x water 5km x poldergrasland >= 200ha binnen 20km)
# Geschikt voor: Zwartkopmeeuw, etc.
# ==============================================================================

bouw_biotoop_kolonievogel_kust_polder <- function(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
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
    lagen <- readRDS(f_cache_rds)
  } else {
    message(paste0("    ⏳ Parquet-lagen extraheren voor ", soort_naam, "..."))
    lagen <- list()
    for (tp in types_in_tabel) {
      r_tp <- extract_bwk_laag(soort_naam, tp, ds_parquet, df_bwk_codes, master_grid)
      lagen[[tp]] <- if (!is.null(r_tp)) r_tp else (master_grid * NA)
    }
    saveRDS(lagen, f_cache_rds)
    message("    💾 [CACHE] Parquet-lagen veilig bewaard.")
  }
  
  r_vp_raw    <- if ("voortplanting1" %in% names(lagen)) lagen[["voortplanting1"]] else (master_grid * NA)
  r_water_raw <- if ("water1"         %in% names(lagen)) lagen[["water1"]]         else (master_grid * NA)
  r_foer_raw  <- if ("foerageer1"     %in% names(lagen)) lagen[["foerageer1"]]     else (master_grid * NA)
  
  # ----------------------------------------------------------------------------
  # 2. STAP 1: BASISCLUSTERING BROEDGEBIED (>= 30 HA) EN WATER (>= 2 HA)
  # ----------------------------------------------------------------------------
  message("    🏝️ Broedgebieden clusteren op >= 30 ha en waterlichamen op >= 2 ha (50m netwerk)...")
  
  # Broedgebied >= 30 ha
  res_vp <- cluster_filter_universeel(
    masker            = r_vp_raw,
    opp_laag          = r_vp_raw,
    min_opp_ha        = 30,
    min_totale_opp_ha = NA,
    afstand_m         = 50,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_vp_30ha <- res_vp$raster
  rm(res_vp, r_vp_raw)
  gc()
  
  ha_vp_30ha <- calc_ha_exact(r_vp_30ha)
  message(sprintf("       -> Goedgekeurd broedgebied (>= 30 ha): %.2f ha", ha_vp_30ha))
  
  # Water >= 2 ha
  res_water <- cluster_filter_universeel(
    masker            = r_water_raw,
    opp_laag          = r_water_raw,
    min_opp_ha        = 2,
    min_totale_opp_ha = NA,
    afstand_m         = 50,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_water_2ha <- res_water$raster
  rm(res_water, r_water_raw)
  gc()
  
  # ----------------------------------------------------------------------------
  # 3. STAP 2: WATERKOPPELING ROND BROEDLOCATIE (5 KM VIA GDAL PROXIMITY)
  # ----------------------------------------------------------------------------
  message("    🌊 5 km waterinteractie toetsen vanaf open water >= 2 ha...")
  
  heeft_vp    <- !all(is.na(suppressWarnings(terra::minmax(r_vp_30ha))))
  heeft_water <- !all(is.na(suppressWarnings(terra::minmax(r_water_2ha))))
  
  if (heeft_vp && heeft_water) {
    f_w2_bin <- tempfile(fileext = ".tif")
    terra::writeRaster(terra::ifel(!is.na(r_water_2ha) & r_water_2ha > 0, 1, NA), f_w2_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 0)
    f_prox_w5k <- hulp_proximity_gdal(f_w2_bin, 5000)
    unlink(f_w2_bin)
    
    p_w5k <- terra::rast(f_prox_w5k)
    mask_water5k <- terra::ifel(!is.na(p_w5k) & p_w5k >= 0 & p_w5k <= 5000, 1, NA)
    unlink(f_prox_w5k)
    rm(p_w5k)
    gc()
    
    # Broedgebied behouden dat binnen 5 km van water ligt
    r_vp_bij_water <- terra::mask(r_vp_30ha, mask_water5k)
    rm(mask_water5k)
  } else {
    r_vp_bij_water <- master_grid * NA
  }
  rm(r_vp_30ha)
  gc()
  
  ha_vp_bij_water <- calc_ha_exact(r_vp_bij_water)
  message(sprintf("       -> Broedgebied binnen 5 km van water: %.2f ha", ha_vp_bij_water))
  
  # ----------------------------------------------------------------------------
  # 4. STAP 3: FOERAGEERCOMPLEX (>= 200 HA) & ACTIERADIUSKOPPELING (20 KM)
  # ----------------------------------------------------------------------------
  message("    🌾 Foerageergebieden (weilanden) clusteren op >= 200 ha...")
  
  res_foer <- cluster_filter_universeel(
    masker            = r_foer_raw,
    opp_laag          = r_foer_raw,
    min_opp_ha        = 200,
    min_totale_opp_ha = NA,
    afstand_m         = 50,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_foer_200ha <- res_foer$raster
  rm(res_foer, r_foer_raw)
  gc()
  
  ha_foer_200ha <- calc_ha_exact(r_foer_200ha)
  message(sprintf("       -> Foerageercomplexen (>= 200 ha): %.2f ha", ha_foer_200ha))
  
  message("    🎯 20 km actieradiuskoppeling berekenen tussen broedgebied en foerageercomplex...")
  
  heeft_vp_w   <- !all(is.na(suppressWarnings(terra::minmax(r_vp_bij_water))))
  heeft_foer20 <- !all(is.na(suppressWarnings(terra::minmax(r_foer_200ha))))
  
  if (heeft_vp_w && heeft_foer20) {
    # 20 km buffer vanaf foerageercomplexen
    f_f200_bin <- tempfile(fileext = ".tif")
    terra::writeRaster(terra::ifel(!is.na(r_foer_200ha) & r_foer_200ha > 0, 1, NA), f_f200_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 0)
    f_prox_f20k <- hulp_proximity_gdal(f_f200_bin, 20000)
    unlink(f_f200_bin)
    
    p_f20k <- terra::rast(f_prox_f20k)
    mask_foer20k <- terra::ifel(!is.na(p_f20k) & p_f20k >= 0 & p_f20k <= 20000, 1, NA)
    unlink(f_prox_f20k)
    rm(p_f20k)
    gc()
    
    # Finaal goedgekeurd broedgebied
    vp_finaal <- terra::mask(r_vp_bij_water, mask_foer20k)
    rm(mask_foer20k)
    
    # Indien broedgebied goedgekeurd is, bepalen we ook de functioneel gekoppelde foerageer- en waterzones
    if (!all(is.na(suppressWarnings(terra::minmax(vp_finaal))))) {
      f_vp_fin_bin <- tempfile(fileext = ".tif")
      terra::writeRaster(terra::ifel(!is.na(vp_finaal) & vp_finaal > 0, 1, NA), f_vp_fin_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 0)
      
      # 20 km proximity vanaf broedgebied voor foerageergebied
      f_prox_vp20k <- hulp_proximity_gdal(f_vp_fin_bin, 20000)
      p_vp20k      <- terra::rast(f_prox_vp20k)
      foer_finaal  <- terra::mask(r_foer_200ha, terra::ifel(!is.na(p_vp20k) & p_vp20k <= 20000, 1, NA))
      unlink(f_prox_vp20k)
      rm(p_vp20k)
      
      # 5 km proximity vanaf broedgebied voor waterpartijen
      f_prox_vp5k <- hulp_proximity_gdal(f_vp_fin_bin, 5000)
      unlink(f_vp_fin_bin)
      p_vp5k      <- terra::rast(f_prox_vp5k)
      water_finaal <- terra::mask(r_water_2ha, terra::ifel(!is.na(p_vp5k) & p_vp5k <= 5000, 1, NA))
      unlink(f_prox_vp5k)
      rm(p_vp5k)
    } else {
      foer_finaal  <- master_grid * NA
      water_finaal <- master_grid * NA
    }
  } else {
    vp_finaal    <- master_grid * NA
    foer_finaal  <- master_grid * NA
    water_finaal <- master_grid * NA
  }
  rm(r_vp_bij_water, r_foer_200ha, r_water_2ha)
  gc()
  
  ha_vp_fin <- calc_ha_exact(vp_finaal)
  message(sprintf("       -> Finaal broedgebied binnen actieradius: %.2f ha", ha_vp_fin))
  
  # ----------------------------------------------------------------------------
  # 5. SAMENSTELLEN COMPOSIET LEEFGEBIED & NETWERK-ID'S
  # ----------------------------------------------------------------------------
  message("    🧩 Composiet leefgebied samenstellen (broedgebied + foerageerweiden + waterpartijen)...")
  
  final_opp <- terra::cover(vp_finaal, terra::cover(foer_finaal, water_finaal))
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
    
    cl_opp <- terra::rast(f_clump)
    terra::set.cats(cl_opp, layer = 1, value = NULL)
    if (terra::has.colors(cl_opp)) {
      try(terra::coltab(cl_opp) <- NULL, silent = TRUE)
    }
    cl_opp <- terra::mask(cl_opp, r_fin_bin)
    
    unlink(c(f_bin, f_clump))
  } else {
    cl_opp <- template_Vlaanderen
  }
  rm(r_fin_bin)
  gc()
  
  return(list(
    final_opp     = final_opp,
    cl_opp        = cl_opp,
    voortplanting = vp_finaal,
    foerageren    = foer_finaal,
    water         = water_finaal
  ))
}
