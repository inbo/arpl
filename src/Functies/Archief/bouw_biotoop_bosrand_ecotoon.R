# ==============================================================================
# bouw_biotoop_bosrand_ecotoon.R
# Archetype: Warmteminnende ecotoonsoorten (bos >= 1ha x ruigte >= 1ha x contact 10m x rand >= 0.25ha x 100m buffer)
# Geschikt voor: Spaanse vlag (Euplagia quadripunctaria), etc.
# ==============================================================================

bouw_biotoop_bosrand_ecotoon <- function(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
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
  
  r_vp_raw  <- if ("voortplanting_bwk" %in% names(lagen)) lagen[["voortplanting_bwk"]] else (master_grid * NA)
  r_bos_raw <- if ("bos_bwk"           %in% names(lagen)) lagen[["bos_bwk"]]           else (master_grid * NA)
  
  # ----------------------------------------------------------------------------
  # 2. KERNDREMPELS VOOR BASISBIOTOPEN (MIN. 1 HA OP 50 M)
  # ----------------------------------------------------------------------------
  message(sprintf("    🌾 Voortplantingsbiotoop clusteren (>= %.1f ha, %dm netwerk)...", min_opp_ha, as.integer(zoekafstand_m)))
  res_vp <- cluster_filter_universeel(
    masker            = r_vp_raw,
    opp_laag          = r_vp_raw,
    min_opp_ha        = min_opp_ha,
    min_totale_opp_ha = NA,
    afstand_m         = zoekafstand_m,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_vp_1ha <- res_vp$raster
  rm(res_vp)
  gc()
  
  ha_vp_1ha <- calc_ha_exact(r_vp_1ha)
  message(sprintf("       -> Goedgekeurde voortplantingsruigtes (>= 1 ha): %.2f ha", ha_vp_1ha))
  
  message(sprintf("    🌳 Bosbiotoop clusteren (>= %.1f ha, %dm netwerk)...", min_opp_ha, as.integer(zoekafstand_m)))
  res_bos <- cluster_filter_universeel(
    masker            = r_bos_raw,
    opp_laag          = r_bos_raw,
    min_opp_ha        = min_opp_ha,
    min_totale_opp_ha = NA,
    afstand_m         = zoekafstand_m,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_bos_1ha <- res_bos$raster
  rm(res_bos)
  gc()
  
  ha_bos_1ha <- calc_ha_exact(r_bos_1ha)
  message(sprintf("       -> Goedgekeurde boskernen (>= 1 ha): %.2f ha", ha_bos_1ha))
  
  # ----------------------------------------------------------------------------
  # 3. CONTACTZONES IDENTIFICEREN (<= 10 M INTERACTIE)
  # ----------------------------------------------------------------------------
  message("    🔍 Scherpe contactzones (<= 10 m grensvlak) bepalen...")
  
  heeft_vp  <- !all(is.na(suppressWarnings(terra::minmax(r_vp_1ha))))
  heeft_bos <- !all(is.na(suppressWarnings(terra::minmax(r_bos_1ha))))
  
  if (heeft_vp && heeft_bos) {
    # Proximity vanaf bos naar ruigte (10 m)
    f_b_bin  <- tempfile(fileext = ".tif")
    terra::writeRaster(terra::ifel(!is.na(r_bos_1ha) & r_bos_1ha > 0, 1, NA), f_b_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    f_prox_b <- hulp_proximity_gdal(f_b_bin, 10)
    unlink(f_b_bin)
    
    p_b <- terra::rast(f_prox_b)
    mask_b10 <- terra::ifel(!is.na(p_b) & p_b >= 0 & p_b <= 10, 1, NA)
    unlink(f_prox_b)
    rm(p_b)
    
    r_rand_vp <- terra::mask(r_vp_1ha, mask_b10)
    rm(mask_b10)
    
    # Proximity vanaf ruigte naar bos (10 m)
    f_v_bin  <- tempfile(fileext = ".tif")
    terra::writeRaster(terra::ifel(!is.na(r_vp_1ha) & r_vp_1ha > 0, 1, NA), f_v_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    f_prox_v <- hulp_proximity_gdal(f_v_bin, 10)
    unlink(f_v_bin)
    
    p_v <- terra::rast(f_prox_v)
    mask_v10 <- terra::ifel(!is.na(p_v) & p_v >= 0 & p_v <= 10, 1, NA)
    unlink(f_prox_v)
    rm(p_v)
    
    r_rand_bos <- terra::mask(r_bos_1ha, mask_v10)
    rm(mask_v10)
    
    # Gecombineerde contactzone (ecotoon)
    r_rand_tot <- terra::cover(r_rand_vp, r_rand_bos)
    r_rand_tot <- terra::clamp(r_rand_tot, upper = 1.0)
    rm(r_rand_vp, r_rand_bos)
    gc()
  } else {
    r_rand_tot <- master_grid * NA
  }
  
  ha_rand_ruw <- calc_ha_exact(r_rand_tot)
  message(sprintf("       -> Totale contactzone vóór lengtetoets: %.2f ha", ha_rand_ruw))
  
  # ----------------------------------------------------------------------------
  # 4. FILTER OP AANEENGESLOTEN BOSRANDLENGTE (>= 0.25 HA / ~100-150 M)
  # ----------------------------------------------------------------------------
  message("    📏 Bosranden filteren op aaneengesloten lengte (>= 0.25 ha, 0m netwerk)...")
  
  res_rand_025 <- cluster_filter_universeel(
    masker            = r_rand_tot,
    opp_laag          = r_rand_tot,
    min_opp_ha        = 0.25,
    min_totale_opp_ha = NA,
    afstand_m         = 0,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_rand_gevalideerd <- res_rand_025$raster
  rm(res_rand_025, r_rand_tot)
  gc()
  
  ha_rand_val <- calc_ha_exact(r_rand_gevalideerd)
  message(sprintf("       -> Goedgekeurde bosrandecotonen (>= 0.25 ha): %.2f ha", ha_rand_val))
  
  # ----------------------------------------------------------------------------
  # 5. FUNCTIONEEL LEEFGEBIED VIA 100 M BOSRANDBUFFER
  # ----------------------------------------------------------------------------
  message("    🎯 Functioneel leefgebied afbakenen binnen 100 m van gevalideerde randen...")
  
  if (!all(is.na(suppressWarnings(terra::minmax(r_rand_gevalideerd))))) {
    f_r_bin  <- tempfile(fileext = ".tif")
    terra::writeRaster(terra::ifel(!is.na(r_rand_gevalideerd) & r_rand_gevalideerd > 0, 1, NA), f_r_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    f_prox_r100 <- hulp_proximity_gdal(f_r_bin, 100)
    unlink(f_r_bin)
    
    p_r100 <- terra::rast(f_prox_r100)
    mask_r100 <- terra::ifel(!is.na(p_r100) & p_r100 >= 0 & p_r100 <= 100, 1, NA)
    unlink(f_prox_r100)
    rm(p_r100)
    
    # Ruigte en bos binnen 100m van de goedgekeurde rand
    vp_finaal  <- terra::mask(r_vp_1ha, mask_r100)
    bos_finaal <- terra::mask(r_bos_1ha, mask_r100)
    rm(mask_r100)
    gc()
    
    # Totaal leefgebied = Ruigte + Bosrand zelf
    final_opp <- terra::cover(vp_finaal, r_rand_gevalideerd)
    final_opp <- terra::clamp(final_opp, upper = 1.0)
  } else {
    vp_finaal  <- master_grid * NA
    bos_finaal <- master_grid * NA
    final_opp  <- master_grid * NA
  }
  rm(r_vp_1ha, r_bos_1ha)
  gc()
  
  ha_finaal <- calc_ha_exact(final_opp)
  message(sprintf("       -> Definitief leefgebied Spaanse vlag: %.2f ha", ha_finaal))
  
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
    final_opp        = final_opp,
    cl_opp           = cl_opp,
    bosrand_ecotoon  = r_rand_gevalideerd,
    voortplanting    = vp_finaal,
    bos_buffer       = bos_finaal
  ))
}
