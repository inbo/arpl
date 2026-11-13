# ==============================================================================
# bouw_biotoop_mozaiek_met_substraat.R - LEEFGEBIEDENMODEL VOOR HEIVLINDER
# Archetype: Open heide/duin (>= 10 ha) x zandsubstraat (>= 0.5 ha) x 50m bosrand x 500m interactie
# Geschikt voor: Heivlinder, etc.
# ==============================================================================

bouw_biotoop_mozaiek_met_substraat <- function(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
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
  
  # ----------------------------------------------------------------------------
  # 2. OPEN BIOTOPEN COMBINEREN (Heide + Duin + Grasland + Groeve/Terril)
  # ----------------------------------------------------------------------------
  message("    🌾 Open biotooplagen (heide, duin, grasland, groeve/terril) combineren...")
  
  r_heide     <- if ("bwk_heide"         %in% names(lagen)) lagen[["bwk_heide"]]         else (master_grid * NA)
  r_duin      <- if ("bwk_duin"          %in% names(lagen)) lagen[["bwk_duin"]]          else (master_grid * NA)
  r_grasland  <- if ("bwk_grasland"      %in% names(lagen)) lagen[["bwk_grasland"]]      else (master_grid * NA)
  r_groeve    <- if ("bwk_groeve_terril" %in% names(lagen)) lagen[["bwk_groeve_terril"]] else (master_grid * NA)
  
  # Combinatie conform Bert (cover met clamp)
  r_open_opp <- terra::cover(r_heide, r_duin)
  r_open_opp <- terra::cover(r_open_opp, r_grasland)
  r_open_opp <- terra::cover(r_open_opp, r_groeve)
  r_open_opp <- terra::clamp(r_open_opp, upper = 1.0)
  
  rm(r_heide, r_duin, r_grasland, r_groeve)
  gc()
  
  ha_open_ruw <- calc_ha_exact(r_open_opp)
  message(sprintf("       -> Open biotoop totaal vóór clustering: %.2f ha", ha_open_ruw))
  
  message("    🏝️ Open biotoop clusteren (>= 10 ha, 50m netwerk)...")
  res_open <- cluster_filter_universeel(
    masker            = r_open_opp,
    opp_laag          = r_open_opp,
    min_opp_ha        = 10,
    min_totale_opp_ha = NA,
    afstand_m         = 50,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_open_goedgekeurd <- res_open$raster
  rm(res_open, r_open_opp)
  gc()
  
  ha_open_ok <- calc_ha_exact(r_open_goedgekeurd)
  message(sprintf("       -> Open biotoop goedgekeurd (>= 10 ha): %.2f ha", ha_open_ok))
  
  # ----------------------------------------------------------------------------
  # 3. ZANDSUBSTRAAT CLUSTEREN (>= 0.5 ha, 100m zoekafstand)
  # ----------------------------------------------------------------------------
  message("    🏖️ Zandsubstraat clusteren (>= 0.5 ha, 100m netwerk)...")
  r_zand_raw  <- if ("zand" %in% names(lagen)) lagen[["zand"]] else (master_grid * NA)
  ha_zand_ruw <- calc_ha_exact(r_zand_raw)
  message(sprintf("       -> Zandsubstraat vóór clustering: %.2f ha", ha_zand_ruw))
  
  res_zand <- cluster_filter_universeel(
    masker            = r_zand_raw,
    opp_laag          = r_zand_raw,
    min_opp_ha        = 0.5,
    min_totale_opp_ha = NA,
    afstand_m         = 100,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_zand_goedgekeurd <- res_zand$raster
  rm(res_zand, r_zand_raw)
  gc()
  
  ha_zand_ok <- calc_ha_exact(r_zand_goedgekeurd)
  message(sprintf("       -> Zandsubstraat goedgekeurd (>= 0.5 ha): %.2f ha", ha_zand_ok))
  
  # ----------------------------------------------------------------------------
  # 4. BESCHUTTING OPBOUWEN (BWK + STRIKT HOOGGROEN KLASSE 1) EN CLUSTEREN (>= 0.04 ha)
  # ----------------------------------------------------------------------------
  message("    🌳 Beschutting inlezen (Groenkaart strikt klasse 1 + BWK)...")
  f_groen_aligned <- hulp_align_gdal(paths$groenkaart, master_grid)
  r_groen_raw     <- terra::rast(f_groen_aligned)
  # Strikt klasse 1 = echte bomen/houtkanten in Groenkaart 2021
  r_groen_1       <- terra::ifel(!is.na(r_groen_raw) & r_groen_raw == 1, 1.0, NA)
  unlink(f_groen_aligned)
  rm(r_groen_raw)
  
  r_besch_bwk <- if ("bwk_beschutting" %in% names(lagen)) lagen[["bwk_beschutting"]] else (master_grid * NA)
  
  # Beschutting combineren
  r_besch_tot <- terra::cover(r_besch_bwk, r_groen_1)
  r_besch_tot <- terra::clamp(r_besch_tot, upper = 1.0)
  rm(r_groen_1, r_besch_bwk)
  gc()
  
  message("    🌲 Beschutting clusteren (>= 0.04 ha = 400m², 100m netwerk)...")
  res_besch <- cluster_filter_universeel(
    masker            = r_besch_tot,
    opp_laag          = r_besch_tot,
    min_opp_ha        = 0.04,
    min_totale_opp_ha = NA,
    afstand_m         = 100,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_besch_goedgekeurd <- res_besch$raster
  rm(res_besch, r_besch_tot)
  gc()
  
  # ----------------------------------------------------------------------------
  # 5. DIRECTE BOSRANDINTEGRATIE (50m randzone rond open biotoop)
  # ----------------------------------------------------------------------------
  message("    🎯 50m bosrand meenemen rondom goedgekeurd open biotoop...")
  r_open_bin <- terra::ifel(!is.na(r_open_goedgekeurd) & r_open_goedgekeurd > 0, 1, NA)
  
  if (!all(is.na(suppressWarnings(terra::minmax(r_open_bin))))) {
    f_open_bin  <- tempfile(fileext = ".tif")
    terra::writeRaster(r_open_bin, f_open_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    f_prox_open <- hulp_proximity_gdal(f_open_bin, 50)
    unlink(f_open_bin)
    
    p_open50    <- terra::rast(f_prox_open)
    mask_rand50 <- terra::ifel(!is.na(p_open50) & p_open50 >= 0 & p_open50 <= 50, 1, NA)
    unlink(f_prox_open)
    rm(p_open50)
    gc()
    
    bosrand_direct_opp <- terra::mask(r_besch_goedgekeurd, mask_rand50)
    rm(mask_rand50)
    gc()
  } else {
    bosrand_direct_opp <- master_grid * NA
  }
  rm(r_open_bin, r_besch_goedgekeurd)
  gc()
  
  # Foerageerbasis = Open Biotoop (>=10ha) + direct aangrenzende bosrand (<=50m)
  foerageer_base_opp     <- terra::cover(r_open_goedgekeurd, bosrand_direct_opp)
  voortplanting_base_opp <- r_zand_goedgekeurd
  
  ha_foer_base <- calc_ha_exact(foerageer_base_opp)
  ha_vpt_base  <- calc_ha_exact(voortplanting_base_opp)
  message(sprintf("       -> Foerageerbasis (open + 50m bosrand): %.2f ha", ha_foer_base))
  message(sprintf("       -> Voortplantingsbasis (zand >= 0.5 ha): %.2f ha", ha_vpt_base))
  
  # ----------------------------------------------------------------------------
  # 6. WEDERZIJDSE KOPPELING OP 500M (VOORTPLANTING X FOERAGEREN)
  # ----------------------------------------------------------------------------
  message("    ↔️ Wederzijdse afstandsfilter (500m) tussen zand en foerageerhabitat...")
  
  heeft_vpt  <- !all(is.na(suppressWarnings(terra::minmax(voortplanting_base_opp))))
  heeft_foer <- !all(is.na(suppressWarnings(terra::minmax(foerageer_base_opp))))
  
  if (heeft_vpt && heeft_foer) {
    r_vpt_bin  <- terra::ifel(!is.na(voortplanting_base_opp) & voortplanting_base_opp > 0, 1, NA)
    r_foer_bin <- terra::ifel(!is.na(foerageer_base_opp) & foerageer_base_opp > 0, 1, NA)
    
    # Zone 500m rond voortplanting (zand)
    f_vpt_bin   <- tempfile(fileext = ".tif")
    terra::writeRaster(r_vpt_bin, f_vpt_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    f_prox_vpt  <- hulp_proximity_gdal(f_vpt_bin, 500)
    unlink(f_vpt_bin)
    
    p_vpt       <- terra::rast(f_prox_vpt)
    mask_vpt500 <- terra::ifel(!is.na(p_vpt) & p_vpt >= 0 & p_vpt <= 500, 1, NA)
    unlink(f_prox_vpt)
    rm(p_vpt)
    
    # Zone 500m rond foerageren (heide/duin/rand)
    f_foer_bin   <- tempfile(fileext = ".tif")
    terra::writeRaster(r_foer_bin, f_foer_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    f_prox_foer  <- hulp_proximity_gdal(f_foer_bin, 500)
    unlink(f_foer_bin)
    
    p_foer       <- terra::rast(f_prox_foer)
    mask_foer500 <- terra::ifel(!is.na(p_foer) & p_foer >= 0 & p_foer <= 500, 1, NA)
    unlink(f_prox_foer)
    rm(p_foer)
    
    foerageer_finaal_opp     <- terra::mask(foerageer_base_opp, mask_vpt500)
    voortplanting_finaal_opp <- terra::mask(voortplanting_base_opp, mask_foer500)
    
    # Aparte deellagen netjes bijsnijden
    open_finaal_opp    <- terra::mask(r_open_goedgekeurd, mask_vpt500)
    bosrand_finaal_opp <- terra::mask(bosrand_direct_opp, mask_vpt500)
    
    rm(mask_vpt500, mask_foer500, r_vpt_bin, r_foer_bin)
    gc()
  } else {
    message("    ⚠️ Voortplanting of foerageerbasis is leeg vóór de wederzijdse toets!")
    foerageer_finaal_opp     <- master_grid * NA
    voortplanting_finaal_opp <- master_grid * NA
    open_finaal_opp          <- master_grid * NA
    bosrand_finaal_opp       <- master_grid * NA
  }
  rm(foerageer_base_opp, voortplanting_base_opp, r_open_goedgekeurd, bosrand_direct_opp, r_zand_goedgekeurd)
  gc()
  
  # ----------------------------------------------------------------------------
  # 7. SAMENVOEGEN TOT DEFINITIEF LEEFGEBIED & NETWERKEN
  # ----------------------------------------------------------------------------
  message("    🧩 Deellagen samenvoegen tot definitief leefgebied...")
  final_opp <- terra::cover(foerageer_finaal_opp, voortplanting_finaal_opp)
  final_opp <- terra::clamp(final_opp, upper = 1.0)
  
  ha_fin <- calc_ha_exact(final_opp)
  message(sprintf("       -> Definitief leefgebied Heivlinder: %.2f ha", ha_fin))
  
  r_fin_bin <- terra::ifel(!is.na(final_opp) & final_opp > 0, 1, NA)
  
  if (!all(is.na(suppressWarnings(terra::minmax(r_fin_bin))))) {
    f_bin   <- tempfile(fileext = ".tif")
    f_clump <- tempfile(fileext = ".tif")
    terra::writeRaster(r_fin_bin, f_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    
    e <- terra::ext(master_grid)
    qgisprocess::qgis_run_algorithm(
      "grass:r.clump",
      input                           = f_bin,
      output                          = f_clump,
      title                           = "cl_heivlinder",
      "-d"                            = FALSE,
      GRASS_REGION_PARAMETER          = sprintf("%f,%f,%f,%f", e$xmin, e$xmax, e$ymin, e$ymax),
      GRASS_REGION_CELLSIZE_PARAMETER = 10,
      .quiet                          = TRUE
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
    final_opp     = final_opp,
    cl_opp        = cl_opp,
    voortplanting = voortplanting_finaal_opp,
    open_biotoop  = open_finaal_opp,
    bosranden     = bosrand_finaal_opp
  ))
}
