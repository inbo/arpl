# ==============================================================================
# bouw_biotoop_mozaiek_met_substraat.R - LEEFGEBIEDENMODEL VOOR HEIVLINDER (CORRECTE COMBINATIE)
# ==============================================================================

bouw_biotoop_mozaiek_met_substraat <- function(params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
  bestandsnaam <- "heivlinder"
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
    # Veilig inpakken vóór wegschrijven naar RDS
    lagen_packed <- lapply(lagen, terra::wrap)
    saveRDS(lagen_packed, f_cache_rds)
    message("    💾 [CACHE] Parquet-lagen veilig bewaard.")
  }
  
  # ----------------------------------------------------------------------------
  # 2. OPEN BIOTOPEN COMBINEREN (Heide + Duin + Grasland + Groeve/Terril)
  # ----------------------------------------------------------------------------
  message("    🌾 Open biotooplagen (heide, duin, grasland, groeve/terril) combineren...")
  
  # Haal de lagen veilig op; vervang NA door 0 voor de sommatie
  r_heide     <- if ("bwk_heide"         %in% names(lagen)) terra::ifel(is.na(lagen[["bwk_heide"]]),         0, lagen[["bwk_heide"]])         else (master_grid * 0)
  r_duin      <- if ("bwk_duin"          %in% names(lagen)) terra::ifel(is.na(lagen[["bwk_duin"]]),          0, lagen[["bwk_duin"]])          else (master_grid * 0)
  r_grasland  <- if ("bwk_grasland"      %in% names(lagen)) terra::ifel(is.na(lagen[["bwk_grasland"]]),      0, lagen[["bwk_grasland"]])      else (master_grid * 0)
  r_groeve    <- if ("bwk_groeve_terril" %in% names(lagen)) terra::ifel(is.na(lagen[["bwk_groeve_terril"]]), 0, lagen[["bwk_groeve_terril"]]) else (master_grid * 0)
  
  # Sommatie en clamp op 1.0 (exact zoals in brondocument)
  r_open_som <- r_heide + r_duin + r_grasland + r_groeve
  r_open_opp <- terra::ifel(r_open_som > 0, terra::clamp(r_open_som, upper = 1.0), NA)
  
  rm(r_heide, r_duin, r_grasland, r_groeve, r_open_som)
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
  ha_open_ok <- calc_ha_exact(r_open_goedgekeurd)
  message(sprintf("       -> Open biotoop goedgekeurd (>= 10 ha): %.2f ha", ha_open_ok))
  
  rm(res_open, r_open_opp)
  gc()
  
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
  ha_zand_ok <- calc_ha_exact(r_zand_goedgekeurd)
  message(sprintf("       -> Zandsubstraat goedgekeurd (>= 0.5 ha): %.2f ha", ha_zand_ok))
  
  rm(res_zand, r_zand_raw)
  gc()
  
  # ----------------------------------------------------------------------------
  # 4. BESCHUTTING OPBOUWEN (BWK + Groenkaart) EN CLUSTEREN (>= 0.04 ha, 100m)
  # ----------------------------------------------------------------------------
  message("    🌳 Beschutting inlezen (Groenkaart + BWK)...")
  f_groen_aligned <- hulp_align_gdal(paths$groenkaart, master_grid)
  r_groen_raw     <- terra::rast(f_groen_aligned)
  r_groen_0       <- terra::ifel(!is.na(r_groen_raw) & r_groen_raw > 0, 1.0, 0)
  
  r_besch_bwk     <- if ("bwk_beschutting" %in% names(lagen)) lagen[["bwk_beschutting"]] else (master_grid * NA)
  r_besch_0       <- terra::ifel(is.na(r_besch_bwk), 0, r_besch_bwk)
  
  r_besch_tot_som <- r_groen_0 + r_besch_0
  r_besch_tot     <- terra::ifel(r_besch_tot_som > 0, terra::clamp(r_besch_tot_som, upper = 1.0), NA)
  
  unlink(f_groen_aligned)
  rm(r_groen_raw, r_groen_0, r_besch_bwk, r_besch_0, r_besch_tot_som)
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
    terra::writeRaster(r_open_bin, f_open_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 0)
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
  w_open <- terra::ifel(is.na(r_open_goedgekeurd), 0, r_open_goedgekeurd)
  w_rand <- terra::ifel(is.na(bosrand_direct_opp), 0, bosrand_direct_opp)
  foerageer_base_opp <- terra::ifel((w_open + w_rand) > 0, terra::clamp(w_open + w_rand, upper = 1.0), NA)
  
  # Voortplantingsbasis = Goedgekeurd zandsubstraat (>= 0.5 ha)
  voortplanting_base_opp <- r_zand_goedgekeurd
  
  ha_foer_base <- calc_ha_exact(foerageer_base_opp)
  ha_vpt_base  <- calc_ha_exact(voortplanting_base_opp)
  message(sprintf("       -> Foerageerbasis (open + 50m bosrand): %.2f ha", ha_foer_base))
  message(sprintf("       -> Voortplantingsbasis (zand >= 0.5 ha): %.2f ha", ha_vpt_base))
  
  rm(w_open, w_rand, r_open_goedgekeurd, r_zand_goedgekeurd)
  gc()
  
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
    terra::writeRaster(r_vpt_bin, f_vpt_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 0)
    f_prox_vpt  <- hulp_proximity_gdal(f_vpt_bin, 500)
    unlink(f_vpt_bin)
    
    p_vpt       <- terra::rast(f_prox_vpt)
    mask_vpt500 <- terra::ifel(!is.na(p_vpt) & p_vpt >= 0 & p_vpt <= 500, 1, NA)
    unlink(f_prox_vpt)
    rm(p_vpt)
    
    # Zone 500m rond foerageren (heide/duin/rand)
    f_foer_bin   <- tempfile(fileext = ".tif")
    terra::writeRaster(r_foer_bin, f_foer_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 0)
    f_prox_foer  <- hulp_proximity_gdal(f_foer_bin, 500)
    unlink(f_foer_bin)
    
    p_foer       <- terra::rast(f_prox_foer)
    mask_foer500 <- terra::ifel(!is.na(p_foer) & p_foer >= 0 & p_foer <= 500, 1, NA)
    unlink(f_prox_foer)
    rm(p_foer)
    
    foerageer_finaal_opp     <- terra::mask(foerageer_base_opp, mask_vpt500)
    voortplanting_finaal_opp <- terra::mask(voortplanting_base_opp, mask_foer500)
    
    unlink(c(mask_vpt500, mask_foer500))
    rm(r_vpt_bin, r_foer_bin)
    gc()
  } else {
    message("    ⚠️ Voortplanting of foerageerbasis is leeg vóór de wederzijdse toets!")
    foerageer_finaal_opp     <- master_grid * NA
    voortplanting_finaal_opp <- master_grid * NA
  }
  rm(foerageer_base_opp, voortplanting_base_opp)
  gc()
  
  # ----------------------------------------------------------------------------
  # 7. SAMENVOEGEN TOT DEFINITIEF LEEFGEBIED & NETWERKEN
  # ----------------------------------------------------------------------------
  message("    🧩 Deellagen samenvoegen tot definitief leefgebied...")
  w_f <- terra::ifel(is.na(foerageer_finaal_opp), 0, foerageer_finaal_opp)
  w_v <- terra::ifel(is.na(voortplanting_finaal_opp), 0, voortplanting_finaal_opp)
  
  final_opp <- terra::ifel((w_f + w_v) > 0, terra::clamp(w_f + w_v, upper = 1.0), NA)
  rm(w_f, w_v)
  gc()
  
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
    
    cl_opp <- terra::rast(f_clump)
    terra::set.cats(cl_opp, layer = 1, value = NULL)
    terra::coltab(cl_opp) <- NULL
    cl_opp <- terra::mask(cl_opp, r_fin_bin)
    
    unlink(c(f_bin, f_clump))
  } else {
    cl_opp <- master_grid * NA
  }
  rm(r_fin_bin)
  gc()
  
  return(list(
    final_opp        = final_opp,
    cl_opp           = cl_opp,
    voortplanting    = voortplanting_finaal_opp,
    foerageren       = foerageer_finaal_opp,
    bosranden_meegen = bosrand_direct_opp
  ))
}
