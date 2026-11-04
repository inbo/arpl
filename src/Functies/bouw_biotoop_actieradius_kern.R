# ==============================================================================
# bouw_biotoop_actieradius_kern.R - KERNBIOTOOP MET MORFOLOGIE & ACTIERADIUS
# Archetype voor soorten met brede broed/rustkernen en uitgestrekt jachtgebied
# ==============================================================================

bouw_biotoop_actieradius_kern <- function(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
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
    # Veilig inpakken vóór wegschrijven naar RDS
    lagen_packed <- lapply(lagen, terra::wrap)
    saveRDS(lagen_packed, f_cache_rds)
    message("    💾 [CACHE] Parquet-lagen veilig bewaard.")
  }
  
  r_vp_bwk   <- if ("voortplanting_bwk" %in% names(lagen)) lagen[["voortplanting_bwk"]] else (master_grid * NA)
  r_foer_raw <- if ("foerageer"         %in% names(lagen)) lagen[["foerageer"]]         else (master_grid * NA)
  
  # ----------------------------------------------------------------------------
  # 2. VOORTPLANTING KOPPELEN MET LANDBOUWGEWASSEN (LGP)
  # ----------------------------------------------------------------------------
  if (file.exists(paths$lgp_landbouw)) {
    message("    🌾 Landbouwgebruikspercelen (LGP graangewassen 311, 321, 331) inlezen...")
    f_lgp_aligned <- hulp_align_gdal(paths$lgp_landbouw, master_grid)
    r_lgp         <- terra::rast(f_lgp_aligned)
    
    # Kiekendief gewascodes: wintertarwe, wintergerst, spelt
    kiek_codes <- c(311, 321, 331)
    r_lgp_kiek <- terra::ifel(!is.na(r_lgp) & r_lgp %in% kiek_codes, 1.0, NA)
    unlink(f_lgp_aligned)
    rm(r_lgp)
    gc()
    
    # Cover LGP over BWK voortplanting
    r_vp_totaal <- terra::cover(r_lgp_kiek, r_vp_bwk)
    rm(r_lgp_kiek)
    gc()
  } else {
    message("    ⚠️ Geen LGP bestand gevonden op paths$lgp_landbouw. Alleen BWK gebruikt.")
    r_vp_totaal <- r_vp_bwk
  }
  
  ha_vp_ruw <- calc_ha_exact(r_vp_totaal)
  message(sprintf("       -> Voortplantingsbiotoop vóór filtering: %.2f ha", ha_vp_ruw))
  
  # ----------------------------------------------------------------------------
  # 3. MORFOLOGISCHE ZUIVERING (OPENING STRAAL 20M: MIN. 40M BREEDTE)
  # ----------------------------------------------------------------------------
  message("    🔬 Morfologische zuivering (min. 40m perceelbreedte / eliminatie smalle rietkragen)...")
  r_vp_bin <- terra::ifel(!is.na(r_vp_totaal) & r_vp_totaal > 0, 1, NA)
  
  if (!all(is.na(suppressWarnings(terra::minmax(r_vp_bin))))) {
    # Snelle raster-opening via $5\times 5$ focal window (2 cellen = 20m) op schijf
    # Stap A: Erosie (alle buren binnen 20m moeten 1 zijn)
    w_20m <- matrix(1, nrow = 5, ncol = 5)
    f_eroded <- tempfile(fileext = ".tif")
    r_eroded <- terra::focal(
      r_vp_bin, 
      w = w_20m, 
      fun = "min", 
      na.rm = FALSE, 
      filename = f_eroded, 
      overwrite = TRUE, 
      gdal = c("COMPRESS=LZW")
    )
    r_eroded_bin <- terra::ifel(!is.na(r_eroded) & r_eroded == 1, 1, NA)
    
    # Stap B: Dilatatie (terug uitbreiden vanaf de behouden kernen)
    f_dilated <- tempfile(fileext = ".tif")
    r_dilated <- terra::focal(
      r_eroded_bin, 
      w = w_20m, 
      fun = "max", 
      na.rm = TRUE, 
      filename = f_dilated, 
      overwrite = TRUE, 
      gdal = c("COMPRESS=LZW")
    )
    mask_breed <- terra::ifel(!is.na(r_dilated) & r_dilated == 1, 1, NA)
    
    # Behoud enkel werkelijke habitatoppervlaktes binnen de brede kernen
    r_vp_gezuiverd <- terra::mask(r_vp_totaal, mask_breed)
    
    unlink(c(f_eroded, f_dilated))
    rm(r_eroded, r_eroded_bin, r_dilated, mask_breed)
    gc()
  } else {
    r_vp_gezuiverd <- master_grid * NA
  }
  rm(r_vp_bin)
  gc()
  
  ha_vp_gezuiverd <- calc_ha_exact(r_vp_gezuiverd)
  message(sprintf("       -> Voortplantingsbiotoop na breedtefilter (>= 40m): %.2f ha", ha_vp_gezuiverd))
  
  # ----------------------------------------------------------------------------
  # 4. CLUSTERING VOORTPLANTING (>= 10 HA, 100M ZOEFAFSTAND)
  # ----------------------------------------------------------------------------
  message("    🏝️ Voortplanting clusteren (>= 10 ha, 100m netwerk)...")
  res_vp <- cluster_filter_universeel(
    masker            = r_vp_gezuiverd,
    opp_laag          = r_vp_gezuiverd,
    min_opp_ha        = 10,
    min_totale_opp_ha = NA,
    afstand_m         = 100,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_vp_10ha <- res_vp$raster
  ha_vp_10ha <- calc_ha_exact(r_vp_10ha)
  message(sprintf("       -> Voortplantingskernen goedgekeurd (>= 10 ha): %.2f ha", ha_vp_10ha))
  rm(res_vp, r_vp_gezuiverd)
  gc()
  
  # ----------------------------------------------------------------------------
  # 5. CLUSTERING FOERAGEERGEBIED (>= 100 HA, 100M ZOEFAFSTAND)
  # ----------------------------------------------------------------------------
  message("    🦅 Foerageergebied clusteren (>= 100 ha, 100m netwerk)...")
  res_foer <- cluster_filter_universeel(
    masker            = r_foer_raw,
    opp_laag          = r_foer_raw,
    min_opp_ha        = 100,
    min_totale_opp_ha = NA,
    afstand_m         = 100,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_foer_100ha <- res_foer$raster
  ha_foer_100ha <- calc_ha_exact(r_foer_100ha)
  message(sprintf("       -> Foerageercomplexen goedgekeurd (>= 100 ha): %.2f ha", ha_foer_100ha))
  rm(res_foer, r_foer_raw)
  gc()
  
  # ----------------------------------------------------------------------------
  # 6. ACTIERADIUSTOETS (5 KM): VOORTPLANTING MOET BINNEN 5 KM VAN FOERAGEERAREAAL LIGGEN
  # ----------------------------------------------------------------------------
  message("    🎯 Omgevingstoetsing: Voortplanting binnen <= 5 km van groot foerageercomplex...")
  
  heeft_vp   <- !all(is.na(suppressWarnings(terra::minmax(r_vp_10ha))))
  heeft_foer <- !all(is.na(suppressWarnings(terra::minmax(r_foer_100ha))))
  
  if (heeft_vp && heeft_foer) {
    r_foer_bin <- terra::ifel(!is.na(r_foer_100ha) & r_foer_100ha > 0, 1, NA)
    
    f_foer_bin  <- tempfile(fileext = ".tif")
    terra::writeRaster(r_foer_bin, f_foer_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 0)
    
    # 5 km proximity vanaf goedgekeurd foerageergebied
    f_prox_5km  <- hulp_proximity_gdal(f_foer_bin, 5000)
    unlink(f_foer_bin)
    rm(r_foer_bin)
    gc()
    
    p_5km       <- terra::rast(f_prox_5km)
    mask_foer_5km <- terra::ifel(!is.na(p_5km) & p_5km >= 0 & p_5km <= 5000, 1, NA)
    unlink(f_prox_5km)
    rm(p_5km)
    gc()
    
    # Voortplanting maskeren: alleen behouden als het binnen 5 km van foerageercomplex ligt
    vp_finaal <- terra::mask(r_vp_10ha, mask_foer_5km)
    rm(mask_foer_5km)
    gc()
  } else {
    message("    ⚠️ Voortplanting of foerageercomplexen ontbreken vóór de 5km-toets!")
    vp_finaal <- master_grid * NA
  }
  rm(r_vp_10ha)
  gc()
  
  ha_vp_finaal <- calc_ha_exact(vp_finaal)
  message(sprintf("       -> Voortplanting finaal goedgekeurd binnen 5 km: %.2f ha", ha_vp_finaal))
  
  # ----------------------------------------------------------------------------
  # 7. SAMENVOEGEN TOT TOTAAL LEEFGEBIED (VOORTPLANTING + FOERAGEREN) & NETWERK-ID'S
  # ----------------------------------------------------------------------------
  message("    🧩 Leefgebied samenvoegen (Voortplanting + Foerageerareaal)...")
  
  # Cover: waar beide aanwezig zijn telt voortplanting eerst, foerageren vult aan
  final_opp <- terra::cover(vp_finaal, r_foer_100ha)
  
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
      title                           = paste0("cl_", bestandsnaam),
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
    final_opp     = final_opp,
    cl_opp        = cl_opp,
    voortplanting = vp_finaal,
    foerageren    = r_foer_100ha
  ))
}
