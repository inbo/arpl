# ==============================================================================
# bouw_biotoop_actieradius_kern.R - KERNBIOTOOP MET MORFOLOGIE & ACTIERADIUS
# Archetype voor soorten met brede broed/rustkernen en uitgestrekt jachtgebied
# Geschikt voor: Bruine kiekendief, Grauwe kiekendief, etc.
# ==============================================================================

bouw_biotoop_actieradius_kern <- function(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
  bestandsnaam <- tolower(gsub(" ", "", soort_naam))
  soort_cache  <- file.path(paths$dir_cache, bestandsnaam)
  dir.create(soort_cache, showWarnings = FALSE, recursive = TRUE)
  
  # Parameters ophalen
  rij_vp   <- params_soort %>% dplyr::filter(grepl("voortplanting", tolower(Type))) %>% dplyr::slice(1)
  rij_foer <- params_soort %>% dplyr::filter(grepl("foerageer", tolower(Type))) %>% dplyr::slice(1)
  
  min_opp_vp    <- as.numeric(ifelse(nrow(rij_vp) > 0 && !is.na(rij_vp$MinOpp_ha), rij_vp$MinOpp_ha, 10.0))
  min_opp_foer  <- as.numeric(ifelse(nrow(rij_foer) > 0 && !is.na(rij_foer$MinOpp_ha), rij_foer$MinOpp_ha, 100.0))
  zoekafstand_m <- as.numeric(ifelse(nrow(rij_vp) > 0 && !is.na(rij_vp$AfstandBiotopen_m), rij_vp$AfstandBiotopen_m, 100.0))
  actieradius_m <- 5000.0 # 5 km actieradius conform Kiekendief-model
  
  # ----------------------------------------------------------------------------
  # 1. PARQUET EXTRACTIE (via modulaire GeoTIFF cache)
  # ----------------------------------------------------------------------------
  types_te_laden <- unique(c("voortplanting_bwk", "foerageer", params_soort$Type))
  lagen <- laad_of_extraheer_parquet_lagen(soort_naam, types_te_laden, ds_parquet, df_bwk_codes, master_grid, soort_cache)
  
  r_vp_bwk   <- if ("voortplanting_bwk" %in% names(lagen)) lagen[["voortplanting_bwk"]] else (master_grid * NA_real_)
  r_foer_raw <- if ("foerageer"         %in% names(lagen)) lagen[["foerageer"]]         else (master_grid * NA_real_)
  
  # ----------------------------------------------------------------------------
  # 2. VOORTPLANTING KOPPELEN MET LANDBOUWGEWASSEN (LGP)
  # ----------------------------------------------------------------------------
  # Kiekendief gewascodes: wintertarwe (311), wintergerst (321), spelt (331)
  r_lgp_kiek <- haal_categorisch_masker(
    pad_raster     = paths$lgp_landbouw, 
    codes_selectie = c(311, 321, 331), 
    master_grid    = master_grid, 
    dir_cache      = paths$dir_cache, 
    naam_prefix    = "lgp_graan"
  )
  
  r_vp_totaal <- terra::clamp(terra::cover(r_lgp_kiek, r_vp_bwk), upper = 1.0)
  rm(r_lgp_kiek, r_vp_bwk)
  gc()
  
  ha_vp_ruw <- calc_ha_exact(r_vp_totaal)
  message(sprintf("       -> Voortplantingsbiotoop vóór filtering: %.2f ha", ha_vp_ruw))
  
  # ----------------------------------------------------------------------------
  # 3. MORFOLOGISCHE ZUIVERING (OPENING: MIN. 40M BREEDTE)
  # ----------------------------------------------------------------------------
  message("    🔬 Morfologische zuivering (min. 40m perceelbreedte / eliminatie smalle linten)...")
  r_vp_gezuiverd <- hulp_morfologische_opening(r_vp_totaal, diameter_m = 40, master_grid = master_grid)
  rm(r_vp_totaal)
  gc()
  
  ha_vp_gezuiverd <- calc_ha_exact(r_vp_gezuiverd)
  message(sprintf("       -> Voortplantingsbiotoop na breedtefilter (>= 40m): %.2f ha", ha_vp_gezuiverd))
  
  # ----------------------------------------------------------------------------
  # 4. CLUSTERING VOORTPLANTING (>= 10 HA, 100M ZOEFAFSTAND)
  # ----------------------------------------------------------------------------
  message(sprintf("    🏝️ Voortplanting clusteren (>= %.1f ha, %dm netwerk)...", min_opp_vp, as.integer(zoekafstand_m)))
  res_vp <- cluster_filter_universeel(
    masker            = r_vp_gezuiverd,
    opp_laag          = r_vp_gezuiverd,
    min_opp_ha        = min_opp_vp,
    min_totale_opp_ha = NA,
    afstand_m         = zoekafstand_m,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_vp_10ha  <- res_vp$raster
  rm(res_vp, r_vp_gezuiverd)
  gc()
  
  ha_vp_10ha <- calc_ha_exact(r_vp_10ha)
  message(sprintf("       -> Voortplantingskernen goedgekeurd (>= %.0f ha): %.2f ha", min_opp_vp, ha_vp_10ha))
  
  # ----------------------------------------------------------------------------
  # 5. CLUSTERING FOERAGEERGEBIED (>= 100 HA, 100M ZOEFAFSTAND)
  # ----------------------------------------------------------------------------
  message(sprintf("    🦅 Foerageergebied clusteren (>= %.1f ha, %dm netwerk)...", min_opp_foer, as.integer(zoekafstand_m)))
  res_foer <- cluster_filter_universeel(
    masker            = r_foer_raw,
    opp_laag          = r_foer_raw,
    min_opp_ha        = min_opp_foer,
    min_totale_opp_ha = NA,
    afstand_m         = zoekafstand_m,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_foer_100ha  <- res_foer$raster
  rm(res_foer, r_foer_raw)
  gc()
  
  ha_foer_100ha <- calc_ha_exact(r_foer_100ha)
  message(sprintf("       -> Foerageercomplexen goedgekeurd (>= %.0f ha): %.2f ha", min_opp_foer, ha_foer_100ha))
  
  # ----------------------------------------------------------------------------
  # 6. ACTIERADIUSTOETS (<= 5 KM VIA GDAL PROXIMITY)
  # ----------------------------------------------------------------------------
  message("    🎯 Omgevingstoetsing: Voortplanting binnen <= 5 km van foerageercomplex...")
  
  heeft_vp   <- !all(is.na(suppressWarnings(terra::minmax(r_vp_10ha))))
  heeft_foer <- !all(is.na(suppressWarnings(terra::minmax(r_foer_100ha))))
  
  if (heeft_vp && heeft_foer) {
    f_foer_bin <- tempfile(fileext = ".tif")
    terra::writeRaster(terra::ifel(!is.na(r_foer_100ha) & r_foer_100ha > 0, 1, NA), f_foer_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    
    f_prox_5km <- hulp_proximity_gdal(f_foer_bin, actieradius_m)
    unlink(f_foer_bin)
    
    p_5km <- terra::rast(f_prox_5km)
    mask_foer_5km <- terra::ifel(!is.na(p_5km) & p_5km >= 0 & p_5km <= actieradius_m, 1, NA)
    unlink(f_prox_5km)
    rm(p_5km)
    
    vp_finaal <- terra::mask(r_vp_10ha, mask_foer_5km)
    rm(mask_foer_5km)
    gc()
  } else {
    message("    ⚠️ Voortplanting of foerageercomplex ontbreekt vóór 5km-toets!")
    vp_finaal <- master_grid * NA_real_
  }
  rm(r_vp_10ha)
  gc()
  
  ha_vp_fin <- calc_ha_exact(vp_finaal)
  message(sprintf("       -> Voortplanting finaal goedgekeurd binnen 5 km: %.2f ha", ha_vp_fin))
  
  # ----------------------------------------------------------------------------
  # 7. SAMENVOEGEN TOT TOTAAL LEEFGEBIED & NETWERK-CLUSTERING
  # ----------------------------------------------------------------------------
  message("    🧩 Leefgebied samenvoegen (Voortplanting + Foerageerareaal)...")
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
    unlink(f_bin)
    
    cl_opp <- terra::rast(f_clump)
    try(terra::set.cats(cl_opp, layer = 1, value = NULL), silent = TRUE)
    try(terra::coltab(cl_opp) <- NULL, silent = TRUE)
    cl_opp <- terra::mask(cl_opp, r_fin_bin)
    unlink(f_clump)
  } else {
    cl_opp <- master_grid * NA_real_
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