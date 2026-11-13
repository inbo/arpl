# ==============================================================================
# bouw_biotoop_bos_landschap_actieradius.R
# Archetype voor bosroofvogels (broedbos x open foerageer x landschap >= 250ha x 5km buffer)
# Geschikt voor: Wespendief, Boomvalk, etc.
# ==============================================================================

bouw_biotoop_bos_landschap_actieradius <- function(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
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
  
  r_bos_raw  <- if ("voortplanting_bwk" %in% names(lagen)) lagen[["voortplanting_bwk"]] else (master_grid * NA)
  r_open_raw <- if ("foerageer_open"     %in% names(lagen)) lagen[["foerageer_open"]]     else (master_grid * NA)
  
  # ----------------------------------------------------------------------------
  # 2. STAP 1: BROEDBOS (VERIFICATIE HOOGGROEN & CLUSTEREN OP >= 50 HA OP 100M)
  # ----------------------------------------------------------------------------
  message("    🌲 Broedbos valideren met Groenkaart (strikt klasse 1) en clusteren op >= 50 ha (100m netwerk)...")
  
  f_groen_al     <- hulp_align_gdal(paths$groenkaart, master_grid)
  r_groen        <- terra::rast(f_groen_al)
  # Strikt klasse 1 = opgaand bos / bomen in Groenkaart 2021
  r_bomen_direct <- terra::ifel(!is.na(r_groen) & r_groen == 1, 1, NA)
  unlink(f_groen_al)
  rm(r_groen)
  
  # Broedbos valideren tegen echte opgaande boomkruinen
  r_bos_opgaand <- terra::mask(r_bos_raw, r_bomen_direct)
  rm(r_bomen_direct, r_bos_raw)
  gc()
  
  res_bos <- cluster_filter_universeel(
    masker            = r_bos_opgaand,
    opp_laag          = r_bos_opgaand,
    min_opp_ha        = 50,
    min_totale_opp_ha = NA,
    afstand_m         = 100,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_bos_50ha <- res_bos$raster
  rm(res_bos, r_bos_opgaand)
  gc()
  
  ha_bos_50ha <- calc_ha_exact(r_bos_50ha)
  message(sprintf("       -> Goedgekeurd broedbos (>= 50 ha): %.2f ha", ha_bos_50ha))
  
  # ----------------------------------------------------------------------------
  # 3. STAP 2: OPEN FOERAGEERHABITAT (CLUSTEREN OP >= 50 HA OP 500M)
  # ----------------------------------------------------------------------------
  message("    🌾 Open foerageerhabitat clusteren op >= 50 ha (500m netwerk)...")
  res_open <- cluster_filter_universeel(
    masker            = r_open_raw,
    opp_laag          = r_open_raw,
    min_opp_ha        = 50,
    min_totale_opp_ha = NA,
    afstand_m         = 500,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_open_50ha <- res_open$raster
  rm(res_open, r_open_raw)
  gc()
  
  ha_open_50ha <- calc_ha_exact(r_open_50ha)
  message(sprintf("       -> Goedgekeurd open foerageerareaal (>= 50 ha): %.2f ha", ha_open_50ha))
  
  # ----------------------------------------------------------------------------
  # 4. STAP 3: COMPOSIET FOERAGEERLANDSCHAP (>= 250 HA OP 500M NETWERK)
  # ----------------------------------------------------------------------------
  message("    🧩 Composiet landschap samenstellen (bos + open) en clusteren op >= 250 ha...")
  
  # Combinatie conform Bert (cover met clamp)
  r_landschap_bruto <- terra::cover(r_bos_50ha, r_open_50ha)
  r_landschap_bruto <- terra::clamp(r_landschap_bruto, upper = 1.0)
  
  res_land250 <- cluster_filter_universeel(
    masker            = r_landschap_bruto,
    opp_laag          = r_landschap_bruto,
    min_opp_ha        = 250,
    min_totale_opp_ha = NA,
    afstand_m         = 500,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_foer_250ha <- res_land250$raster
  rm(res_land250, r_landschap_bruto)
  gc()
  
  ha_foer_250ha <- calc_ha_exact(r_foer_250ha)
  message(sprintf("       -> Totaal foerageerlandschap (>= 250 ha): %.2f ha", ha_foer_250ha))
  
  # ----------------------------------------------------------------------------
  # 5. STAP 4: ACTIERADIUS-KOPPELING (5 KM) TUSSEN BROEDBOS EN FOERAGEERLANDSCHAP
  # ----------------------------------------------------------------------------
  message("    🎯 5 km actieradiustoets uitvoeren via GDAL proximity...")
  
  heeft_bos  <- !all(is.na(suppressWarnings(terra::minmax(r_bos_50ha))))
  heeft_foer <- !all(is.na(suppressWarnings(terra::minmax(r_foer_250ha))))
  
  if (heeft_bos && heeft_foer) {
    r_foer_bin <- terra::ifel(!is.na(r_foer_250ha) & r_foer_250ha > 0, 1, NA)
    
    f_foer_bin  <- tempfile(fileext = ".tif")
    terra::writeRaster(r_foer_bin, f_foer_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    f_prox_5km  <- hulp_proximity_gdal(f_foer_bin, 5000)
    unlink(f_foer_bin)
    rm(r_foer_bin)
    gc()
    
    p_5km       <- terra::rast(f_prox_5km)
    mask_foer5k <- terra::ifel(!is.na(p_5km) & p_5km >= 0 & p_5km <= 5000, 1, NA)
    unlink(f_prox_5km)
    rm(p_5km)
    gc()
    
    # Broedbossen filteren: moeten binnen 5 km van het grote foerageerlandschap liggen
    bos_finaal <- terra::mask(r_bos_50ha, mask_foer5k)
    rm(mask_foer5k)
    gc()
    
    # Foerageergebied finaal: open biotoop dat deel uitmaakt van het 250ha complex
    foer_finaal <- terra::mask(r_open_50ha, r_foer_250ha)
  } else {
    bos_finaal  <- master_grid * NA
    foer_finaal <- master_grid * NA
  }
  rm(r_bos_50ha, r_open_50ha, r_foer_250ha)
  gc()
  
  ha_bos_fin  <- calc_ha_exact(bos_finaal)
  ha_foer_fin <- calc_ha_exact(foer_finaal)
  message(sprintf("       -> Broedbos goedgekeurd binnen 5 km: %.2f ha", ha_bos_fin))
  message(sprintf("       -> Open foerageergebied binnen goedgekeurd complex: %.2f ha", ha_foer_fin))
  
  # ----------------------------------------------------------------------------
  # 6. SAMENSTELLEN COMPOSIET LEEFGEBIED (BROEDBOS + FOERAGEERAREAAL) & NETWERK-ID'S
  # ----------------------------------------------------------------------------
  message("    🧩 Totaal leefgebiedscomplex samenstellen (broedbos + foerageerareaal)...")
  final_opp <- terra::cover(bos_finaal, foer_finaal)
  final_opp <- terra::clamp(final_opp, upper = 1.0)
  
  ha_tot <- calc_ha_exact(final_opp)
  message(sprintf("       -> Totaal leefgebied Wespendief: %.2f ha", ha_tot))
  
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
    final_opp     = final_opp,
    cl_opp        = cl_opp,
    voortplanting = bos_finaal,
    foerageren    = foer_finaal
  ))
}
