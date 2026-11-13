# ==============================================================================
# bouw_biotoop_open_zand_stuifduin.R
# Archetype voor pioniers van open zand & heide (complex >= 75ha x stuifduinkern >= 5ha binnen 10m)
# Geschikt voor: Duinpieper, etc.
# ==============================================================================

bouw_biotoop_open_zand_stuifduin <- function(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
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
  
  r_complex_raw  <- if ("voortplanting" %in% names(lagen)) lagen[["voortplanting"]] else (master_grid * NA)
  r_stuif_raw    <- if ("stuifduin"     %in% names(lagen)) lagen[["stuifduin"]]     else (master_grid * NA)
  
  # ----------------------------------------------------------------------------
  # 2. STAP 1: STUIFDUIN- EN STUIFZANDKERNEN (>= 5 HA OP 20M NETWERK)
  # ----------------------------------------------------------------------------
  message("    🏜️ Stuifduinkernen clusteren op >= 5 ha (20m netwerk)...")
  
  res_stuif <- cluster_filter_universeel(
    masker            = r_stuif_raw,
    opp_laag          = r_stuif_raw,
    min_opp_ha        = 5,
    min_totale_opp_ha = NA,
    afstand_m         = 20,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_stuif_5ha <- res_stuif$raster
  rm(res_stuif, r_stuif_raw)
  gc()
  
  ha_stuif_5ha <- calc_ha_exact(r_stuif_5ha)
  message(sprintf("       -> Goedgekeurde stuifduinkernen (>= 5 ha): %.2f ha", ha_stuif_5ha))
  
  # ----------------------------------------------------------------------------
  # 3. STAP 2: LANDSCHAPPELIJK LEEFGEBIEDCOMPLEX (>= 75 HA OP 50M NETWERK)
  # ----------------------------------------------------------------------------
  message("    🌾 Totaal heide- en duincomplex clusteren op >= 75 ha (50m netwerk)...")
  
  res_complex <- cluster_filter_universeel(
    masker            = r_complex_raw,
    opp_laag          = r_complex_raw,
    min_opp_ha        = 75,
    min_totale_opp_ha = NA,
    afstand_m         = 50,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_complex_75ha <- res_complex$raster
  rm(res_complex, r_complex_raw)
  gc()
  
  ha_complex_75ha <- calc_ha_exact(r_complex_75ha)
  message(sprintf("       -> Goedgekeurde heide-/duincomplexen (>= 75 ha): %.2f ha", ha_complex_75ha))
  
  # ----------------------------------------------------------------------------
  # 4. STAP 3: KOPPELING EN INTERACTIEZONE (10M BUFFER ROND STUIFDUINKERNEN)
  # ----------------------------------------------------------------------------
  message("    🎯 Complex koppelen aan stuifduinkernen (<= 10m afstand via GDAL proximity)...")
  
  heeft_complex <- !all(is.na(suppressWarnings(terra::minmax(r_complex_75ha))))
  heeft_stuif   <- !all(is.na(suppressWarnings(terra::minmax(r_stuif_5ha))))
  
  if (heeft_complex && heeft_stuif) {
    r_stuif_bin <- terra::ifel(!is.na(r_stuif_5ha) & r_stuif_5ha > 0, 1, NA)
    
    f_s_bin    <- tempfile(fileext = ".tif")
    terra::writeRaster(r_stuif_bin, f_s_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 0)
    f_prox_s10 <- hulp_proximity_gdal(f_s_bin, 10)
    unlink(f_s_bin)
    rm(r_stuif_bin)
    
    p_s10 <- terra::rast(f_prox_s10)
    mask_stuif10m <- terra::ifel(!is.na(p_s10) & p_s10 >= 0 & p_s10 <= 10, 1, NA)
    unlink(f_prox_s10)
    rm(p_s10)
    gc()
    
    # Complexen selecteren die in contact staan met een 5ha stuifduinkern
    # 1. Bepaal welke 75ha clusters contact maken
    # Om het hele 75ha complex te behouden dat contact maakt met de stuifduin:
    cl_complex <- terra::patches(r_complex_75ha, directions = 8, zeroAsNA = TRUE)
    
    cl_contact <- terra::mask(cl_complex, mask_stuif10m)
    contact_ids <- unique(na.omit(terra::values(cl_contact, mat = FALSE)))
    rm(cl_contact, mask_stuif10m)
    gc()
    
    if (length(contact_ids) > 0) {
      final_opp <- terra::ifel(cl_complex %in% as.integer(contact_ids), r_complex_75ha, NA)
      stuif_uitvoer <- terra::mask(r_stuif_5ha, final_opp)
    } else {
      final_opp     <- master_grid * NA
      stuif_uitvoer <- master_grid * NA
    }
    rm(cl_complex)
    gc()
  } else {
    final_opp     <- master_grid * NA
    stuif_uitvoer <- master_grid * NA
  }
  rm(r_complex_75ha, r_stuif_5ha)
  gc()
  
  ha_finaal <- calc_ha_exact(final_opp)
  message(sprintf("       -> Finaal leefgebied na stuifduinkoppeling: %.2f ha", ha_finaal))
  
  # ----------------------------------------------------------------------------
  # 5. NETWERK-ID'S GENEREREN
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
    
    cl_opp <- terra::rast(f_clump)
    terra::set.cats(cl_opp, layer = 1, value = NULL)
    if (terra::has.colors(cl_opp)) {
      try(terra::coltab(cl_opp) <- NULL, silent = TRUE)
    }
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
    stuifduinkern = stuif_uitvoer
  ))
}
