# ==============================================================================
# bouw_biotoop_grootschalige_weidevogel.R
# Archetype voor grootschalige open-landschap steltlopers (valleigrasland x heide/slik x ontbossing x >= 100ha)
# Geschikt voor: Wulp, etc.
# ==============================================================================

bouw_biotoop_grootschalige_weidevogel <- function(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
  bestandsnaam <- tolower(gsub(" ", "", soort_naam))
  soort_cache  <- file.path(paths$dir_cache, bestandsnaam)
  dir.create(soort_cache, showWarnings = FALSE, recursive = TRUE)
  
  f_cache_rds <- file.path(soort_cache, "01_raw_parquet_layers.rds")
  
  # Parameters ophalen
  rij_params    <- params_soort %>% dplyr::slice(1)
  min_opp_ha    <- as.numeric(ifelse(is.na(rij_params$MinOpp_ha), 100, rij_params$MinOpp_ha))
  zoekafstand_m <- as.numeric(ifelse(is.na(rij_params$AfstandBiotopen_m), 100, rij_params$AfstandBiotopen_m))
  
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
  
  r_gras_raw   <- if ("biotoop_grasland_nat" %in% names(lagen)) lagen[["biotoop_grasland_nat"]] else (master_grid * NA)
  r_overig_raw <- if ("biotoop_overig"       %in% names(lagen)) lagen[["biotoop_overig"]]       else (master_grid * NA)
  
  # Akkers veilig uitsluiten via centrale Parquet helper (^b codes)
  r_akker_bin <- vul_raster_vanuit_parquet(
    ds_parquet  = ds_parquet,
    codes_bevat = "b",
    master_grid = master_grid,
    waarde_vast = 1
  )
  
  # Grasland ontdoen van akkers
  if (!all(is.na(suppressWarnings(terra::minmax(r_akker_bin))))) {
    r_gras_zuiver <- terra::mask(r_gras_raw, r_akker_bin, inverse = TRUE)
  } else {
    r_gras_zuiver <- r_gras_raw
  }
  rm(r_akker_bin, r_gras_raw)
  gc()
  
  ha_gras_zuiver <- calc_ha_exact(r_gras_zuiver)
  message(sprintf("       -> Graslandbiotoop na akker-uitsluiting: %.2f ha", ha_gras_zuiver))
  
  # ----------------------------------------------------------------------------
  # 2. VALLEISTURING OP HET GRASLANDBIOTOOP
  # ----------------------------------------------------------------------------
  message("    🌊 Valleisturing toepassen op nat grasland (drainage, profiel & overstroming)...")
  
  # A. Drainage via centrale helper (gecached)
  mask_drainage <- master_grid * NA
  if (exists("haal_drainage_masker", mode = "function")) {
    mask_drainage <- haal_drainage_masker(paths$drainage, master_grid, paths$dir_cache)
  } else if (!is.null(paths$drainage) && file.exists(paths$drainage)) {
    f_drain_al <- hulp_align_gdal(paths$drainage, master_grid)
    r_drain    <- terra::rast(f_drain_al)
    leg_drain  <- terra::cats(r_drain)[[1]]
    drain_codes_ok <- c("d", "e", "f", "g", "h", "i", "e-f", "h-i", "e-i")
    ok_d_ids <- integer(0)
    if (!is.null(leg_drain)) {
      kol_label_d <- colnames(leg_drain)[tolower(colnames(leg_drain)) %in% c("label", "omschrijving", "drainage")][1]
      kol_val_d   <- colnames(leg_drain)[tolower(colnames(leg_drain)) %in% c("value", "id")][1]
      if (is.na(kol_val_d)) kol_val_d <- 1
      if (!is.na(kol_label_d)) {
        ok_d_ids <- leg_drain[[kol_val_d]][tolower(trimws(as.character(leg_drain[[kol_label_d]]))) %in% drain_codes_ok]
      }
    }
    if (length(ok_d_ids) == 0) ok_d_ids <- c(4, 5, 6, 7, 8, 9, 45, 56, 89)
    mask_drainage <- terra::ifel(!is.na(r_drain) & (r_drain %in% ok_d_ids), 1, NA)
    unlink(f_drain_al)
    rm(r_drain)
  }
  
  # B. Bodemprofiel (p, p+x) gecached
  f_cache_prof <- file.path(paths$dir_cache, "masker_profiel_p_px.tif")
  if (file.exists(f_cache_prof)) {
    mask_profiel <- terra::rast(f_cache_prof)
  } else if (!is.null(paths$profiel) && file.exists(paths$profiel)) {
    r_prof_orig <- terra::rast(paths$profiel)
    leg_prof    <- terra::cats(r_prof_orig)[[1]]
    ok_prof_ids <- integer(0)
    if (!is.null(leg_prof)) {
      kol_label <- colnames(leg_prof)[tolower(colnames(leg_prof)) %in% c("label", "omschrijving", "naam", "profiel")][1]
      kol_val   <- colnames(leg_prof)[tolower(colnames(leg_prof)) %in% c("value", "id")][1]
      if (is.na(kol_val)) kol_val <- 1
      if (!is.na(kol_label)) {
        ok_prof_ids <- leg_prof[[kol_val]][tolower(trimws(as.character(leg_prof[[kol_label]]))) %in% c("p", "p+x")]
      }
    }
    rm(r_prof_orig)
    
    f_prof_al    <- hulp_align_gdal(paths$profiel, master_grid)
    r_prof       <- terra::rast(f_prof_al)
    mask_profiel <- if (length(ok_prof_ids) > 0) terra::ifel(!is.na(r_prof) & r_prof %in% ok_prof_ids, 1, NA) else (master_grid * NA)
    unlink(f_prof_al)
    rm(r_prof)
    terra::writeRaster(mask_profiel, f_cache_prof, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
  } else {
    mask_profiel <- master_grid * NA
  }
  
  # C. Overstroming (1, 2) gecached
  f_cache_ovstr <- file.path(paths$dir_cache, "masker_overstroming_1_2.tif")
  if (file.exists(f_cache_ovstr)) {
    mask_overstroming <- terra::rast(f_cache_ovstr)
  } else if (!is.null(paths$overstroming) && file.exists(paths$overstroming)) {
    f_ovstr_al        <- hulp_align_gdal(paths$overstroming, master_grid)
    r_ovstr           <- terra::rast(f_ovstr_al)
    mask_overstroming <- terra::ifel(!is.na(r_ovstr) & (r_ovstr %in% c(1, 2)), 1, NA)
    unlink(f_ovstr_al)
    rm(r_ovstr)
    terra::writeRaster(mask_overstroming, f_cache_ovstr, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
  } else {
    mask_overstroming <- master_grid * NA
  }
  
  vallei_masker <- terra::cover(mask_overstroming, terra::cover(mask_drainage, mask_profiel))
  vallei_masker <- terra::ifel(!is.na(vallei_masker) & vallei_masker > 0, 1, NA)
  rm(mask_overstroming, mask_drainage, mask_profiel)
  gc()
  
  # Grasland maskeren met de vallei-contour
  r_gras_vallei <- terra::mask(r_gras_zuiver, vallei_masker)
  rm(r_gras_zuiver, vallei_masker)
  gc()
  
  # Samenvoegen met overig biotoop (heide/veen/slikken zijn vrijgesteld van valleitest)
  r_bio_gecombineerd <- terra::cover(r_gras_vallei, r_overig_raw)
  r_bio_gecombineerd <- terra::clamp(r_bio_gecombineerd, upper = 1.0)
  
  ha_bio_comb <- calc_ha_exact(r_bio_gecombineerd)
  message(sprintf("       -> Biotoop na valleitest (grasland in vallei + overig biotoop): %.2f ha", ha_bio_comb))
  
  # ----------------------------------------------------------------------------
  # 3. VERBOSSINGSFILTER (BOSCOMPLEXEN >= 20 ARE UITSLUITEN)
  # ----------------------------------------------------------------------------
  message("    🌲 Verbossingsfilter: Hooggroencomplexen >= 20 are integraal uitsluiten...")
  
  f_cache_verbossing <- file.path(paths$dir_cache, "masker_verbossing_20are.tif")
  
  if (file.exists(f_cache_verbossing)) {
    message("       -> [CACHE] Verbossingsmasker ingeladen van schijf...")
    mask_verbossing <- terra::rast(f_cache_verbossing)
  } else {
    f_groen_al  <- hulp_align_gdal(paths$groenkaart, master_grid)
    r_groen     <- terra::rast(f_groen_al)
    # Strikt klasse 1 = bomen/bos in Groenkaart_2021
    r_groen_bin <- terra::ifel(!is.na(r_groen) & r_groen == 1, 1, NA)
    unlink(f_groen_al)
    rm(r_groen)
    
    # Filter bosclusters >= 20 cellen (20 are = 2000 m2) via GRASS r.clump
    f_g_bin  <- tempfile(fileext = ".tif")
    f_cl_g   <- tempfile(fileext = ".tif")
    terra::writeRaster(r_groen_bin, f_g_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    rm(r_groen_bin)
    gc()
    
    e <- terra::ext(master_grid)
    qgisprocess::qgis_run_algorithm(
      "grass:r.clump", input = f_g_bin, output = f_cl_g, title = "cl_bos20are_wulp",
      "-d" = FALSE, GRASS_REGION_PARAMETER = sprintf("%f,%f,%f,%f", e$xmin, e$xmax, e$ymin, e$ymax),
      GRASS_REGION_CELLSIZE_PARAMETER = 10, .quiet = TRUE
    )
    unlink(f_g_bin)
    
    cl_bos <- terra::rast(f_cl_g)
    try(terra::set.cats(cl_bos, layer = 1, value = NULL), silent = TRUE)
    try(terra::coltab(cl_bos) <- NULL, silent = TRUE)
    
    dt_bos <- data.table::as.data.table(terra::freq(cl_bos))
    setnames(dt_bos, c("layer", "ID", "count"))
    grote_bos_ids <- dt_bos[!is.na(ID) & ID > 0 & count >= 20, ID]
    rm(dt_bos)
    
    if (length(grote_bos_ids) > 0) {
      mask_verbossing <- terra::ifel(cl_bos %in% as.integer(grote_bos_ids), 1, NA)
    } else {
      mask_verbossing <- master_grid * NA
    }
    unlink(f_cl_g)
    rm(cl_bos)
    gc()
    
    terra::writeRaster(mask_verbossing, f_cache_verbossing, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
  }
  
  # Uitsluiten van boscomplexen
  if (!all(is.na(suppressWarnings(terra::minmax(mask_verbossing))))) {
    r_bio_ontbost <- terra::mask(r_bio_gecombineerd, mask_verbossing, inverse = TRUE)
  } else {
    r_bio_ontbost <- r_bio_gecombineerd
  }
  rm(r_bio_gecombineerd, mask_verbossing)
  gc()
  
  ha_ontbost <- calc_ha_exact(r_bio_ontbost)
  message(sprintf("       -> Biotoop na vallei- en verbossingsfilter: %.2f ha", ha_ontbost))
  
  # ----------------------------------------------------------------------------
  # 4. LANDSCHAPPELIJKE CLUSTERING (>= 100 HA OP 100M NETWERK)
  # ----------------------------------------------------------------------------
  message(sprintf("    🏝️ Landschappelijke clustering (>= %.1f ha, %dm overbrugging)...", min_opp_ha, as.integer(zoekafstand_m)))
  
  res_100ha <- cluster_filter_universeel(
    masker            = r_bio_ontbost,
    opp_laag          = r_bio_ontbost,
    min_opp_ha        = min_opp_ha,
    min_totale_opp_ha = NA,
    afstand_m         = zoekafstand_m,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  final_opp <- res_100ha$raster
  cl_opp    <- res_100ha$clusters
  rm(res_100ha, r_bio_ontbost)
  gc()
  
  ha_finaal <- calc_ha_exact(final_opp)
  message(sprintf("       -> Definitief leefgebied Wulp: %.2f ha", ha_finaal))
  
  # Deellagen isoleren binnen de goedgekeurde complexen
  if (!all(is.na(suppressWarnings(terra::minmax(cl_opp))))) {
    gras_uitvoer   <- terra::mask(r_gras_vallei, cl_opp)
    overig_uitvoer <- terra::mask(r_overig_raw,  cl_opp)
  } else {
    gras_uitvoer   <- master_grid * NA
    overig_uitvoer <- master_grid * NA
  }
  rm(r_gras_vallei, r_overig_raw)
  gc()
  
  return(list(
    final_opp        = final_opp,
    cl_opp           = cl_opp,
    grasland_vallei  = gras_uitvoer,
    overig           = overig_uitvoer
  ))
}
