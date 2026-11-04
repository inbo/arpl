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
  
  r_gras_raw   <- if ("biotoop_grasland_nat" %in% names(lagen)) lagen[["biotoop_grasland_nat"]] else (master_grid * NA)
  r_overig_raw <- if ("biotoop_overig"       %in% names(lagen)) lagen[["biotoop_overig"]]       else (master_grid * NA)
  
  # Akkers uitsluiten uit Parquet
  t_akker <- ds_parquet %>%
    dplyr::filter(stringr::str_detect(tolower(CODE), "^b")) %>%
    dplyr::collect() %>%
    data.table::as.data.table()
  
  r_akker_bin <- master_grid * NA
  if (nrow(t_akker) > 0) {
    t_akker[, CODE := tolower(trimws(CODE))]
    t_akker_cel <- unique(t_akker, by = c("cel_id", "CODE"))[, .(Val = 1), by = cel_id]
    r_akker_bin[t_akker_cel$cel_id] <- 1
    rm(t_akker_cel)
  }
  rm(t_akker)
  gc()
  
  # Grasland ontdoen van akkers
  r_gras_zuiver <- terra::mask(r_gras_raw, r_akker_bin, inverse = TRUE)
  rm(r_akker_bin, r_gras_raw)
  gc()
  
  # ----------------------------------------------------------------------------
  # 2. VALLEIGASTURING OP HET GRASLANDBIOTOOP
  # ----------------------------------------------------------------------------
  message("    🌊 Valleisturing toepassen op nat grasland (drainage, profiel & overstroming)...")
  
  # Drainage via centrale helper (gecached)
  mask_drainage <- haal_drainage_masker(paths$drainage, master_grid, paths$dir_cache)
  
  # Bodemprofiel (p, p+x) gecached
  f_cache_prof <- file.path(paths$dir_cache, "masker_profiel_p_px.tif")
  if (file.exists(f_cache_prof)) {
    mask_profiel <- terra::rast(f_cache_prof)
  } else if (file.exists(paths$profiel)) {
    r_prof_orig <- terra::rast(paths$profiel)
    leg_prof    <- terra::cats(r_prof_orig)[[1]]
    val_kol_p   <- if (!is.null(leg_prof) && "Value" %in% colnames(leg_prof)) "Value" else 1
    ok_prof_ids <- if (!is.null(leg_prof) && "Label" %in% colnames(leg_prof)) {
      leg_prof[[val_kol_p]][tolower(trimws(as.character(leg_prof$Label))) %in% c("p", "p+x")]
    } else integer(0)
    rm(r_prof_orig)
    
    f_prof_al <- hulp_align_gdal(paths$profiel, master_grid)
    r_prof    <- terra::rast(f_prof_al)
    mask_profiel <- if (length(ok_prof_ids) > 0) terra::ifel(!is.na(r_prof) & r_prof %in% ok_prof_ids, 1, NA) else (master_grid * NA)
    unlink(f_prof_al)
    rm(r_prof)
    terra::writeRaster(mask_profiel, f_cache_prof, overwrite = TRUE, datatype = "INT1U", gdal = "COMPRESS=LZW")
  } else {
    mask_profiel <- master_grid * NA
  }
  
  # Overstroming (1, 2) gecached
  f_cache_ovstr <- file.path(paths$dir_cache, "masker_overstroming_1_2.tif")
  if (file.exists(f_cache_ovstr)) {
    mask_overstroming <- terra::rast(f_cache_ovstr)
  } else if (file.exists(paths$overstroming)) {
    f_ovstr_al <- hulp_align_gdal(paths$overstroming, master_grid)
    r_ovstr    <- terra::rast(f_ovstr_al)
    mask_overstroming <- terra::ifel(!is.na(r_ovstr) & r_ovstr %in% c(1, 2), 1, NA)
    unlink(f_ovstr_al)
    rm(r_ovstr)
    terra::writeRaster(mask_overstroming, f_cache_ovstr, overwrite = TRUE, datatype = "INT1U", gdal = "COMPRESS=LZW")
  } else {
    mask_overstroming <- master_grid * NA
  }
  
  vallei_masker <- terra::cover(mask_overstroming, terra::cover(mask_drainage, mask_profiel))
  rm(mask_overstroming, mask_drainage, mask_profiel)
  gc()
  
  # Grasland maskeren met de vallei-contour
  r_gras_vallei <- terra::mask(r_gras_zuiver, vallei_masker)
  rm(r_gras_zuiver, vallei_masker)
  gc()
  
  # Samenvoegen met overig biotoop (heide/veen/slikken zijn vrijgesteld van valleitest)
  g_c <- terra::ifel(is.na(r_gras_vallei), 0, r_gras_vallei)
  o_c <- terra::ifel(is.na(r_overig_raw),  0, r_overig_raw)
  r_bio_gecombineerd <- terra::ifel((g_c + o_c) > 0, terra::clamp(g_c + o_c, upper = 1.0), NA)
  rm(g_c, o_c)
  gc()
  
  # ----------------------------------------------------------------------------
  # 3. VERBOSSINGSFILTER (BOSCOMPLEXEN >= 20 ARE UITSLUITEN)
  # ----------------------------------------------------------------------------
  message("    🌲 Verbossingsfilter: Hooggroencomplexen >= 20 are integraal uitsluiten...")
  
  f_groen_al <- hulp_align_gdal(paths$groenkaart, master_grid)
  r_groen    <- terra::rast(f_groen_al)
  r_groen_bin <- terra::ifel(!is.na(r_groen) & (r_groen == 1 | r_groen > 0), 1, NA)
  unlink(f_groen_al)
  rm(r_groen)
  
  # Filter bosclusters >= 20 cellen (20 are = 2000 m2) via GRASS
  f_g_bin   <- tempfile(fileext = ".tif")
  f_cl_g    <- tempfile(fileext = ".tif")
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
  terra::set.cats(cl_bos, layer = 1, value = NULL)
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
  
  r_bio_ontbost <- terra::mask(r_bio_gecombineerd, mask_verbossing, inverse = TRUE)
  rm(r_bio_gecombineerd, mask_verbossing)
  gc()
  
  ha_ontbost <- calc_ha_exact(r_bio_ontbost)
  message(sprintf("       -> Biotoop na vallei- en verbossingsfilter: %.2f ha", ha_ontbost))
  
  # ----------------------------------------------------------------------------
  # 4. LANDSCHAPPELIJKE CLUSTERING (>= 100 HA OP 100M NETWERK)
  # ----------------------------------------------------------------------------
  message("    🏝️ Landschappelijke clustering (>= 100 ha, 100m overbrugging)...")
  res_100ha <- cluster_filter_universeel(
    masker            = r_bio_ontbost,
    opp_laag          = r_bio_ontbost,
    min_opp_ha        = 100,
    min_totale_opp_ha = NA,
    afstand_m         = 100,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  final_opp <- res_100ha$raster
  cl_opp    <- res_100ha$clusters
  rm(res_100ha, r_bio_ontbost)
  gc()
  
  # Deellagen isoleren binnen de goedgekeurde 100 ha complexen
  gras_uitvoer   <- terra::mask(r_gras_vallei, cl_opp)
  overig_uitvoer <- terra::mask(r_overig_raw,  cl_opp)
  rm(r_gras_vallei, r_overig_raw)
  gc()
  
  return(list(
    final_opp        = final_opp,
    cl_opp           = cl_opp,
    grasland_vallei  = gras_uitvoer,
    overig           = overig_uitvoer
  ))
}
