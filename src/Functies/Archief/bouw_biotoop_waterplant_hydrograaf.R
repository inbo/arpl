# ==============================================================================
# bouw_biotoop_waterplant_hydrograaf.R
# Archetype voor hydrofieten en waterplanten (BWK stilstaand water x Huet beken/kanalen x 20m overbrugging)
# Geschikt voor: Drijvende waterweegbree (Luronium natans), etc.
# ==============================================================================

bouw_biotoop_waterplant_hydrograaf <- function(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
  bestandsnaam <- tolower(gsub(" ", "", soort_naam))
  soort_cache  <- file.path(paths$dir_cache, bestandsnaam)
  dir.create(soort_cache, showWarnings = FALSE, recursive = TRUE)
  
  f_cache_rds <- file.path(soort_cache, "01_raw_parquet_layers.rds")
  
  # Parameters ophalen uit Excel
  rij_params    <- params_soort %>% dplyr::slice(1)
  min_opp_ha    <- as.numeric(ifelse(is.na(rij_params$MinOpp_ha), 0.16, rij_params$MinOpp_ha))
  zoekafstand_m <- as.numeric(ifelse(is.na(rij_params$AfstandBiotopen_m), 20, rij_params$AfstandBiotopen_m))
  
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
  
  r_bwk_raw <- if ("voortplanting_bwk" %in% names(lagen)) lagen[["voortplanting_bwk"]] else (master_grid * NA)
  
  # ----------------------------------------------------------------------------
  # 2. HYDROGRAFISCHE HUET-ZONES INLADEN & SYNCHRONISEREN
  # ----------------------------------------------------------------------------
  message("    🌊 Geschikte Huetzones (kanalen, vaarten, beken) filteren en cachen...")
  
  f_cache_huet_weegbree <- file.path(paths$dir_cache, "masker_huetzon_weegbree.tif")
  
  if (file.exists(f_cache_huet_weegbree)) {
    message("       -> [CACHE] Huetzone-laag voor waterweegbree ingeladen...")
    r_huet_ok <- terra::rast(f_cache_huet_weegbree)
  } else if (!is.null(paths$huetzon) && file.exists(paths$huetzon)) {
    f_huet_al <- hulp_align_gdal(paths$huetzon, master_grid)
    r_huet    <- terra::rast(f_huet_al)
    leg_huet  <- terra::cats(r_huet)[[1]]
    
    # Doellabels uit Bert's script: beken, kanalen en gekanaliseerde waterlopen
    target_labels <- c("bkk1", "bkk2", "kabpkunst", "kan", "kanpkunst", "kunst")
    ok_huet_ids   <- integer(0)
    
    if (!is.null(leg_huet)) {
      kol_val   <- colnames(leg_huet)[tolower(colnames(leg_huet)) %in% c("value", "id")][1]
      kol_label <- colnames(leg_huet)[tolower(colnames(leg_huet)) %in% c("label", "omschrijving", "huet")][1]
      if (is.na(kol_val)) kol_val <- 1
      if (!is.na(kol_label)) {
        huet_labels_clean <- tolower(trimws(as.character(leg_huet[[kol_label]])))
        ok_huet_ids <- leg_huet[[kol_val]][huet_labels_clean %in% target_labels]
      }
    }
    
    # Veilige numerieke fallback conform Bert's script indien TIFF geen tekstcats heeft
    if (length(ok_huet_ids) == 0) ok_huet_ids <- 5:10
    
    r_huet_ok <- terra::ifel(!is.na(r_huet) & (r_huet %in% ok_huet_ids), 0.01, NA)
    unlink(f_huet_al)
    rm(r_huet)
    terra::writeRaster(r_huet_ok, f_cache_huet_weegbree, overwrite = TRUE, datatype = "FLT4S", gdal = "COMPRESS=LZW")
  } else {
    r_huet_ok <- master_grid * NA
  }
  
  # ----------------------------------------------------------------------------
  # 3. FUSIE WATERHABITAT (BWK + NOMINALE HUET FRACTIE 0.01)
  # ----------------------------------------------------------------------------
  message("    💧 BWK-water en Huetzones samenvoegen...")
  
  # Fusie via cover: BWK-fracties primeren, Huet vult aan met 0.01
  r_water_tot <- terra::cover(r_bwk_raw, r_huet_ok)
  r_water_tot <- terra::clamp(r_water_tot, upper = 1.0)
  
  ha_bwk_tot  <- calc_ha_exact(r_bwk_raw)
  ha_huet_tot <- calc_ha_exact(r_huet_ok)
  ha_water_tot <- calc_ha_exact(r_water_tot)
  message(sprintf("       -> BWK stilstaand water: %.2f ha", ha_bwk_tot))
  message(sprintf("       -> Huet stromend water & kanalen: %.2f ha", ha_huet_tot))
  message(sprintf("       -> Totaal waterhabitat vóór clustering: %.2f ha", ha_water_tot))
  
  # ----------------------------------------------------------------------------
  # 4. RUIMTELIJKE CLUSTERING & OPPERVLAKTEFILTER (20M OVERBRUGGING & 0.16 HA)
  # ----------------------------------------------------------------------------
  message(sprintf("    🏝️ Duikers/bruggen overbruggen (%dm) en filteren op >= %.2f ha...", as.integer(zoekafstand_m), min_opp_ha))
  
  res_cluster <- cluster_filter_universeel(
    masker            = r_water_tot,
    opp_laag          = r_water_tot,
    min_opp_ha        = min_opp_ha,
    min_totale_opp_ha = NA,
    afstand_m         = zoekafstand_m,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  
  final_opp <- res_cluster$raster
  cl_opp    <- res_cluster$clusters
  rm(res_cluster, r_water_tot)
  gc()
  
  ha_finaal <- calc_ha_exact(final_opp)
  message(sprintf("       -> Definitief leefgebied Drijvende waterweegbree: %.2f ha", ha_finaal))
  
  # Deellagen isoleren binnen de goedgekeurde watercomplexen
  if (!all(is.na(suppressWarnings(terra::minmax(cl_opp))))) {
    bwk_uitvoer  <- terra::mask(r_bwk_raw, cl_opp)
    huet_uitvoer <- terra::mask(r_huet_ok, cl_opp)
  } else {
    bwk_uitvoer  <- master_grid * NA
    huet_uitvoer <- master_grid * NA
  }
  rm(r_bwk_raw, r_huet_ok)
  gc()
  
  return(list(
    final_opp = final_opp,
    cl_opp    = cl_opp,
    bwk_water = bwk_uitvoer,
    huetzones = huet_uitvoer
  ))
}