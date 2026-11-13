# ==============================================================================
# bouw_biotoop_saproxyl_kle_park.R
# Archetype voor holle-bomen kevers (parken/knotbomen kj/kp x Groenkaart bomenrijen x >= 0.5 ha op 50m)
# Geschikt voor: Juchtleerkever (Osmoderma eremita), etc.
# ==============================================================================

bouw_biotoop_saproxyl_kle_park <- function(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
  bestandsnaam <- tolower(gsub(" ", "", soort_naam))
  soort_cache  <- file.path(paths$dir_cache, bestandsnaam)
  dir.create(soort_cache, showWarnings = FALSE, recursive = TRUE)
  
  f_cache_rds <- file.path(soort_cache, "01_raw_parquet_layers.rds")
  
  # Parameters ophalen uit Excel
  rij_params    <- params_soort %>% dplyr::slice(1)
  min_opp_ha    <- as.numeric(ifelse(is.na(rij_params$MinOpp_ha), 0.5, rij_params$MinOpp_ha))
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
  
  r_bwk_raw <- if ("bwk" %in% names(lagen)) lagen[["bwk"]] else (master_grid * NA)
  
  # ----------------------------------------------------------------------------
  # 2. GROENKAART OPGAANDE BREDERE KLE'S & RANDZONES (KLASSE 1 GRENSEND AAN OPEN)
  # ----------------------------------------------------------------------------
  message("    🌳 Groenkaart-bomenrijen en zonnige randen isoleren (edge detectie)...")
  
  f_cache_groen_kle <- file.path(paths$dir_cache, "masker_groenkaart_kle_bomen.tif")
  
  if (file.exists(f_cache_groen_kle)) {
    message("       -> [CACHE] Groenkaart-KLE laag ingeladen van schijf...")
    r_kle_groen <- terra::rast(f_cache_groen_kle)
  } else if (!is.null(paths$groenkaart) && file.exists(paths$groenkaart)) {
    f_groen_al <- hulp_align_gdal(paths$groenkaart, master_grid)
    r_groen    <- terra::rast(f_groen_al)
    
    # Klasse 1 = bomen, Klassen 2, 3, 4 = open / laag groen
    r_bomen_bin <- terra::ifel(!is.na(r_groen) & r_groen == 1, 1, NA)
    r_open_bin  <- terra::ifel(!is.na(r_groen) & (r_groen %in% c(2, 3, 4)), 1, NA)
    unlink(f_groen_al)
    rm(r_groen)
    
    # Bomen die binnen 10m van open terrein liggen (dreven, houtkanten, randen)
    if (!all(is.na(suppressWarnings(terra::minmax(r_bomen_bin))))) {
      f_op_bin  <- tempfile(fileext = ".tif")
      terra::writeRaster(r_open_bin, f_op_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
      f_prox_op <- hulp_proximity_gdal(f_op_bin, 10)
      unlink(f_op_bin)
      rm(r_open_bin)
      
      p_op <- terra::rast(f_prox_op)
      mask_op10 <- terra::ifel(!is.na(p_op) & p_op >= 0 & p_op <= 10, 1, NA)
      unlink(f_prox_op)
      rm(p_op)
      
      r_kle_groen <- terra::mask(r_bomen_bin, mask_op10)
      rm(mask_op10, r_bomen_bin)
    } else {
      r_kle_groen <- master_grid * NA
      rm(r_bomen_bin, r_open_bin)
    }
    
    terra::writeRaster(r_kle_groen, f_cache_groen_kle, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
  } else {
    r_kle_groen <- master_grid * NA
  }
  
  # ----------------------------------------------------------------------------
  # 3. FUSIE BASISBIOTOOP (BWK PARKEN/KNOTBOMEN + GROENKAART BREDERE DREVEN)
  # ----------------------------------------------------------------------------
  message("    🧩 BWK parken en Groenkaart bomenrijen samenvoegen...")
  
  # Cover: BWK-fractie primeert, Groenkaart vult solitaire dreven aan
  r_kle_tot <- terra::cover(r_bwk_raw, r_kle_groen)
  r_kle_tot <- terra::clamp(r_kle_tot, upper = 1.0)
  
  ha_bwk_ruw  <- calc_ha_exact(r_bwk_raw)
  ha_groen_ruw <- calc_ha_exact(r_kle_groen)
  ha_kle_tot  <- calc_ha_exact(r_kle_tot)
  
  message(sprintf("       -> BWK parken/knotbomen (kp, kj): %.2f ha", ha_bwk_ruw))
  message(sprintf("       -> Groenkaart bomenrijen/randen: %.2f ha", ha_groen_ruw))
  message(sprintf("       -> Totaal geschikt bomenhabitat vóór clustering: %.2f ha", ha_kle_tot))
  
  # ----------------------------------------------------------------------------
  # 4. RUIMTELIJKE CLUSTERING & OPPERVLAKTEDREMPEL (>= 0.5 HA OP 50 M)
  # ----------------------------------------------------------------------------
  message(sprintf("    🏝️ Dreven en parken clusteren (%dm netwerk) en filteren op >= %.2f ha...", as.integer(zoekafstand_m), min_opp_ha))
  
  res_cluster <- cluster_filter_universeel(
    masker            = r_kle_tot,
    opp_laag          = r_kle_tot,
    min_opp_ha        = min_opp_ha,
    min_totale_opp_ha = NA,
    afstand_m         = zoekafstand_m,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  
  final_opp <- res_cluster$raster
  cl_opp    <- res_cluster$clusters
  rm(res_cluster, r_kle_tot)
  gc()
  
  ha_finaal <- calc_ha_exact(final_opp)
  message(sprintf("       -> Definitief leefgebied Juchtleerkever: %.2f ha", ha_finaal))
  
  # Deellagen isoleren binnen goedgekeurde netwerken
  if (!all(is.na(suppressWarnings(terra::minmax(cl_opp))))) {
    bwk_uitvoer   <- terra::mask(r_bwk_raw, cl_opp)
    groen_uitvoer <- terra::mask(r_kle_groen, cl_opp)
  } else {
    bwk_uitvoer   <- master_grid * NA
    groen_uitvoer <- master_grid * NA
  }
  rm(r_bwk_raw, r_kle_groen)
  gc()
  
  return(list(
    final_opp             = final_opp,
    cl_opp                = cl_opp,
    bwk_parken_knotbomen  = bwk_uitvoer,
    groenkaart_bomenrijen = groen_uitvoer
  ))
}