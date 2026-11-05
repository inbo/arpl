# ==============================================================================
# bouw_biotoop_commuter_watermoeras.R
# Archetype voor reigerachtigen/moerasvogels (broedbos x moeras/water x 5km commuter)
# Geschikt voor: Kwak, Kleine zilverreiger, etc.
# ==============================================================================

bouw_biotoop_commuter_watermoeras <- function(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
  bestandsnaam <- tolower(gsub(" ", "", soort_naam))
  soort_cache  <- file.path(paths$dir_cache, bestandsnaam)
  dir.create(soort_cache, showWarnings = FALSE, recursive = TRUE)
  
  f_cache_rds <- file.path(soort_cache, "01_raw_parquet_layers.rds")
  
  # Parameters flexibel uitlezen uit de Excel-rijen
  rij_vp   <- params_soort %>% dplyr::filter(tolower(Type) %in% c("voortplanting", "voortplanting0", "broedbos")) %>% dplyr::slice(1)
  rij_foer <- params_soort %>% dplyr::filter(tolower(Type) %in% c("foerageer", "moeras", "water")) %>% dplyr::slice(1)
  
  min_ha_vp     <- as.numeric(ifelse(nrow(rij_vp) > 0 && !is.na(rij_vp$MinOpp_ha), rij_vp$MinOpp_ha, 5))
  dist_vp_m     <- as.numeric(ifelse(nrow(rij_vp) > 0 && !is.na(rij_vp$AfstandBiotopen_m), rij_vp$AfstandBiotopen_m, 50))
  
  min_ha_foer   <- as.numeric(ifelse(nrow(rij_foer) > 0 && !is.na(rij_foer$MinOpp_ha), rij_foer$MinOpp_ha, 10))
  dist_foer_m   <- as.numeric(ifelse(nrow(rij_foer) > 0 && !is.na(rij_foer$AfstandBiotopen_m), rij_foer$AfstandBiotopen_m, 50))
  pct_water_min <- as.numeric(ifelse(nrow(rij_foer) > 0 && !is.na(rij_foer$Min_Water_Pct), rij_foer$Min_Water_Pct, 20)) / 100
  
  commuter_dist_m <- as.numeric(ifelse(nrow(rij_vp) > 0 && !is.na(rij_vp$Wederzijdse_Afstand_m), rij_vp$Wederzijdse_Afstand_m, 5000))
  
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
  
  # Flexibele toewijzing van lagen ongeacht type-naamgeving
  r_vp_raw <- if ("voortplanting" %in% names(lagen)) {
    lagen[["voortplanting"]]
  } else if ("voortplanting0" %in% names(lagen)) {
    lagen[["voortplanting0"]]
  } else {
    master_grid * NA
  }
  
  r_not_eik <- if ("not_eikenbos" %in% names(lagen)) {
    lagen[["not_eikenbos"]]
  } else if ("eikenbos" %in% names(lagen)) {
    lagen[["eikenbos"]]
  } else {
    master_grid * NA
  }
  
  r_moeras_raw <- if ("foerageer" %in% names(lagen)) {
    lagen[["foerageer"]]
  } else if ("moeras" %in% names(lagen)) {
    lagen[["moeras"]]
  } else {
    master_grid * NA
  }
  
  r_water_raw <- if ("openwater" %in% names(lagen)) {
    lagen[["openwater"]]
  } else if ("water" %in% names(lagen)) {
    lagen[["water"]]
  } else {
    master_grid * NA
  }
  
  # ----------------------------------------------------------------------------
  # 2. SPOOR A: BROEDBIOTOOP (UITSLUITING DROGE BOSSEN + CLUSTERING)
  # ----------------------------------------------------------------------------
  message(sprintf("    🌳 Spoor A: Broedbiotoop zuiveren en clusteren (>= %.1f ha op %dm)...", min_ha_vp, as.integer(dist_vp_m)))
  
  if (!all(is.na(suppressWarnings(terra::minmax(r_not_eik))))) {
    r_vp_gezuiverd <- terra::mask(r_vp_raw, r_not_eik > 0, inverse = TRUE)
  } else {
    r_vp_gezuiverd <- r_vp_raw
  }
  rm(r_vp_raw, r_not_eik)
  gc()
  
  res_vp <- cluster_filter_universeel(
    masker            = r_vp_gezuiverd,
    opp_laag          = r_vp_gezuiverd,
    min_opp_ha        = min_ha_vp,
    min_totale_opp_ha = NA,
    afstand_m         = dist_vp_m,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_vp_ok <- res_vp$raster
  rm(res_vp, r_vp_gezuiverd)
  gc()
  
  ha_vp_ok <- calc_ha_exact(r_vp_ok)
  message(sprintf("       -> Broedkernen goedgekeurd (>= %.1f ha): %.2f ha", min_ha_vp, ha_vp_ok))
  
  # ----------------------------------------------------------------------------
  # 3. SPOOR B: FOERAGEERBIOTOOP (MOERAS + WATERLOPEN + WATERAANDEEL TOETS)
  # ----------------------------------------------------------------------------
  message(sprintf("    🐟 Spoor B: Foerageergebied clusteren (>= %.1f ha op %dm) en wateraandeel (>= %d%%) toetsen...", 
                  min_ha_foer, as.integer(dist_foer_m), as.integer(pct_water_min * 100)))
  
  # Huet-zones integreren indien van toepassing
  r_huet_fractie <- master_grid * NA
  if (nrow(rij_foer) > 0 && isTRUE(as.logical(rij_foer$Gebruik_Huetzones)) && file.exists(paths$huetzon)) {
    f_hz_al <- hulp_align_gdal(paths$huetzon, master_grid)
    r_hz    <- terra::rast(f_hz_al)
    
    huet_codes <- c(1:6, 9, 11:19)
    if (!is.na(rij_foer$Huet_Labels) && nzchar(trimws(rij_foer$Huet_Labels))) {
      huet_codes <- as.numeric(trimws(unlist(strsplit(as.character(rij_foer$Huet_Labels), ","))))
    }
    r_huet_fractie <- terra::ifel(!is.na(r_hz) & r_hz %in% huet_codes, 0.01, NA)
    unlink(f_hz_al)
    rm(r_hz)
  }
  
  # Moeras / Foerageerbasis samenvoegen met waterlopen
  m_clean <- terra::ifel(is.na(r_moeras_raw), 0, r_moeras_raw)
  h_clean <- terra::ifel(is.na(r_huet_fractie), 0, r_huet_fractie)
  r_foer_totaal <- terra::ifel((m_clean + h_clean) > 0, terra::clamp(m_clean + h_clean, upper = 1.0), NA)
  rm(m_clean, h_clean, r_moeras_raw)
  
  # Zuiver open watergrid voor waterfractietoets
  w_clean <- terra::ifel(is.na(r_water_raw), 0, r_water_raw)
  hw_clean <- terra::ifel(is.na(r_huet_fractie), 0, r_huet_fractie)
  r_water_puur <- terra::ifel((w_clean + hw_clean) > 0, terra::clamp(w_clean + hw_clean, upper = 1.0), NA)
  rm(w_clean, hw_clean, r_water_raw, r_huet_fractie)
  gc()
  
  # Clusteren van foerageercomplex
  res_foer <- cluster_filter_universeel(
    masker            = r_foer_totaal,
    opp_laag          = r_foer_totaal,
    min_opp_ha        = min_ha_foer,
    min_totale_opp_ha = NA,
    afstand_m         = dist_foer_m,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_foer_drempel <- res_foer$raster
  cl_foer        <- res_foer$clusters
  rm(res_foer, r_foer_totaal)
  gc()
  
  # Wateraandeel-toets (indien pct_water_min > 0)
  if (!all(is.na(suppressWarnings(terra::minmax(cl_foer)))) && pct_water_min > 0) {
    cells_f <- terra::cells(!is.na(cl_foer) & cl_foer > 0)
    dt_water_toets <- data.table::data.table(
      ID         = as.integer(cl_foer[cells_f][[1]]),
      Moeras_Val = as.numeric(r_foer_drempel[cells_f][[1]]),
      Water_Val  = as.numeric(terra::ifel(is.na(r_water_puur[cells_f][[1]]), 0, r_water_puur[cells_f][[1]]))
    )
    rm(cells_f)
    gc()
    
    dt_agg <- dt_water_toets[, .(
      Moeras_ha = sum(Moeras_Val, na.rm = TRUE) * 0.01,
      Water_ha  = sum(Water_Val,  na.rm = TRUE) * 0.01
    ), by = ID]
    dt_agg[, Water_Aandeel := Water_ha / Moeras_ha]
    
    ok_foer_ids <- dt_agg[!is.na(ID) & (Water_Aandeel >= pct_water_min | Water_ha >= (min_ha_foer * pct_water_min)), ID]
    rm(dt_water_toets, dt_agg)
    gc()
    
    if (length(ok_foer_ids) > 0) {
      cl_foer_ok   <- terra::ifel(cl_foer %in% as.integer(ok_foer_ids), cl_foer, NA)
      r_foer_geflt <- terra::mask(r_foer_drempel, cl_foer_ok)
      rm(cl_foer_ok)
    } else {
      r_foer_geflt <- master_grid * NA
    }
  } else {
    r_foer_geflt <- r_foer_drempel
  }
  rm(cl_foer, r_foer_drempel, r_water_puur)
  gc()
  
  ha_foer_ok <- calc_ha_exact(r_foer_geflt)
  message(sprintf("       -> Foerageercomplex met voldoende open water: %.2f ha", ha_foer_ok))
  
  # ----------------------------------------------------------------------------
  # 4. WEDERZIJDSE COMMUTER-KOPPELING (GDAL PROXIMITY)
  # ----------------------------------------------------------------------------
  message(sprintf("    ↔️ Wederzijdse %dm commuter-toets uitvoeren...", as.integer(commuter_dist_m)))
  
  heeft_vp   <- !all(is.na(suppressWarnings(terra::minmax(r_vp_ok))))
  heeft_foer <- !all(is.na(suppressWarnings(terra::minmax(r_foer_geflt))))
  
  if (heeft_vp && heeft_foer) {
    r_vp_bin   <- terra::ifel(!is.na(r_vp_ok) & r_vp_ok > 0, 1, NA)
    r_foer_bin <- terra::ifel(!is.na(r_foer_geflt) & r_foer_geflt > 0, 1, NA)
    
    # Proximity vanaf foerageercomplex
    f_foer_bin  <- tempfile(fileext = ".tif")
    terra::writeRaster(r_foer_bin, f_foer_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    f_prox_foer <- hulp_proximity_gdal(f_foer_bin, commuter_dist_m)
    unlink(f_foer_bin)
    
    p_foer      <- terra::rast(f_prox_foer)
    mask_foer5k <- terra::ifel(!is.na(p_foer) & p_foer >= 0 & p_foer <= commuter_dist_m, 1, NA)
    unlink(f_prox_foer)
    rm(p_foer)
    
    # Proximity vanaf broedplaats
    f_vp_bin  <- tempfile(fileext = ".tif")
    terra::writeRaster(r_vp_bin, f_vp_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    f_prox_vp <- hulp_proximity_gdal(f_vp_bin, commuter_dist_m)
    unlink(f_vp_bin)
    
    p_vp      <- terra::rast(f_prox_vp)
    mask_vp5k <- terra::ifel(!is.na(p_vp) & p_vp >= 0 & p_vp <= commuter_dist_m, 1, NA)
    unlink(f_prox_vp)
    rm(p_vp)
    
    # Wederzijds maskeren
    vp_gekoppeld   <- terra::mask(r_vp_ok, mask_foer5k)
    foer_gekoppeld <- terra::mask(r_foer_geflt, mask_vp5k)
    
    unlink(c(mask_foer5k, mask_vp5k))
    rm(r_vp_bin, r_foer_bin)
    gc()
  } else {
    vp_gekoppeld   <- master_grid * NA
    foer_gekoppeld <- master_grid * NA
  }
  rm(r_vp_ok, r_foer_geflt)
  gc()
  
  # ----------------------------------------------------------------------------
  # 5. SAMENVOEGEN TOT DEFINITIEF LEEFGEBIED & OVERKOEPELENDE CLUSTERING
  # ----------------------------------------------------------------------------
  min_ha_finaal <- max(min_ha_vp, min_ha_foer)
  message(sprintf("    🧩 Leefgebied samenvoegen en finaal clusteren (>= %.1f ha, 50m netwerk)...", min_ha_finaal))
  
  w_v <- terra::ifel(is.na(vp_gekoppeld), 0, vp_gekoppeld)
  w_f <- terra::ifel(is.na(foer_gekoppeld), 0, foer_gekoppeld)
  r_leef_bruto <- terra::ifel((w_v + w_f) > 0, terra::clamp(w_v + w_f, upper = 1.0), NA)
  rm(w_v, w_f)
  gc()
  
  res_finaal <- cluster_filter_universeel(
    masker            = r_leef_bruto,
    opp_laag          = r_leef_bruto,
    min_opp_ha        = min_ha_finaal,
    min_totale_opp_ha = NA,
    afstand_m         = 50,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  final_opp <- res_finaal$raster
  cl_opp    <- res_finaal$clusters
  rm(res_finaal, r_leef_bruto)
  gc()
  
  vp_uitvoer   <- terra::mask(vp_gekoppeld, cl_opp)
  foer_uitvoer <- terra::mask(foer_gekoppeld, cl_opp)
  rm(vp_gekoppeld, foer_gekoppeld)
  gc()
  
  return(list(
    final_opp     = final_opp,
    cl_opp        = cl_opp,
    voortplanting = vp_uitvoer,
    foerageren    = foer_uitvoer
  ))
}
