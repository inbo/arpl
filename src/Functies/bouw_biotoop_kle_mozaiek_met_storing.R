# ==============================================================================
# bouw_biotoop_kle_mozaiek_met_storing.R
# Archetype: Halfopen cultuurlandschappen met storingsbuffers, KLE-eis en matrixfilter
# ==============================================================================

bouw_biotoop_kle_mozaiek_met_storing <- function(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
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
  
  r_vp_raw   <- if ("voortplanting_bwk" %in% names(lagen)) lagen[["voortplanting_bwk"]] else (master_grid * NA)
  r_foer_raw <- if ("foerageerbiotoop"  %in% names(lagen)) lagen[["foerageerbiotoop"]]  else (master_grid * NA)
  
  # ----------------------------------------------------------------------------
  # 2. NEGATIEVE STOORFACTOREN (BEBOUWING >= 10 HA & GROTE WEGEN: 100M BUFFER)
  # ----------------------------------------------------------------------------
  message("    🚫 Negatieve storingsbuffers (bebouwing & grote wegen 100m) berekenen...")
  
  # A. Bebouwing >= 10 ha via Parquet
  f_cache_bebouw <- file.path(paths$dir_cache, "masker_bebouwing_groot_100m.tif")
  
  if (file.exists(f_cache_bebouw)) {
    message("       -> [CACHE] Bebouwingsbuffer ingeladen van schijf...")
    masker_bebouw_stoor <- terra::rast(f_cache_bebouw)
  } else {
    message("       -> GDX-bebouwingscodes ophalen en clusteren op >= 10 ha...")
    gdx_codes <- c("u", "ua", "ua-", "uc", "ud", "ui", "ui-", "un", "un-", "ur", "ur-", "ur+", "uv", "uv-", "uv+")
    t_bebouw  <- ds_parquet %>%
      dplyr::filter(tolower(CODE) %in% gdx_codes) %>%
      dplyr::collect() %>%
      data.table::as.data.table()
    
    if (nrow(t_bebouw) > 0) {
      t_bebouw[, CODE := tolower(trimws(CODE))]
      t_bebouw_cel <- unique(t_bebouw, by = c("cel_id", "CODE"))[, .(Val = 1), by = cel_id]
      r_bebouw_bin <- master_grid * NA
      r_bebouw_bin[t_bebouw_cel$cel_id] <- 1
      
      # Clusteren op 10 ha (25m zoekafstand om dorpskernen te sluiten)
      res_bebouw <- cluster_filter_universeel(
        masker            = r_bebouw_bin,
        opp_laag          = r_bebouw_bin,
        min_opp_ha        = 10,
        min_totale_opp_ha = NA,
        afstand_m         = 25,
        dispersie_m       = 0,
        werkelijk         = FALSE
      )
      
      r_bebouw_groot <- res_bebouw$raster
      if (!all(is.na(suppressWarnings(terra::minmax(r_bebouw_groot))))) {
        f_b_bin  <- tempfile(fileext = ".tif")
        terra::writeRaster(r_bebouw_groot, f_b_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 0)
        f_prox_b <- hulp_proximity_gdal(f_b_bin, 100)
        unlink(f_b_bin)
        
        p_b <- terra::rast(f_prox_b)
        masker_bebouw_stoor <- terra::ifel(!is.na(p_b) & p_b >= 0 & p_b <= 100, 1, NA)
        unlink(f_prox_b)
        rm(p_b)
      } else {
        masker_bebouw_stoor <- master_grid * NA
      }
      rm(res_bebouw, r_bebouw_groot, r_bebouw_bin)
    } else {
      masker_bebouw_stoor <- master_grid * NA
    }
    terra::writeRaster(masker_bebouw_stoor, f_cache_bebouw, overwrite = TRUE, datatype = "INT1U", gdal = "COMPRESS=LZW")
    rm(t_bebouw)
    gc()
  }
  
  # B. Wegenstoring (100m)
  masker_weg_stoor <- master_grid * NA
  if (file.exists(paths$grb_wegen)) {
    f_cache_wegen <- file.path(paths$dir_cache, "masker_wegen_stoor_100m.tif")
    if (file.exists(f_cache_wegen)) {
      message("       -> [CACHE] Wegenbuffer ingeladen van schijf...")
      masker_weg_stoor <- terra::rast(f_cache_wegen)
    } else {
      message("       -> Snelwegen en gescheiden rijbanen bufferen met 100m...")
      f_wegen_aligned <- hulp_align_gdal(paths$grb_wegen, master_grid)
      r_wegen         <- terra::rast(f_wegen_aligned)
      
      weg_cats        <- terra::cats(r_wegen)[[1]]
      gewenste_labels <- c("autosnelweg", "weg met gescheiden rijbanen die geen autosnelweg is")
      ids_wegen       <- weg_cats$value[tolower(trimws(weg_cats$Label)) %in% gewenste_labels]
      if (length(ids_wegen) == 0) {
        ids_wegen <- weg_cats$value[grepl("autosnelweg|gescheiden rijbanen", tolower(weg_cats$Label))]
      }
      
      r_wegen_bin <- terra::ifel(!is.na(r_wegen) & r_wegen %in% ids_wegen, 1, NA)
      unlink(f_wegen_aligned)
      rm(r_wegen)
      
      if (!all(is.na(suppressWarnings(terra::minmax(r_wegen_bin))))) {
        f_w_bin  <- tempfile(fileext = ".tif")
        terra::writeRaster(r_wegen_bin, f_w_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 0)
        f_prox_w <- hulp_proximity_gdal(f_w_bin, 100)
        unlink(f_w_bin)
        
        p_w <- terra::rast(f_prox_w)
        masker_weg_stoor <- terra::ifel(!is.na(p_w) & p_w >= 0 & p_w <= 100, 1, NA)
        unlink(f_prox_w)
        rm(p_w)
      }
      terra::writeRaster(masker_weg_stoor, f_cache_wegen, overwrite = TRUE, datatype = "INT1U", gdal = "COMPRESS=LZW")
      rm(r_wegen_bin)
      gc()
    }
  }
  
  # Stoorfactoren toepassen op voortplantingsbiotoop
  stoor_comb <- terra::cover(masker_bebouw_stoor, masker_weg_stoor)
  r_vp_gezuiverd <- terra::mask(r_vp_raw, stoor_comb, inverse = TRUE)
  rm(masker_bebouw_stoor, masker_weg_stoor, stoor_comb)
  gc()
  
  ha_vp_gezuiverd <- calc_ha_exact(r_vp_gezuiverd)
  message(sprintf("       -> Voortplantingsbiotoop na stoorfactoren: %.2f ha", ha_vp_gezuiverd))
  
  # ----------------------------------------------------------------------------
  # 3. KERNBIOTOOP CLUSTEREN (>= 10 HA, 50M) & LINEAIRE KLE-TOETS (>= 500M)
  # ----------------------------------------------------------------------------
  message("    🌳 Kernbiotoop clusteren (>= 10 ha, 50m netwerk)...")
  res_kern <- cluster_filter_universeel(
    masker            = r_vp_gezuiverd,
    opp_laag          = r_vp_gezuiverd,
    min_opp_ha        = 10,
    min_totale_opp_ha = NA,
    afstand_m         = 50,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_kern_10ha <- res_kern$raster
  cl_kern     <- res_kern$clusters
  rm(res_kern, r_vp_gezuiverd)
  gc()
  
  # KLE-eis via Groenkaart 2021 (waarde 1 en 2 = bomenrijen, houtkanten, heggen)
  message("    🌿 Lineaire KLE-eis toetsen (minstens 500m / 50 cellen bomen/struwelen per kern)...")
  if (!all(is.na(suppressWarnings(terra::minmax(cl_kern))))) {
    f_groen_aligned <- hulp_align_gdal(paths$groenkaart, master_grid)
    r_groen         <- terra::rast(f_groen_aligned)
    r_kle_bin       <- terra::ifel(!is.na(r_groen) & (r_groen == 1 | r_groen == 2), 1, 0)
    unlink(f_groen_aligned)
    rm(r_groen)
    
    # Tel KLE cellen per cluster ID via data.table
    cells_k <- terra::cells(!is.na(cl_kern) & cl_kern > 0)
    dt_kle  <- data.table::data.table(
      ID  = as.integer(cl_kern[cells_k][[1]]),
      KLE = as.integer(r_kle_bin[cells_k][[1]])
    )
    rm(cells_k, r_kle_bin)
    gc()
    
    ok_kle_ids <- dt_kle[!is.na(ID) & ID > 0, .(KLE_pixels = sum(KLE, na.rm = TRUE)), by = ID][KLE_pixels >= 50, ID]
    rm(dt_kle)
    gc()
    
    if (length(ok_kle_ids) > 0) {
      cl_kern_ok   <- terra::ifel(cl_kern %in% as.integer(ok_kle_ids), cl_kern, NA)
      r_kern_finaal <- terra::mask(r_kern_10ha, cl_kern_ok)
      rm(cl_kern_ok)
    } else {
      r_kern_finaal <- master_grid * NA
    }
  } else {
    r_kern_finaal <- master_grid * NA
  }
  rm(cl_kern, r_kern_10ha)
  gc()
  
  ha_kern_finaal <- calc_ha_exact(r_kern_finaal)
  message(sprintf("       -> Kernbiotoop met voldoende KLE (>= 500m): %.2f ha", ha_kern_finaal))
  
  # ----------------------------------------------------------------------------
  # 4. FOERAGEERKOPPELING (100M BUFFER) & TOTAAL LEFGEBIED (>= 15 HA)
  # ----------------------------------------------------------------------------
  message("    🌾 Foerageerhabitat koppelen binnen 100m van goedgekeurde kernen...")
  r_kern_bin <- terra::ifel(!is.na(r_kern_finaal) & r_kern_finaal > 0, 1, NA)
  
  if (!all(is.na(suppressWarnings(terra::minmax(r_kern_bin))))) {
    f_k_bin   <- tempfile(fileext = ".tif")
    terra::writeRaster(r_kern_bin, f_k_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 0)
    f_prox_k  <- hulp_proximity_gdal(f_k_bin, 100)
    unlink(f_k_bin)
    
    p_k        <- terra::rast(f_prox_k)
    mask_k100  <- terra::ifel(!is.na(p_k) & p_k >= 0 & p_k <= 100, 1, NA)
    unlink(f_prox_k)
    rm(p_k)
    gc()
    
    r_foer_gekoppeld <- terra::mask(r_foer_raw, mask_k100)
    rm(mask_k100)
    gc()
  } else {
    r_foer_gekoppeld <- master_grid * NA
  }
  rm(r_kern_bin, r_foer_raw)
  gc()
  
  # Totaal bruto leefgebied = Kern + Gekoppeld foerageerhabitat
  w_kern <- terra::ifel(is.na(r_kern_finaal), 0, r_kern_finaal)
  w_foer <- terra::ifel(is.na(r_foer_gekoppeld), 0, r_foer_gekoppeld)
  r_totaal_bruto <- terra::ifel((w_kern + w_foer) > 0, terra::clamp(w_kern + w_foer, upper = 1.0), NA)
  rm(w_kern, w_foer)
  gc()
  
  # Clusteren totaal complex op >= 15 ha (50m netwerk)
  message("    🧩 Totaal leefgebiedscomplex clusteren (>= 15 ha, 50m netwerk)...")
  res_totaal <- cluster_filter_universeel(
    masker            = r_totaal_bruto,
    opp_laag          = r_totaal_bruto,
    min_opp_ha        = 15,
    min_totale_opp_ha = NA,
    afstand_m         = 50,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_lg_15ha <- res_totaal$raster
  cl_lg     <- res_totaal$clusters
  rm(res_totaal, r_totaal_bruto)
  gc()
  
  # ----------------------------------------------------------------------------
  # 5. GRASLANDSTRUCTUURTEST (< 90% INTENSIEF PRODUCTIEGRASLAND HP)
  # ----------------------------------------------------------------------------
  message("    📊 Matrixfilter: Aandeel intensief cultuurgrasland (hp) toetsen (< 90%)...")
  if (!all(is.na(suppressWarnings(terra::minmax(cl_lg))))) {
    t_hp <- ds_parquet %>%
      dplyr::filter(stringr::str_detect(tolower(CODE), "^hp")) %>%
      dplyr::collect() %>%
      data.table::as.data.table()
    
    r_hp <- master_grid * NA
    if (nrow(t_hp) > 0) {
      t_hp[, CODE := tolower(trimws(CODE))]
      t_hp_cel <- unique(t_hp, by = c("cel_id", "CODE"))[, .(Opp_HP = pmin(sum(BWK_FRAC, na.rm = TRUE), 1.0)), by = cel_id]
      r_hp[t_hp_cel$cel_id] <- t_hp_cel$Opp_HP
      rm(t_hp_cel)
    }
    rm(t_hp)
    gc()
    
    cells_lg <- terra::cells(!is.na(cl_lg) & cl_lg > 0)
    dt_hp_check <- data.table::data.table(
      ID     = as.integer(cl_lg[cells_lg][[1]]),
      Tot_ha = as.numeric(r_lg_15ha[cells_lg][[1]]) * 0.01,
      HP_ha  = as.numeric(terra::ifel(is.na(r_hp[cells_lg][[1]]), 0, r_hp[cells_lg][[1]])) * 0.01
    )
    rm(cells_lg, r_hp)
    gc()
    
    dt_agg <- dt_hp_check[, .(
      Tot_ha = sum(Tot_ha, na.rm = TRUE),
      HP_ha  = sum(HP_ha, na.rm = TRUE)
    ), by = ID]
    dt_agg[, Fractie_HP := HP_ha / Tot_ha]
    
    ok_complex_ids <- dt_agg[!is.na(ID) & Fractie_HP < 0.90, ID]
    rm(dt_hp_check, dt_agg)
    gc()
    
    if (length(ok_complex_ids) > 0) {
      cl_finaal  <- terra::ifel(cl_lg %in% as.integer(ok_complex_ids), cl_lg, NA)
      final_opp  <- terra::mask(r_lg_15ha, cl_finaal)
      vp_uitvoer <- terra::mask(r_kern_finaal, cl_finaal)
      foer_uitvoer <- terra::mask(r_foer_gekoppeld, cl_finaal)
    } else {
      final_opp    <- master_grid * NA
      cl_finaal    <- master_grid * NA
      vp_uitvoer   <- master_grid * NA
      foer_uitvoer <- master_grid * NA
    }
  } else {
    final_opp    <- master_grid * NA
    cl_finaal    <- master_grid * NA
    vp_uitvoer   <- master_grid * NA
    foer_uitvoer <- master_grid * NA
  }
  rm(cl_lg, r_lg_15ha, r_kern_finaal, r_foer_gekoppeld)
  gc()
  
  return(list(
    final_opp     = final_opp,
    cl_opp        = cl_finaal,
    voortplanting = vp_uitvoer,
    foerageren    = foer_uitvoer
  ))
}
