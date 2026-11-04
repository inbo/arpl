# ==============================================================================
# bouw_biotoop_kleinschalig_cultuurland_driehoek.R
# Archetype voor akkervogels & kleinschalig cultuurland (neststruweel x 300m voedsel x 500m water)
# Geschikt voor: Zomertortel, Grauwe klauwier, etc.
# ==============================================================================

bouw_biotoop_kleinschalig_cultuurland_driehoek <- function(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
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
  
  r_vp_raw    <- if ("voortplanting" %in% names(lagen)) lagen[["voortplanting"]] else (master_grid * NA)
  r_foer_raw  <- if ("foerageer"     %in% names(lagen)) lagen[["foerageer"]]     else (master_grid * NA)
  r_water_raw <- if ("water"         %in% names(lagen)) lagen[["water"]]         else (master_grid * NA)
  
  # ----------------------------------------------------------------------------
  # 2. SAMENSTELLEN EXTERNE COMPONENTEN (LGP LANDBOUW & WATER)
  # ----------------------------------------------------------------------------
  message("    🌾 LGP-landbouwpercelen en waterbronnen integreren...")
  
  # A. Foerageerhabitat: BWK + LGP teelten (gewascodes 201, 202, 311, 321, 601, 800)
  if (file.exists(paths$lgp_landbouw)) {
    f_lgp_al <- hulp_align_gdal(paths$lgp_landbouw, master_grid)
    r_lgp    <- terra::rast(f_lgp_al)
    zomer_codes <- c(201, 202, 311, 321, 601, 800)
    r_lgp_fractie <- terra::ifel(!is.na(r_lgp) & r_lgp %in% zomer_codes, 0.01, NA)
    unlink(f_lgp_al)
    rm(r_lgp)
  } else {
    r_lgp_fractie <- master_grid * NA
  }
  
  f_bwk_c <- terra::ifel(is.na(r_foer_raw), 0, r_foer_raw)
  f_lgp_c <- terra::ifel(is.na(r_lgp_fractie), 0, r_lgp_fractie)
  r_foer_totaal <- terra::ifel((f_bwk_c + f_lgp_c) > 0, terra::clamp(f_bwk_c + f_lgp_c, upper = 1.0), NA)
  rm(f_bwk_c, f_lgp_c, r_lgp_fractie, r_foer_raw)
  gc()
  
  # B. Drinkwater: BWK + Watervlakken 2024 + Huet-zones
  r_watervlak_fractie <- master_grid * NA
  if (file.exists(paths$watervlakken_raster)) {
    f_wv_al <- hulp_align_gdal(paths$watervlakken_raster, master_grid)
    r_wv    <- terra::rast(f_wv_al)
    r_watervlak_fractie <- terra::ifel(!is.na(r_wv) & r_wv > 0, 0.01, NA)
    unlink(f_wv_al)
    rm(r_wv)
  }
  
  r_huet_fractie <- master_grid * NA
  if (file.exists(paths$huetzon)) {
    f_hz_al <- hulp_align_gdal(paths$huetzon, master_grid)
    r_hz    <- terra::rast(f_hz_al)
    r_huet_fractie <- terra::ifel(!is.na(r_hz) & r_hz > 0, 0.01, NA)
    unlink(f_hz_al)
    rm(r_hz)
  }
  
  w_bwk_c  <- terra::ifel(is.na(r_water_raw), 0, r_water_raw)
  w_vlak_c <- terra::ifel(is.na(r_watervlak_fractie), 0, r_watervlak_fractie)
  w_huet_c <- terra::ifel(is.na(r_huet_fractie), 0, r_huet_fractie)
  r_water_totaal <- terra::ifel((w_bwk_c + w_vlak_c + w_huet_c) > 0, terra::clamp(w_bwk_c + w_vlak_c + w_huet_c, upper = 1.0), NA)
  rm(w_bwk_c, w_vlak_c, w_huet_c, r_watervlak_fractie, r_huet_fractie, r_water_raw)
  gc()
  
  # ----------------------------------------------------------------------------
  # 3. NEGATIEVE FILTERS OP BROEDBIOTOOP (GESLOTEN BOS >= 20 ARE & BEBOUWING >= 10 HA)
  # ----------------------------------------------------------------------------
  message("    🚫 Negatieve filters toepassen op neststruweel (bos > 20 are en bebouwing > 10 ha)...")
  
  # A. Gesloten boscluster >= 20 are (5x5 focal window met >= 12 cellen)
  f_groen_al  <- hulp_align_gdal(paths$groenkaart, master_grid)
  r_groen_raw <- terra::rast(f_groen_al)
  r_groen_bin <- terra::ifel(!is.na(r_groen_raw) & (r_groen_raw == 1 | r_groen_raw > 0), 1, NA)
  unlink(f_groen_al)
  rm(r_groen_raw)
  
  w_5x5 <- matrix(1, nrow = 5, ncol = 5)
  f_focal_g <- tempfile(fileext = ".tif")
  r_groen_som <- terra::focal(r_groen_bin, w = w_5x5, fun = "sum", na.rm = TRUE, filename = f_focal_g, overwrite = TRUE)
  r_bos_20are <- terra::ifel(!is.na(r_groen_bin) & r_groen_som >= 12, 1, NA)
  unlink(f_focal_g)
  rm(r_groen_bin, r_groen_som)
  gc()
  
  # B. Bebouwingsbuffer (reeds gecached of berekenen)
  f_cache_bebouw <- file.path(paths$dir_cache, "masker_bebouwing_groot_100m.tif")
  if (file.exists(f_cache_bebouw)) {
    masker_bebouw_100m <- terra::rast(f_cache_bebouw)
  } else {
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
      
      res_bebouw <- cluster_filter_universeel(
        masker = r_bebouw_bin, opp_laag = r_bebouw_bin, min_opp_ha = 10,
        min_totale_opp_ha = NA, afstand_m = 25, dispersie_m = 0, werkelijk = FALSE
      )
      r_bebouw_groot <- res_bebouw$raster
      if (!all(is.na(suppressWarnings(terra::minmax(r_bebouw_groot))))) {
        f_b_bin  <- tempfile(fileext = ".tif")
        terra::writeRaster(r_bebouw_groot, f_b_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 0)
        f_prox_b <- hulp_proximity_gdal(f_b_bin, 100)
        unlink(f_b_bin)
        
        p_b <- terra::rast(f_prox_b)
        masker_bebouw_100m <- terra::ifel(!is.na(p_b) & p_b >= 0 & p_b <= 100, 1, NA)
        unlink(f_prox_b)
        rm(p_b)
      } else {
        masker_bebouw_100m <- master_grid * NA
      }
      rm(res_bebouw, r_bebouw_groot, r_bebouw_bin)
    } else {
      masker_bebouw_100m <- master_grid * NA
    }
    terra::writeRaster(masker_bebouw_100m, f_cache_bebouw, overwrite = TRUE, datatype = "INT1U", gdal = "COMPRESS=LZW")
    rm(t_bebouw)
    gc()
  }
  
  # Neststruweel schonen
  stoor_nest <- terra::cover(r_bos_20are, masker_bebouw_100m)
  r_vp_gezuiverd <- terra::mask(r_vp_raw, stoor_nest, inverse = TRUE)
  rm(r_bos_20are, masker_bebouw_100m, stoor_nest, r_vp_raw)
  gc()
  
  # ----------------------------------------------------------------------------
  # 4. CLUSTERING BASISBIOTOPEN
  # ----------------------------------------------------------------------------
  message("    🏝️ Basisbiotopen filteren op minimale drempels...")
  
  # A. Neststruweel: >= 300 m2 (0.03 ha) op 50m
  res_vp <- cluster_filter_universeel(
    masker            = r_vp_gezuiverd,
    opp_laag          = r_vp_gezuiverd,
    min_opp_ha        = 0.03,
    min_totale_opp_ha = NA,
    afstand_m         = 50,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_vp_300m2 <- res_vp$raster
  rm(res_vp, r_vp_gezuiverd)
  gc()
  
  # B. Foerageerhabitat: >= 1 ha op 50m
  res_foer <- cluster_filter_universeel(
    masker            = r_foer_totaal,
    opp_laag          = r_foer_totaal,
    min_opp_ha        = 1.0,
    min_totale_opp_ha = NA,
    afstand_m         = 50,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_foer_1ha <- res_foer$raster
  rm(res_foer, r_foer_totaal)
  gc()
  
  # ----------------------------------------------------------------------------
  # 5. RUIMTELIJKE DRIEHOEK: 300M VOEDSEL EN 500M DRINKWATER
  # ----------------------------------------------------------------------------
  message("    📐 Ruimtelijke driehoek berekenen (300m naar voedsel, 500m naar drinkwater)...")
  
  heeft_vp    <- !all(is.na(suppressWarnings(terra::minmax(r_vp_300m2))))
  heeft_foer  <- !all(is.na(suppressWarnings(terra::minmax(r_foer_1ha))))
  heeft_water <- !all(is.na(suppressWarnings(terra::minmax(r_water_totaal))))
  
  if (heeft_vp && heeft_foer && heeft_water) {
    # 300m proximity vanaf voedselgronden
    f_f_bin <- tempfile(fileext = ".tif")
    terra::writeRaster(terra::ifel(!is.na(r_foer_1ha) & r_foer_1ha > 0, 1, NA), f_f_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 0)
    f_prox_f300 <- hulp_proximity_gdal(f_f_bin, 300)
    unlink(f_f_bin)
    
    p_f300 <- terra::rast(f_prox_f300)
    mask_foer300m <- terra::ifel(!is.na(p_f300) & p_f300 >= 0 & p_f300 <= 300, 1, NA)
    unlink(f_prox_f300)
    rm(p_f300)
    
    # 500m proximity vanaf drinkwater
    f_w_bin <- tempfile(fileext = ".tif")
    terra::writeRaster(terra::ifel(!is.na(r_water_totaal) & r_water_totaal > 0, 1, NA), f_w_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 0)
    f_prox_w500 <- hulp_proximity_gdal(f_w_bin, 500)
    unlink(f_w_bin)
    
    p_w500 <- terra::rast(f_prox_w500)
    mask_water500m <- terra::ifel(!is.na(p_w500) & p_w500 >= 0 & p_w500 <= 500, 1, NA)
    unlink(f_prox_w500)
    rm(p_w500)
    
    # De interactiezone waar neststruweel mag liggen
    mask_driehoek <- terra::mask(mask_foer300m, mask_water500m)
    
    # Goedgekeurd neststruweel
    vp_finaal <- terra::mask(r_vp_300m2, mask_driehoek)
    
    # Reikwijdte vanuit goedgekeurd neststruweel naar voedsel en water
    if (!all(is.na(suppressWarnings(terra::minmax(vp_finaal))))) {
      f_vp_ok_bin <- tempfile(fileext = ".tif")
      terra::writeRaster(terra::ifel(!is.na(vp_finaal) & vp_finaal > 0, 1, NA), f_vp_ok_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 0)
      
      # Buffer 300m rond neststruweel voor foerageerhabitat
      f_prox_vp300 <- hulp_proximity_gdal(f_vp_ok_bin, 300)
      p_vp300      <- terra::rast(f_prox_vp300)
      foer_finaal  <- terra::mask(r_foer_1ha, terra::ifel(!is.na(p_vp300) & p_vp300 <= 300, 1, NA))
      unlink(f_prox_vp300)
      rm(p_vp300)
      
      # Buffer 500m rond neststruweel voor drinkwater
      f_prox_vp500 <- hulp_proximity_gdal(f_vp_ok_bin, 500)
      unlink(f_vp_ok_bin)
      p_vp500      <- terra::rast(f_prox_vp500)
      water_finaal <- terra::mask(r_water_totaal, terra::ifel(!is.na(p_vp500) & p_vp500 <= 500, 1, NA))
      unlink(f_prox_vp500)
      rm(p_vp500)
    } else {
      foer_finaal  <- master_grid * NA
      water_finaal <- master_grid * NA
    }
    
    unlink(c(mask_foer300m, mask_water500m, mask_driehoek))
    gc()
  } else {
    vp_finaal    <- master_grid * NA
    foer_finaal  <- master_grid * NA
    water_finaal <- master_grid * NA
  }
  rm(r_vp_300m2, r_foer_1ha, r_water_totaal)
  gc()
  
  ha_vp_fin <- calc_ha_exact(vp_finaal)
  message(sprintf("       -> Goedgekeurd neststruweel in interactiezone: %.2f ha", ha_vp_fin))
  
  # ----------------------------------------------------------------------------
  # 6. SAMENSTELLEN DEFINITIEF COMPOSIET LEEFGEBIED & NETWERK-ID'S
  # ----------------------------------------------------------------------------
  message("    🧩 Composiet leefgebied samenstellen (neststruweel + foerageer + drinkwater)...")
  
  final_opp <- terra::cover(vp_finaal, terra::cover(foer_finaal, water_finaal))
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
    cl_opp <- template_Vlaanderen
  }
  rm(r_fin_bin)
  gc()
  
  return(list(
    final_opp     = final_opp,
    cl_opp        = cl_opp,
    voortplanting = vp_finaal,
    foerageren    = foer_finaal,
    drinkwater    = water_finaal
  ))
}
