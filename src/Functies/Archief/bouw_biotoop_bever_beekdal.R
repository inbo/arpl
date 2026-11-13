# ==============================================================================
# bouw_biotoop_bever_beekdal.R
# Archetype voor Europese bever: 50 ha beekdal/water x 5 ha houtige opslag (50m koppeling)
# Conform officiële INBO-standaard ANB2019_EuropeseBever2018.asc
# ==============================================================================

bouw_biotoop_bever_beekdal <- function(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
  bestandsnaam <- tolower(gsub(" ", "", soort_naam))
  soort_cache  <- file.path(paths$dir_cache, bestandsnaam)
  dir.create(soort_cache, showWarnings = FALSE, recursive = TRUE)
  
  # ----------------------------------------------------------------------------
  # 0. PARAMETERS OPHALEN UIT EXCEL
  # ----------------------------------------------------------------------------
  rij_alg  <- params_soort %>% dplyr::filter(tolower(trimws(Type)) == "leefgebied_bwk") %>% dplyr::slice(1)
  rij_hout <- params_soort %>% dplyr::filter(tolower(trimws(Type)) == "houtig") %>% dplyr::slice(1)
  
  min_opp_alg   <- as.numeric(ifelse(nrow(rij_alg) > 0 && !is.na(rij_alg$MinOpp_ha), rij_alg$MinOpp_ha, 50.0))
  min_opp_hout  <- as.numeric(ifelse(nrow(rij_hout) > 0 && !is.na(rij_hout$MinOpp_ha), rij_hout$MinOpp_ha, 5.0))
  zoekafstand_m <- as.numeric(ifelse(nrow(rij_alg) > 0 && !is.na(rij_alg$AfstandBiotopen_m), rij_alg$AfstandBiotopen_m, 50.0))
  koppel_dist_m <- 50.0 # Vaste wederzijdse koppeling van 50 m conform 2018-model
  
  message(sprintf("    ⚙️ Parameters: Beekdalcomplex >= %.1f ha | Hout >= %.1f ha | Zoekafstand = %d m | Koppeling = %d m",
                  min_opp_alg, min_opp_hout, as.integer(zoekafstand_m), as.integer(koppel_dist_m)))
  
  types_in_tabel <- unique(params_soort$Type)
  
  # ----------------------------------------------------------------------------
  # 1. PARQUET EXTRACTIE (Garandeer leefgebied_bwk en houtig voor europesebever)
  # ----------------------------------------------------------------------------
  types_te_laden <- c("leefgebied_bwk", "houtig")
  
  # Zorg dat de zoeknaam in de crosswalk matcht met "europesebever"
  cw_naam <- if (any(grepl("^europesebever$", tolower(df_bwk_codes$Soort)))) "europesebever" else soort_naam
  
  lagen <- list()
  for (tp in types_te_laden) {
    f_tif_cache <- file.path(soort_cache, paste0("laag_", tp, ".tif"))
    if (file.exists(f_tif_cache)) {
      message(paste0("    📦 [CACHE] Laag '", tp, "' ingeladen van schijf..."))
      lagen[[tp]] <- terra::rast(f_tif_cache)
    } else {
      message(paste0("    ⏳ Parquet-laag '", tp, "' extraheren voor ", cw_naam, "..."))
      r_tp <- extract_bwk_laag(cw_naam, tp, ds_parquet, df_bwk_codes, master_grid)
      if (is.null(r_tp)) r_tp <- master_grid * NA_real_
      terra::writeRaster(r_tp, f_tif_cache, overwrite = TRUE, datatype = "FLT4S", gdal = "COMPRESS=LZW")
      lagen[[tp]] <- terra::rast(f_tif_cache)
    }
  }
  
  r_alg_bwk  <- if ("leefgebied_bwk" %in% names(lagen)) lagen[["leefgebied_bwk"]] else (master_grid * NA_real_)
  r_hout_bwk <- if ("houtig"         %in% names(lagen)) lagen[["houtig"]]         else (master_grid * NA_real_)
  
  # ----------------------------------------------------------------------------
  # 2. HYDROGRAFISCHE WATERLOPEN EN MAAS KOPPELEN MET SPOOR A
  # ----------------------------------------------------------------------------
  message("    🌊 Hydrografische waterlopen (Huetzon: WGb, WKb, WGr, WKr, WZgr) toevoegen...")
  
  r_huet_ok <- master_grid * NA_real_
  if (!is.null(paths$huetzon) && file.exists(paths$huetzon)) {
    f_huet_al <- hulp_align_gdal(paths$huetzon, master_grid)
    r_huet    <- terra::rast(f_huet_al)
    
    # Exacte IDs conform de vertaalsleutel:
    # 1: Bg  (WGb)            | 2: BgK  (WGb Kempen)
    # 3: Bk1 (WKb)            | 4: Bk2  (WKb)
    # 5: BkK1 (WKb Kempen)    | 6: BkK2 (WKb Kempen)
    # 17: Rg (WGr)            | 18: Rk  (WKr)        | 19: Rzg (WZgr)
    valide_huet_ids <- c(1, 2, 3, 4, 5, 6, 17, 18, 19)
    
    # Extra controle via attribuuttabel indien aanwezig
    leg_huet <- terra::cats(r_huet)[[1]]
    if (!is.null(leg_huet)) {
      kol_val   <- colnames(leg_huet)[tolower(colnames(leg_huet)) %in% c("value", "id")][1]
      kol_label <- colnames(leg_huet)[tolower(colnames(leg_huet)) %in% c("label", "omschrijving", "huet", "naam")][1]
      if (is.na(kol_val)) kol_val <- 1
      
      if (!is.na(kol_label)) {
        labels_raw <- trimws(as.character(leg_huet[[kol_label]]))
        target_labels <- c("Bg", "BgK", "Bk1", "Bk2", "BkK1", "BkK2", "Rg", "Rk", "Rzg")
        ids_via_label <- leg_huet[[kol_val]][labels_raw %in% target_labels]
        if (length(ids_via_label) > 0) valide_huet_ids <- ids_via_label
      }
    }
    
    r_huet_ok <- terra::ifel(!is.na(r_huet) & (r_huet %in% valide_huet_ids), 1.0, NA_real_)
    unlink(f_huet_al)
    rm(r_huet)
  }
  
  # Optionele extra Maas-laag
  r_maas_ok <- master_grid * NA_real_
  if (!is.null(paths$maas) && file.exists(paths$maas)) {
    f_maas_al <- hulp_align_gdal(paths$maas, master_grid)
    r_maas    <- terra::rast(f_maas_al)
    r_maas_ok <- terra::ifel(!is.na(r_maas) & r_maas >= 1, 1.0, NA_real_)
    unlink(f_maas_al)
    rm(r_maas)
  }
  
  # Stilstaand water (plassen en vennen)
  r_wv_ok <- master_grid * NA_real_
  if (!is.null(paths$watervlakken_raster) && file.exists(paths$watervlakken_raster)) {
    f_wv_al <- hulp_align_gdal(paths$watervlakken_raster, master_grid)
    r_wv    <- terra::rast(f_wv_al)
    r_wv_ok <- terra::ifel(!is.na(r_wv) & r_wv > 0, 1.0, NA_real_)
    unlink(f_wv_al)
    rm(r_wv)
  }
  
  # Spoor A samenstellen: BWK beekdal + Huet beken/rivieren + Maas + Watervlakken
  r_alg_totaal <- terra::cover(r_alg_bwk, r_huet_ok)
  r_alg_totaal <- terra::cover(r_alg_totaal, r_maas_ok)
  r_alg_totaal <- terra::cover(r_alg_totaal, r_wv_ok)
  r_alg_totaal <- terra::clamp(r_alg_totaal, upper = 1.0)
  rm(r_huet_ok, r_maas_ok, r_wv_ok)
  gc()
  
  ha_alg_ruw  <- calc_ha_exact(r_alg_totaal)
  ha_hout_ruw <- calc_ha_exact(r_hout_bwk)
  message(sprintf("       -> Algemeen beekdal-/watercomplex vóór clustering: %.2f ha", ha_alg_ruw))
  message(sprintf("       -> Houtige opslag vóór clustering: %.2f ha", ha_hout_ruw))
  
  # ----------------------------------------------------------------------------
  # 3. ZELFSTANDIGE CLUSTERING (CONFORM 2018: BRUTO COMPLEXGROOTTE)
  # ----------------------------------------------------------------------------
  message(sprintf("    🏝️ Spoor A clusteren op >= %.1f ha (%dm zoekafstand)...", min_opp_alg, as.integer(zoekafstand_m)))
  res_alg <- cluster_filter_universeel(
    masker            = r_alg_totaal,
    opp_laag          = r_alg_totaal,
    min_opp_ha        = min_opp_alg,
    min_totale_opp_ha = NA,
    afstand_m         = zoekafstand_m,
    dispersie_m       = 0,
    werkelijk         = TRUE  # 2018-model toetst de totale mozaïekomvang
  )
  r_alg_50ha <- res_alg$raster
  rm(res_alg, r_alg_totaal)
  gc()
  
  message(sprintf("    🪵 Spoor B clusteren op >= %.1f ha (%dm zoekafstand)...", min_opp_hout, as.integer(zoekafstand_m)))
  res_hout <- cluster_filter_universeel(
    masker            = r_hout_bwk,
    opp_laag          = r_hout_bwk,
    min_opp_ha        = min_opp_hout,
    min_totale_opp_ha = NA,
    afstand_m         = zoekafstand_m,
    dispersie_m       = 0,
    werkelijk         = TRUE  # 2018-model toetst de houtige clusteromvang
  )
  r_hout_5ha <- res_hout$raster
  rm(res_hout, r_hout_bwk)
  gc()
  
  ha_alg_cl  <- calc_ha_exact(r_alg_50ha)
  ha_hout_cl <- calc_ha_exact(r_hout_5ha)
  message(sprintf("       -> Geclusterd water-/dalcomplex (>= %.0f ha): %.2f ha", min_opp_alg, ha_alg_cl))
  message(sprintf("       -> Geclusterde houtige opslag (>= %.0f ha): %.2f ha", min_opp_hout, ha_hout_cl))
  
  # ----------------------------------------------------------------------------
  # 4. WEDERZIJDSE KOPPELING WATERCOMPLEX <-> HOUTIGE OPSLAG (<= 50 METER)
  # ----------------------------------------------------------------------------
  message(sprintf("    🔗 Wederzijdse oeverkoppeling berekenen (<= %d m via GDAL proximity)...", as.integer(koppel_dist_m)))
  
  heeft_alg  <- isTRUE(ha_alg_cl > 0)
  heeft_hout <- isTRUE(ha_hout_cl > 0)
  
  if (heeft_alg && heeft_hout) {
    # A. Proximity tot Hout berekenen (max. 50 m)
    f_h_bin  <- tempfile(fileext = ".tif")
    terra::writeRaster(terra::ifel(!is.na(r_hout_5ha) & r_hout_5ha > 0, 1, NA), f_h_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    f_prox_h <- hulp_proximity_gdal(f_h_bin, koppel_dist_m)
    unlink(f_h_bin)
    
    p_h <- terra::rast(f_prox_h)
    mask_alg_ok <- terra::ifel(!is.na(p_h) & p_h >= 0 & p_h <= koppel_dist_m, 1, NA)
    unlink(f_prox_h)
    rm(p_h)
    
    r_alg_gekoppeld <- terra::mask(r_alg_50ha, mask_alg_ok)
    rm(mask_alg_ok, r_alg_50ha)
    gc()
    
    # B. Proximity tot Watercomplex berekenen (max. 50 m)
    f_a_bin  <- tempfile(fileext = ".tif")
    terra::writeRaster(terra::ifel(!is.na(r_alg_gekoppeld) & r_alg_gekoppeld > 0, 1, NA), f_a_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    f_prox_a <- hulp_proximity_gdal(f_a_bin, koppel_dist_m)
    unlink(f_a_bin)
    
    p_a <- terra::rast(f_prox_a)
    mask_hout_ok <- terra::ifel(!is.na(p_a) & p_a >= 0 & p_a <= koppel_dist_m, 1, NA)
    unlink(f_prox_a)
    rm(p_a)
    
    r_hout_gekoppeld <- terra::mask(r_hout_5ha, mask_hout_ok)
    rm(mask_hout_ok, r_hout_5ha)
    gc()
    
    # C. Samensmelten tot definitief leefgebied
    final_opp <- terra::cover(r_alg_gekoppeld, r_hout_gekoppeld)
    final_opp <- terra::clamp(final_opp, upper = 1.0)
    
    # D. Netwerk clustering via GRASS
    f_fin_bin <- tempfile(fileext = ".tif")
    f_cl_fin  <- tempfile(fileext = ".tif")
    terra::writeRaster(terra::ifel(!is.na(final_opp) & final_opp > 0, 1, NA), f_fin_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    
    e <- terra::ext(master_grid)
    qgisprocess::qgis_run_algorithm(
      "grass:r.clump", input = f_fin_bin, output = f_cl_fin, title = "cl_bever",
      "-d" = FALSE,
      GRASS_REGION_PARAMETER = sprintf("%f,%f,%f,%f", e$xmin, e$xmax, e$ymin, e$ymax),
      GRASS_REGION_CELLSIZE_PARAMETER = 10, .quiet = TRUE
    )
    unlink(f_fin_bin)
    
    cl_opp <- terra::rast(f_cl_fin)
    try(terra::set.cats(cl_opp, layer = 1, value = NULL), silent = TRUE)
    try(terra::coltab(cl_opp) <- NULL, silent = TRUE)
    cl_opp <- terra::mask(cl_opp, final_opp)
    unlink(f_cl_fin)
    
  } else {
    message("    ⚠️ Geen overlappend beekdal- of houtcomplex na koppeling.")
    final_opp <- master_grid * NA_real_
    cl_opp    <- master_grid * NA_real_
  }
  
  ha_finaal <- calc_ha_exact(final_opp)
  message(sprintf("       -> Definitief leefgebied Europese bever: %.2f ha", ha_finaal))
  
  return(list(
    final_opp = final_opp,
    cl_opp    = cl_opp
  ))
}
