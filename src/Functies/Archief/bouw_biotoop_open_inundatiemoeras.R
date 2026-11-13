# ==============================================================================
# bouw_biotoop_open_inundatiemoeras.R
# Archetype: Inundatiegraslanden & zeggenmoerassen x drainage x hooggroen 100m x meertrapscluster
# Geschikt voor: Porseleinhoen, Kleinst waterhoen, etc.
# ==============================================================================

bouw_biotoop_open_inundatiemoeras <- function(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
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
  
  r_gras_raw   <- if ("grasland" %in% names(lagen)) lagen[["grasland"]] else (master_grid * NA)
  r_moeras_raw <- if ("moeras"   %in% names(lagen)) lagen[["moeras"]]   else (master_grid * NA)
  
  # Terrils veilig uitsluiten via centrale Parquet helper (^kg codes)
  r_terril_bin <- vul_raster_vanuit_parquet(
    ds_parquet  = ds_parquet,
    codes_bevat = "kg",
    master_grid = master_grid,
    waarde_vast = 1
  )
  
  # ----------------------------------------------------------------------------
  # 2. SAMENVOEGING & DIRECTE UITSLUITING (HOOGGROEN & TERRILS)
  # ----------------------------------------------------------------------------
  message("    🌾 Grasland en moeras combineren, directe hooggroen- en terril-uitsluiting...")
  
  # Combinatie conform Bert (cover met clamp)
  r_bio_bruto <- terra::cover(r_gras_raw, r_moeras_raw)
  r_bio_bruto <- terra::clamp(r_bio_bruto, upper = 1.0)
  
  # Directe uitsluiting: uitsluitend klasse 1 van de Groenkaart (bomen/bos) en terrils
  f_groen_al     <- hulp_align_gdal(paths$groenkaart, master_grid)
  r_groen        <- terra::rast(f_groen_al)
  r_bomen_direct <- terra::ifel(!is.na(r_groen) & r_groen == 1, 1, NA)
  unlink(f_groen_al)
  rm(r_groen)
  
  # Maskeren met directe bomen en terrils
  r_bio_geschoond <- terra::mask(r_bio_bruto, r_bomen_direct, inverse = TRUE)
  if (!all(is.na(suppressWarnings(terra::minmax(r_terril_bin))))) {
    r_bio_geschoond <- terra::mask(r_bio_geschoond, r_terril_bin, inverse = TRUE)
  }
  rm(r_bio_bruto, r_terril_bin)
  gc()
  
  ha_geschoond <- calc_ha_exact(r_bio_geschoond)
  message(sprintf("       -> Biotoop na directe bomen- en terrilfilter: %.2f ha", ha_geschoond))
  
  # ----------------------------------------------------------------------------
  # 3. NATTE GRONDEN FILTER (BODEMDRAINAGEKAART)
  # ----------------------------------------------------------------------------
  message("    💧 Beperken tot natte drainageklassen (d, e, f, h, i, g, e-f, h-i, e-i)...")
  
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
  
  if (!all(is.na(suppressWarnings(terra::minmax(mask_drainage))))) {
    r_bio_nat <- terra::mask(r_bio_geschoond, mask_drainage)
  } else {
    r_bio_nat <- r_bio_geschoond
  }
  rm(r_bio_geschoond, mask_drainage)
  gc()
  
  ha_nat <- calc_ha_exact(r_bio_nat)
  message(sprintf("       -> Biotoop na drainagefilter: %.2f ha", ha_nat))
  
  # ----------------------------------------------------------------------------
  # 4. EERSTE CLUSTERING & FILTERING (>= 3 HA OP 50M NETWERK)
  # ----------------------------------------------------------------------------
  message("    🏝️ Eerste clustertoets uitvoeren (>= 3 ha, 50m netwerk)...")
  res_3ha <- cluster_filter_universeel(
    masker            = r_bio_nat,
    opp_laag          = r_bio_nat,
    min_opp_ha        = 3,
    min_totale_opp_ha = NA,
    afstand_m         = 50,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_bio_3ha <- res_3ha$raster
  rm(res_3ha, r_bio_nat)
  gc()
  
  ha_3ha <- calc_ha_exact(r_bio_3ha)
  message(sprintf("       -> Biotoop in clusters >= 3 ha: %.2f ha", ha_3ha))
  
  # ----------------------------------------------------------------------------
  # 5. OPEN ZICHT FILTER (100M EXCLUSIE ROND HOOGGROEN >= 20 ARE & HERKEURING 3 HA)
  # ----------------------------------------------------------------------------
  message("    🌲 Open zicht filter: 100m buffer rond hooggroenclusters (>= 20 are)...")
  
  f_cache_groen100 <- file.path(paths$dir_cache, "masker_groen_100m_cluster.tif")
  if (file.exists(f_cache_groen100)) {
    message("       -> [CACHE] Hooggroen-100m storingsbuffer ingeladen van schijf...")
    masker_groen_100m <- terra::rast(f_cache_groen100)
  } else {
    # 5x5 focal window op bomen (strikt klasse 1)
    r_groen_bin <- terra::ifel(!is.na(r_bomen_direct), 1, 0)
    w_5x5       <- matrix(1, nrow = 5, ncol = 5)
    f_focal_g   <- tempfile(fileext = ".tif")
    r_groen_som <- terra::focal(r_groen_bin, w = w_5x5, fun = "sum", na.rm = TRUE, filename = f_focal_g, overwrite = TRUE)
    
    # Drempel van Bert: >= 12 cellen in 5x5 venster en zelf boom zijn
    r_groen_20are <- terra::ifel(!is.na(r_groen_som) & r_groen_som >= 12 & !is.na(r_bomen_direct), 1, NA)
    unlink(f_focal_g)
    rm(r_groen_bin, r_groen_som)
    
    if (!all(is.na(suppressWarnings(terra::minmax(r_groen_20are))))) {
      f_g20_bin <- tempfile(fileext = ".tif")
      terra::writeRaster(r_groen_20are, f_g20_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
      f_prox_g  <- hulp_proximity_gdal(f_g20_bin, 100)
      unlink(f_g20_bin)
      
      p_g <- terra::rast(f_prox_g)
      masker_groen_100m <- terra::ifel(!is.na(p_g) & p_g >= 0 & p_g <= 100, 1, NA)
      unlink(f_prox_g)
      rm(p_g)
    } else {
      masker_groen_100m <- master_grid * NA
    }
    terra::writeRaster(masker_groen_100m, f_cache_groen100, overwrite = TRUE, datatype = "INT1U", gdal = "COMPRESS=LZW")
    rm(r_groen_20are)
    gc()
  }
  
  # Aftrek van de 100m buffer
  r_bio_open <- terra::mask(r_bio_3ha, masker_groen_100m, inverse = TRUE)
  rm(r_bio_3ha, masker_groen_100m, r_bomen_direct)
  gc()
  
  # Directe herkeuring op de 3 ha drempel na aftrek van de bosrandstoring
  res_3ha_post <- cluster_filter_universeel(
    masker            = r_bio_open,
    opp_laag          = r_bio_open,
    min_opp_ha        = 3,
    min_totale_opp_ha = NA,
    afstand_m         = 50,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_bio_open_3ha <- res_3ha_post$raster
  rm(res_3ha_post, r_bio_open)
  gc()
  
  ha_open_3ha <- calc_ha_exact(r_bio_open_3ha)
  message(sprintf("       -> Biotoop buiten 100m bosrandbuffer (>= 3 ha): %.2f ha", ha_open_3ha))
  
  # ----------------------------------------------------------------------------
  # 6. FINALE CLUSTERING & EINDDREMPEL (>= 10 HA OP 50M NETWERK)
  # ----------------------------------------------------------------------------
  message("    🧩 Finale clustering op functioneel complex (>= 10 ha, 50m netwerk)...")
  res_10ha <- cluster_filter_universeel(
    masker            = r_bio_open_3ha,
    opp_laag          = r_bio_open_3ha,
    min_opp_ha        = 10,
    min_totale_opp_ha = NA,
    afstand_m         = 50,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  final_opp <- res_10ha$raster
  cl_opp    <- res_10ha$clusters
  rm(res_10ha, r_bio_open_3ha)
  gc()
  
  ha_finaal <- calc_ha_exact(final_opp)
  message(sprintf("       -> Definitief leefgebied Porseleinhoen: %.2f ha", ha_finaal))
  
  # Deellagen isoleren binnen de definitieve clusters
  if (!all(is.na(suppressWarnings(terra::minmax(cl_opp))))) {
    gras_uitvoer   <- terra::mask(r_gras_raw, cl_opp)
    moeras_uitvoer <- terra::mask(r_moeras_raw, cl_opp)
  } else {
    gras_uitvoer   <- master_grid * NA
    moeras_uitvoer <- master_grid * NA
  }
  rm(r_gras_raw, r_moeras_raw)
  gc()
  
  return(list(
    final_opp = final_opp,
    cl_opp    = cl_opp,
    grasland  = gras_uitvoer,
    moeras    = moeras_uitvoer
  ))
}
