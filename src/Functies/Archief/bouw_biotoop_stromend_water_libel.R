# ==============================================================================
# bouw_biotoop_stromend_water_libel.R
# Archetype voor stromend-water libellen (beken >= 0.4ha x bos >= 5ha in matrix >= 50ha x 500m interactie)
# Geschikt voor: Beekrombout (Gomphus vulgatissimus), etc.
# ==============================================================================

bouw_biotoop_stromend_water_libel <- function(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
  bestandsnaam <- tolower(gsub(" ", "", soort_naam))
  soort_cache  <- file.path(paths$dir_cache, bestandsnaam)
  dir.create(soort_cache, showWarnings = FALSE, recursive = TRUE)
  
  f_cache_rds <- file.path(soort_cache, "01_raw_parquet_layers.rds")
  
  # Parameters ophalen uit Excel
  rij_wat <- params_soort %>% dplyr::filter(grepl("water", tolower(Type))) %>% dplyr::slice(1)
  min_opp_water <- as.numeric(ifelse(nrow(rij_wat) > 0 && !is.na(rij_wat$MinOpp_ha), rij_wat$MinOpp_ha, 0.4))
  
  rij_bos <- params_soort %>% dplyr::filter(grepl("bos", tolower(Type))) %>% dplyr::slice(1)
  min_opp_bos <- as.numeric(ifelse(nrow(rij_bos) > 0 && !is.na(rij_bos$MinOpp_ha), rij_bos$MinOpp_ha, 5.0))
  
  rij_rest <- params_soort %>% dplyr::filter(grepl("rest", tolower(Type))) %>% dplyr::slice(1)
  min_opp_rest <- as.numeric(ifelse(nrow(rij_rest) > 0 && !is.na(rij_rest$MinOpp_ha), rij_rest$MinOpp_ha, 50.0))
  
  koppel_m <- as.numeric(ifelse(nrow(rij_wat) > 0 && !is.na(rij_wat$Dispersiecap_m), rij_wat$Dispersiecap_m, 500))
  
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
  
  r_wat_raw  <- if ("waterlopen1_bwk"     %in% names(lagen)) lagen[["waterlopen1_bwk"]]     else (master_grid * NA)
  r_bos_raw  <- if ("landbiotoop_bos_bwk"  %in% names(lagen)) lagen[["landbiotoop_bos_bwk"]]  else (master_grid * NA)
  r_rest_raw <- if ("landbiotoop_rest_bwk" %in% names(lagen)) lagen[["landbiotoop_rest_bwk"]] else (master_grid * NA)
  
  # ----------------------------------------------------------------------------
  # 2. HUETZONES WATERLOPEN INTEGREREN (Bg, BgK, Bk1, Bk2, BkK1, BkK2)
  # ----------------------------------------------------------------------------
  message("    🌊 Hydrografische beekzones (Bg, Bk) koppelen met BWK-water...")
  
  f_cache_huet_rombout <- file.path(paths$dir_cache, "masker_huetzon_beekrombout.tif")
  
  if (file.exists(f_cache_huet_rombout)) {
    r_huet_ok <- terra::rast(f_cache_huet_rombout)
  } else if (!is.null(paths$huetzon) && file.exists(paths$huetzon)) {
    f_huet_al <- hulp_align_gdal(paths$huetzon, master_grid)
    r_huet    <- terra::rast(f_huet_al)
    leg_huet  <- terra::cats(r_huet)[[1]]
    
    target_labels <- c("bg", "bgk", "bk1", "bk2", "bkk1", "bkk2")
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
    
    if (length(ok_huet_ids) == 0) ok_huet_ids <- c(1, 2, 3, 4, 11, 12)
    
    r_huet_ok <- terra::ifel(!is.na(r_huet) & (r_huet %in% ok_huet_ids), 0.01, NA)
    unlink(f_huet_al)
    rm(r_huet)
    terra::writeRaster(r_huet_ok, f_cache_huet_rombout, overwrite = TRUE, datatype = "FLT4S", gdal = "COMPRESS=LZW")
  } else {
    r_huet_ok <- master_grid * NA
  }
  
  # Waterbasis samenstellen en clusteren op 10 m
  r_water_tot <- terra::cover(r_wat_raw, r_huet_ok)
  r_water_tot <- terra::clamp(r_water_tot, upper = 1.0)
  rm(r_wat_raw, r_huet_ok)
  gc()
  
  message(sprintf("    💧 Beken clusteren (10m duiker-overbrugging, >= %.2f ha)...", min_opp_water))
  res_wat <- cluster_filter_universeel(
    masker            = r_water_tot,
    opp_laag          = r_water_tot,
    min_opp_ha        = min_opp_water,
    min_totale_opp_ha = NA,
    afstand_m         = 10,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_water_valide <- res_wat$raster
  rm(res_wat, r_water_tot)
  gc()
  
  ha_wat_val <- calc_ha_exact(r_water_valide)
  message(sprintf("       -> Goedgekeurde beektrajecten (>= 0.4 ha): %.2f ha", ha_wat_val))
  
  # ----------------------------------------------------------------------------
  # 3. LANDBIOTOOP: BOS (>= 5 HA) EN MATRIX (>= 50 HA) OP 50 M NETWERK
  # ----------------------------------------------------------------------------
  message(sprintf("    🌳 Boskernen clusteren (>= %.1f ha, 50m netwerk)...", min_opp_bos))
  res_bos <- cluster_filter_universeel(
    masker            = r_bos_raw,
    opp_laag          = r_bos_raw,
    min_opp_ha        = min_opp_bos,
    min_totale_opp_ha = NA,
    afstand_m         = 50,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_bos_5ha  <- res_bos$raster
  cl_bos_ids <- res_bos$clusters
  rm(res_bos, r_bos_raw)
  gc()
  
  message(sprintf("    🌾 Dragende landschapsmatrix clusteren (>= %.1f ha, 50m netwerk)...", min_opp_rest))
  res_rest <- cluster_filter_universeel(
    masker            = r_rest_raw,
    opp_laag          = r_rest_raw,
    min_opp_ha        = min_opp_rest,
    min_totale_opp_ha = NA,
    afstand_m         = 50,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_rest_50ha <- res_rest$raster
  rm(res_rest, r_rest_raw)
  gc()
  
  # Inbeddingstoets: Bos (>= 5 ha) moet fysiek grenzen aan / liggen in matrix (>= 50 ha) met 10m speling
  message("    🔍 Inbeddingstoets: controleren of boskernen in de extensieve matrix liggen...")
  
  heeft_bos  <- !all(is.na(suppressWarnings(terra::minmax(r_bos_5ha))))
  heeft_rest <- !all(is.na(suppressWarnings(terra::minmax(r_rest_50ha))))
  
  if (heeft_bos && heeft_rest) {
    # 10m buffer rond matrix via 3x3 focal window
    w_mat <- matrix(1, nrow = 3, ncol = 3)
    rest_zone10 <- terra::focal(!is.na(r_rest_50ha) & r_rest_50ha > 0, w = w_mat, fun = "max", na.rm = TRUE)
    
    # Zonale statistiek: welke bosclusters raken de matrix?
    df_overlap <- terra::zonal(rest_zone10, cl_bos_ids, fun = "max", na.rm = TRUE)
    colnames(df_overlap) <- c("ID", "Raakt_Matrix")
    ids_bos_ok <- df_overlap$ID[!is.na(df_overlap$Raakt_Matrix) & df_overlap$Raakt_Matrix > 0]
    rm(rest_zone10, df_overlap)
    
    if (length(ids_bos_ok) > 0) {
      r_bos_ingebed <- terra::mask(r_bos_5ha, cl_bos_ids %in% ids_bos_ok)
    } else {
      r_bos_ingebed <- master_grid * NA
    }
  } else {
    r_bos_ingebed <- master_grid * NA
  }
  rm(r_bos_5ha, cl_bos_ids)
  gc()
  
  ha_bos_ingebed <- calc_ha_exact(r_bos_ingebed)
  message(sprintf("       -> Boskernen goedgekeurd in landschapsmatrix: %.2f ha", ha_bos_ingebed))
  
  # ----------------------------------------------------------------------------
  # 4. WEDERZIJDSE 500 M AFSTANDSFILTER TUSSEN WATER EN RIJPINGSBOS
  # ----------------------------------------------------------------------------
  message(sprintf("    🎯 Wederzijdse interactie (%dm) berekenen via GDAL proximity...", as.integer(koppel_m)))
  
  heeft_water_ok <- !all(is.na(suppressWarnings(terra::minmax(r_water_valide))))
  heeft_bos_ok   <- !all(is.na(suppressWarnings(terra::minmax(r_bos_ingebed))))
  
  if (heeft_water_ok && heeft_bos_ok) {
    # 1. Proximity vanaf bos naar water (water moet binnen 500m van bos liggen)
    f_b_bin  <- tempfile(fileext = ".tif")
    terra::writeRaster(terra::ifel(!is.na(r_bos_ingebed) & r_bos_ingebed > 0, 1, NA), f_b_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    f_prox_b <- hulp_proximity_gdal(f_b_bin, koppel_m)
    unlink(f_b_bin)
    
    p_b <- terra::rast(f_prox_b)
    mask_b500 <- terra::ifel(!is.na(p_b) & p_b >= 0 & p_b <= koppel_m, 1, NA)
    unlink(f_prox_b)
    rm(p_b)
    
    water_finaal <- terra::mask(r_water_valide, mask_b500)
    rm(mask_b500)
    
    # 2. Proximity vanaf het geselecteerde water naar bos (bos afgekapt op 500m van water)
    if (!all(is.na(suppressWarnings(terra::minmax(water_finaal))))) {
      f_w_bin  <- tempfile(fileext = ".tif")
      terra::writeRaster(terra::ifel(!is.na(water_finaal) & water_finaal > 0, 1, NA), f_w_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
      f_prox_w <- hulp_proximity_gdal(f_w_bin, koppel_m)
      unlink(f_w_bin)
      
      p_w <- terra::rast(f_prox_w)
      mask_w500 <- terra::ifel(!is.na(p_w) & p_w >= 0 & p_w <= koppel_m, 1, NA)
      unlink(f_prox_w)
      rm(p_w)
      
      bos_finaal <- terra::mask(r_bos_ingebed, mask_w500)
      rm(mask_w500)
    } else {
      bos_finaal <- master_grid * NA
    }
  } else {
    water_finaal <- master_grid * NA
    bos_finaal   <- master_grid * NA
  }
  rm(r_water_valide, r_bos_ingebed)
  gc()
  
  ha_wat_fin <- calc_ha_exact(water_finaal)
  ha_bos_fin <- calc_ha_exact(bos_finaal)
  message(sprintf("       -> Voortplantingsbeken binnen 500m van rijpingsbos: %.2f ha", ha_wat_fin))
  message(sprintf("       -> Rijpingsbos binnen 500m van goedgekeurde beken: %.2f ha", ha_bos_fin))
  
  # ----------------------------------------------------------------------------
  # 5. TOTAAL LEEFGEBIED & NETWERK CLUSTERING
  # ----------------------------------------------------------------------------
  message("    🧩 Totaal leefgebiedscomplex samenstellen (beek + rijpingsbos)...")
  final_opp <- terra::cover(water_finaal, bos_finaal)
  final_opp <- terra::clamp(final_opp, upper = 1.0)
  
  ha_tot <- calc_ha_exact(final_opp)
  message(sprintf("       -> Definitief leefgebied Beekrombout: %.2f ha", ha_tot))
  
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
    unlink(f_bin)
    
    cl_opp <- terra::rast(f_clump)
    try(terra::set.cats(cl_opp, layer = 1, value = NULL), silent = TRUE)
    try(terra::coltab(cl_opp) <- NULL, silent = TRUE)
    cl_opp <- terra::mask(cl_opp, r_fin_bin)
    unlink(f_clump)
  } else {
    cl_opp <- master_grid * NA
  }
  rm(r_fin_bin)
  gc()
  
  return(list(
    final_opp        = final_opp,
    cl_opp           = cl_opp,
    water            = water_finaal,
    rijpingsbos      = bos_finaal,
    matrix_landschap = r_rest_50ha
  ))
}
