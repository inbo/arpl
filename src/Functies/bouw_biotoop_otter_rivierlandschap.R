# ==============================================================================
# bouw_biotoop_otter_rivierlandschap.R
# Archetype voor Europese otter: 
# - Vlaanderen algemeen: >= 1.000 ha landmatrix
# - Grensmaas (Ecodistricten 17 & 36): verlaagde drempel >= 400 ha landmatrix
# - Wederzijdse water-landkoppeling binnen <= 1.000 m
# Conform officiële INBO-standaard ANB2019_EuropeseOtter2018.asc
# ==============================================================================

bouw_biotoop_otter_rivierlandschap <- function(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
  bestandsnaam <- tolower(gsub(" ", "", soort_naam))
  soort_cache  <- file.path(paths$dir_cache, bestandsnaam)
  dir.create(soort_cache, showWarnings = FALSE, recursive = TRUE)
  
  f_cache_rds <- file.path(soort_cache, "01_raw_parquet_layers.rds")
  
  # Parameters ophalen
  rij_land      <- params_soort %>% dplyr::filter(tolower(Type) == "landbiotoop") %>% dplyr::slice(1)
  min_land_ha   <- as.numeric(ifelse(nrow(rij_land) > 0 && !is.na(rij_land$MinOpp_ha), rij_land$MinOpp_ha, 1000))
  min_maas_ha   <- as.numeric(ifelse(nrow(rij_land) > 0 && !is.na(rij_land$Kern_MinOpp_ha), rij_land$Kern_MinOpp_ha, 400))
  zoekafstand_m <- as.numeric(ifelse(nrow(rij_land) > 0 && !is.na(rij_land$AfstandBiotopen_m), rij_land$AfstandBiotopen_m, 100))
  koppel_dist_m <- as.numeric(ifelse(nrow(rij_land) > 0 && !is.na(rij_land$Wederzijdse_Afstand_m), rij_land$Wederzijdse_Afstand_m, 1000))
  
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
  
  r_land_raw <- if ("landbiotoop" %in% names(lagen)) lagen[["landbiotoop"]] else (master_grid * NA)
  
  ha_land_ruw <- calc_ha_exact(r_land_raw)
  message(sprintf("       -> Totale landmatrix vóór filtering: %.2f ha", ha_land_ruw))
  
  # ----------------------------------------------------------------------------
  # 2. INTEGRAAL WATERNETWERK (HUET-ZONES WAARDE >= 1 + WATERVLAKKEN)
  # ----------------------------------------------------------------------------
  message("    🌊 Integraal waternetwerk samenstellen (waterlopen + watervlakken)...")
  
  f_cache_water <- file.path(paths$dir_cache, "otter_water_integraal.tif")
  
  if (file.exists(f_cache_water)) {
    r_water_tot <- terra::rast(f_cache_water)
  } else {
    # A. Huet-zones: Primaire waterloopvakken (waarde >= 1 conform script)
    if (file.exists(paths$huetzon)) {
      f_huet_al   <- hulp_align_gdal(paths$huetzon, master_grid)
      r_huet      <- terra::rast(f_huet_al)
      r_waterloop <- terra::ifel(!is.na(r_huet) & r_huet >= 1, 1, NA)
      unlink(f_huet_al)
      rm(r_huet)
    } else {
      r_waterloop <- master_grid * NA
    }
    
    # B. Watervlakken (plassen, meren, grindgaten)
    if (file.exists(paths$watervlakken)) {
      f_wat_al <- hulp_align_gdal(paths$watervlakken, master_grid)
      r_wat    <- terra::rast(f_wat_al)
      r_plas   <- terra::ifel(!is.na(r_wat) & r_wat >= 1, 1, NA)
      unlink(f_wat_al)
      rm(r_wat)
    } else {
      r_plas <- master_grid * NA
    }
    
    # C. Optionele extra Maas-laag (indien aanwezig)
    if (!is.null(paths$maas) && file.exists(paths$maas)) {
      message("       🌊 Extra Grensmaas-laag integreren in waternetwerk...")
      f_maas_al <- hulp_align_gdal(paths$maas, master_grid)
      r_maas    <- terra::rast(f_maas_al)
      r_maas_ok <- terra::ifel(!is.na(r_maas) & r_maas >= 1, 1, NA)
      unlink(f_maas_al)
      r_water_tot <- terra::cover(r_water_tot, r_maas_ok)
      rm(r_maas, r_maas_ok)
    }
    
    terra::writeRaster(r_water_tot, f_cache_water, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
  }
  
  ha_water_tot <- calc_ha_exact(r_water_tot)
  message(sprintf("       -> Integraal waternetwerk: %.2f ha", ha_water_tot))
  
  # ----------------------------------------------------------------------------
  # 3. LANDMATRIX CLUSTEREN (100M ZOEFAFSTAND) & DUBBELE DREMPELTOETS
  # ----------------------------------------------------------------------------
  message(sprintf("    🌳 Landmatrix clusteren op %dm zoekafstand...", as.integer(zoekafstand_m)))
  
  r_bin_land <- terra::ifel(!is.na(r_land_raw) & r_land_raw > 0, 1, NA)
  
  if (all(is.na(suppressWarnings(terra::minmax(r_bin_land))))) {
    return(list(final_opp = master_grid * NA, cl_opp = master_grid * NA))
  }
  
  f_bin_l   <- tempfile(fileext = ".tif")
  f_clump_l <- tempfile(fileext = ".tif")
  
  # Buffer 50m voor 100m overbrugging
  terra::writeRaster(r_bin_land, f_bin_l, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
  f_prox_l100 <- hulp_proximity_gdal(f_bin_l, zoekafstand_m / 2)
  unlink(f_bin_l)
  
  p_l100   <- terra::rast(f_prox_l100)
  r_buf100 <- terra::ifel(!is.na(p_l100) & p_l100 >= 0 & p_l100 <= (zoekafstand_m / 2), 1, NA)
  unlink(f_prox_l100)
  rm(p_l100)
  
  f_buf100 <- tempfile(fileext = ".tif")
  terra::writeRaster(r_buf100, f_buf100, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
  rm(r_buf100)
  
  e <- terra::ext(master_grid)
  qgisprocess::qgis_run_algorithm(
    "grass:r.clump", input = f_buf100, output = f_clump_l, title = "cl_land_otter",
    "-d" = FALSE, GRASS_REGION_PARAMETER = sprintf("%f,%f,%f,%f", e$xmin, e$xmax, e$ymin, e$ymax),
    GRASS_REGION_CELLSIZE_PARAMETER = 10, .quiet = TRUE
  )
  unlink(f_buf100)
  
  cl_l_raw <- terra::rast(f_clump_l)
  terra::set.cats(cl_l_raw, layer = 1, value = NULL)
  if (terra::has.colors(cl_l_raw)) try(terra::coltab(cl_l_raw) <- NULL, silent = TRUE)
  cl_land <- terra::mask(cl_l_raw, r_bin_land)
  unlink(f_clump_l)
  rm(cl_l_raw, r_bin_land)
  gc()
  
  # Bereken oppervlakte per landcluster via zonale statistiek
  df_land_stats <- terra::zonal(r_land_raw, cl_land, fun = "sum", na.rm = TRUE)
  colnames(df_land_stats) <- c("Cluster_ID", "Som_Val")
  dt_land <- data.table::as.data.table(df_land_stats)
  dt_land[, Opp_ha := Som_Val * 0.01]
  rm(df_land_stats)
  
  # A. Standaard Vlaamse clusters >= 1.000 ha
  ids_1000ha <- dt_land[!is.na(Cluster_ID) & Opp_ha >= min_land_ha, Cluster_ID]
  
  # B. Grensmaas clusters: 400 ha <= Opp < 1.000 ha getoetst aan Maasdistrict (17, 36)
  ids_kandidaat_maas <- dt_land[!is.na(Cluster_ID) & Opp_ha >= min_maas_ha & Opp_ha < min_land_ha, Cluster_ID]
  
  ids_maas_ok <- integer(0)
  if (length(ids_kandidaat_maas) > 0 && file.exists(paths$ecoregio)) {
    message("    🌊 Grensmaas-toets uitvoeren voor clusters tussen 400 en 1.000 ha...")
    f_eco_al <- hulp_align_gdal(paths$ecoregio, master_grid)
    r_eco    <- terra::rast(f_eco_al)
    
    # Ecodistricten 17 en 36
    mask_maas_eco <- terra::ifel(!is.na(r_eco) & (r_eco %in% c(17, 36)), 1, NA)
    unlink(f_eco_al)
    rm(r_eco)
    
    # Controleer overlap van kandidaat-clusters met Maasdistrict
    cl_kandidaat <- terra::ifel(cl_land %in% as.integer(ids_kandidaat_maas), cl_land, NA)
    df_maas_over <- terra::zonal(mask_maas_eco, cl_kandidaat, fun = "max", na.rm = TRUE)
    colnames(df_maas_over) <- c("Cluster_ID", "In_Maas")
    
    ids_maas_ok <- df_maas_over$Cluster_ID[!is.na(df_maas_over$In_Maas) & df_maas_over$In_Maas == 1]
    rm(mask_maas_eco, cl_kandidaat, df_maas_over)
    gc()
    message(sprintf("       -> Extra Grensmaas-clusters goedgekeurd (>= 400 ha): %d", length(ids_maas_ok)))
  }
  
  alle_ok_land_ids <- unique(c(ids_1000ha, ids_maas_ok))
  rm(dt_land)
  
  if (length(alle_ok_land_ids) > 0) {
    cl_land_valide <- terra::ifel(cl_land %in% as.integer(alle_ok_land_ids), 1, NA)
    r_land_gefilterd <- terra::mask(r_land_raw, cl_land_valide)
    rm(cl_land_valide)
  } else {
    r_land_gefilterd <- master_grid * NA
  }
  rm(cl_land, r_land_raw)
  gc()
  
  ha_land_ok <- calc_ha_exact(r_land_gefilterd)
  message(sprintf("       -> Landmatrix na drempeltoets (>= 1.000 ha / >= 400 ha Maas): %.2f ha", ha_land_ok))
  
  # ----------------------------------------------------------------------------
  # 4. WEDERZIJDSE KOPPELING WATER <-> LAND (<= 1.000 METER)
  # ----------------------------------------------------------------------------
  message(sprintf("    🔗 Wederzijdse actieradiuskoppeling berekenen (<= %d m via GDAL proximity)...", as.integer(koppel_dist_m)))
  
  heeft_land  <- (ha_land_ok > 0)
  heeft_water <- (ha_water_tot > 0)
  
  if (heeft_land && heeft_water) {
    # A. Proximity tot water berekenen (max. 1.000 m)
    f_w_bin <- tempfile(fileext = ".tif")
    terra::writeRaster(r_water_tot, f_w_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    f_prox_w <- hulp_proximity_gdal(f_w_bin, koppel_dist_m)
    unlink(f_w_bin)
    
    p_w <- terra::rast(f_prox_w)
    mask_land_ok <- terra::ifel(!is.na(p_w) & p_w >= 0 & p_w <= koppel_dist_m, 1, NA)
    unlink(f_prox_w)
    rm(p_w)
    
    # Filter het land: landmatrix binnen 1 km van water
    r_land_gekoppeld <- terra::mask(r_land_gefilterd, mask_land_ok)
    rm(mask_land_ok, r_land_gefilterd)
    gc()
    
    # B. Proximity tot land berekenen (max. 1.000 m)
    r_land_bin <- terra::ifel(!is.na(r_land_gekoppeld) & r_land_gekoppeld > 0, 1, NA)
    f_l_bin <- tempfile(fileext = ".tif")
    terra::writeRaster(r_land_bin, f_l_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    rm(r_land_bin)
    
    f_prox_l <- hulp_proximity_gdal(f_l_bin, koppel_dist_m)
    unlink(f_l_bin)
    
    p_l <- terra::rast(f_prox_l)
    mask_wat_ok <- terra::ifel(!is.na(p_l) & p_l >= 0 & p_l <= koppel_dist_m, 1, NA)
    unlink(f_prox_l)
    rm(p_l)
    
    # Filter het water: watergangen/plassen binnen 1 km van landmatrix
    r_water_gekoppeld <- terra::mask(r_water_tot, mask_wat_ok)
    rm(mask_wat_ok, r_water_tot)
    gc()
    
    # C. Samensmelten tot definitief leefgebied
    l_c <- terra::ifel(is.na(r_land_gekoppeld), 0, r_land_gekoppeld)
    w_c <- terra::ifel(is.na(r_water_gekoppeld), 0, r_water_gekoppeld)
    
    r_tot_raw <- terra::ifel((l_c + w_c) > 0, terra::clamp(l_c + w_c, upper = 1.0), NA)
    final_opp <- r_tot_raw
    rm(l_c, w_c, r_tot_raw)
    gc()
    
    # D. Unieke netwerk ID's via GRASS r.clump (8-richtingen)
    f_fin_bin <- tempfile(fileext = ".tif")
    f_cl_fin  <- tempfile(fileext = ".tif")
    terra::writeRaster(terra::ifel(!is.na(final_opp), 1, NA), f_fin_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    
    e <- terra::ext(master_grid)
    qgisprocess::qgis_run_algorithm(
      "grass:r.clump", input = f_fin_bin, output = f_cl_fin, title = "cl_otter",
      "-d" = TRUE, # directions = 8
      GRASS_REGION_PARAMETER = sprintf("%f,%f,%f,%f", e$xmin, e$xmax, e$ymin, e$ymax),
      GRASS_REGION_CELLSIZE_PARAMETER = 10, .quiet = TRUE
    )
    unlink(f_fin_bin)
    
    cl_opp <- terra::rast(f_cl_fin)
    terra::set.cats(cl_opp, layer = 1, value = NULL)
    if (terra::has.colors(cl_opp)) try(terra::coltab(cl_opp) <- NULL, silent = TRUE)
    cl_opp <- terra::mask(cl_opp, final_opp)
    unlink(f_cl_fin)
    
  } else {
    message("    ⚠️ Geen geldige landmatrix of watergangen na filtering.")
    final_opp <- master_grid * NA
    cl_opp    <- master_grid * NA
  }
  
  ha_finaal <- calc_ha_exact(final_opp)
  message(sprintf("       -> Definitief leefgebied Europese otter: %.2f ha", ha_finaal))
  
  return(list(
    final_opp = final_opp,
    cl_opp    = cl_opp
  ))
}
