# ==============================================================================
# bouw_biotoop_hoogveen_libel.R
# Archetype voor hoogveenlibellen (veenmoswater >= 0.01ha x land >= 5ha met 60% bos op 500m x ruim land >= 50ha)
# Geschikt voor: Hoogveenglanslibel (Somatochlora arctica), etc.
# ==============================================================================

bouw_biotoop_hoogveen_libel <- function(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
  bestandsnaam <- tolower(gsub(" ", "", soort_naam))
  soort_cache  <- file.path(paths$dir_cache, bestandsnaam)
  dir.create(soort_cache, showWarnings = FALSE, recursive = TRUE)
  
  f_cache_rds <- file.path(soort_cache, "01_raw_parquet_layers.rds")
  
  # Parameters ophalen uit Excel
  rij_wat <- params_soort %>% dplyr::filter(grepl("waterbiotoop_bwk1", tolower(Type))) %>% dplyr::slice(1)
  min_opp_water <- as.numeric(ifelse(nrow(rij_wat) > 0 && !is.na(rij_wat$MinOpp_ha), rij_wat$MinOpp_ha, 0.01))
  
  rij_land <- params_soort %>% dplyr::filter(grepl("landbiotoop_bwk1", tolower(Type))) %>% dplyr::slice(1)
  min_opp_land <- as.numeric(ifelse(nrow(rij_land) > 0 && !is.na(rij_land$MinOpp_ha), rij_land$MinOpp_ha, 5.0))
  
  rij_ruim <- params_soort %>% dplyr::filter(grepl("landbiotoop_bwk_ruim", tolower(Type))) %>% dplyr::slice(1)
  min_opp_ruim <- as.numeric(ifelse(nrow(rij_ruim) > 0 && !is.na(rij_ruim$MinOpp_ha), rij_ruim$MinOpp_ha, 50.0))
  
  koppel_m <- 500 # 500 m vlieg- en foerageerstraal rond water
  
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
  
  r_wat1_raw      <- if ("waterbiotoop_bwk1"     %in% names(lagen)) lagen[["waterbiotoop_bwk1"]]     else (master_grid * NA)
  r_not_wat_raw   <- if ("waterbiotoop_notbwk"   %in% names(lagen)) lagen[["waterbiotoop_notbwk"]]   else (master_grid * NA)
  r_land1_raw     <- if ("landbiotoop_bwk1"      %in% names(lagen)) lagen[["landbiotoop_bwk1"]]      else (master_grid * NA)
  r_not_land_raw  <- if ("landbiotoop_notbwk"    %in% names(lagen)) lagen[["landbiotoop_notbwk"]]    else (master_grid * NA)
  r_bos_raw       <- if ("landbiotoop_bos"       %in% names(lagen)) lagen[["landbiotoop_bos"]]       else (master_grid * NA)
  r_ruim_raw      <- if ("landbiotoop_bwk_ruim"  %in% names(lagen)) lagen[["landbiotoop_bwk_ruim"]]  else (master_grid * NA)
  
  # ----------------------------------------------------------------------------
  # 2. WATERBIOTOOP SCHONEN EN CLUSTEREN (>= 0.01 HA OP 50 M)
  # ----------------------------------------------------------------------------
  message(sprintf("    💧 Veenmos- en hoogveenwater schonen en clusteren (>= %.2f ha op 50m)...", min_opp_water))
  
  # Negatieve aftrek van eutroof/ongeschikt water
  r_not_wat_clean <- terra::ifel(is.na(r_not_wat_raw), 0, r_not_wat_raw)
  r_wat_geschoond <- terra::clamp(r_wat1_raw - r_not_wat_clean, lower = 0, upper = 1.0)
  r_wat_geschoond <- terra::ifel(r_wat_geschoond > 0, r_wat_geschoond, NA)
  rm(r_not_wat_clean, r_wat1_raw, r_not_wat_raw)
  gc()
  
  res_wat <- cluster_filter_universeel(
    masker            = r_wat_geschoond,
    opp_laag          = r_wat_geschoond,
    min_opp_ha        = min_opp_water,
    min_totale_opp_ha = NA,
    afstand_m         = 50,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_wat_geclusterd <- res_wat$raster
  rm(res_wat, r_wat_geschoond)
  gc()
  
  ha_wat_gecl <- calc_ha_exact(r_wat_geclusterd)
  message(sprintf("       -> Goedgekeurde hoogveenslenken/water (>= 0.01 ha): %.2f ha", ha_wat_gecl))
  
  # ----------------------------------------------------------------------------
  # 3. RUIM LANDBIOTOOP CLUSTEREN (>= 50 HA OP 50 M)
  # ----------------------------------------------------------------------------
  message(sprintf("    🌾 Ruim hoogveenlandschap clusteren (>= %.1f ha, 50m netwerk)...", min_opp_ruim))
  res_ruim <- cluster_filter_universeel(
    masker            = r_ruim_raw,
    opp_laag          = r_ruim_raw,
    min_opp_ha        = min_opp_ruim,
    min_totale_opp_ha = NA,
    afstand_m         = 50,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_ruim_50ha <- res_ruim$raster
  rm(res_ruim, r_ruim_raw)
  gc()
  
  ha_ruim <- calc_ha_exact(r_ruim_50ha)
  message(sprintf("       -> Ruim landschapscomplex (>= 50 ha): %.2f ha", ha_ruim))
  
  # ----------------------------------------------------------------------------
  # 4. AANPALEND LANDBIOTOOP FILTEREN (>= 5 HA OP 10 M, >= 60% BOS)
  # ----------------------------------------------------------------------------
  message(sprintf("    🌳 Aanpalend land schonen, clusteren (>= %.1f ha) en filteren op >= 60%% bos...", min_opp_land))
  
  # Negatieve aftrek van droge heide/stuifzand
  r_not_land_clean <- terra::ifel(is.na(r_not_land_raw), 0, r_not_land_raw)
  r_land_geschoond <- terra::clamp(r_land1_raw - r_not_land_clean, lower = 0, upper = 1.0)
  r_land_geschoond <- terra::ifel(r_land_geschoond > 0, r_land_geschoond, NA)
  rm(r_not_land_clean, r_land1_raw, r_not_land_raw)
  gc()
  
  res_land5 <- cluster_filter_universeel(
    masker            = r_land_geschoond,
    opp_laag          = r_land_geschoond,
    min_opp_ha        = min_opp_land,
    min_totale_opp_ha = NA,
    afstand_m         = 10,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_land_5ha  <- res_land5$raster
  cl_land_5ha <- res_land5$clusters
  rm(res_land5, r_land_geschoond)
  gc()
  
  # Toetsen op minstens 60% bosdekking
  if (!all(is.na(suppressWarnings(terra::minmax(cl_land_5ha))))) {
    r_bos_clean <- terra::ifel(!is.na(r_bos_raw) & r_bos_raw > 0, r_bos_raw, 0)
    
    df_land_tot <- terra::zonal(r_land_5ha, cl_land_5ha, fun = "sum", na.rm = TRUE)
    colnames(df_land_tot) <- c("ID", "Totaal_Sum")
    
    df_bos_tot <- terra::zonal(r_bos_clean, cl_land_5ha, fun = "sum", na.rm = TRUE)
    colnames(df_bos_tot) <- c("ID", "Bos_Sum")
    
    dt_check <- merge(data.table::as.data.table(df_land_tot), data.table::as.data.table(df_bos_tot), by = "ID", all.x = TRUE)
    dt_check[is.na(Bos_Sum), Bos_Sum := 0]
    dt_check[, Pct_Bos := (Bos_Sum / Totaal_Sum) * 100]
    
    valide_land_ids <- dt_check[!is.na(ID) & Pct_Bos >= 60, ID]
    rm(df_land_tot, df_bos_tot, dt_check, r_bos_clean)
    
    if (length(valide_land_ids) > 0) {
      r_land_bos60 <- terra::mask(r_land_5ha, cl_land_5ha %in% as.integer(valide_land_ids))
    } else {
      r_land_bos60 <- master_grid * NA
    }
  } else {
    r_land_bos60 <- master_grid * NA
  }
  rm(r_land_5ha, cl_land_5ha, r_bos_raw)
  gc()
  
  ha_bos60 <- calc_ha_exact(r_land_bos60)
  message(sprintf("       -> Landbiotoop met >= 60%% bosdekking: %.2f ha", ha_bos60))
  
  # ----------------------------------------------------------------------------
  # 5. WATER-FIRST AFBAKENING & RUIME OMGEVINGSFACTOR (500 M AFSTANDSTOETS)
  # ----------------------------------------------------------------------------
  message("    🎯 500m koppeling berekenen: Land afsnijden rond water & water toetsen aan 50ha landschap...")
  
  heeft_water <- !all(is.na(suppressWarnings(terra::minmax(r_wat_geclusterd))))
  heeft_land  <- !all(is.na(suppressWarnings(terra::minmax(r_land_bos60))))
  heeft_ruim  <- !all(is.na(suppressWarnings(terra::minmax(r_ruim_50ha))))
  
  if (heeft_water && heeft_ruim) {
    # 1. Water toetsen aan 50 ha ruim landcomplex (moet binnen <= 500m liggen)
    f_rm_bin  <- tempfile(fileext = ".tif")
    terra::writeRaster(terra::ifel(!is.na(r_ruim_50ha) & r_ruim_50ha > 0, 1, NA), f_rm_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    f_prox_rm <- hulp_proximity_gdal(f_rm_bin, koppel_m)
    unlink(f_rm_bin)
    
    p_rm <- terra::rast(f_prox_rm)
    mask_rm500 <- terra::ifel(!is.na(p_rm) & p_rm >= 0 & p_rm <= koppel_m, 1, NA)
    unlink(f_prox_rm)
    rm(p_rm)
    
    water_finaal <- terra::mask(r_wat_geclusterd, mask_rm500)
    rm(mask_rm500)
  } else {
    water_finaal <- master_grid * NA
  }
  rm(r_wat_geclusterd)
  gc()
  
  # 2. Land (>= 60% bos) afsnijden op 500 m rond goedgekeurd water
  if (!all(is.na(suppressWarnings(terra::minmax(water_finaal)))) && heeft_land) {
    f_wf_bin  <- tempfile(fileext = ".tif")
    terra::writeRaster(terra::ifel(!is.na(water_finaal) & water_finaal > 0, 1, NA), f_wf_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    f_prox_wf <- hulp_proximity_gdal(f_wf_bin, koppel_m)
    unlink(f_wf_bin)
    
    p_wf <- terra::rast(f_prox_wf)
    mask_wf500 <- terra::ifel(!is.na(p_wf) & p_wf >= 0 & p_wf <= koppel_m, 1, NA)
    unlink(f_prox_wf)
    rm(p_wf)
    
    land_finaal <- terra::mask(r_land_bos60, mask_wf500)
    rm(mask_wf500)
  } else {
    land_finaal <- master_grid * NA
  }
  rm(r_land_bos60)
  gc()
  
  ha_wat_fin  <- calc_ha_exact(water_finaal)
  ha_land_fin <- calc_ha_exact(land_finaal)
  message(sprintf("       -> Waterbiotoop binnen 500m van ruim 50 ha complex: %.2f ha", ha_wat_fin))
  message(sprintf("       -> Jachtbos (>= 60%% bos) binnen 500m van geschikt water: %.2f ha", ha_land_fin))
  
  # ----------------------------------------------------------------------------
  # 6. SAMENVOEGING TOT DEFINITIEF LEEFGEBIED & NETWERKEN
  # ----------------------------------------------------------------------------
  message("    🧩 Totaal leefgebiedscomplex samenstellen (hoogveenwater + jachtbos)...")
  final_opp <- terra::cover(water_finaal, land_finaal)
  final_opp <- terra::clamp(final_opp, upper = 1.0)
  
  ha_tot <- calc_ha_exact(final_opp)
  message(sprintf("       -> Definitief leefgebied Hoogveenglanslibel: %.2f ha", ha_tot))
  
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
    final_opp           = final_opp,
    cl_opp              = cl_opp,
    water               = water_finaal,
    jachtbos_500m       = land_finaal,
    ruim_landschap_50ha = r_ruim_50ha
  ))
}
