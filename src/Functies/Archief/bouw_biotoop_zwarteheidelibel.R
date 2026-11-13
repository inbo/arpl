# ==============================================================================
# bouw_biotoop_zwarteheidelibel.R
# Archetype voor kleine heidelibellen (water >= 0.01ha x land >= 5ha op 20m x ruim land >= 40ha met 40% bos op 500m x land afgesneden op 500m)
# Geschikt voor: Zwarte heidelibel (Sympetrum danae), etc.
# ==============================================================================

bouw_biotoop_zwarteheidelibel <- function(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
  bestandsnaam <- tolower(gsub(" ", "", soort_naam))
  soort_cache  <- file.path(paths$dir_cache, bestandsnaam)
  dir.create(soort_cache, showWarnings = FALSE, recursive = TRUE)
  
  f_cache_rds <- file.path(soort_cache, "01_raw_parquet_layers.rds")
  
  # Parameters ophalen uit Excel
  rij_wat <- params_soort %>% dplyr::filter(grepl("waterbiotoop_bwk1", tolower(Type))) %>% dplyr::slice(1)
  min_opp_water <- as.numeric(ifelse(nrow(rij_wat) > 0 && !is.na(rij_wat$MinOpp_ha), rij_wat$MinOpp_ha, 0.01))
  zoekafstand_wat_m <- as.numeric(ifelse(nrow(rij_wat) > 0 && !is.na(rij_wat$AfstandBiotopen_m), rij_wat$AfstandBiotopen_m, 50))
  
  rij_land <- params_soort %>% dplyr::filter(grepl("landbiotoop_bwk$", tolower(Type))) %>% dplyr::slice(1)
  min_opp_land <- as.numeric(ifelse(nrow(rij_land) > 0 && !is.na(rij_land$MinOpp_ha), rij_land$MinOpp_ha, 5.0))
  zoekafstand_land_m <- as.numeric(ifelse(nrow(rij_land) > 0 && !is.na(rij_land$AfstandBiotopen_m), rij_land$AfstandBiotopen_m, 50))
  
  rij_ruim <- params_soort %>% dplyr::filter(grepl("landbiotoop_bwk_ruim$", tolower(Type))) %>% dplyr::slice(1)
  min_opp_ruim <- as.numeric(ifelse(nrow(rij_ruim) > 0 && !is.na(rij_ruim$MinOpp_ha), rij_ruim$MinOpp_ha, 40.0))
  zoekafstand_ruim_m <- as.numeric(ifelse(nrow(rij_ruim) > 0 && !is.na(rij_ruim$AfstandBiotopen_m), rij_ruim$AfstandBiotopen_m, 50))
  
  oever_koppel_m <- 20 # 20 m oeverinteractie tussen water en land
  koppel_m       <- 500 # 500 m vlieg- en interactiezone rond water
  
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
  
  r_wat1_raw     <- if ("waterbiotoop_bwk1"        %in% names(lagen)) lagen[["waterbiotoop_bwk1"]]        else (master_grid * NA)
  r_not_wat_raw  <- if ("waterbiotoop_notbwk"      %in% names(lagen)) lagen[["waterbiotoop_notbwk"]]      else (master_grid * NA)
  r_land1_raw    <- if ("landbiotoop_bwk"          %in% names(lagen)) lagen[["landbiotoop_bwk"]]          else (master_grid * NA)
  r_ruim_raw     <- if ("landbiotoop_bwk_ruim"     %in% names(lagen)) lagen[["landbiotoop_bwk_ruim"]]     else (master_grid * NA)
  r_bos_ruim_raw <- if ("landbiotoop_bwk_ruim_bos" %in% names(lagen)) lagen[["landbiotoop_bwk_ruim_bos"]] else (master_grid * NA)
  
  # ----------------------------------------------------------------------------
  # 2. WATERBIOTOOP SCHONEN EN CLUSTEREN (>= 0.01 HA OP 50 M)
  # ----------------------------------------------------------------------------
  message(sprintf("    💧 Dystrofe vennen schonen en clusteren (>= %.2f ha op %dm)...", min_opp_water, as.integer(zoekafstand_wat_m)))
  
  r_not_wat_clean <- terra::ifel(is.na(r_not_wat_raw), 0, r_not_wat_raw)
  r_wat_geschoond <- terra::clamp(r_wat1_raw - r_not_wat_clean, lower = 0, upper = 1.0)
  r_wat_geschoond <- terra::ifel(r_wat_geschoond > 0, r_wat_geschoond, NA)
  rm(r_not_wat_clean, r_not_wat_raw, r_wat1_raw)
  gc()
  
  res_wat <- cluster_filter_universeel(
    masker            = r_wat_geschoond,
    opp_laag          = r_wat_geschoond,
    min_opp_ha        = min_opp_water,
    min_totale_opp_ha = NA,
    afstand_m         = zoekafstand_wat_m,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_wat_geclusterd <- res_wat$raster
  rm(res_wat, r_wat_geschoond)
  gc()
  
  ha_wat_gecl <- calc_ha_exact(r_wat_geclusterd)
  message(sprintf("       -> Goedgekeurde vennen (>= 0.01 ha): %.2f ha", ha_wat_gecl))
  
  # ----------------------------------------------------------------------------
  # 3. AANPALEND LANDBIOTOOP CLUSTEREN (>= 5 HA OP 50 M)
  # ----------------------------------------------------------------------------
  message(sprintf("    🌳 Aanpalend land clusteren (>= %.1f ha op %dm)...", min_opp_land, as.integer(zoekafstand_land_m)))
  
  res_land <- cluster_filter_universeel(
    masker            = r_land1_raw,
    opp_laag          = r_land1_raw,
    min_opp_ha        = min_opp_land,
    min_totale_opp_ha = NA,
    afstand_m         = zoekafstand_land_m,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_land_geclusterd <- res_land$raster
  rm(res_land, r_land1_raw)
  gc()
  
  ha_land_gecl <- calc_ha_exact(r_land_geclusterd)
  message(sprintf("       -> Goedgekeurd aanpalend land (>= 5 ha): %.2f ha", ha_land_gecl))
  
  # ----------------------------------------------------------------------------
  # 4. RUIM LANDBIOTOOP (>= 40 HA OP 50 M) & 40% BOSEIS
  # ----------------------------------------------------------------------------
  message(sprintf("    🌾 Ruim landschap clusteren (>= %.1f ha op %dm) en filteren op >= 40%% bos...", min_opp_ruim, as.integer(zoekafstand_ruim_m)))
  res_ruim <- cluster_filter_universeel(
    masker            = r_ruim_raw,
    opp_laag          = r_ruim_raw,
    min_opp_ha        = min_opp_ruim,
    min_totale_opp_ha = NA,
    afstand_m         = zoekafstand_ruim_m,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_ruim_40ha  <- res_ruim$raster
  cl_ruim_40ha <- res_ruim$clusters
  rm(res_ruim, r_ruim_raw)
  gc()
  
  if (!all(is.na(suppressWarnings(terra::minmax(cl_ruim_40ha))))) {
    r_bos_ruim_clean <- terra::ifel(!is.na(r_bos_ruim_raw) & r_bos_ruim_raw > 0, r_bos_ruim_raw, 0)
    
    df_ruim_tot <- terra::zonal(r_ruim_40ha, cl_ruim_40ha, fun = "sum", na.rm = TRUE)
    colnames(df_ruim_tot) <- c("ID", "Totaal_Sum")
    
    df_bos_ruim <- terra::zonal(r_bos_ruim_clean, cl_ruim_40ha, fun = "sum", na.rm = TRUE)
    colnames(df_bos_ruim) <- c("ID", "Bos_Sum")
    
    dt_ruim <- merge(data.table::as.data.table(df_ruim_tot), data.table::as.data.table(df_bos_ruim), by = "ID", all.x = TRUE)
    dt_ruim[is.na(Bos_Sum), Bos_Sum := 0]
    dt_ruim[, Pct_Bos := (Bos_Sum / Totaal_Sum) * 100]
    
    valide_ruim_ids <- dt_ruim[!is.na(ID) & Pct_Bos >= 40, ID]
    rm(df_ruim_tot, df_bos_ruim, dt_ruim, r_bos_ruim_clean)
    
    if (length(valide_ruim_ids) > 0) {
      r_ruim_valide <- terra::mask(r_ruim_40ha, cl_ruim_40ha %in% as.integer(valide_ruim_ids))
    } else {
      r_ruim_valide <- master_grid * NA
    }
  } else {
    r_ruim_valide <- master_grid * NA
  }
  rm(r_ruim_40ha, cl_ruim_40ha, r_bos_ruim_raw)
  gc()
  
  ha_ruim_val <- calc_ha_exact(r_ruim_valide)
  message(sprintf("       -> Ruim heide/boslandschap (>= 40 ha, >= 40%% bos): %.2f ha", ha_ruim_val))
  
  # ----------------------------------------------------------------------------
  # 5. KOPPELINGEN: 20 M OEVERCONTACT & 500 M INBEDDING IN RUIM COMPLEX
  # ----------------------------------------------------------------------------
  message("    🎯 Wederzijdse 20 m oeverkoppeling en 500 m landschapsinbedding berekenen...")
  
  heeft_wat  <- !all(is.na(suppressWarnings(terra::minmax(r_wat_geclusterd))))
  heeft_land <- !all(is.na(suppressWarnings(terra::minmax(r_land_geclusterd))))
  heeft_ruim <- !all(is.na(suppressWarnings(terra::minmax(r_ruim_valide))))
  
  if (heeft_wat && heeft_land && heeft_ruim) {
    # A. 20 m oeverkoppeling tussen water en land
    f_w_bin   <- tempfile(fileext = ".tif")
    terra::writeRaster(terra::ifel(!is.na(r_wat_geclusterd) & r_wat_geclusterd > 0, 1, NA), f_w_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    f_prox_w20 <- hulp_proximity_gdal(f_w_bin, oever_koppel_m)
    unlink(f_w_bin)
    
    p_w20 <- terra::rast(f_prox_w20)
    mask_w20 <- terra::ifel(!is.na(p_w20) & p_w20 >= 0 & p_w20 <= oever_koppel_m, 1, NA)
    unlink(f_prox_w20)
    rm(p_w20)
    
    f_l_bin   <- tempfile(fileext = ".tif")
    terra::writeRaster(terra::ifel(!is.na(r_land_geclusterd) & r_land_geclusterd > 0, 1, NA), f_l_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    f_prox_l20 <- hulp_proximity_gdal(f_l_bin, oever_koppel_m)
    unlink(f_l_bin)
    
    p_l20 <- terra::rast(f_prox_l20)
    mask_l20 <- terra::ifel(!is.na(p_l20) & p_l20 >= 0 & p_l20 <= oever_koppel_m, 1, NA)
    unlink(f_prox_l20)
    rm(p_l20)
    
    r_wat_20m  <- terra::mask(r_wat_geclusterd, mask_l20)
    r_land_20m <- terra::mask(r_land_geclusterd, mask_w20)
    rm(mask_w20, mask_l20)
    
    # B. 500 m inbedding in ruim complex
    f_rm_bin  <- tempfile(fileext = ".tif")
    terra::writeRaster(terra::ifel(!is.na(r_ruim_valide) & r_ruim_valide > 0, 1, NA), f_rm_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    f_prox_rm500 <- hulp_proximity_gdal(f_rm_bin, koppel_m)
    unlink(f_rm_bin)
    
    p_rm <- terra::rast(f_prox_rm500)
    mask_rm500 <- terra::ifel(!is.na(p_rm) & p_rm >= 0 & p_rm <= koppel_m, 1, NA)
    unlink(f_prox_rm500)
    rm(p_rm)
    
    water_ingebed <- terra::mask(r_wat_20m, mask_rm500)
    land_ingebed  <- terra::mask(r_land_20m, mask_rm500)
    rm(mask_rm500, r_wat_20m, r_land_20m)
    
    # C. Water-first afsnijding: land afkappen op 500 m rondom het goedgekeurde water
    if (!all(is.na(suppressWarnings(terra::minmax(water_ingebed))))) {
      f_wi_bin  <- tempfile(fileext = ".tif")
      terra::writeRaster(terra::ifel(!is.na(water_ingebed) & water_ingebed > 0, 1, NA), f_wi_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
      f_prox_wi500 <- hulp_proximity_gdal(f_wi_bin, koppel_m)
      unlink(f_wi_bin)
      
      p_wi <- terra::rast(f_prox_wi500)
      mask_wi500 <- terra::ifel(!is.na(p_wi) & p_wi >= 0 & p_wi <= koppel_m, 1, NA)
      unlink(f_prox_wi500)
      rm(p_wi)
      
      water_finaal <- water_ingebed
      land_finaal  <- terra::mask(land_ingebed, mask_wi500)
      rm(mask_wi500)
    } else {
      water_finaal <- master_grid * NA
      land_finaal  <- master_grid * NA
    }
    rm(water_ingebed, land_ingebed)
  } else {
    water_finaal <- master_grid * NA
    land_finaal  <- master_grid * NA
  }
  rm(r_wat_geclusterd, r_land_geclusterd)
  gc()
  
  ha_wat_fin  <- calc_ha_exact(water_finaal)
  ha_land_fin <- calc_ha_exact(land_finaal)
  message(sprintf("       -> Voortplantingswater na koppeling en inbedding: %.2f ha", ha_wat_fin))
  message(sprintf("       -> Aanpalend land binnen 500m van goedgekeurd water: %.2f ha", ha_land_fin))
  
  # ----------------------------------------------------------------------------
  # 6. SAMENVOEGING TOT DEFINITIEF LEEFGEBIED & NETWERKEN
  # ----------------------------------------------------------------------------
  message("    🧩 Totaal leefgebiedscomplex samenstellen (vennen + jachtland)...")
  final_opp <- terra::cover(water_finaal, land_finaal)
  final_opp <- terra::clamp(final_opp, upper = 1.0)
  
  ha_tot <- calc_ha_exact(final_opp)
  message(sprintf("       -> Definitief leefgebied Zwarte heidelibel: %.2f ha", ha_tot))
  
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
    land_500m           = land_finaal,
    ruim_landschap_40ha = r_ruim_valide
  ))
}
