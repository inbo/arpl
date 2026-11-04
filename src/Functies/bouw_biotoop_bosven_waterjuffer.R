# ==============================================================================
# bouw_biotoop_bosven_waterjuffer.R
# Archetype voor bosven-waterjuffers (vennen >= 0.02ha x land >= 10ha met 40% bos op 500m x ruim land >= 50ha met 50% bos)
# Geschikt voor: Maanwaterjuffer (Coenagrion lunulatum), etc.
# ==============================================================================

bouw_biotoop_bosven_waterjuffer <- function(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
  bestandsnaam <- tolower(gsub(" ", "", soort_naam))
  soort_cache  <- file.path(paths$dir_cache, bestandsnaam)
  dir.create(soort_cache, showWarnings = FALSE, recursive = TRUE)
  
  f_cache_rds <- file.path(soort_cache, "01_raw_parquet_layers.rds")
  
  # Parameters ophalen uit Excel
  rij_wat <- params_soort %>% dplyr::filter(grepl("waterbiotoop_bwk1", tolower(Type))) %>% dplyr::slice(1)
  min_opp_water <- as.numeric(ifelse(nrow(rij_wat) > 0 && !is.na(rij_wat$MinOpp_ha), rij_wat$MinOpp_ha, 0.02))
  zoekafstand_wat_m <- as.numeric(ifelse(nrow(rij_wat) > 0 && !is.na(rij_wat$AfstandBiotopen_m), rij_wat$AfstandBiotopen_m, 10))
  
  rij_land <- params_soort %>% dplyr::filter(grepl("landbiotoop_bwk1", tolower(Type))) %>% dplyr::slice(1)
  min_opp_land <- as.numeric(ifelse(nrow(rij_land) > 0 && !is.na(rij_land$MinOpp_ha), rij_land$MinOpp_ha, 10.0))
  zoekafstand_land_m <- as.numeric(ifelse(nrow(rij_land) > 0 && !is.na(rij_land$AfstandBiotopen_m), rij_land$AfstandBiotopen_m, 20))
  
  rij_ruim <- params_soort %>% dplyr::filter(grepl("landbiotoop_bwk_ruim", tolower(Type))) %>% dplyr::slice(1)
  min_opp_ruim <- as.numeric(ifelse(nrow(rij_ruim) > 0 && !is.na(rij_ruim$MinOpp_ha), rij_ruim$MinOpp_ha, 50.0))
  zoekafstand_ruim_m <- as.numeric(ifelse(nrow(rij_ruim) > 0 && !is.na(rij_ruim$AfstandBiotopen_m), rij_ruim$AfstandBiotopen_m, 50))
  
  koppel_m <- 500 # 500 m vlieg- en foerageeractieradius vanaf water
  
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
  
  r_wat1_raw       <- if ("waterbiotoop_bwk1"      %in% names(lagen)) lagen[["waterbiotoop_bwk1"]]      else (master_grid * NA)
  r_not_wat_raw    <- if ("waterbiotoop_notbwk"    %in% names(lagen)) lagen[["waterbiotoop_notbwk"]]    else (master_grid * NA)
  r_land1_raw      <- if ("landbiotoop_bwk1"       %in% names(lagen)) lagen[["landbiotoop_bwk1"]]       else (master_grid * NA)
  r_not_land_raw   <- if ("landbiotoop_notbwk"     %in% names(lagen)) lagen[["landbiotoop_notbwk"]]     else (master_grid * NA)
  r_bos_nabij_raw  <- if ("landbiotoop_nabij_bos"  %in% names(lagen)) lagen[["landbiotoop_nabij_bos"]]  else (master_grid * NA)
  r_ruim_raw       <- if ("landbiotoop_bwk_ruim"   %in% names(lagen)) lagen[["landbiotoop_bwk_ruim"]]   else (master_grid * NA)
  r_bos_ruim_raw   <- if ("landbiotoop_ruim_bos"   %in% names(lagen)) lagen[["landbiotoop_ruim_bos"]]   else (master_grid * NA)
  
  # ----------------------------------------------------------------------------
  # 2. WATERBIOTOOP SCHONEN, KWALITEITSDREMPEL (>= 20%) & CLUSTEREN (>= 0.02 HA)
  # ----------------------------------------------------------------------------
  message(sprintf("    💧 Vennen schonen en filteren (fractie >= 0.20, drempel >= %.2f ha)...", min_opp_water))
  
  # Kwaliteitsfilter op initiële fractie >= 0.20
  r_wat1_gefilterd <- terra::ifel(!is.na(r_wat1_raw) & r_wat1_raw >= 0.20, r_wat1_raw, NA)
  
  # Negatieve aftrek van eutroof water/visvijvers
  r_not_wat_clean <- terra::ifel(is.na(r_not_wat_raw), 0, r_not_wat_raw)
  r_wat_geschoond <- terra::clamp(r_wat1_gefilterd - r_not_wat_clean, lower = 0, upper = 1.0)
  r_wat_geschoond <- terra::ifel(r_wat_geschoond > 0, r_wat_geschoond, NA)
  rm(r_not_wat_clean, r_wat1_gefilterd, r_wat1_raw, r_not_wat_raw)
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
  message(sprintf("       -> Goedgekeurde voortplantingsvennen (>= 0.02 ha): %.2f ha", ha_wat_gecl))
  
  # ----------------------------------------------------------------------------
  # 3. RUIM LANDBIOTOOP CLUSTEREN (>= 50 HA OP 50 M) & BOSSTRUCTUUR (> 50%)
  # ----------------------------------------------------------------------------
  message(sprintf("    🌾 Ruim landschap clusteren (>= %.1f ha) en filteren op > 50%% bos...", min_opp_ruim))
  res_ruim <- cluster_filter_universeel(
    masker            = r_ruim_raw,
    opp_laag          = r_ruim_raw,
    min_opp_ha        = min_opp_ruim,
    min_totale_opp_ha = NA,
    afstand_m         = zoekafstand_ruim_m,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_ruim_50ha  <- res_ruim$raster
  cl_ruim_50ha <- res_ruim$clusters
  rm(res_ruim, r_ruim_raw)
  gc()
  
  # Toetsen op > 50% bos in ruime schil
  if (!all(is.na(suppressWarnings(terra::minmax(cl_ruim_50ha))))) {
    r_bos_ruim_clean <- terra::ifel(!is.na(r_bos_ruim_raw) & r_bos_ruim_raw > 0, r_bos_ruim_raw, 0)
    
    df_ruim_tot <- terra::zonal(r_ruim_50ha, cl_ruim_50ha, fun = "sum", na.rm = TRUE)
    colnames(df_ruim_tot) <- c("ID", "Totaal_Sum")
    
    df_bos_ruim <- terra::zonal(r_bos_ruim_clean, cl_ruim_50ha, fun = "sum", na.rm = TRUE)
    colnames(df_bos_ruim) <- c("ID", "Bos_Sum")
    
    dt_ruim <- merge(data.table::as.data.table(df_ruim_tot), data.table::as.data.table(df_bos_ruim), by = "ID", all.x = TRUE)
    dt_ruim[is.na(Bos_Sum), Bos_Sum := 0]
    dt_ruim[, Pct_Bos := (Bos_Sum / Totaal_Sum) * 100]
    
    valide_ruim_ids <- dt_ruim[!is.na(ID) & Pct_Bos > 50, ID]
    rm(df_ruim_tot, df_bos_ruim, dt_ruim, r_bos_ruim_clean)
    
    if (length(valide_ruim_ids) > 0) {
      r_ruim_valide <- terra::mask(r_ruim_50ha, cl_ruim_50ha %in% as.integer(valide_ruim_ids))
    } else {
      r_ruim_valide <- master_grid * NA
    }
  } else {
    r_ruim_valide <- master_grid * NA
  }
  rm(r_ruim_50ha, cl_ruim_50ha, r_bos_ruim_raw)
  gc()
  
  ha_ruim_val <- calc_ha_exact(r_ruim_valide)
  message(sprintf("       -> Ruim landschapscomplex (> 50%% bos): %.2f ha", ha_ruim_val))
  
  # ----------------------------------------------------------------------------
  # 4. WATER-FIRST FILTERING: LAND AFSNIJDEN BINNEN 500 M VAN WATER
  # ----------------------------------------------------------------------------
  message(sprintf("    🎯 Land afsnijden binnen %dm van vennen (>= %.1f ha, >= 40%% bos)...", as.integer(koppel_m), min_opp_land))
  
  # Zuiveren landbiotoop van droge heide/stuifzand
  r_not_land_clean <- terra::ifel(is.na(r_not_land_raw), 0, r_not_land_raw)
  r_land_geschoond <- terra::clamp(r_land1_raw - r_not_land_clean, lower = 0, upper = 1.0)
  r_land_geschoond <- terra::ifel(r_land_geschoond > 0, r_land_geschoond, NA)
  rm(r_not_land_clean, r_land1_raw, r_not_land_raw)
  gc()
  
  heeft_water <- !all(is.na(suppressWarnings(terra::minmax(r_wat_geclusterd))))
  heeft_land  <- !all(is.na(suppressWarnings(terra::minmax(r_land_geschoond))))
  
  if (heeft_water && heeft_land) {
    # 500m buffer rondom goedgekeurde vennen via GDAL proximity
    f_w_bin  <- tempfile(fileext = ".tif")
    terra::writeRaster(terra::ifel(!is.na(r_wat_geclusterd) & r_wat_geclusterd > 0, 1, NA), f_w_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    f_prox_w <- hulp_proximity_gdal(f_w_bin, koppel_m)
    unlink(f_w_bin)
    
    p_w <- terra::rast(f_prox_w)
    mask_w500 <- terra::ifel(!is.na(p_w) & p_w >= 0 & p_w <= koppel_m, 1, NA)
    unlink(f_prox_w)
    rm(p_w)
    
    # Land afsnijden op 500 m
    r_land_500m <- terra::mask(r_land_geschoond, mask_w500)
    rm(mask_w500)
    
    # Clusteren van land binnen de buffer (>= 10 ha op 20 m)
    res_land10 <- cluster_filter_universeel(
      masker            = r_land_500m,
      opp_laag          = r_land_500m,
      min_opp_ha        = min_opp_land,
      min_totale_opp_ha = NA,
      afstand_m         = zoekafstand_land_m,
      dispersie_m       = 0,
      werkelijk         = TRUE
    )
    r_land_10ha  <- res_land10$raster
    cl_land_10ha <- res_land10$clusters
    rm(res_land10, r_land_500m)
    gc()
    
    # Toetsen op minstens 40% bos in het 10 ha landcluster
    if (!all(is.na(suppressWarnings(terra::minmax(cl_land_10ha))))) {
      r_bos_nabij_clean <- terra::ifel(!is.na(r_bos_nabij_raw) & r_bos_nabij_raw > 0, r_bos_nabij_raw, 0)
      
      df_land_tot <- terra::zonal(r_land_10ha, cl_land_10ha, fun = "sum", na.rm = TRUE)
      colnames(df_land_tot) <- c("ID", "Totaal_Sum")
      
      df_bos_tot <- terra::zonal(r_bos_nabij_clean, cl_land_10ha, fun = "sum", na.rm = TRUE)
      colnames(df_bos_tot) <- c("ID", "Bos_Sum")
      
      dt_land <- merge(data.table::as.data.table(df_land_tot), data.table::as.data.table(df_bos_tot), by = "ID", all.x = TRUE)
      dt_land[is.na(Bos_Sum), Bos_Sum := 0]
      dt_land[, Pct_Bos := (Bos_Sum / Totaal_Sum) * 100]
      
      valide_land_ids <- dt_land[!is.na(ID) & Pct_Bos >= 40, ID]
      rm(df_land_tot, df_bos_tot, dt_land, r_bos_nabij_clean)
      
      if (length(valide_land_ids) > 0) {
        r_land_finaal <- terra::mask(r_land_10ha, cl_land_10ha %in% as.integer(valide_land_ids))
      } else {
        r_land_finaal <- master_grid * NA
      }
    } else {
      r_land_finaal <- master_grid * NA
    }
    rm(r_land_10ha, cl_land_10ha)
    gc()
  } else {
    r_land_finaal <- master_grid * NA
  }
  rm(r_land_geschoond, r_bos_nabij_raw)
  gc()
  
  ha_land_fin <- calc_ha_exact(r_land_finaal)
  message(sprintf("       -> Aanpalend beschut land (500m, >= 10 ha, >= 40%% bos): %.2f ha", ha_land_fin))
  
  # ----------------------------------------------------------------------------
  # 5. SAMENVOEGING TOT DEFINITIEF LEEFGEBIED & NETWERKEN
  # ----------------------------------------------------------------------------
  message("    🧩 Totaal leefgebiedscomplex samenstellen (vennen + beschut land)...")
  final_opp <- terra::cover(r_wat_geclusterd, r_land_finaal)
  final_opp <- terra::clamp(final_opp, upper = 1.0)
  
  ha_tot <- calc_ha_exact(final_opp)
  message(sprintf("       -> Definitief leefgebied Maanwaterjuffer: %.2f ha", ha_tot))
  
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
    water               = r_wat_geclusterd,
    land_500m           = r_land_finaal,
    ruim_landschap_50ha = r_ruim_valide
  ))
}
