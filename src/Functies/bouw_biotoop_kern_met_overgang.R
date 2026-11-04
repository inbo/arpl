# ==============================================================================
# bouw_biotoop_kern_met_overgang.R - LEEFGEBIEDENMODEL VOOR KOMMAVLINDER
# ==============================================================================

bouw_biotoop_kern_met_overgang <- function(params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
  bestandsnaam <- "kommavlinder"
  soort_cache  <- file.path(paths$dir_cache, bestandsnaam)
  dir.create(soort_cache, showWarnings = FALSE, recursive = TRUE)
  
  f_cache_rds <- file.path(soort_cache, "01_raw_parquet_layers.rds")
  
  # 1. PARQUET EXTRACTIE (met persistente cache)
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
  
  r_bwk <- if ("bwk" %in% names(lagen)) lagen[["bwk"]] else (master_grid * NA)
  r_ku  <- if ("ku"  %in% names(lagen)) lagen[["ku"]]  else (master_grid * NA)
  
  ha_bwk_ruw <- calc_ha_exact(r_bwk)
  ha_ku_ruw  <- calc_ha_exact(r_ku)
  message(sprintf("       -> Primair kernbiotoop (bwk): %.2f ha", ha_bwk_ruw))
  message(sprintf("       -> Secundair biotoop (ku) vóór filter: %.2f ha", ha_ku_ruw))
  
  # 2. 'KU' FILTEREN OP BASIS VAN AFSTAND TOT 'BWK' (50M) VIA GDAL
  message("    🎯 'ku'-cellen filteren (enkel binnen <= 50m van primair bwk-biotoop)...")
  r_bwk_bin <- terra::ifel(!is.na(r_bwk) & r_bwk > 0, 1, NA)
  
  if (!all(is.na(suppressWarnings(terra::minmax(r_bwk_bin))))) {
    f_bwk_bin  <- tempfile(fileext = ".tif")
    terra::writeRaster(r_bwk_bin, f_bwk_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 0)
    
    # Gebruik de beproefde hulp_proximity_gdal
    f_prox_50m <- hulp_proximity_gdal(f_bwk_bin, 50)
    unlink(f_bwk_bin)
    
    p50       <- terra::rast(f_prox_50m)
    mask_bwk50 <- terra::ifel(!is.na(p50) & p50 >= 0 & p50 <= 50, 1, NA)
    unlink(f_prox_50m)
    rm(p50)
    gc()
    
    r_ku_gefilterd <- terra::mask(r_ku, mask_bwk50)
    rm(mask_bwk50)
    gc()
  } else {
    r_ku_gefilterd <- master_grid * NA
  }
  rm(r_bwk_bin)
  gc()
  
  ha_ku_ok <- calc_ha_exact(r_ku_gefilterd)
  message(sprintf("       -> Secundair biotoop (ku) binnen 50m van bwk: %.2f ha", ha_ku_ok))
  
  # 3. SOMMATIEF SAMENVOEGEN TOT COMBINATIEBIOTOOP
  message("    ➕ Kernbiotoop en goedgekeurd 'ku'-biotoop sommatief samenvoegen...")
  w_bwk <- terra::ifel(is.na(r_bwk), 0, r_bwk)
  w_ku  <- terra::ifel(is.na(r_ku_gefilterd), 0, r_ku_gefilterd)
  
  r_som <- w_bwk + w_ku
  r_gecombineerd <- terra::ifel(r_som > 0, terra::clamp(r_som, upper = 1.0), NA)
  rm(w_bwk, w_ku, r_som)
  gc()
  
  ha_comb <- calc_ha_exact(r_gecombineerd)
  message(sprintf("       -> Totaal gecombineerd biotoop vóór clustering: %.2f ha", ha_comb))
  
  # 4. CLUSTEREN EN OPPERVLAKTEFILTERING (>= 5 ha, 10m zoekafstand)
  rij_params <- params_soort %>% dplyr::slice(1)
  min_opp_ha <- as.numeric(ifelse(is.na(rij_params$MinOpp_ha), 5, rij_params$MinOpp_ha))
  zoekafstand_m <- as.numeric(ifelse(is.na(rij_params$AfstandBiotopen_m), 10, rij_params$AfstandBiotopen_m))
  
  message(sprintf("    🏝️ Clusteren op minimumoppervlakte (>= %.1f ha, zoekafstand %dm)...", min_opp_ha, zoekafstand_m))
  res_cluster <- cluster_filter_universeel(
    masker            = r_gecombineerd,
    opp_laag          = r_gecombineerd,
    min_opp_ha        = min_opp_ha,
    min_totale_opp_ha = NA,
    afstand_m         = zoekafstand_m,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  
  final_opp <- res_cluster$raster
  cl_opp    <- res_cluster$clusters
  rm(res_cluster, r_gecombineerd)
  gc()
  
  # Deellagen isoleren binnen de goedgekeurde clusters
  bwk_finaal <- terra::mask(r_bwk, cl_opp)
  ku_finaal  <- terra::mask(r_ku_gefilterd, cl_opp)
  
  return(list(
    final_opp  = final_opp,
    cl_opp     = cl_opp,
    deellaag_a = bwk_finaal,
    deellaag_b = ku_finaal
  ))
}