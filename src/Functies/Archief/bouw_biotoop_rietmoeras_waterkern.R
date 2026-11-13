# ==============================================================================
# bouw_biotoop_rietmoeras_waterkern.R
# Archetype voor kleine moerasvogels (moeras >= 10ha x populierenfilter x water >= 1ha binnen 100m)
# Geschikt voor: Woudaap, etc.
# ==============================================================================

bouw_biotoop_rietmoeras_waterkern <- function(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
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
  
  r_vp_raw    <- if ("voortplanting_bwk" %in% names(lagen)) lagen[["voortplanting_bwk"]] else (master_grid * NA)
  r_not_pop   <- if ("not_populier"       %in% names(lagen)) lagen[["not_populier"]]       else (master_grid * NA)
  r_water_raw <- if ("openwater"          %in% names(lagen)) lagen[["openwater"]]          else (master_grid * NA)
  
  # ----------------------------------------------------------------------------
  # 2. STAP 1: MOERASHABITAT SCHONEN VAN POPULIEREN & CLUSTEREN OP >= 10 HA (50M)
  # ----------------------------------------------------------------------------
  message("    🌾 Moerashabitat schonen van populierenaanplanten en clusteren op >= 10 ha (50m netwerk)...")
  
  if (!all(is.na(suppressWarnings(terra::minmax(r_not_pop))))) {
    r_vp_gezuiverd <- terra::mask(r_vp_raw, r_not_pop > 0, inverse = TRUE)
  } else {
    r_vp_gezuiverd <- r_vp_raw
  }
  rm(r_vp_raw, r_not_pop)
  gc()
  
  res_vp <- cluster_filter_universeel(
    masker            = r_vp_gezuiverd,
    opp_laag          = r_vp_gezuiverd,
    min_opp_ha        = 10,
    min_totale_opp_ha = NA,
    afstand_m         = 50,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_vp_10ha <- res_vp$raster
  rm(res_vp, r_vp_gezuiverd)
  gc()
  
  ha_vp_10ha <- calc_ha_exact(r_vp_10ha)
  message(sprintf("       -> Goedgekeurd moerashabitat (>= 10 ha): %.2f ha", ha_vp_10ha))
  
  # ----------------------------------------------------------------------------
  # 3. STAP 2: OPEN WATER FILTEREN OP >= 1 HA (10M NETWERK) & 100M SCHIL
  # ----------------------------------------------------------------------------
  message("    🌊 Open waterpartijen clusteren op >= 1 ha (10m) en koppelen binnen 100m...")
  
  res_water <- cluster_filter_universeel(
    masker            = r_water_raw,
    opp_laag          = r_water_raw,
    min_opp_ha        = 1,
    min_totale_opp_ha = NA,
    afstand_m         = 10,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_water_1ha <- res_water$raster
  rm(res_water, r_water_raw)
  gc()
  
  heeft_vp    <- !all(is.na(suppressWarnings(terra::minmax(r_vp_10ha))))
  heeft_water <- !all(is.na(suppressWarnings(terra::minmax(r_water_1ha))))
  
  if (heeft_vp && heeft_water) {
    r_w1_bin  <- terra::ifel(!is.na(r_water_1ha) & r_water_1ha > 0, 1, NA)
    f_w1_bin  <- tempfile(fileext = ".tif")
    terra::writeRaster(r_w1_bin, f_w1_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 0)
    f_prox_w1 <- hulp_proximity_gdal(f_w1_bin, 100)
    unlink(f_w1_bin)
    rm(r_w1_bin)
    
    p_w1      <- terra::rast(f_prox_w1)
    mask_w100 <- terra::ifel(!is.na(p_w1) & p_w1 >= 0 & p_w1 <= 100, 1, NA)
    unlink(f_prox_w1)
    rm(p_w1)
    gc()
    
    # Moeras behouden dat binnen 100m van >= 1 ha water ligt
    r_vp_bij_water <- terra::mask(r_vp_10ha, mask_w100)
    rm(mask_w100)
    gc()
  } else {
    r_vp_bij_water <- master_grid * NA
  }
  rm(r_vp_10ha, r_water_1ha)
  gc()
  
  ha_vp_water <- calc_ha_exact(r_vp_bij_water)
  message(sprintf("       -> Moeras binnen 100m van open water: %.2f ha", ha_vp_water))
  
  # ----------------------------------------------------------------------------
  # 4. STAP 3: FINALE NETWERKCLUSTERING (>= 10 HA OP 100M NETWERK)
  # ----------------------------------------------------------------------------
  message("    🧩 Finale netwerkclustering (>= 10 ha, 100m netwerkafstand)...")
  res_finaal <- cluster_filter_universeel(
    masker            = r_vp_bij_water,
    opp_laag          = r_vp_bij_water,
    min_opp_ha        = 10,
    min_totale_opp_ha = NA,
    afstand_m         = 100,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  final_opp <- res_finaal$raster
  cl_opp    <- res_finaal$clusters
  rm(res_finaal, r_vp_bij_water)
  gc()
  
  return(list(
    final_opp = final_opp,
    cl_opp    = cl_opp
  ))
}