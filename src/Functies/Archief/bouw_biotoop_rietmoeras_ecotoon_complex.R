# ==============================================================================
# bouw_biotoop_rietmoeras_ecotoon_complex.R
# Archetype voor grote rietmoerasvogels (3 schillen: water >= 10ha, ecotoon >= 0.5ha, 
# nestkern >= 0.5ha x metapopulatie >= 25ha op 500m)
# Geschikt voor: Roerdomp, etc.
# ==============================================================================

bouw_biotoop_rietmoeras_ecotoon_complex <- function(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
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
    lagen_packed <- lapply(lagen, terra::wrap)
    saveRDS(lagen_packed, f_cache_rds)
    message("    💾 [CACHE] Parquet-lagen veilig bewaard.")
  }
  
  r_water_raw    <- if ("openwater" %in% names(lagen)) lagen[["openwater"]] else (master_grid * NA)
  r_riet_raw     <- if ("rietland"  %in% names(lagen)) lagen[["rietland"]]  else (master_grid * NA)
  r_zegge_raw    <- if ("zegge"     %in% names(lagen)) lagen[["zegge"]]     else (master_grid * NA)
  r_moeras_raw   <- if ("moeras"    %in% names(lagen)) lagen[["moeras"]]    else (master_grid * NA)
  r_gras_raw     <- if ("grasland"  %in% names(lagen)) lagen[["grasland"]]  else (master_grid * NA)
  r_struweel_raw <- if ("struweel"  %in% names(lagen)) lagen[["struweel"]]  else (master_grid * NA)
  
  # ----------------------------------------------------------------------------
  # 2. NEGATIEVE STOORLAGEN (HOOGGROEN & ESTUARIUM 1130)
  # ----------------------------------------------------------------------------
  message("    🚫 Stoorlagen isoleren (hooggroen en estuaria 1130)...")
  
  # A. Hooggroen (strikt klasse 1 in Groenkaart 2021)
  f_groen_al     <- hulp_align_gdal(paths$groenkaart, master_grid)
  r_groen        <- terra::rast(f_groen_al)
  r_bomen_direct <- terra::ifel(!is.na(r_groen) & r_groen == 1, 1, NA)
  unlink(f_groen_al)
  rm(r_groen)
  
  # B. Estuarium (code 1130 veilig via Parquet helper)
  r_est_bin <- vul_raster_vanuit_parquet(
    ds_parquet  = ds_parquet,
    codes_exact = "1130",
    master_grid = master_grid,
    waarde_vast = 1
  )
  
  stoorlaag <- terra::cover(r_bomen_direct, r_est_bin)
  rm(r_bomen_direct, r_est_bin)
  gc()
  
  # ----------------------------------------------------------------------------
  # 3. BASISHABITAT SAMENSTELLEN EN SCHONEN
  # ----------------------------------------------------------------------------
  message("    🌾 Basishabitat (6 componenten) samenvoegen en schonen van stoorlagen...")
  
  # Samenvoegen basishabitat via cover conform Bert
  r_bio_bruto <- terra::cover(r_water_raw, r_riet_raw)
  r_bio_bruto <- terra::cover(r_bio_bruto, r_zegge_raw)
  r_bio_bruto <- terra::cover(r_bio_bruto, r_moeras_raw)
  r_bio_bruto <- terra::cover(r_bio_bruto, r_gras_raw)
  r_bio_bruto <- terra::cover(r_bio_bruto, r_struweel_raw)
  r_bio_bruto <- terra::clamp(r_bio_bruto, upper = 1.0)
  
  rm(r_zegge_raw, r_moeras_raw, r_gras_raw, r_struweel_raw)
  gc()
  
  # Filteren met de stoorlaag
  r_leefgebied1 <- terra::mask(r_bio_bruto, stoorlaag, inverse = TRUE)
  r_riet_clean  <- terra::mask(r_riet_raw,  stoorlaag, inverse = TRUE)
  r_water_clean <- terra::mask(r_water_raw, stoorlaag, inverse = TRUE)
  rm(r_bio_bruto, stoorlaag)
  gc()
  
  ha_leef1 <- calc_ha_exact(r_leefgebied1)
  message(sprintf("       -> Basishabitat na stoorfilter: %.2f ha", ha_leef1))
  
  # ----------------------------------------------------------------------------
  # 4. SCHIL 1: GROTE OPEN WATERPARTIJEN (>= 10 HA) & 100M SCHIL
  # ----------------------------------------------------------------------------
  message("    🌊 Schil 1: Grote waterplassen clusteren (>= 10 ha) en 100m schil toepassen...")
  
  res_water10 <- cluster_filter_universeel(
    masker            = r_water_clean,
    opp_laag          = r_water_clean,
    min_opp_ha        = 10,
    min_totale_opp_ha = NA,
    afstand_m         = 50,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_water_10ha <- res_water10$raster
  rm(res_water10)
  gc()
  
  if (!all(is.na(suppressWarnings(terra::minmax(r_water_10ha))))) {
    f_w10_bin <- tempfile(fileext = ".tif")
    terra::writeRaster(terra::ifel(!is.na(r_water_10ha) & r_water_10ha > 0, 1, NA), f_w10_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    f_prox_w10 <- hulp_proximity_gdal(f_w10_bin, 100)
    unlink(f_w10_bin)
    
    p_w10      <- terra::rast(f_prox_w10)
    mask_w100m <- terra::ifel(!is.na(p_w10) & p_w10 >= 0 & p_w10 <= 100, 1, NA)
    unlink(f_prox_w10)
    rm(p_w10)
    gc()
    
    r_leefgebied2 <- terra::mask(r_leefgebied1, mask_w100m)
    rm(mask_w100m)
  } else {
    r_leefgebied2 <- master_grid * NA
  }
  rm(r_water_10ha, r_leefgebied1)
  gc()
  
  # ----------------------------------------------------------------------------
  # 5. SCHIL 2: ROBUUSTE RIET-WATERKANTEN (>= 0.5 HA WATER OP 10M BUFFER) & 100M SCHIL
  # ----------------------------------------------------------------------------
  message("    📏 Schil 2: Riet-waterkanten (>= 0.5 ha) clusteren en 100m schil toepassen...")
  
  r_riet_bin <- terra::ifel(!is.na(r_riet_clean) & r_riet_clean > 0, 1, NA)
  
  if (!all(is.na(suppressWarnings(terra::minmax(r_riet_bin)))) && !all(is.na(suppressWarnings(terra::minmax(r_water_clean))))) {
    # 10m buffer rond riet via GDAL proximity
    f_r_bin    <- tempfile(fileext = ".tif")
    terra::writeRaster(r_riet_bin, f_r_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    f_prox_r10 <- hulp_proximity_gdal(f_r_bin, 10)
    unlink(f_r_bin)
    
    p_r10        <- terra::rast(f_prox_r10)
    mask_riet10m <- terra::ifel(!is.na(p_r10) & p_r10 >= 0 & p_r10 <= 10, 1, NA)
    unlink(f_prox_r10)
    rm(p_r10)
    gc()
    
    # Waterkant = open water binnen 10m van rietland
    r_waterkant <- terra::mask(r_water_clean, mask_riet10m)
    rm(mask_riet10m)
    gc()
    
    # Clusteren waterkant op >= 0.5 ha (50m netwerk)
    res_wk <- cluster_filter_universeel(
      masker            = r_waterkant,
      opp_laag          = r_waterkant,
      min_opp_ha        = 0.5,
      min_totale_opp_ha = NA,
      afstand_m         = 50,
      dispersie_m       = 0,
      werkelijk         = TRUE
    )
    r_wk_05ha <- res_wk$raster
    rm(res_wk, r_waterkant)
    gc()
    
    if (!all(is.na(suppressWarnings(terra::minmax(r_wk_05ha))))) {
      f_wk_bin  <- tempfile(fileext = ".tif")
      terra::writeRaster(terra::ifel(!is.na(r_wk_05ha) & r_wk_05ha > 0, 1, NA), f_wk_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
      f_prox_wk <- hulp_proximity_gdal(f_wk_bin, 100)
      unlink(f_wk_bin)
      
      p_wk        <- terra::rast(f_prox_wk)
      mask_wk100m <- terra::ifel(!is.na(p_wk) & p_wk >= 0 & p_wk <= 100, 1, NA)
      unlink(f_prox_wk)
      rm(p_wk)
      gc()
      
      r_leefgebied3 <- terra::mask(r_leefgebied2, mask_wk100m)
      rm(mask_wk100m)
    } else {
      r_leefgebied3 <- master_grid * NA
    }
    rm(r_wk_05ha)
  } else {
    r_leefgebied3 <- master_grid * NA
  }
  rm(r_water_clean, r_leefgebied2)
  gc()
  
  # ----------------------------------------------------------------------------
  # 6. SCHIL 3: NESTKERNEN IN ZUIVER RIETLAND (>= 0.5 HA) & 100M SCHIL
  # ----------------------------------------------------------------------------
  message("    🪺 Schil 3: Zuivere rietnestkernen (>= 0.5 ha) clusteren en 100m schil toepassen...")
  
  res_nest <- cluster_filter_universeel(
    masker            = r_riet_clean,
    opp_laag          = r_riet_clean,
    min_opp_ha        = 0.5,
    min_totale_opp_ha = NA,
    afstand_m         = 50,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_nest_05ha <- res_nest$raster
  rm(res_nest, r_riet_clean, r_riet_bin)
  gc()
  
  if (!all(is.na(suppressWarnings(terra::minmax(r_nest_05ha))))) {
    f_nest_bin  <- tempfile(fileext = ".tif")
    terra::writeRaster(terra::ifel(!is.na(r_nest_05ha) & r_nest_05ha > 0, 1, NA), f_nest_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    f_prox_nest <- hulp_proximity_gdal(f_nest_bin, 100)
    unlink(f_nest_bin)
    
    p_nest        <- terra::rast(f_prox_nest)
    mask_nest100m <- terra::ifel(!is.na(p_nest) & p_nest >= 0 & p_nest <= 100, 1, NA)
    unlink(f_prox_nest)
    rm(p_nest)
    gc()
    
    r_leefgebied4 <- terra::mask(r_leefgebied3, mask_nest100m)
    rm(mask_nest100m)
  } else {
    r_leefgebied4 <- master_grid * NA
  }
  rm(r_nest_05ha, r_leefgebied3)
  gc()
  
  ha_leef4 <- calc_ha_exact(r_leefgebied4)
  message(sprintf("       -> Habitat na de 3 schillen (water, ecotoon, nestkern): %.2f ha", ha_leef4))
  
  # ----------------------------------------------------------------------------
  # 7. FINALE METAPOPULATIE-CLUSTERING (>= 25 HA OP 500M NETWERK)
  # ----------------------------------------------------------------------------
  message("    🧩 Finale metapopulatie-clustering (>= 25 ha, 500m netwerkafstand)...")
  res_meta <- cluster_filter_universeel(
    masker            = r_leefgebied4,
    opp_laag          = r_leefgebied4,
    min_opp_ha        = 25,
    min_totale_opp_ha = NA,
    afstand_m         = 500,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  final_opp <- res_meta$raster
  cl_opp    <- res_meta$clusters
  rm(res_meta, r_leefgebied4)
  gc()
  
  ha_finaal <- calc_ha_exact(final_opp)
  message(sprintf("       -> Definitief metapopulatie-leefgebied Roerdomp: %.2f ha", ha_finaal))
  
  # Deellagen isoleren binnen het definitieve netwerkcomplex
  if (!all(is.na(suppressWarnings(terra::minmax(cl_opp))))) {
    riet_uitvoer  <- terra::mask(r_riet_raw, cl_opp)
    water_uitvoer <- terra::mask(r_water_raw, cl_opp)
  } else {
    riet_uitvoer  <- master_grid * NA
    water_uitvoer <- master_grid * NA
  }
  rm(r_riet_raw, r_water_raw)
  gc()
  
  return(list(
    final_opp = final_opp,
    cl_opp    = cl_opp,
    rietland  = riet_uitvoer,
    openwater = water_uitvoer
  ))
}
