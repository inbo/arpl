# ==============================================================================
# bouw_biotoop_bos_metapopulatie_netwerk.R
# Archetype: Boskernen (>= 300 ha) x satellietbossen (50-300 ha) x 4 km interactie
# Geschikt voor: Zwarte specht, etc.
# ==============================================================================

bouw_biotoop_bos_metapopulatie_netwerk <- function(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
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
    lagen <- readRDS(f_cache_rds)
  } else {
    message(paste0("    ⏳ Parquet-lagen extraheren voor ", soort_naam, "..."))
    lagen <- list()
    for (tp in types_in_tabel) {
      r_tp <- extract_bwk_laag(soort_naam, tp, ds_parquet, df_bwk_codes, master_grid)
      lagen[[tp]] <- if (!is.null(r_tp)) r_tp else (master_grid * NA)
    }
    saveRDS(lagen, f_cache_rds)
    message("    💾 [CACHE] Parquet-lagen veilig bewaard.")
  }
  
  r_bwk_raw <- if ("bwk" %in% names(lagen)) lagen[["bwk"]] else (master_grid * NA)
  
  # ----------------------------------------------------------------------------
  # 2. FUSIE BWK-BOS EN GROENKAART (HOOGGROEN = 1.0)
  # ----------------------------------------------------------------------------
  message("    🌲 BWK-bos combineren met Groenkaart (hooggroen = 1.0)...")
  
  f_groen_al <- hulp_align_gdal(paths$groenkaart, master_grid)
  r_groen    <- terra::rast(f_groen_al)
  r_groen_opp <- terra::ifel(!is.na(r_groen) & (r_groen == 1 | r_groen > 0), 1.0, 0)
  unlink(f_groen_al)
  rm(r_groen)
  
  b_clean <- terra::ifel(is.na(r_bwk_raw), 0, r_bwk_raw)
  r_bos_som <- b_clean + r_groen_opp
  r_bos_basis <- terra::ifel(r_bos_som > 0, terra::clamp(r_bos_som, upper = 1.0), NA)
  rm(b_clean, r_groen_opp, r_bos_som, r_bwk_raw)
  gc()
  
  ha_bos_basis <- calc_ha_exact(r_bos_basis)
  message(sprintf("       -> Totaal bosareaal vóór clustering: %.2f ha", ha_bos_basis))
  
  # ----------------------------------------------------------------------------
  # 3. PRIMAIRE CLUSTERING OP 100M NETWERK (DREMPEL = 0 HA OM ALLE PATCHES TE LABELEN)
  # ----------------------------------------------------------------------------
  message("    🏝️ Primaire 100m clustering uitvoeren voor bosgrootte-analyse...")
  
  res_alle_patches <- cluster_filter_universeel(
    masker            = r_bos_basis,
    opp_laag          = r_bos_basis,
    min_opp_ha        = 0,
    min_totale_opp_ha = NA,
    afstand_m         = 100,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  cl_patches <- res_alle_patches$clusters
  rm(res_alle_patches)
  gc()
  
  # ----------------------------------------------------------------------------
  # 4. CLASSIFICATIE BOSGROOTTE (>= 300 HA VS. 50-300 HA)
  # ----------------------------------------------------------------------------
  message("    📊 Bosclusters classificeren: Grote kernen (>= 300 ha) en satellieten (50-300 ha)...")
  
  if (!all(is.na(suppressWarnings(terra::minmax(cl_patches))))) {
    cells_p <- terra::cells(!is.na(cl_patches) & cl_patches > 0)
    dt_opp <- data.table::data.table(
      ID       = as.integer(cl_patches[cells_p][[1]]),
      Opp_Val  = as.numeric(r_bos_basis[cells_p][[1]])
    )
    rm(cells_p)
    gc()
    
    dt_agg <- dt_opp[, .(Opp_ha = sum(Opp_Val, na.rm = TRUE) * 0.01), by = ID]
    rm(dt_opp)
    
    ids_300ha    <- dt_agg[!is.na(ID) & Opp_ha >= 300, ID]
    ids_50_300ha <- dt_agg[!is.na(ID) & Opp_ha >= 50 & Opp_ha < 300, ID]
    rm(dt_agg)
    
    if (length(ids_300ha) > 0) {
      r_bos_groot <- terra::mask(r_bos_basis, cl_patches %in% as.integer(ids_300ha))
    } else {
      r_bos_groot <- master_grid * NA
    }
    
    if (length(ids_50_300ha) > 0) {
      r_bos_klein <- terra::mask(r_bos_basis, cl_patches %in% as.integer(ids_50_300ha))
    } else {
      r_bos_klein <- master_grid * NA
    }
    rm(cl_patches)
    gc()
  } else {
    r_bos_groot <- master_grid * NA
    r_bos_klein <- master_grid * NA
  }
  
  ha_groot <- calc_ha_exact(r_bos_groot)
  ha_klein <- calc_ha_exact(r_bos_klein)
  message(sprintf("       -> Grote zelfstandige kernen (>= 300 ha): %.2f ha", ha_groot))
  message(sprintf("       -> Satellietbossen (50-300 ha): %.2f ha", ha_klein))
  
  # ----------------------------------------------------------------------------
  # 5. RECONSTRUCTIE NETWERK VIA 4 KM BUFFER (GDAL PROXIMITY)
  # ----------------------------------------------------------------------------
  message("    🎯 4 km netwerkinteractie berekenen voor satellietbossen...")
  
  heeft_groot <- !all(is.na(suppressWarnings(terra::minmax(r_bos_groot))))
  heeft_klein <- !all(is.na(suppressWarnings(terra::minmax(r_bos_klein))))
  
  if (heeft_klein) {
    # Proximity vanaf kleine bossen (om interactie tussen satellieten onderling te vinden)
    f_k_bin <- tempfile(fileext = ".tif")
    terra::writeRaster(terra::ifel(!is.na(r_bos_klein) & r_bos_klein > 0, 1, NA), f_k_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 0)
    f_prox_k <- hulp_proximity_gdal(f_k_bin, 4000)
    unlink(f_k_bin)
    
    p_k <- terra::rast(f_prox_k)
    mask_klein_4km <- terra::ifel(!is.na(p_k) & p_k >= 0 & p_k <= 4000, 1, NA)
    unlink(f_prox_k)
    rm(p_k)
    
    if (heeft_groot) {
      # Proximity vanaf grote bossen
      f_g_bin <- tempfile(fileext = ".tif")
      terra::writeRaster(terra::ifel(!is.na(r_bos_groot) & r_bos_groot > 0, 1, NA), f_g_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 0)
      f_prox_g <- hulp_proximity_gdal(f_g_bin, 4000)
      unlink(f_g_bin)
      
      p_g <- terra::rast(f_prox_g)
      mask_groot_4km <- terra::ifel(!is.na(p_g) & p_g >= 0 & p_g <= 4000, 1, NA)
      unlink(f_prox_g)
      rm(p_g)
      
      netwerk_masker <- terra::cover(mask_groot_4km, mask_klein_4km)
      rm(mask_groot_4km, mask_klein_4km)
    } else {
      netwerk_masker <- mask_klein_4km
      rm(mask_klein_4km)
    }
    
    r_bos_klein_ok <- terra::mask(r_bos_klein, netwerk_masker)
    rm(netwerk_masker)
  } else {
    r_bos_klein_ok <- master_grid * NA
  }
  rm(r_bos_klein)
  gc()
  
  # Samenvoegen van goedgekeurde grote kernen en goedgekeurde satellieten
  g_c <- terra::ifel(is.na(r_bos_groot),    0, r_bos_groot)
  k_c <- terra::ifel(is.na(r_bos_klein_ok), 0, r_bos_klein_ok)
  r_leef_bruto <- terra::ifel((g_c + k_c) > 0, terra::clamp(g_c + k_c, upper = 1.0), NA)
  rm(g_c, k_c, r_bos_groot, r_bos_klein_ok, r_bos_basis)
  gc()
  
  # ----------------------------------------------------------------------------
  # 6. FINALE CLUSTERING OP 50 HA (100M NETWERK)
  # ----------------------------------------------------------------------------
  message("    🧩 Finale clustering (>= 50 ha, 100m overbruggingsafstand)...")
  res_finaal <- cluster_filter_universeel(
    masker            = r_leef_bruto,
    opp_laag          = r_leef_bruto,
    min_opp_ha        = 50,
    min_totale_opp_ha = NA,
    afstand_m         = 100,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  final_opp <- res_finaal$raster
  cl_opp    <- res_finaal$clusters
  rm(res_finaal, r_leef_bruto)
  gc()
  
  return(list(
    final_opp = final_opp,
    cl_opp    = cl_opp
  ))
}