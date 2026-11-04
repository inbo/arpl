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
  
  r_bwk_raw <- if ("bwk" %in% names(lagen)) lagen[["bwk"]] else (master_grid * NA)
  
  # ----------------------------------------------------------------------------
  # 2. FUSIE BWK-BOS EN GROENKAART (STRIKT HOOGGROEN KLASSE 1)
  # ----------------------------------------------------------------------------
  message("    🌲 BWK-bos combineren met Groenkaart (strikt klasse 1 bomen/bos)...")
  
  f_groen_al     <- hulp_align_gdal(paths$groenkaart, master_grid)
  r_groen        <- terra::rast(f_groen_al)
  # Strikt klasse 1 = opgaand bos / bomen in Groenkaart 2021 (geen landbouw 3 of gras 4/5)
  r_bomen_direct <- terra::ifel(!is.na(r_groen) & r_groen == 1, 1.0, NA)
  unlink(f_groen_al)
  rm(r_groen)
  
  # Fusie conform Bert (cover met clamp)
  r_bos_basis <- terra::cover(r_bwk_raw, r_bomen_direct)
  r_bos_basis <- terra::clamp(r_bos_basis, upper = 1.0)
  rm(r_bomen_direct, r_bwk_raw)
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
  # 4. CLASSIFICATIE BOSGROOTTE (>= 300 HA VS. 50-300 HA VIA ZONAL)
  # ----------------------------------------------------------------------------
  message("    📊 Bosclusters classificeren: Grote kernen (>= 300 ha) en satellieten (50-300 ha)...")
  
  if (!all(is.na(suppressWarnings(terra::minmax(cl_patches))))) {
    # Zonale statistiek via terra (geen geheugencrashes)
    df_opp <- terra::zonal(r_bos_basis, cl_patches, fun = "sum", na.rm = TRUE)
    colnames(df_opp) <- c("ID", "Opp_Sum")
    dt_agg <- data.table::as.data.table(df_opp)
    dt_agg[, Opp_ha := Opp_Sum * 0.01]
    rm(df_opp)
    
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
  # 5. RECONSTRUCTIE NETWERK VIA 4 KM INTERACTIE
  # ----------------------------------------------------------------------------
  message("    🎯 4 km netwerkinteractie berekenen voor satellietbossen...")
  
  heeft_groot <- !all(is.na(suppressWarnings(terra::minmax(r_bos_groot))))
  heeft_klein <- !all(is.na(suppressWarnings(terra::minmax(r_bos_klein))))
  
  if (heeft_klein) {
    # 1. Satellieten die binnen 4 km van een groot bos (>= 300 ha) liggen
    if (heeft_groot) {
      f_g_bin  <- tempfile(fileext = ".tif")
      terra::writeRaster(terra::ifel(!is.na(r_bos_groot) & r_bos_groot > 0, 1, NA), f_g_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
      f_prox_g <- hulp_proximity_gdal(f_g_bin, 4000)
      unlink(f_g_bin)
      
      p_g <- terra::rast(f_prox_g)
      mask_groot_4km <- terra::ifel(!is.na(p_g) & p_g >= 0 & p_g <= 4000, 1, NA)
      unlink(f_prox_g)
      rm(p_g)
      
      satelliet_bij_groot <- terra::mask(r_bos_klein, mask_groot_4km)
      rm(mask_groot_4km)
    } else {
      satelliet_bij_groot <- master_grid * NA
    }
    
    # 2. Satellieten die binnen 4 km van een ANDER satellietbos liggen
    # Twee satellieten van elk >= 50 ha binnen 4 km vormen samen een cluster >= 100 ha
    res_sat_cluster <- cluster_filter_universeel(
      masker            = r_bos_klein,
      opp_laag          = r_bos_klein,
      min_opp_ha        = 100,
      min_totale_opp_ha = NA,
      afstand_m         = 4000,
      dispersie_m       = 0,
      werkelijk         = TRUE
    )
    satelliet_bij_satelliet <- res_sat_cluster$raster
    rm(res_sat_cluster)
    
    # Goedgekeurde satellieten = bij groot bos OF bij andere satelliet
    r_bos_klein_ok <- terra::cover(satelliet_bij_groot, satelliet_bij_satelliet)
    rm(satelliet_bij_groot, satelliet_bij_satelliet)
    gc()
  } else {
    r_bos_klein_ok <- master_grid * NA
  }
  rm(r_bos_klein)
  gc()
  
  ha_sat_ok <- calc_ha_exact(r_bos_klein_ok)
  message(sprintf("       -> Satellietbossen goedgekeurd binnen 4 km netwerk: %.2f ha", ha_sat_ok))
  
  # Samenvoegen van grote kernen en goedgekeurde satellieten
  r_leef_bruto <- terra::cover(r_bos_groot, r_bos_klein_ok)
  r_leef_bruto <- terra::clamp(r_leef_bruto, upper = 1.0)
  rm(r_bos_basis)
  gc()
  
  # ----------------------------------------------------------------------------
  # 6. FINALE CLUSTERING OP 50 HA (100M NETWERK) & NETWERK-ID'S
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
  
  ha_finaal <- calc_ha_exact(final_opp)
  message(sprintf("       -> Definitief leefgebied Zwarte specht: %.2f ha", ha_finaal))
  
  # Deellagen isoleren binnen goedgekeurde netwerken
  if (!all(is.na(suppressWarnings(terra::minmax(cl_opp))))) {
    groot_uitvoer <- terra::mask(r_bos_groot, cl_opp)
    klein_uitvoer <- terra::mask(r_bos_klein_ok, cl_opp)
  } else {
    groot_uitvoer <- master_grid * NA
    klein_uitvoer <- master_grid * NA
  }
  rm(r_bos_groot, r_bos_klein_ok)
  gc()
  
  return(list(
    final_opp    = final_opp,
    cl_opp       = cl_opp,
    grote_kernen = groot_uitvoer,
    satellieten  = klein_uitvoer
  ))
}
