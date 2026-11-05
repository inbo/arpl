# ==============================================================================
# bouw_biotoop_drasland_steltloper.R
# Archetype voor plas-dras steltlopers (drainage x matrixkwaliteit x terrestrische kern)
# Geschikt voor: Watersnip, etc.
# ==============================================================================

bouw_biotoop_drasland_steltloper <- function(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
  bestandsnaam <- tolower(gsub(" ", "", soort_naam))
  soort_cache  <- file.path(paths$dir_cache, bestandsnaam)
  dir.create(soort_cache, showWarnings = FALSE, recursive = TRUE)
  
  f_cache_rds <- file.path(soort_cache, "01_raw_parquet_layers.rds")
  
  # Parameters ophalen
  rij_params    <- params_soort %>% dplyr::slice(1)
  min_opp_ha    <- as.numeric(ifelse(is.na(rij_params$MinOpp_ha), 10, rij_params$MinOpp_ha))
  zoekafstand_m <- as.numeric(ifelse(is.na(rij_params$AfstandBiotopen_m), 100, rij_params$AfstandBiotopen_m))
  
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
  
  r_bio_raw <- if ("biotopen" %in% names(lagen)) lagen[["biotopen"]] else (master_grid * NA)
  
  ha_bio_ruw <- calc_ha_exact(r_bio_raw)
  message(sprintf("       -> Doelbiotoop vóór drainagefilter: %.2f ha", ha_bio_ruw))
  
  # ----------------------------------------------------------------------------
  # 2. BODEMDRAINAGEFILTER OP BASISOPPERVLAKTE (d, e, f, h, i, g, e-f, h-i, e-i)
  # ----------------------------------------------------------------------------
  message("    💧 Bodemdrainage filteren (natte klassen)...")
  drain_let <- c("d", "e", "f", "h", "i", "g", "e-f", "h-i", "e-i")
  
  if (!is.null(paths$drainage) && file.exists(paths$drainage)) {
    r_drain_orig <- terra::rast(paths$drainage)
    leg_drain    <- terra::cats(r_drain_orig)[[1]]
    
    ok_drain_ids <- integer(0)
    if (!is.null(leg_drain)) {
      kol_val   <- colnames(leg_drain)[tolower(colnames(leg_drain)) %in% c("value", "id")][1]
      kol_label <- colnames(leg_drain)[tolower(colnames(leg_drain)) %in% c("label", "omschrijving", "drainage")][1]
      
      if (!is.na(kol_val) && !is.na(kol_label)) {
        ok_drain_ids <- leg_drain[[kol_val]][tolower(trimws(as.character(leg_drain[[kol_label]]))) %in% tolower(drain_let)]
      }
    }
    if (length(ok_drain_ids) == 0) ok_drain_ids <- c(4, 5, 6, 7, 8, 9, 45, 56, 89)
    rm(r_drain_orig)
    
    f_drain_al <- hulp_align_gdal(paths$drainage, master_grid)
    r_drain    <- terra::rast(f_drain_al)
    
    mask_drain <- terra::ifel(!is.na(r_drain) & (r_drain %in% ok_drain_ids), 1, NA)
    r_bio_nat  <- terra::mask(r_bio_raw, mask_drain)
    unlink(f_drain_al)
    rm(r_drain, mask_drain)
  } else {
    r_bio_nat <- r_bio_raw
  }
  rm(r_bio_raw)
  gc()
  
  ha_nat <- calc_ha_exact(r_bio_nat)
  message(sprintf("       -> Biotoop na drainagefilter: %.2f ha", ha_nat))
  
  # ----------------------------------------------------------------------------
  # 3. PRIMAIRE LANDSCHAPSCLUSTERING (>= 10 HA OP 100M NETWERK)
  # ----------------------------------------------------------------------------
  message(sprintf("    🏝️ Primaire clustering op landschapsniveau (>= %.1f ha, %dm overbrugging)...", min_opp_ha, as.integer(zoekafstand_m)))
  
  res_10ha <- cluster_filter_universeel(
    masker            = r_bio_nat,
    opp_laag          = r_bio_nat,
    min_opp_ha        = min_opp_ha,
    min_totale_opp_ha = NA,
    afstand_m         = zoekafstand_m,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_10ha_netwerk <- res_10ha$raster
  cl_10ha        <- res_10ha$clusters
  rm(res_10ha)
  gc()
  
  ha_10ha <- calc_ha_exact(r_10ha_netwerk)
  message(sprintf("       -> Biotoop in clusters >= %.1f ha: %.2f ha", min_opp_ha, ha_10ha))
  
  # ----------------------------------------------------------------------------
  # 4. LANDSCHAPSKWALITEITSTOETS (< 75% BOS & < 90% WATER PER CLUSTER)
  # ----------------------------------------------------------------------------
  message("    🌳 Landschapskwaliteit toetsen (< 75% bos en < 90% open water per cluster via zonale streaming)...")
  
  # Waterlaag klaarmaken voor kwaliteitscheck en terrestrische kern
  if (!is.null(paths$grb_water) && file.exists(paths$grb_water)) {
    f_grb_al  <- hulp_align_gdal(paths$grb_water, master_grid)
    r_grb     <- terra::rast(f_grb_al)
    r_wat_bin <- terra::ifel(!is.na(r_grb) & r_grb > 0, 1.0, NA)
    unlink(f_grb_al)
    rm(r_grb)
  } else {
    r_wat_bin <- master_grid * NA
  }
  
  ids_kwaliteit_ok <- integer(0)
  
  if (!all(is.na(suppressWarnings(terra::minmax(cl_10ha))))) {
    # A. Totale werkelijke oppervlakte per cluster
    df_tot <- terra::zonal(r_bio_nat, cl_10ha, fun = "sum", na.rm = TRUE)
    colnames(df_tot) <- c("Cluster_ID", "Biotoop_Som")
    dt_tot <- data.table::as.data.table(df_tot)
    dt_tot[, Totaal_ha := Biotoop_Som * 0.01]
    rm(df_tot)
    
    # B. Aandeel Hooggroen (strikt klasse 1 = bomen/bos in Groenkaart_2021)
    f_groen_al <- hulp_align_gdal(paths$groenkaart, master_grid)
    r_groen    <- terra::rast(f_groen_al)
    r_bos_bin  <- terra::ifel(!is.na(r_groen) & r_groen == 1, 1.0, NA)
    unlink(f_groen_al)
    rm(r_groen)
    
    bos_in_zone <- terra::mask(r_bos_bin, cl_10ha)
    rm(r_bos_bin)
    df_bos <- terra::zonal(bos_in_zone, cl_10ha, fun = "sum", na.rm = TRUE)
    colnames(df_bos) <- c("Cluster_ID", "Bos_Som")
    dt_bos <- data.table::as.data.table(df_bos)
    dt_bos[, Bos_ha := Bos_Som * 0.01]
    rm(bos_in_zone, df_bos)
    
    # C. Aandeel Open Water
    wat_in_zone <- terra::mask(r_wat_bin, cl_10ha)
    df_wat <- terra::zonal(wat_in_zone, cl_10ha, fun = "sum", na.rm = TRUE)
    colnames(df_wat) <- c("Cluster_ID", "Water_Som")
    dt_wat <- data.table::as.data.table(df_wat)
    dt_wat[, Water_ha := Water_Som * 0.01]
    rm(wat_in_zone, df_wat)
    
    # Samenvoegen en toetsen
    dt_merged <- merge(dt_tot, dt_bos, by = "Cluster_ID", all.x = TRUE)
    dt_merged <- merge(dt_merged, dt_wat, by = "Cluster_ID", all.x = TRUE)
    dt_merged[is.na(Bos_ha), Bos_ha := 0]
    dt_merged[is.na(Water_ha), Water_ha := 0]
    
    dt_merged[, Pct_Bos   := (Bos_ha / Totaal_ha) * 100]
    dt_merged[, Pct_Water := (Water_ha / Totaal_ha) * 100]
    
    ids_kwaliteit_ok <- dt_merged[!is.na(Cluster_ID) & Pct_Bos < 75 & Pct_Water < 90, Cluster_ID]
    message(sprintf("       -> Clusters goedgekeurd op openheid (<75%% bos, <90%% water): %d van de %d clusters", 
                    length(ids_kwaliteit_ok), nrow(dt_merged)))
    
    rm(dt_tot, dt_bos, dt_wat, dt_merged)
    gc()
  }
  
  # ----------------------------------------------------------------------------
  # 5. TERRESTRISCHE KERN CHECK (>= 10 HA LAND OP 0M CLUSTERING)
  # ----------------------------------------------------------------------------
  message("    🌾 Terrestrische kerncheck: Minstens 10 ha zuiver landbiotoop (0m clustering)...")
  
  if (length(ids_kwaliteit_ok) > 0) {
    # Sluit open water uit van het biotoop voor de zuiver terrestrische kern
    r_terr_biotoop <- terra::mask(r_bio_nat, r_wat_bin, inverse = TRUE)
    rm(r_wat_bin)
    gc()
    
    # 0m clustering (aaneengesloten terrestrische kern >= 10 ha)
    res_terr <- cluster_filter_universeel(
      masker            = r_terr_biotoop,
      opp_laag          = r_terr_biotoop,
      min_opp_ha        = 10,
      min_totale_opp_ha = NA,
      afstand_m         = 0,
      dispersie_m       = 0,
      werkelijk         = TRUE
    )
    cl_terr_10ha <- res_terr$clusters
    rm(res_terr, r_terr_biotoop)
    gc()
    
    # Controleer overlap tussen goedgekeurde netwerken en de 10 ha landkernen via crosstab
    if (!all(is.na(suppressWarnings(terra::minmax(cl_terr_10ha))))) {
      ct <- terra::crosstab(c(cl_10ha, cl_terr_10ha), long = TRUE)
      colnames(ct) <- c("Cluster_ID", "Terr_ID", "n")
      dt_ct <- data.table::as.data.table(ct)[!is.na(Cluster_ID) & !is.na(Terr_ID) & n > 0]
      rm(ct, cl_terr_10ha)
      
      valide_cluster_ids <- unique(dt_ct[Cluster_ID %in% ids_kwaliteit_ok, Cluster_ID])
      rm(dt_ct)
    } else {
      valide_cluster_ids <- integer(0)
      rm(cl_terr_10ha)
    }
    
    if (length(valide_cluster_ids) > 0) {
      cl_valide_mask <- terra::ifel(cl_10ha %in% as.integer(valide_cluster_ids), 1, NA)
      final_opp      <- terra::mask(r_10ha_netwerk, cl_valide_mask)
      cl_opp         <- terra::mask(cl_10ha, cl_valide_mask)
      rm(cl_valide_mask)
    } else {
      final_opp <- master_grid * NA
      cl_opp    <- master_grid * NA
    }
  } else {
    rm(r_wat_bin)
    final_opp <- master_grid * NA
    cl_opp    <- master_grid * NA
  }
  rm(r_bio_nat, r_10ha_netwerk, cl_10ha)
  gc()
  
  ha_finaal <- calc_ha_exact(final_opp)
  message(sprintf("       -> Definitief leefgebied Watersnip: %.2f ha", ha_finaal))
  
  return(list(
    final_opp = final_opp,
    cl_opp    = cl_opp
  ))
}
