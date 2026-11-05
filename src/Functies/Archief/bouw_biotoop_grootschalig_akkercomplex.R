# ==============================================================================
# bouw_biotoop_grootschalig_akkercomplex.R
# Archetype: Grootschalig open akkercomplex (BWK x LGP teelten x >= 3000 ha op 200m)
# Geschikt voor: Grauwe kiekendief, etc.
# ==============================================================================

bouw_biotoop_grootschalig_akkercomplex <- function(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
  bestandsnaam <- tolower(gsub(" ", "", soort_naam))
  soort_cache  <- file.path(paths$dir_cache, bestandsnaam)
  dir.create(soort_cache, showWarnings = FALSE, recursive = TRUE)
  
  f_cache_rds <- file.path(soort_cache, "01_raw_parquet_layers.rds")
  
  # Parameters ophalen
  rij_params    <- params_soort %>% dplyr::slice(1)
  min_opp_ha    <- as.numeric(ifelse(is.na(rij_params$MinOpp_ha), 3000, rij_params$MinOpp_ha))
  zoekafstand_m <- as.numeric(ifelse(is.na(rij_params$AfstandBiotopen_m), 200, rij_params$AfstandBiotopen_m))
  
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
  # 2. FILTERING OP LGP-LANDBOUWGEWASSEN (AND-OPERATIE)
  # ----------------------------------------------------------------------------
  message("    🌾 Doelbiotoop filteren op geschikte landbouwteelten (LGP 2025)...")
  
  if (file.exists(paths$lgp_landbouw)) {
    f_lgp_al <- hulp_align_gdal(paths$lgp_landbouw, master_grid)
    r_lgp    <- terra::rast(f_lgp_al)
    
    # Historische EPR-codes: 6, 8, 38, 45, 60, 138, 174
    gewascodes <- c(6, 8, 38, 45, 60, 138, 174)
    if (!is.na(rij_params$LGP_Gewascodes) && nzchar(trimws(rij_params$LGP_Gewascodes))) {
      gewascodes <- as.numeric(trimws(unlist(strsplit(as.character(rij_params$LGP_Gewascodes), ","))))
    }
    
    mask_teelten <- terra::ifel(!is.na(r_lgp) & r_lgp %in% gewascodes, 1, NA)
    unlink(f_lgp_al)
    rm(r_lgp)
    
    r_akker_geschikt <- terra::mask(r_bwk_raw, mask_teelten)
    rm(mask_teelten)
  } else {
    message("    ⚠️ Waarschuwing: Geen LGP-raster gevonden. BWK-biotoop wordt ongefilterd gebruikt.")
    r_akker_geschikt <- r_bwk_raw
  }
  rm(r_bwk_raw)
  gc()
  
  ha_geschikt <- calc_ha_exact(r_akker_geschikt)
  message(sprintf("       -> Geschikt akkerbiotoop vóór clustering: %.2f ha", ha_geschikt))
  
  # ----------------------------------------------------------------------------
  # 3. LANDSCHAPPELIJKE CLUSTERING (>= 3000 HA OP 200M ZOEFAFSTAND)
  # ----------------------------------------------------------------------------
  message(sprintf("    🏝️ Grootschalig akkercomplex clusteren (>= %.1f ha, zoekafstand %dm)...", min_opp_ha, as.integer(zoekafstand_m)))
  
  res_cluster <- cluster_filter_universeel(
    masker            = r_akker_geschikt,
    opp_laag          = r_akker_geschikt,
    min_opp_ha        = min_opp_ha,
    min_totale_opp_ha = NA,
    afstand_m         = zoekafstand_m,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  
  final_opp <- res_cluster$raster
  cl_opp    <- res_cluster$clusters
  rm(res_cluster, r_akker_geschikt)
  gc()
  
  return(list(
    final_opp = final_opp,
    cl_opp    = cl_opp
  ))
}