# ==============================================================================
# bouw_biotoop_kleinschalig_struweel_loofbos.R
# Archetype voor struweel- en boszoogdieren (enkelvoudig biotoop met overbrugging)
# Geschikt voor: Eikelmuis, Hazelmuis, etc.
# ==============================================================================

bouw_biotoop_kleinschalig_struweel_loofbos <- function(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
  bestandsnaam <- tolower(gsub(" ", "", soort_naam))
  soort_cache  <- file.path(paths$dir_cache, bestandsnaam)
  dir.create(soort_cache, showWarnings = FALSE, recursive = TRUE)
  
  f_cache_rds <- file.path(soort_cache, "01_raw_parquet_layers.rds")
  
  # Parameters uitlezen
  rij_params    <- params_soort %>% dplyr::slice(1)
  min_opp_ha    <- as.numeric(ifelse(is.na(rij_params$MinOpp_ha), 5, rij_params$MinOpp_ha))
  zoekafstand_m <- as.numeric(ifelse(is.na(rij_params$AfstandBiotopen_m), 20, rij_params$AfstandBiotopen_m))
  
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
  
  r_bio_raw <- if ("bwk" %in% names(lagen)) lagen[["bwk"]] else (master_grid * NA)
  
  ha_bio_ruw <- calc_ha_exact(r_bio_raw)
  message(sprintf("       -> Doelbiotoop vóór clustering: %.2f ha", ha_bio_ruw))
  
  # ----------------------------------------------------------------------------
  # 2. CLUSTERING EN MINIMUMOPPERVLAKTE (>= 5 HA OP 20M ZOEFAFSTAND)
  # ----------------------------------------------------------------------------
  message(sprintf("    🏝️ Clusteren op minimumoppervlakte (>= %.1f ha, zoekafstand %dm)...", min_opp_ha, as.integer(zoekafstand_m)))
  
  res_cluster <- cluster_filter_universeel(
    masker            = r_bio_raw,
    opp_laag          = r_bio_raw,
    min_opp_ha        = min_opp_ha,
    min_totale_opp_ha = NA,
    afstand_m         = zoekafstand_m,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  
  final_opp <- res_cluster$raster
  cl_opp    <- res_cluster$clusters
  rm(res_cluster, r_bio_raw)
  gc()
  
  ha_finaal <- calc_ha_exact(final_opp)
  message(sprintf("       -> Goedgekeurd leefgebied na clustering: %.2f ha", ha_finaal))
  
  return(list(
    final_opp = final_opp,
    cl_opp    = cl_opp
  ))
}