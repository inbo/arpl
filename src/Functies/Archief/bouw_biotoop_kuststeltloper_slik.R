# ==============================================================================
# bouw_biotoop_kuststeltloper_slik.R
# Archetype voor kust-/sliksteltlopers (broedlocatie >= 5ha x slikfoerageer >= 5ha binnen 2km)
# Geschikt voor: Kluut, etc.
# ==============================================================================

bouw_biotoop_kuststeltloper_slik <- function(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
  bestandsnaam <- tolower(gsub(" ", "", soort_naam))
  soort_cache  <- file.path(paths$dir_cache, bestandsnaam)
  dir.create(soort_cache, showWarnings = FALSE, recursive = TRUE)
  
  f_cache_rds <- file.path(soort_cache, "01_raw_parquet_layers.rds")
  
  # Parameters uitlezen
  rij_vp   <- params_soort %>% dplyr::filter(tolower(Type) == "voortplanting") %>% dplyr::slice(1)
  rij_foer <- params_soort %>% dplyr::filter(tolower(Type) == "foerageer") %>% dplyr::slice(1)
  
  min_ha_vp   <- as.numeric(ifelse(nrow(rij_vp) > 0 && !is.na(rij_vp$MinOpp_ha), rij_vp$MinOpp_ha, 5))
  dist_vp_m   <- as.numeric(ifelse(nrow(rij_vp) > 0 && !is.na(rij_vp$AfstandBiotopen_m), rij_vp$AfstandBiotopen_m, 50))
  
  min_ha_foer <- as.numeric(ifelse(nrow(rij_foer) > 0 && !is.na(rij_foer$MinOpp_ha), rij_foer$MinOpp_ha, 5))
  dist_foer_m <- as.numeric(ifelse(nrow(rij_foer) > 0 && !is.na(rij_foer$AfstandBiotopen_m), rij_foer$AfstandBiotopen_m, 10))
  
  koppel_dist_m <- as.numeric(ifelse(nrow(rij_vp) > 0 && !is.na(rij_vp$Wederzijdse_Afstand_m), rij_vp$Wederzijdse_Afstand_m, 2000))
  
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
  
  r_vp_raw    <- if ("voortplanting" %in% names(lagen)) lagen[["voortplanting"]] else (master_grid * NA)
  r_foer_raw  <- if ("foerageer"     %in% names(lagen)) lagen[["foerageer"]]     else (master_grid * NA)
  r_water_raw <- if ("openwater"     %in% names(lagen)) lagen[["openwater"]]     else (master_grid * NA)
  
  # ----------------------------------------------------------------------------
  # 2. BASISLAGEN FILTEREN EN CLUSTEREN
  # ----------------------------------------------------------------------------
  message(sprintf("    🏝️ Broedbiotoop zuiveren en clusteren (>= %.1f ha op %dm)...", min_ha_vp, as.integer(dist_vp_m)))
  
  # Broedbiotoop zuiveren van direct diep open water / moeras indien openwaterlaag aanwezig is
  if (!all(is.na(suppressWarnings(terra::minmax(r_water_raw))))) {
    r_vp_gezuiverd <- terra::mask(r_vp_raw, r_water_raw > 0, inverse = TRUE)
  } else {
    r_vp_gezuiverd <- r_vp_raw
  }
  rm(r_vp_raw, r_water_raw)
  gc()
  
  res_vp <- cluster_filter_universeel(
    masker            = r_vp_gezuiverd,
    opp_laag          = r_vp_gezuiverd,
    min_opp_ha        = min_ha_vp,
    min_totale_opp_ha = NA,
    afstand_m         = dist_vp_m,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_vp_ok <- res_vp$raster
  rm(res_vp, r_vp_gezuiverd)
  gc()
  
  message(sprintf("    🦀 Foerageerslikken clusteren (>= %.1f ha op %dm)...", min_ha_foer, as.integer(dist_foer_m)))
  res_foer <- cluster_filter_universeel(
    masker            = r_foer_raw,
    opp_laag          = r_foer_raw,
    min_opp_ha        = min_ha_foer,
    min_totale_opp_ha = NA,
    afstand_m         = dist_foer_m,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_foer_ok <- res_foer$raster
  rm(res_foer, r_foer_raw)
  gc()
  
  ha_vp_ok   <- calc_ha_exact(r_vp_ok)
  ha_foer_ok <- calc_ha_exact(r_foer_ok)
  message(sprintf("       -> Broedgebieden (>= %.1f ha): %.2f ha", min_ha_vp, ha_vp_ok))
  message(sprintf("       -> Foerageerslikken (>= %.1f ha): %.2f ha", min_ha_foer, ha_foer_ok))
  
  # ----------------------------------------------------------------------------
  # 3. WEDERZIJDSE ACTIERADIUSKOPPELING (GDAL PROXIMITY)
  # ----------------------------------------------------------------------------
  message(sprintf("    ↔️ Wederzijdse %dm actieradiuskoppeling tussen nest- en voedselgebieden...", as.integer(koppel_dist_m)))
  
  heeft_vp   <- !all(is.na(suppressWarnings(terra::minmax(r_vp_ok))))
  heeft_foer <- !all(is.na(suppressWarnings(terra::minmax(r_foer_ok))))
  
  if (heeft_vp && heeft_foer) {
    r_vp_bin   <- terra::ifel(!is.na(r_vp_ok) & r_vp_ok > 0, 1, NA)
    r_foer_bin <- terra::ifel(!is.na(r_foer_ok) & r_foer_ok > 0, 1, NA)
    
    # Proximity vanaf foerageergebied
    f_foer_bin  <- tempfile(fileext = ".tif")
    terra::writeRaster(r_foer_bin, f_foer_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    f_prox_foer <- hulp_proximity_gdal(f_foer_bin, koppel_dist_m)
    unlink(f_foer_bin)
    
    p_foer      <- terra::rast(f_prox_foer)
    mask_foer2k <- terra::ifel(!is.na(p_foer) & p_foer >= 0 & p_foer <= koppel_dist_m, 1, NA)
    unlink(f_prox_foer)
    rm(p_foer)
    
    # Proximity vanaf broedgebied
    f_vp_bin  <- tempfile(fileext = ".tif")
    terra::writeRaster(r_vp_bin, f_vp_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    f_prox_vp <- hulp_proximity_gdal(f_vp_bin, koppel_dist_m)
    unlink(f_vp_bin)
    
    p_vp      <- terra::rast(f_prox_vp)
    mask_vp2k <- terra::ifel(!is.na(p_vp) & p_vp >= 0 & p_vp <= koppel_dist_m, 1, NA)
    unlink(f_prox_vp)
    rm(p_vp)
    
    # Wederzijds maskeren
    vp_gekoppeld   <- terra::mask(r_vp_ok, mask_foer2k)
    foer_gekoppeld <- terra::mask(r_foer_ok, mask_vp2k)
    
    unlink(c(mask_foer2k, mask_vp2k))
    rm(r_vp_bin, r_foer_bin)
    gc()
  } else {
    vp_gekoppeld   <- master_grid * NA
    foer_gekoppeld <- master_grid * NA
  }
  rm(r_vp_ok, r_foer_ok)
  gc()
  
  # ----------------------------------------------------------------------------
  # 4. SAMENVOEGEN COMPOSIET LEEFGEBIED & OVERKOEPELENDE CLUSTERING
  # ----------------------------------------------------------------------------
  min_ha_finaal <- max(min_ha_vp, min_ha_foer)
  message(sprintf("    🧩 Leefgebied samenvoegen en finaal clusteren (>= %.1f ha, 50m netwerk)...", min_ha_finaal))
  
  w_v <- terra::ifel(is.na(vp_gekoppeld), 0, vp_gekoppeld)
  w_f <- terra::ifel(is.na(foer_gekoppeld), 0, foer_gekoppeld)
  r_leef_bruto <- terra::ifel((w_v + w_f) > 0, terra::clamp(w_v + w_f, upper = 1.0), NA)
  rm(w_v, w_f)
  gc()
  
  res_finaal <- cluster_filter_universeel(
    masker            = r_leef_bruto,
    opp_laag          = r_leef_bruto,
    min_opp_ha        = min_ha_finaal,
    min_totale_opp_ha = NA,
    afstand_m         = 50,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  final_opp <- res_finaal$raster
  cl_opp    <- res_finaal$clusters
  rm(res_finaal, r_leef_bruto)
  gc()
  
  vp_uitvoer   <- terra::mask(vp_gekoppeld, cl_opp)
  foer_uitvoer <- terra::mask(foer_gekoppeld, cl_opp)
  rm(vp_gekoppeld, foer_gekoppeld)
  gc()
  
  return(list(
    final_opp     = final_opp,
    cl_opp        = cl_opp,
    voortplanting = vp_uitvoer,
    foerageren    = foer_uitvoer
  ))
}
