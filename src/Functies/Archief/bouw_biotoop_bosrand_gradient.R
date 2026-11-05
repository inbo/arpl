# ==============================================================================
# # ==============================================================================
# bouw_biotoop_bosrand_gradient.R - LEEFGEBIEDENMODEL VOOR BOOMPIEPER
# ==============================================================================

bouw_biotoop_bosrand_gradient <- function(params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
  bestandsnaam <- "boompieper"
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
  
  # Haal de 4 benodigde lagen op
  r_open_raw   <- if ("bwk"        %in% names(lagen)) lagen[["bwk"]]        else (master_grid * NA)
  r_bos_raw    <- if ("bwk_bossen" %in% names(lagen)) lagen[["bwk_bossen"]] else (master_grid * NA)
  r_heide_raw  <- if ("heide"      %in% names(lagen)) lagen[["heide"]]      else (master_grid * NA)
  r_naald_raw  <- if ("naaldbos"   %in% names(lagen)) lagen[["naaldbos"]]   else (master_grid * NA)
  
  # ----------------------------------------------------------------------------
  # 2. SPOOR A: GROOT LEEFGEBIED (>= 5 HA MET 100 M BOSRAND EN CONTACT)
  # ----------------------------------------------------------------------------
  message("    🌳 Spoor A: Open biotoop en bos onafhankelijk clusteren op 5 ha...")
  
  # A. Open biotoop clusteren (>= 5 ha, 50m zoekafstand)
  res_open <- cluster_filter_universeel(
    masker            = r_open_raw,
    opp_laag          = r_open_raw,
    min_opp_ha        = 5,
    min_totale_opp_ha = NA,
    afstand_m         = 50,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_open_5ha <- res_open$raster
  ha_open_5ha <- calc_ha_exact(r_open_5ha)
  message(sprintf("       -> Open biotoop goedgekeurd (>= 5 ha): %.2f ha", ha_open_5ha))
  rm(res_open)
  gc()
  
  # B. Bos clusteren (>= 5 ha, 50m zoekafstand)
  res_bos <- cluster_filter_universeel(
    masker            = r_bos_raw,
    opp_laag          = r_bos_raw,
    min_opp_ha        = 5,
    min_totale_opp_ha = NA,
    afstand_m         = 50,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_bos_5ha <- res_bos$raster
  ha_bos_5ha <- calc_ha_exact(r_bos_5ha)
  message(sprintf("       -> Boskernen goedgekeurd (>= 5 ha): %.2f ha", ha_bos_5ha))
  rm(res_bos)
  gc()
  
  # C. Buitenste 100m bosrand bepalen via GDAL proximity
  # De rand van het bos bevindt zich waar niet-bos grenst aan bos
  message("    🌲 Buitenste 100m bosrand uitsnijden van grote boscomplexen via GDAL...")
  r_niet_bos_bin <- terra::ifel(is.na(r_bos_5ha), 1, NA)
  
  f_niet_bos_bin <- tempfile(fileext = ".tif")
  terra::writeRaster(r_niet_bos_bin, f_niet_bos_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 0)
  f_prox_bosrand <- hulp_proximity_gdal(f_niet_bos_bin, 100)
  unlink(f_niet_bos_bin)
  rm(r_niet_bos_bin)
  gc()
  
  p_bosrand <- terra::rast(f_prox_bosrand)
  # Cellen binnen het bos die <= 100m van niet-bos liggen
  mask_bos_100m <- terra::ifel(!is.na(p_bosrand) & p_bosrand >= 0 & p_bosrand <= 100, 1, NA)
  unlink(f_prox_bosrand)
  rm(p_bosrand)
  gc()
  
  r_bos_100m <- terra::mask(r_bos_5ha, mask_bos_100m)
  rm(mask_bos_100m, r_bos_5ha)
  gc()
  
  # D. Contactzone: Alleen bosranden behouden die binnen <= 20m van het goedgekeurde open biotoop liggen
  message("    🎯 Bosranden koppelen aan goedgekeurd open biotoop (<= 20m contact)...")
  r_open_bin <- terra::ifel(!is.na(r_open_5ha) & r_open_5ha > 0, 1, NA)
  
  if (!all(is.na(suppressWarnings(terra::minmax(r_open_bin))))) {
    f_open_bin  <- tempfile(fileext = ".tif")
    terra::writeRaster(r_open_bin, f_open_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 0)
    f_prox_open <- hulp_proximity_gdal(f_open_bin, 20)
    unlink(f_open_bin)
    
    p_open20   <- terra::rast(f_prox_open)
    mask_open20 <- terra::ifel(!is.na(p_open20) & p_open20 >= 0 & p_open20 <= 20, 1, NA)
    unlink(f_prox_open)
    rm(p_open20)
    gc()
    
    bosrand_gekoppeld <- terra::mask(r_bos_100m, mask_open20)
    rm(mask_open20)
    gc()
  } else {
    bosrand_gekoppeld <- master_grid * NA
  }
  rm(r_open_bin, r_bos_100m)
  gc()
  
  # Spoor A = Goedgekeurd open biotoop (>= 5 ha) + Gekoppelde bosrand (100m diep)
  spoor_a_opp <- terra::cover(r_open_5ha, bosrand_gekoppeld)
  ha_spoor_a  <- calc_ha_exact(spoor_a_opp)
  message(sprintf("       -> Spoor A totaal (Groot leefgebied): %.2f ha", ha_spoor_a))
  
  # ----------------------------------------------------------------------------
  # 3. SPOOR B: KLEIN LEEFGEBIED (HEIDE-NAALDBOS OVERGANGSZONES >= 2 HA)
  # ----------------------------------------------------------------------------
  message("    🌾 Spoor B: Heide grenzend aan naaldbos (<= 20m) clusteren op 2 ha...")
  r_naald_bin <- terra::ifel(!is.na(r_naald_raw) & r_naald_raw > 0, 1, NA)
  r_heide_bin <- terra::ifel(!is.na(r_heide_raw) & r_heide_raw > 0, 1, NA)
  
  heeft_naald <- !all(is.na(suppressWarnings(terra::minmax(r_naald_bin))))
  heeft_heide <- !all(is.na(suppressWarnings(terra::minmax(r_heide_bin))))
  
  if (heeft_naald && heeft_heide) {
    f_naald_bin  <- tempfile(fileext = ".tif")
    terra::writeRaster(r_naald_bin, f_naald_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 0)
    f_prox_naald <- hulp_proximity_gdal(f_naald_bin, 20)
    unlink(f_naald_bin)
    
    p_naald20    <- terra::rast(f_prox_naald)
    mask_naald20 <- terra::ifel(!is.na(p_naald20) & p_naald20 >= 0 & p_naald20 <= 20, 1, NA)
    unlink(f_prox_naald)
    rm(p_naald20)
    gc()
    
    heide_bij_naald <- terra::mask(r_heide_raw, mask_naald20)
    rm(mask_naald20)
    gc()
    
    # Clusteren op 2 ha (50m zoekafstand)
    res_klein <- cluster_filter_universeel(
      masker            = heide_bij_naald,
      opp_laag          = heide_bij_naald,
      min_opp_ha        = 2,
      min_totale_opp_ha = NA,
      afstand_m         = 50,
      dispersie_m       = 0,
      werkelijk         = TRUE
    )
    spoor_b_opp <- res_klein$raster
    rm(res_klein, heide_bij_naald)
    gc()
  } else {
    spoor_b_opp <- master_grid * NA
  }
  rm(r_naald_bin, r_heide_bin, r_naald_raw, r_heide_raw)
  gc()
  
  ha_spoor_b <- calc_ha_exact(spoor_b_opp)
  message(sprintf("       -> Spoor B totaal (Klein leefgebied): %.2f ha", ha_spoor_b))
  
  # ----------------------------------------------------------------------------
  # 4. SAMENVOEGEN TOT DEFINITIEF LEEFGEBIED & NETWERKEN
  # ----------------------------------------------------------------------------
  message("    🧩 Spoor A en Spoor B combineren (cover)...")
  final_opp <- terra::cover(spoor_a_opp, spoor_b_opp)
  
  # Cluster-ID netwerken via GRASS
  r_fin_bin <- terra::ifel(!is.na(final_opp) & final_opp > 0, 1, NA)
  
  if (!all(is.na(suppressWarnings(terra::minmax(r_fin_bin))))) {
    f_bin   <- tempfile(fileext = ".tif")
    f_clump <- tempfile(fileext = ".tif")
    terra::writeRaster(r_fin_bin, f_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    
    e <- terra::ext(master_grid)
    qgisprocess::qgis_run_algorithm(
      "grass:r.clump",
      input                           = f_bin,
      output                          = f_clump,
      title                           = "cl_boompieper",
      "-d"                            = FALSE,
      GRASS_REGION_PARAMETER          = sprintf("%f,%f,%f,%f", e$xmin, e$xmax, e$ymin, e$ymax),
      GRASS_REGION_CELLSIZE_PARAMETER = 10,
      .quiet                          = TRUE
    )
    
    cl_opp <- terra::rast(f_clump)
    terra::set.cats(cl_opp, layer = 1, value = NULL)
    terra::coltab(cl_opp) <- NULL
    cl_opp <- terra::mask(cl_opp, r_fin_bin)
    
    unlink(c(f_bin, f_clump))
  } else {
    cl_opp <- master_grid * NA
  }
  rm(r_fin_bin)
  gc()
  
  return(list(
    final_opp        = final_opp,
    cl_opp           = cl_opp,
    groot_open       = r_open_5ha,
    groot_bosrand    = bosrand_gekoppeld,
    klein_heide      = spoor_b_opp
  ))
}
