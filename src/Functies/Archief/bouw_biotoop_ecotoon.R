# ==============================================================================
# bouw_biotoop_ecotoon.R - GENERIEKE ECOTONALE ANALYSE (DAGVLINDERS & BOSRANDEN)
# ==============================================================================

# Helper: Generieke ecotonale bosrandtoets
evalueer_bosrand_contact_universeel <- function(r_hab_opp, r_bos_opp, drempel_ha = 5, min_raaklijn_m = 50, master_grid) {
  
  if (all(is.na(suppressWarnings(terra::minmax(r_hab_opp))))) {
    return(list(opp = master_grid * NA, clusters = master_grid * NA, hab_deel = master_grid * NA, bos_zoom = master_grid * NA))
  }
  if (all(is.na(suppressWarnings(terra::minmax(r_bos_opp))))) {
    return(list(opp = master_grid * NA, clusters = master_grid * NA, hab_deel = master_grid * NA, bos_zoom = master_grid * NA))
  }
  
  r_hab_bin <- terra::ifel(!is.na(r_hab_opp) & r_hab_opp > 0, 1, NA)
  if (all(is.na(suppressWarnings(terra::minmax(r_hab_bin))))) {
    return(list(opp = master_grid * NA, clusters = master_grid * NA, hab_deel = master_grid * NA, bos_zoom = master_grid * NA))
  }
  
  e <- terra::ext(master_grid)
  ext_str <- sprintf("%f,%f,%f,%f", e$xmin, e$xmax, e$ymin, e$ymax)
  
  # A. Clump habitat
  f_hab_bin <- tempfile(fileext = ".tif")
  f_cl_hab  <- tempfile(fileext = ".tif")
  terra::writeRaster(r_hab_bin, f_hab_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
  
  qgisprocess::qgis_run_algorithm(
    "grass:r.clump", input = f_hab_bin, output = f_cl_hab, title = "cl_hab",
    "-d" = FALSE, GRASS_REGION_PARAMETER = ext_str, GRASS_REGION_CELLSIZE_PARAMETER = 10, .quiet = TRUE
  )
  unlink(f_hab_bin)
  
  if (!file.exists(f_cl_hab)) {
    return(list(opp = master_grid * NA, clusters = master_grid * NA, hab_deel = master_grid * NA, bos_zoom = master_grid * NA))
  }
  
  cl_hab <- terra::rast(f_cl_hab)
  terra::set.cats(cl_hab, layer = 1, value = NULL)
  if (terra::has.colors(cl_hab)) try(terra::coltab(cl_hab) <- NULL, silent = TRUE)
  
  cells_h <- terra::cells(r_hab_bin)
  dt_h    <- data.table::data.table(ID = as.integer(cl_hab[cells_h][[1]]), Waarde = as.numeric(r_hab_opp[cells_h][[1]]))
  ok_ids  <- dt_h[!is.na(ID) & ID > 0, .(ha = sum(Waarde, na.rm = TRUE) * 0.01), by = ID][ha >= drempel_ha, ID]
  rm(cells_h, dt_h)
  gc()
  
  if (length(ok_ids) == 0) {
    unlink(f_cl_hab)
    return(list(opp = master_grid * NA, clusters = master_grid * NA, hab_deel = master_grid * NA, bos_zoom = master_grid * NA))
  }
  
  f_hab_groot <- tempfile(fileext = ".tif")
  r_hab_groot_bin <- terra::ifel(cl_hab %in% as.integer(ok_ids), 1, NA)
  terra::writeRaster(r_hab_groot_bin, f_hab_groot, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
  rm(r_hab_groot_bin)
  gc()
  
  # B. Zoek contactlijnen met bos (binnen 10m via GDAL)
  f_prox_10m  <- hulp_proximity_gdal(f_hab_groot, 10)
  unlink(f_hab_groot)
  
  p10       <- terra::rast(f_prox_10m)
  r_bos_bin <- terra::ifel(!is.na(r_bos_opp) & r_bos_opp > 0, 1, NA)
  r_contact <- terra::ifel(!is.na(p10) & p10 >= 0 & p10 <= 10 & !is.na(r_bos_bin), 1, NA)
  
  unlink(f_prox_10m)
  rm(p10, r_bos_bin)
  gc()
  
  if (all(is.na(suppressWarnings(terra::minmax(r_contact))))) {
    unlink(f_cl_hab)
    return(list(opp = master_grid * NA, clusters = master_grid * NA, hab_deel = master_grid * NA, bos_zoom = master_grid * NA))
  }
  
  # C. Clump contactpixels
  f_contact_bin <- tempfile(fileext = ".tif")
  f_cl_contact  <- tempfile(fileext = ".tif")
  terra::writeRaster(r_contact, f_contact_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
  rm(r_contact)
  gc()
  
  qgisprocess::qgis_run_algorithm(
    "grass:r.clump", input = f_contact_bin, output = f_cl_contact, title = "cl_contact",
    "-d" = FALSE, GRASS_REGION_PARAMETER = ext_str, GRASS_REGION_CELLSIZE_PARAMETER = 10, .quiet = TRUE
  )
  unlink(f_contact_bin)
  
  if (!file.exists(f_cl_contact)) {
    unlink(f_cl_hab)
    return(list(opp = master_grid * NA, clusters = master_grid * NA, hab_deel = master_grid * NA, bos_zoom = master_grid * NA))
  }
  
  cl_contact <- terra::rast(f_cl_contact)
  terra::set.cats(cl_contact, layer = 1, value = NULL)
  if (terra::has.colors(cl_contact)) try(terra::coltab(cl_contact) <- NULL, silent = TRUE)
  
  dt_contact <- data.table::as.data.table(terra::freq(cl_contact))
  setnames(dt_contact, c("layer", "ID", "count"))
  
  min_cellen_lijn <- ceiling(min_raaklijn_m / 10)
  ok_contact_ids  <- dt_contact[!is.na(ID) & ID > 0 & count >= min_cellen_lijn, ID]
  rm(dt_contact)
  gc()
  
  if (length(ok_contact_ids) == 0) {
    unlink(c(f_cl_hab, f_cl_contact))
    return(list(opp = master_grid * NA, clusters = master_grid * NA, hab_deel = master_grid * NA, bos_zoom = master_grid * NA))
  }
  
  # D. Koppel goedgekeurde contactlijnen terug via GDAL proximity 15m
  f_raaklijn    <- tempfile(fileext = ".tif")
  r_raaklijn_ok <- terra::ifel(cl_contact %in% as.integer(ok_contact_ids), 1, NA)
  terra::writeRaster(r_raaklijn_ok, f_raaklijn, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
  unlink(f_cl_contact)
  rm(r_raaklijn_ok, cl_contact)
  gc()
  
  f_prox_15m <- hulp_proximity_gdal(f_raaklijn, 15)
  unlink(f_raaklijn)
  
  p15         <- terra::rast(f_prox_15m)
  r_touch_bin <- terra::ifel(!is.na(p15) & p15 >= 0 & p15 <= 15, 1, NA)
  cells_touch <- terra::cells(r_touch_bin)
  unlink(f_prox_15m)
  rm(p15, r_touch_bin)
  gc()
  
  if (length(cells_touch) == 0) {
    unlink(f_cl_hab)
    return(list(opp = master_grid * NA, clusters = master_grid * NA, hab_deel = master_grid * NA, bos_zoom = master_grid * NA))
  }
  
  voldoet_bos_ids <- unique(na.omit(as.integer(cl_hab[cells_touch][[1]])))
  rm(cells_touch)
  gc()
  
  if (length(voldoet_bos_ids) == 0) {
    unlink(f_cl_hab)
    return(list(opp = master_grid * NA, clusters = master_grid * NA, hab_deel = master_grid * NA, bos_zoom = master_grid * NA))
  }
  
  # E. Goedgekeurd habitat en 20m boszoom
  r_hab_goedgekeurd_opp <- terra::mask(r_hab_opp, cl_hab %in% as.integer(voldoet_bos_ids))
  r_hab_goedgekeurd_bin <- terra::ifel(!is.na(r_hab_goedgekeurd_opp), 1, NA)
  
  f_hab_finaal <- tempfile(fileext = ".tif")
  terra::writeRaster(r_hab_goedgekeurd_bin, f_hab_finaal, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
  rm(r_hab_goedgekeurd_bin)
  gc()
  
  f_prox_20m   <- hulp_proximity_gdal(f_hab_finaal, 20)
  unlink(f_hab_finaal)
  
  p20          <- terra::rast(f_prox_20m)
  bos_zoom_opp <- terra::mask(r_bos_opp, terra::ifel(!is.na(p20) & p20 >= 0 & p20 <= 20, 1, NA))
  unlink(f_prox_20m)
  rm(p20)
  gc()
  
  w_hab  <- terra::ifel(is.na(r_hab_goedgekeurd_opp), 0, r_hab_goedgekeurd_opp)
  w_bos  <- terra::ifel(is.na(bos_zoom_opp), 0, bos_zoom_opp)
  r_comb <- terra::ifel((w_hab + w_bos) > 0, terra::clamp(w_hab + w_bos, upper = 1.0), NA)
  
  cl_finaal <- terra::mask(cl_hab, r_hab_goedgekeurd_opp)
  unlink(f_cl_hab)
  
  return(list(
    opp      = r_comb,
    clusters = cl_finaal,
    hab_deel = r_hab_goedgekeurd_opp,
    bos_zoom = bos_zoom_opp
  ))
}

# ==============================================================================
# HOOFDFUNCTIE VOOR ALLE ECOTOON- EN MANTELVEGETATIESOORTEN
# ==============================================================================
bouw_biotoop_ecotoon <- function(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
  bestandsnaam <- tolower(gsub(" ", "", soort_naam))
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
    lagen_packed <- lapply(lagen, terra::wrap)
    saveRDS(lagen_packed, f_cache_rds)
    message("    💾 [CACHE] Parquet-lagen veilig bewaard.")
  }
  
  # 2. SAMENSTELLEN SPOOR A (Voortplanting / Eiafzet)
  spoor_a_types <- params_soort %>% 
    dplyr::filter(Spoor == "A" & !tolower(Type) %in% c("bosranden", "water_bwk")) %>% 
    dplyr::pull(Type)
  
  if (length(spoor_a_types) == 1) {
    r_spoor_a <- lagen[[spoor_a_types[1]]]
  } else if (length(spoor_a_types) > 1) {
    message("    ➕ Meerdere Spoor A types sommatief samenvoegen...")
    r_spoor_a <- lagen[[spoor_a_types[1]]]
    for (tp in spoor_a_types[-1]) {
      r1 <- terra::ifel(is.na(r_spoor_a), 0, r_spoor_a)
      r2 <- terra::ifel(is.na(lagen[[tp]]), 0, lagen[[tp]])
      r_som <- r1 + r2
      r_spoor_a <- terra::ifel(r_som > 0, terra::clamp(r_som, upper = 1.0), NA)
    }
  } else {
    r_spoor_a <- master_grid * NA
  }
  
  # 3. SAMENSTELLEN SPOOR B (Foerageren / Nectar)
  spoor_b_types <- params_soort %>% dplyr::filter(Spoor == "B") %>% dplyr::pull(Type)
  if (length(spoor_b_types) >= 1) {
    r_spoor_b <- lagen[[spoor_b_types[1]]]
    if (length(spoor_b_types) > 1) {
      for (tp in spoor_b_types[-1]) {
        r_spoor_b <- terra::cover(r_spoor_b, lagen[[tp]])
      }
    }
  } else {
    r_spoor_b <- NULL
  }
  
  # 4. OPTIONELE ABIOTISCHE FILTERS VOOR SPOOR A (Drainage en/of Water)
  rij_a <- params_soort %>% dplyr::filter(Spoor == "A" & !is.na(MinOpp_ha)) %>% dplyr::slice(1)
  
  # A. Waterfilter
  if (isTRUE(as.logical(rij_a$Gebruik_GRB_Water))) {
    max_w_dist <- as.numeric(rij_a$Water_Afstand_m)
    message(sprintf("    💧 Waterlaag uitlijnen en %dm nabijheid berekenen...", max_w_dist))
    
    f_grb_aligned <- hulp_align_gdal(paths$grb_water, master_grid)
    r_grb         <- terra::rast(f_grb_aligned)
    r_grb_bin     <- terra::ifel(!is.na(r_grb) & r_grb == 1, 1, NA)
    r_water_bwk   <- if ("water_bwk" %in% names(lagen)) lagen[["water_bwk"]] else (master_grid * NA)
    r_water_tot   <- terra::cover(r_grb_bin, r_water_bwk)
    
    f_water_bin     <- tempfile(fileext = ".tif")
    terra::writeRaster(terra::ifel(!is.na(r_water_tot), 1, NA), f_water_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    f_water_prox    <- hulp_proximity_gdal(f_water_bin, max_w_dist)
    p_water         <- terra::rast(f_water_prox)
    mask_water      <- terra::ifel(!is.na(p_water) & p_water >= 0 & p_water <= max_w_dist, 1, NA)
    
    unlink(c(f_grb_aligned, f_water_bin, f_water_prox))
    rm(r_grb, r_grb_bin, r_water_tot, p_water)
    gc()
    
    r_spoor_a <- terra::mask(r_spoor_a, mask_water)
    rm(mask_water)
    gc()
  }
  
  # B. Bodemdrainagefilter
  if (isTRUE(as.logical(rij_a$Gebruik_Drainage)) && !is.na(rij_a$Drainage_Codes)) {
    message("    🌱 Bodemdrainage filter toepassen...")
    codes_str <- unlist(strsplit(as.character(rij_a$Drainage_Codes), ","))
    codes_str <- trimws(codes_str)
    
    r_drain_raw   <- terra::rast(paths$drainage)[[1]]
    drain_cats    <- terra::cats(r_drain_raw)[[1]]
    drain_ids_ok  <- drain_cats$value[as.character(drain_cats$Label) %in% codes_str]
    rm(r_drain_raw, drain_cats)
    
    f_drain_aligned <- hulp_align_gdal(paths$drainage, master_grid)
    r_drain         <- terra::rast(f_drain_aligned)
    mask_drain      <- terra::ifel(!is.na(r_drain) & r_drain %in% drain_ids_ok, 1, NA)
    unlink(f_drain_aligned)
    rm(r_drain)
    gc()
    
    r_spoor_a <- terra::mask(r_spoor_a, mask_drain)
    rm(mask_drain)
    gc()
  }
  
  # 5. ECOTONALE BOSRANDTOETS
  min_opp_a  <- as.numeric(rij_a$MinOpp_ha)
  raaklijn_m <- as.numeric(ifelse(is.na(rij_a$AfstandBiotopen_m), 50, rij_a$AfstandBiotopen_m))
  r_bos      <- if ("bosranden" %in% names(lagen)) lagen[["bosranden"]] else (master_grid * NA)
  
  message(sprintf("    🦋 Spoor A toetsen aan %dm bosrand-eis (>= %.1f ha)...", raaklijn_m, min_opp_a))
  res_a <- evalueer_bosrand_contact_universeel(
    r_hab_opp      = r_spoor_a,
    r_bos_opp      = r_bos,
    drempel_ha     = min_opp_a,
    min_raaklijn_m = raaklijn_m,
    master_grid    = master_grid
  )
  
  res_b <- NULL
  if (!is.null(r_spoor_b)) {
    rij_b     <- params_soort %>% dplyr::filter(Spoor == "B") %>% dplyr::slice(1)
    min_opp_b <- as.numeric(ifelse(is.na(rij_b$MinOpp_ha), min_opp_a, rij_b$MinOpp_ha))
    
    message(sprintf("    🌼 Spoor B toetsen aan %dm bosrand-eis (>= %.1f ha)...", raaklijn_m, min_opp_b))
    res_b <- evalueer_bosrand_contact_universeel(
      r_hab_opp      = r_spoor_b,
      r_bos_opp      = r_bos,
      drempel_ha     = min_opp_b,
      min_raaklijn_m = raaklijn_m,
      master_grid    = master_grid
    )
  }
  
  # 6. ASYMMETRISCHE WEDERZIJDSE AFSTANDSFILTER (Bruine eikenpage)
  if (bestandsnaam == "bruineeikenpage" && !is.null(res_b)) {
    message("    ↔️ Asymmetrische wederzijdse afstandsfilter (250m / 500m)...")
    
    heeft_a <- !all(is.na(suppressWarnings(terra::minmax(res_a$opp))))
    heeft_b <- !all(is.na(suppressWarnings(terra::minmax(res_b$opp))))
    
    if (heeft_a && heeft_b) {
      hab_a_bin <- terra::ifel(!is.na(res_a$opp) & res_a$opp > 0, 1, NA)
      hab_b_bin <- terra::ifel(!is.na(res_b$opp) & res_b$opp > 0, 1, NA)
      
      # 1. Nectar buffer 250m rondom Eiafzet
      f_b_bin <- tempfile(fileext = ".tif")
      terra::writeRaster(hab_b_bin, f_b_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
      f_prox_b250 <- hulp_proximity_gdal(f_b_bin, 250)
      unlink(f_b_bin)
      
      p_b250 <- terra::rast(f_prox_b250)
      mask_b250 <- terra::ifel(!is.na(p_b250) & p_b250 >= 0 & p_b250 <= 250, 1, NA)
      
      # 2. Eiafzet buffer 500m rondom Nectar
      f_a_bin <- tempfile(fileext = ".tif")
      terra::writeRaster(hab_a_bin, f_a_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
      f_prox_a500 <- hulp_proximity_gdal(f_a_bin, 500)
      unlink(f_a_bin)
      
      p_a500 <- terra::rast(f_prox_a500)
      mask_a500 <- terra::ifel(!is.na(p_a500) & p_a500 >= 0 & p_a500 <= 500, 1, NA)
      
      # Maskeer beide deellagen
      res_a$opp      <- terra::mask(res_a$opp, mask_b250)
      res_a$clusters <- terra::mask(res_a$clusters, mask_b250)
      res_a$hab_deel <- terra::mask(res_a$hab_deel, mask_b250)
      
      res_b$opp      <- terra::mask(res_b$opp, mask_a500)
      res_b$clusters <- terra::mask(res_b$clusters, mask_a500)
      res_b$hab_deel <- terra::mask(res_b$hab_deel, mask_a500)
      
      unlink(c(f_prox_b250, f_prox_a500))
      rm(mask_b250, mask_a500, p_b250, p_a500, hab_a_bin, hab_b_bin)
      gc()
    } else {
      message("    ⚠️ Minstens één van de twee sporen (eiafzet of nectar) heeft 0 ha na bosrandtoets.")
      res_a$opp <- master_grid * NA
      res_b$opp <- master_grid * NA
    }
  }
  
  # 7. SAMENVOEGEN TOT DEFINITIEF LEEFGEBIED
  message("    🧩 Deellagen samenvoegen tot definitief leefgebied...")
  w_a <- terra::ifel(is.na(res_a$opp), 0, res_a$opp)
  
  if (!is.null(res_b)) {
    w_b       <- terra::ifel(is.na(res_b$opp), 0, res_b$opp)
    final_opp <- terra::ifel((w_a + w_b) > 0, terra::clamp(w_a + w_b, upper = 1.0), NA)
    
    cl_a <- res_a$clusters
    cl_b <- res_b$clusters
    max_a_id <- suppressWarnings(terra::global(cl_a, "max", na.rm = TRUE)[[1]])
    if (is.na(max_a_id)) max_a_id <- 0
    
    cl_b_uniek <- cl_b + (max_a_id + 1000)
    cl_totaal  <- terra::cover(cl_a, cl_b_uniek)
    bos_zoom   <- terra::cover(res_a$bos_zoom, res_b$bos_zoom)
    hab_b_deel <- res_b$hab_deel
  } else {
    final_opp  <- terra::ifel(w_a > 0, w_a, NA)
    cl_totaal  <- res_a$clusters
    bos_zoom   <- res_a$bos_zoom
    hab_b_deel <- master_grid * NA
  }
  
  return(list(
    final_opp        = final_opp,
    cl_opp           = cl_totaal,
    deellaag_a       = res_a$hab_deel,
    deellaag_b       = hab_b_deel,
    bosranden_meegen = bos_zoom
  ))
}
