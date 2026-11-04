# ==============================================================================
# bouw_biotoop_complex.R - MODULAIRE AMFIBIEËN- EN COMPLEXE SOORTEN PIPELINE
# ==============================================================================

# Snelle proximity helper via GDAL
hulp_snelle_proximity <- function(r_bin, max_dist_m) {
  if (all(is.na(suppressWarnings(terra::minmax(r_bin))))) return(r_bin * NA)
  
  f_in   <- tempfile(fileext = ".tif")
  f_prox <- tempfile(fileext = ".tif")
  
  terra::writeRaster(r_bin, f_in, overwrite = TRUE, datatype = "INT1U", NAflag = 0)
  
  qgisprocess::qgis_run_algorithm(
    "gdal:proximity",
    INPUT        = f_in,
    BAND         = 1,
    VALUES       = "1",
    UNITS        = 0,
    MAX_DISTANCE = as.numeric(max_dist_m),
    NODATA       = 9999,
    DATA_TYPE    = 2,
    OPTIONS      = "COMPRESS=LZW",
    OUTPUT       = f_prox,
    .quiet       = TRUE
  )
  
  p   <- terra::rast(f_prox)
  out <- terra::ifel(!is.na(p) & p >= 0 & p <= as.numeric(max_dist_m), 1, NA)
  
  unlink(c(f_in, f_prox))
  return(out)
}

bouw_biotoop_complex <- function(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
  bestandsnaam <- tolower(gsub(" ", "", soort_naam))
  f_cache_soort <- file.path(tempdir(), paste0("raw_layers_", bestandsnaam, ".rds"))
  
  # 1. Parquet lagen extraheren voor alle types van deze soort (met caching)
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
  
  # 2. Externe watervlakken toevoegen aan waterbiotoop (<= 4000 m2)
  water_rij <- params_soort %>% 
    dplyr::filter(as.logical(Gebruik_Watervlakken_Shp) == TRUE | grepl("water", tolower(Type))) %>% 
    dplyr::slice(1)
  
  r_water_totaal <- if ("waterbiotoop_bwk" %in% names(lagen)) lagen[["waterbiotoop_bwk"]] else (master_grid * NA)
  kleine_plassen_v <- NULL
  
  if (nrow(water_rij) > 0 && as.logical(water_rij$Gebruik_Watervlakken_Shp) && file.exists(paths$watervlakken)) {
    message("    💧 Watervlakken inlezen en filteren op <= 4000 m²...")
    watervlakken_clean <- terra::vect(paths$watervlakken) %>%
      terra::project(terra::crs(master_grid)) %>%
      terra::makeValid()
    
    watervlakken_clean$opp_m2 <- terra::expanse(watervlakken_clean)
    kleine_plassen_v <- watervlakken_clean[watervlakken_clean$opp_m2 <= 4000, ]
    
    if (nrow(kleine_plassen_v) > 0) {
      r_kleinwater <- terra::rasterize(kleine_plassen_v, master_grid, field = 1, background = NA)
      r_water_totaal <- terra::cover(r_water_totaal, r_kleinwater)
    }
  }
  
  # 3. Landbiotoop opbouwen: Wederzijdse nabijheid (500m) tussen HMS en K
  rijen_wederzijds <- params_soort %>% dplyr::filter(as.logical(Gebruik_Wederzijds) == TRUE)
  land_totaal <- NULL
  
  if (nrow(rijen_wederzijds) >= 2) {
    dist_w <- as.numeric(rijen_wederzijds$Wederzijdse_Afstand_m[1])
    t1 <- rijen_wederzijds$Type[1]
    t2 <- rijen_wederzijds$Type[2]
    
    if (t1 %in% names(lagen) && t2 %in% names(lagen)) {
      message(sprintf("    🌿 Wederzijdse nabijheid (%d m) tussen %s en %s...", dist_w, t1, t2))
      r_bin1 <- terra::ifel(!is.na(lagen[[t1]]) & lagen[[t1]] > 0, 1, NA)
      r_bin2 <- terra::ifel(!is.na(lagen[[t2]]) & lagen[[t2]] > 0, 1, NA)
      
      prox1 <- hulp_snelle_proximity(r_bin1, dist_w)
      prox2 <- hulp_snelle_proximity(r_bin2, dist_w)
      
      land_1 <- terra::mask(lagen[[t1]], prox2)
      land_2 <- terra::mask(lagen[[t2]], prox1)
      land_totaal <- terra::cover(land_1, land_2)
    }
  } else {
    land_types <- setdiff(names(lagen), c("waterbiotoop_bwk", "akkers"))
    if (length(land_types) > 0) {
      land_totaal <- lagen[[land_types[1]]]
      if (length(land_types) > 1) {
        for (lt in land_types[-1]) land_totaal <- terra::cover(land_totaal, lagen[[lt]])
      }
    }
  }
  
  # 4. Akkers uitsluiten van landbiotoop
  if (any(as.logical(params_soort$Filter_Akker), na.rm = TRUE) && "akkers" %in% names(lagen) && !is.null(land_totaal)) {
    message("    🚜 Akkers uitsluiten van landbiotoop...")
    r_akkers_bin <- terra::ifel(!is.na(lagen[["akkers"]]) & lagen[["akkers"]] > 0, 1, NA)
    land_totaal  <- terra::mask(land_totaal, r_akkers_bin, inverse = TRUE)
  }
  
  # 5. Eilandfilter (MinOpp_ha) op het LANDBIOTOOP ALLEEN (bv. >= 5 ha)
  rij_land <- params_soort %>% dplyr::filter(grepl("land", tolower(Type)) & !is.na(MinOpp_ha)) %>% dplyr::slice(1)
  min_land_ha <- if (nrow(rij_land) > 0) as.numeric(rij_land$MinOpp_ha) else 0
  
  if (min_land_ha > 0 && !is.null(land_totaal) && !all(is.na(suppressWarnings(terra::minmax(land_totaal))))) {
    message(sprintf("    🏝️ Landbiotoop filteren op minimale patchgrootte (>= %.1f ha)...", min_land_ha))
    
    r_bin_land <- terra::ifel(!is.na(land_totaal) & land_totaal > 0, 1, NA)
    
    f_in_l <- tempfile(fileext = ".tif")
    f_cl_l <- tempfile(fileext = ".tif")
    terra::writeRaster(r_bin_land, f_in_l, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    
    e <- terra::ext(master_grid)
    extent_str <- sprintf("%f,%f,%f,%f", e$xmin, e$xmax, e$ymin, e$ymax)
    
    qgisprocess::qgis_run_algorithm(
      "grass:r.clump", 
      input                           = f_in_l, 
      output                          = f_cl_l,
      title                           = "clumps_land",
      "-d"                            = FALSE,
      GRASS_REGION_PARAMETER          = extent_str, 
      GRASS_REGION_CELLSIZE_PARAMETER = 10, 
      .quiet                          = TRUE
    )
    
    cl_land <- terra::rast(f_cl_l)
    terra::set.cats(cl_land, layer = 1, value = NULL)
    terra::coltab(cl_land) <- NULL
    
    cells_l <- terra::cells(r_bin_land)
    dt_l <- data.table::data.table(
      ID     = as.integer(cl_land[cells_l][[1]]),
      Waarde = as.numeric(land_totaal[cells_l][[1]])
    )
    
    # Filter op werkelijke oppervlakte >= min_land_ha
    ok_land_ids <- dt_l[!is.na(ID) & ID > 0, .(ha = sum(Waarde, na.rm = TRUE) * 0.01), by = ID][ha >= min_land_ha, ID]
    
    if (length(ok_land_ids) > 0) {
      land_totaal <- terra::mask(land_totaal, terra::ifel(cl_land %in% ok_land_ids, 1, NA))
    } else {
      land_totaal <- master_grid * NA
    }
    unlink(c(f_in_l, f_cl_l))
    rm(dt_l, cl_land, cells_l)
  }
  
  # 6. Samenvoegen tot gecombineerd biotoop
  w_clean <- terra::ifel(is.na(r_water_totaal), 0, r_water_totaal)
  l_clean <- terra::ifel(is.na(land_totaal), 0, land_totaal)
  
  r_som <- terra::clamp(w_clean + l_clean, upper = 1.0)
  r_gecombineerd <- terra::ifel(r_som > 0, r_som, NA)
  
  return(list(
    r_gecombineerd   = r_gecombineerd,
    r_land_deel      = land_totaal,
    r_water_deel     = r_water_totaal,
    kleine_plassen_v = kleine_plassen_v
  ))
}
