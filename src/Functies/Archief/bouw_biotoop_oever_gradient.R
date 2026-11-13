# ==============================================================================
# bouw_biotoop_oever_gradient.R
# Archetype voor viseters en oevergebonden soorten (water x oeverlengte >= 1km x 20m zoom)
# ==============================================================================

bouw_biotoop_oever_gradient <- function(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
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
    # Veilig inpakken vóór wegschrijven naar RDS
    lagen_packed <- lapply(lagen, terra::wrap)
    saveRDS(lagen_packed, f_cache_rds)
    message("    💾 [CACHE] Parquet-lagen veilig bewaard.")
  }
  
  r_bwk_water <- if ("bwk" %in% names(lagen)) lagen[["bwk"]] else (master_grid * NA)
  
  # ----------------------------------------------------------------------------
  # 2. BWK-WATER CLUSTEREN OP 5 HA (50M NETWERK)
  # ----------------------------------------------------------------------------
  message("    🌊 Vlakvormig BWK-water clusteren (>= 5 ha, 50m netwerk)...")
  res_bwk <- cluster_filter_universeel(
    masker            = r_bwk_water,
    opp_laag          = r_bwk_water,
    min_opp_ha        = 5,
    min_totale_opp_ha = NA,
    afstand_m         = 50,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_bwk_5ha <- res_bwk$raster
  rm(res_bwk, r_bwk_water)
  gc()
  
  # ----------------------------------------------------------------------------
  # 3. HYDROGRAFISCHE DATA (HUET-ZONES WATERLOPEN)
  # ----------------------------------------------------------------------------
  message("    🐟 Huet-zones waterlopen inlezen (beken, rivieren, kanalen)...")
  if (file.exists(paths$huetzon)) {
    f_huet_aligned <- hulp_align_gdal(paths$huetzon, master_grid)
    r_huet         <- terra::rast(f_huet_aligned)
    
    leg_huet      <- terra::cats(r_huet)[[1]]
    target_labels <- c("Bg", "BgK", "Bk1", "Bk2", "BkK1", "BkK2", "KabPKunst", "Kan", 
                       "KanPKunst", "Kunst", "Mlz", "O1brak", "O1o", "O2zout", "P1", "P2", "Rg", "Rk", "Rzg")
    
    if (!is.null(leg_huet) && "Label" %in% colnames(leg_huet)) {
      valide_huet_ids <- leg_huet$value[leg_huet$Label %in% target_labels]
    } else {
      valide_huet_ids <- 5:10
    }
    
    r_huet_water <- terra::ifel(!is.na(r_huet) & r_huet %in% valide_huet_ids, 0.01, NA)
    unlink(f_huet_aligned)
    rm(r_huet, leg_huet)
  } else {
    r_huet_water <- master_grid * NA
  }
  
  # Fusie BWK water (>= 5 ha) + Huet-waterlopen
  w_bwk  <- terra::ifel(is.na(r_bwk_5ha), 0, r_bwk_5ha)
  w_huet <- terra::ifel(is.na(r_huet_water), 0, r_huet_water)
  r_water_totaal <- terra::ifel((w_bwk + w_huet) > 0, terra::clamp(w_bwk + w_huet, upper = 1.0), NA)
  rm(w_bwk, w_huet, r_bwk_5ha, r_huet_water)
  gc()
  
  ha_water_tot <- calc_ha_exact(r_water_totaal)
  message(sprintf("       -> Totaal waterhabitat vóór oeverfilter: %.2f ha", ha_water_tot))
  
  # ----------------------------------------------------------------------------
  # 4. OEVERLIJNEN ISOLEREN & LENGTETOETS (>= 1000 M AANEENGESLOTEN)
  # ----------------------------------------------------------------------------
  message("    📏 Oeveranalyse: Water-land grenzen isoleren en toetsen op >= 1 km lengte...")
  
  f_cache_oevers <- file.path(paths$dir_cache, "oevers_gevalideerd_1km.tif")
  
  if (file.exists(f_cache_oevers)) {
    message("       -> [CACHE] Goedgekeurde oevers (>= 1 km) ingeladen van schijf...")
    r_oever_1km <- terra::rast(f_cache_oevers)
  } else {
    if (file.exists(paths$watervlakken_raster)) {
      f_water_aligned <- hulp_align_gdal(paths$watervlakken_raster, master_grid)
      r_water_poly    <- terra::rast(f_water_aligned)
      
      # Grenzen isoleren (binnen en buiten oeverlijn)
      r_is_water <- terra::ifel(!is.na(r_water_poly) & r_water_poly == 1, 1, 0)
      r_is_land  <- terra::ifel(!is.na(r_water_poly) & r_water_poly == 0, 1, 0)
      unlink(f_water_aligned)
      rm(r_water_poly)
      
      f_bnd1 <- tempfile(fileext = ".tif")
      f_bnd2 <- tempfile(fileext = ".tif")
      bnd_water <- terra::boundaries(r_is_water, classes = FALSE, inner = TRUE, filename = f_bnd1, overwrite = TRUE)
      bnd_land  <- terra::boundaries(r_is_land,  classes = FALSE, inner = TRUE, filename = f_bnd2, overwrite = TRUE)
      
      r_oever_raw <- terra::ifel((bnd_water == 1 | bnd_land == 1), 1, NA)
      unlink(c(f_bnd1, f_bnd2))
      rm(r_is_water, r_is_land, bnd_water, bnd_land)
      gc()
      
      # Buffer 25m om kleine onderbrekingen te overbruggen
      f_oever_bin <- tempfile(fileext = ".tif")
      terra::writeRaster(r_oever_raw, f_oever_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 0)
      f_prox_oev  <- hulp_proximity_gdal(f_oever_bin, 25)
      unlink(f_oever_bin)
      
      p_oev      <- terra::rast(f_prox_oev)
      mask_oev25 <- terra::ifel(!is.na(p_oev) & p_oev >= 0 & p_oev <= 25, 1, NA)
      unlink(f_prox_oev)
      rm(p_oev)
      gc()
      
      # Clump oevertrajecten via GRASS
      f_buf_bin <- tempfile(fileext = ".tif")
      f_cl_oev  <- tempfile(fileext = ".tif")
      terra::writeRaster(mask_oev25, f_buf_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
      rm(mask_oev25)
      gc()
      
      e <- terra::ext(master_grid)
      qgisprocess::qgis_run_algorithm(
        "grass:r.clump", input = f_buf_bin, output = f_cl_oev, title = "cl_oevers",
        "-d" = FALSE, GRASS_REGION_PARAMETER = sprintf("%f,%f,%f,%f", e$xmin, e$xmax, e$ymin, e$ymax),
        GRASS_REGION_CELLSIZE_PARAMETER = 10, .quiet = TRUE
      )
      unlink(f_buf_bin)
      
      cl_oevers <- terra::rast(f_cl_oev)
      terra::set.cats(cl_oevers, layer = 1, value = NULL)
      terra::coltab(cl_oevers) <- NULL
      
      # Maskeer terug naar zuivere oevercellen
      cl_oever_only <- terra::mask(cl_oevers, r_oever_raw)
      unlink(f_cl_oev)
      rm(r_oever_raw, cl_oevers)
      gc()
      
      # 1 km oeverlengte = minstens 100 oevercellen (10m per cel)
      dt_oev <- data.table::as.data.table(terra::freq(cl_oever_only))
      setnames(dt_oev, c("layer", "ID", "count"))
      valide_oev_ids <- dt_oev[!is.na(ID) & ID > 0 & count >= 100, ID]
      rm(dt_oev)
      
      if (length(valide_oev_ids) > 0) {
        r_oever_1km <- terra::ifel(cl_oever_only %in% as.integer(valide_oev_ids), 1, NA)
      } else {
        r_oever_1km <- master_grid * NA
      }
      rm(cl_oever_only)
      gc()
    } else {
      r_oever_1km <- master_grid * NA
    }
    terra::writeRaster(r_oever_1km, f_cache_oevers, overwrite = TRUE, datatype = "INT1U", gdal = "COMPRESS=LZW")
  }
  
  # ----------------------------------------------------------------------------
  # 5. FUNCTIONELE 20M OEVERSTROOK TOEPASSEN OP HET WATERHABITAT
  # ----------------------------------------------------------------------------
  message("    🎯 20m functionele oeverzone bufferen en toepassen op water...")
  
  if (!all(is.na(suppressWarnings(terra::minmax(r_oever_1km))))) {
    f_oev1_bin  <- tempfile(fileext = ".tif")
    terra::writeRaster(r_oever_1km, f_oev1_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 0)
    f_prox_20m  <- hulp_proximity_gdal(f_oev1_bin, 20)
    unlink(f_oev1_bin)
    
    p_20m       <- terra::rast(f_prox_20m)
    mask_zone20 <- terra::ifel(!is.na(p_20m) & p_20m >= 0 & p_20m <= 20, 1, NA)
    unlink(f_prox_20m)
    rm(p_20m)
    gc()
    
    final_opp <- terra::mask(r_water_totaal, mask_zone20)
    rm(mask_zone20)
    gc()
  } else {
    final_opp <- master_grid * NA
  }
  rm(r_water_totaal, r_oever_1km)
  gc()
  
  # ----------------------------------------------------------------------------
  # 6. NETWERK-ID'S GENEREREN VIA GRASS
  # ----------------------------------------------------------------------------
  r_fin_bin <- terra::ifel(!is.na(final_opp) & final_opp > 0, 1, NA)
  
  if (!all(is.na(suppressWarnings(terra::minmax(r_fin_bin))))) {
    f_bin   <- tempfile(fileext = ".tif")
    f_clump <- tempfile(fileext = ".tif")
    terra::writeRaster(r_fin_bin, f_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    
    e <- terra::ext(master_grid)
    qgisprocess::qgis_run_algorithm(
      "grass:r.clump", input = f_bin, output = f_clump, title = paste0("cl_", bestandsnaam),
      "-d" = FALSE, GRASS_REGION_PARAMETER = sprintf("%f,%f,%f,%f", e$xmin, e$xmax, e$ymin, e$ymax),
      GRASS_REGION_CELLSIZE_PARAMETER = 10, .quiet = TRUE
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
    final_opp = final_opp,
    cl_opp    = cl_opp
  ))
}
