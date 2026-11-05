# ==============================================================================
# bouw_biotoop_vennencomplex_libel.R
# Archetype voor vennen- en petgatenlibellen (water 0.02-1ha met >= 4 poelen op 200m x land 5ha met 30% bos op 200m x ruim 30ha)
# Geschikt voor: Gevlekte witsnuitlibel (Leucorrhinia pectoralis), etc.
# ==============================================================================

bouw_biotoop_vennencomplex_libel <- function(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
  bestandsnaam <- tolower(gsub(" ", "", soort_naam))
  soort_cache  <- file.path(paths$dir_cache, bestandsnaam)
  dir.create(soort_cache, showWarnings = FALSE, recursive = TRUE)
  
  f_cache_rds <- file.path(soort_cache, "01_raw_parquet_layers.rds")
  
  # Parameters ophalen uit Excel
  rij_wat <- params_soort %>% dplyr::filter(grepl("voortplanting", tolower(Type))) %>% dplyr::slice(1)
  min_opp_water <- as.numeric(ifelse(nrow(rij_wat) > 0 && !is.na(rij_wat$MinOpp_ha), rij_wat$MinOpp_ha, 0.02))
  max_opp_water <- 1.0 # 1 ha bovengrens voor kleine/middelgrote vennen
  
  rij_nabij <- params_soort %>% dplyr::filter(grepl("land_nabij", tolower(Type))) %>% dplyr::slice(1)
  min_opp_nabij <- as.numeric(ifelse(nrow(rij_nabij) > 0 && !is.na(rij_nabij$MinOpp_ha), rij_nabij$MinOpp_ha, 5.0))
  
  rij_ruim <- params_soort %>% dplyr::filter(grepl("land_ruim", tolower(Type))) %>% dplyr::slice(1)
  min_opp_ruim <- as.numeric(ifelse(nrow(rij_ruim) > 0 && !is.na(rij_ruim$MinOpp_ha), rij_ruim$MinOpp_ha, 30.0))
  
  straal_water_m <- 200 # 200m interactiezone rond water
  
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
  
  r_wat_raw    <- if ("voortplanting_bwk" %in% names(lagen)) lagen[["voortplanting_bwk"]] else (master_grid * NA)
  r_nabij_raw  <- if ("land_nabij"        %in% names(lagen)) lagen[["land_nabij"]]        else (master_grid * NA)
  r_ruim_raw   <- if ("land_ruim"         %in% names(lagen)) lagen[["land_ruim"]]         else (master_grid * NA)
  r_struw_raw  <- if ("struweel1"         %in% names(lagen)) lagen[["struweel1"]]         else (master_grid * NA)
  
  # ----------------------------------------------------------------------------
  # 2. KLEINE WATERVLAKKEN INTEGREREN & OPPERVLAKTEBAND TOEPASSEN (0.02 - 1.0 HA)
  # ----------------------------------------------------------------------------
  message("    💧 Watervlakken koppelen en oppervlakteband (0.02 - 1.0 ha) filteren...")
  
  # Extra kleine watervlakken meenemen indien beschikbaar
  if (!is.null(paths$watervlakken_raster) && file.exists(paths$watervlakken_raster)) {
    f_wv_al <- hulp_align_gdal(paths$watervlakken_raster, master_grid)
    r_wv    <- terra::rast(f_wv_al)
    r_wv_bin <- terra::ifel(!is.na(r_wv) & r_wv > 0, 1.0, NA)
    unlink(f_wv_al)
    rm(r_wv)
    r_wat_tot <- terra::cover(r_wat_raw, r_wv_bin)
    rm(r_wv_bin)
  } else {
    r_wat_tot <- r_wat_raw
  }
  r_wat_tot <- terra::clamp(r_wat_tot, upper = 1.0)
  rm(r_wat_raw)
  gc()
  
  # Clusteren van water op 10 m
  res_wat10 <- cluster_filter_universeel(
    masker            = r_wat_tot,
    opp_laag          = r_wat_tot,
    min_opp_ha        = min_opp_water,
    min_totale_opp_ha = NA,
    afstand_m         = 10,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_wat_cand  <- res_wat10$raster
  cl_wat_cand <- res_wat10$clusters
  rm(res_wat10, r_wat_tot)
  gc()
  
  # Maximale oppervlaktefilter: poelen > 1.0 ha (10.000 m²) uitsluiten
  if (!all(is.na(suppressWarnings(terra::minmax(cl_wat_cand))))) {
    df_opp_wat <- terra::zonal(r_wat_cand, cl_wat_cand, fun = "sum", na.rm = TRUE)
    colnames(df_opp_wat) <- c("ID", "Opp_Sum")
    dt_opp_wat <- data.table::as.data.table(df_opp_wat)
    dt_opp_wat[, Opp_ha := Opp_Sum * 0.01]
    rm(df_opp_wat)
    
    valide_wat_ids <- dt_opp_wat[!is.na(ID) & Opp_ha >= min_opp_water & Opp_ha <= max_opp_water, ID]
    rm(dt_opp_wat)
    
    if (length(valide_wat_ids) > 0) {
      r_wat_geband <- terra::mask(r_wat_cand, cl_wat_cand %in% as.integer(valide_wat_ids))
      cl_wat_geband <- terra::mask(cl_wat_cand, cl_wat_cand %in% as.integer(valide_wat_ids))
    } else {
      r_wat_geband  <- master_grid * NA
      cl_wat_geband <- master_grid * NA
    }
  } else {
    r_wat_geband  <- master_grid * NA
    cl_wat_geband <- master_grid * NA
  }
  rm(r_wat_cand, cl_wat_cand)
  gc()
  
  ha_wat_geband <- calc_ha_exact(r_wat_geband)
  message(sprintf("       -> Waterpartijen binnen 0.02 - 1.0 ha: %.2f ha", ha_wat_geband))
  
  # ----------------------------------------------------------------------------
  # 3. DICHTHEIDSEIS: VENNENCOMPLEX (MINSTENS 4 POELEN BINNEN 200 M)
  # ----------------------------------------------------------------------------
  message("    🔎 Dichtheidseis toetsen: minstens 4 poelen binnen 200 m...")
  
  if (!all(is.na(suppressWarnings(terra::minmax(cl_wat_geband))))) {
    p_wat_poly <- terra::as.polygons(cl_wat_geband, dissolve = TRUE)
    
    if (!is.null(p_wat_poly) && nrow(p_wat_poly) > 0) {
      pts_wat <- terra::centroids(p_wat_poly, inside = TRUE)
      buf_200m <- terra::buffer(pts_wat, width = 200)
      rel_mat <- terra::is.related(buf_200m, pts_wat, "intersects")
      
      pts_wat$buren_count <- rowSums(rel_mat)
      ok_pts <- pts_wat[pts_wat$buren_count >= 4, ] # zichzelf + minstens 3 buren
      
      if (nrow(ok_pts) > 0) {
        id_col_name <- names(p_wat_poly)[1]
        valide_complex_ids <- ok_pts[[id_col_name]]
        r_water_complex <- terra::mask(r_wat_geband, cl_wat_geband %in% as.integer(valide_complex_ids))
      } else {
        r_water_complex <- master_grid * NA
      }
      rm(pts_wat, buf_200m, rel_mat, p_wat_poly)
    } else {
      r_water_complex <- master_grid * NA
    }
  } else {
    r_water_complex <- master_grid * NA
  }
  rm(r_wat_geband, cl_wat_geband)
  gc()
  
  ha_complex <- calc_ha_exact(r_water_complex)
  message(sprintf("       -> Goedgekeurde vennencomplexen (>= 4 poelen / 200m): %.2f ha", ha_complex))
  
  # ----------------------------------------------------------------------------
  # 4. RUIM LANDBIOTOOP CLUSTEREN (>= 30 HA OP 50 M)
  # ----------------------------------------------------------------------------
  message(sprintf("    🌾 Ruim landschapscomplex clusteren (>= %.1f ha, 50m netwerk)...", min_opp_ruim))
  res_ruim <- cluster_filter_universeel(
    masker            = r_ruim_raw,
    opp_laag          = r_ruim_raw,
    min_opp_ha        = min_opp_ruim,
    min_totale_opp_ha = NA,
    afstand_m         = 50,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  r_ruim_30ha <- res_ruim$raster
  rm(res_ruim, r_ruim_raw)
  gc()
  
  # ----------------------------------------------------------------------------
  # 5. WATER-FIRST LANDBIOTOOP (200 M BUFFER, >= 5 HA, >= 30% STRUWEEL)
  # ----------------------------------------------------------------------------
  message("    🌳 Landbiotoop binnen 200m van vennen filteren (>= 5 ha, >= 30% struweel)...")
  
  heeft_complex <- !all(is.na(suppressWarnings(terra::minmax(r_water_complex))))
  heeft_nabij   <- !all(is.na(suppressWarnings(terra::minmax(r_nabij_raw))))
  
  if (heeft_complex && heeft_nabij) {
    # 200m buffer rondom goedgekeurde watercomplexen via GDAL proximity
    f_wc_bin  <- tempfile(fileext = ".tif")
    terra::writeRaster(terra::ifel(!is.na(r_water_complex) & r_water_complex > 0, 1, NA), f_wc_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
    f_prox_wc <- hulp_proximity_gdal(f_wc_bin, straal_water_m)
    unlink(f_wc_bin)
    
    p_wc <- terra::rast(f_prox_wc)
    mask_wc200 <- terra::ifel(!is.na(p_wc) & p_wc >= 0 & p_wc <= straal_water_m, 1, NA)
    unlink(f_prox_wc)
    rm(p_wc)
    
    # Land afsnijden op 200m
    r_nabij_200m <- terra::mask(r_nabij_raw, mask_wc200)
    rm(mask_wc200)
    
    # Clusteren op 10 m, minimum 5 ha
    res_land5 <- cluster_filter_universeel(
      masker            = r_nabij_200m,
      opp_laag          = r_nabij_200m,
      min_opp_ha        = min_opp_nabij,
      min_totale_opp_ha = NA,
      afstand_m         = 10,
      dispersie_m       = 0,
      werkelijk         = TRUE
    )
    r_land_5ha  <- res_land5$raster
    cl_land_5ha <- res_land5$clusters
    rm(res_land5, r_nabij_200m)
    gc()
    
    # Toetsen op minstens 30% struweeldekking
    if (!all(is.na(suppressWarnings(terra::minmax(cl_land_5ha))))) {
      r_struw_clean <- terra::ifel(!is.na(r_struw_raw) & r_struw_raw > 0, r_struw_raw, 0)
      
      df_land_tot <- terra::zonal(r_land_5ha, cl_land_5ha, fun = "sum", na.rm = TRUE)
      colnames(df_land_tot) <- c("ID", "Totaal_Sum")
      
      df_struw_tot <- terra::zonal(r_struw_clean, cl_land_5ha, fun = "sum", na.rm = TRUE)
      colnames(df_struw_tot) <- c("ID", "Struweel_Sum")
      
      dt_check <- merge(data.table::as.data.table(df_land_tot), data.table::as.data.table(df_struw_tot), by = "ID", all.x = TRUE)
      dt_check[is.na(Struweel_Sum), Struweel_Sum := 0]
      dt_check[, Pct_Bos := (Struweel_Sum / Totaal_Sum) * 100]
      
      valide_land_ids <- dt_check[!is.na(ID) & Pct_Bos >= 30, ID]
      rm(df_land_tot, df_struw_tot, dt_check, r_struw_clean)
      
      if (length(valide_land_ids) > 0) {
        r_land_struw_ok <- terra::mask(r_land_5ha, cl_land_5ha %in% as.integer(valide_land_ids))
      } else {
        r_land_struw_ok <- master_grid * NA
      }
    } else {
      r_land_struw_ok <- master_grid * NA
    }
    rm(r_land_5ha, cl_land_5ha)
    gc()
    
    # Inbedding in ruim landschap (30 ha)
    if (!all(is.na(suppressWarnings(terra::minmax(r_land_struw_ok)))) && !all(is.na(suppressWarnings(terra::minmax(r_ruim_30ha))))) {
      r_land_finaal <- terra::mask(r_land_struw_ok, r_ruim_30ha)
    } else {
      r_land_finaal <- r_land_struw_ok
    }
    rm(r_land_struw_ok)
  } else {
    r_land_finaal <- master_grid * NA
  }
  rm(r_nabij_raw, r_struw_raw)
  gc()
  
  ha_land_fin <- calc_ha_exact(r_land_finaal)
  message(sprintf("       -> Goedgekeurd beschut landbiotoop (200m, >= 30%% bos): %.2f ha", ha_land_fin))
  
  # ----------------------------------------------------------------------------
  # 6. SAMENVOEGING TOT DEFINITIEF LEEFGEBIED & NETWERKEN
  # ----------------------------------------------------------------------------
  message("    🧩 Totaal leefgebiedscomplex samenstellen (vennencomplex + beschut land)...")
  final_opp <- terra::cover(r_water_complex, r_land_finaal)
  final_opp <- terra::clamp(final_opp, upper = 1.0)
  
  ha_tot <- calc_ha_exact(final_opp)
  message(sprintf("       -> Definitief leefgebied Gevlekte witsnuitlibel: %.2f ha", ha_tot))
  
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
    unlink(f_bin)
    
    cl_opp <- terra::rast(f_clump)
    try(terra::set.cats(cl_opp, layer = 1, value = NULL), silent = TRUE)
    try(terra::coltab(cl_opp) <- NULL, silent = TRUE)
    cl_opp <- terra::mask(cl_opp, r_fin_bin)
    unlink(f_clump)
  } else {
    cl_opp <- master_grid * NA
  }
  rm(r_fin_bin)
  gc()
  
  return(list(
    final_opp      = final_opp,
    cl_opp         = cl_opp,
    water_complex  = r_water_complex,
    land_200m      = r_land_finaal,
    land_ruim_30ha = r_ruim_30ha
  ))
}
