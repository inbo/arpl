# ==============================================================================
# bouw_biotoop_universeel.R
# Universele data-gedreven motor voor alle ecologische archetypes (Geoptimaliseerd)
# ==============================================================================

bouw_biotoop_universeel <- function(soort_naam, 
                                    params_soort, 
                                    ds_parquet, 
                                    df_bwk_codes, 
                                    master_grid, 
                                    paths, 
                                    force_overwrite = FALSE, 
                                    ...) {
  
  bestandsnaam <- tolower(gsub(" ", "", soort_naam))
  soort_cache  <- file.path(paths$dir_cache, bestandsnaam)
  dir.create(soort_cache, showWarnings = FALSE, recursive = TRUE)
  
  # ----------------------------------------------------------------------------
  # STAP 1: PARQUET EXTRACTIE VOOR ALLE TYPES VAN DE SOORT
  # ----------------------------------------------------------------------------
  types_te_laden <- unique(params_soort$Type)
  lagen_ruw <- list()
  
  for (tp in types_te_laden) {
    f_tif <- file.path(soort_cache, paste0("laag_", tp, ".tif"))
    if (file.exists(f_tif) && !force_overwrite) {
      message(sprintf("    📦 [CACHE] Laag '%s' ingeladen...", tp))
      lagen_ruw[[tp]] <- terra::rast(f_tif)
    } else {
      message(sprintf("    ⏳ Parquet-laag '%s' extraheren voor %s...", tp, soort_naam))
      r_tp <- extract_bwk_laag(
        soort_naam      = soort_naam, 
        tp              = tp, 
        ds_parquet      = ds_parquet, 
        df_bwk_codes    = df_bwk_codes, 
        master_grid     = master_grid,
        paths           = paths,
        force_overwrite = force_overwrite
      )
      if (is.null(r_tp)) r_tp <- master_grid * NA_real_
      terra::writeRaster(r_tp, f_tif, overwrite = TRUE, datatype = "FLT4S", gdal = "COMPRESS=LZW")
      lagen_ruw[[tp]] <- terra::rast(f_tif)
    }
    message(sprintf("       📊 Ruwe oppervlakte '%s': %.2f ha", tp, calc_ha_exact(lagen_ruw[[tp]])))
  }
  
  # Synchroniseer master_grid altijd op exact de dimensies/extent van lagen_ruw[[1]]
  if (length(lagen_ruw) > 0 && inherits(lagen_ruw[[1]], "SpatRaster")) {
    if (terra::nrow(master_grid) != terra::nrow(lagen_ruw[[1]]) ||
        terra::ncol(master_grid) != terra::ncol(lagen_ruw[[1]]) ||
        !isTRUE(all.equal(as.vector(terra::ext(master_grid)), as.vector(terra::ext(lagen_ruw[[1]]))))) {
      message("    🔍 Subgebied gedetecteerd in lagen_ruw: master_grid synchroniseren met uitsnede...")
      master_grid <- lagen_ruw[[1]] * NA_real_
    }
  }
  
  verwerkte_sporen <- list()
  cl_sporen        <- list()
  
  # ----------------------------------------------------------------------------
  # STAP 2: VERWERKING PER TYPE (RIJ IN EXCEL)
  # ----------------------------------------------------------------------------
  rij_volgorde <- order(
    params_soort$Spoor == "C",
    params_soort$Spoor == "A",
    params_soort$Spoor == "B",
    decreasing = TRUE
  )
  params_gesorteerd <- params_soort[rij_volgorde, ]
  
  for (i in 1:nrow(params_gesorteerd)) {
    rij          <- params_gesorteerd[i, ]
    tp_naam      <- trimws(rij$Type)
    huidig_spoor <- get_param(rij, "Spoor", "A")
    
    min_patch_val <- get_param(rij, "MinPatch_ha", NA)
    if (grepl("^not_", tolower(tp_naam)) || tolower(tp_naam) == "bwkbos" || 
        grepl("_bos$", tolower(tp_naam)) || is.na(min_patch_val)) {
      next
    }
    
    # --------------------------------------------------------------------------
    # 2-weidevogel: Kritische weidevogels (Grutto) - volledig zelfstandig spoor
    # --------------------------------------------------------------------------
    is_weidevogel <- tolower(bestandsnaam) %in% c("grutto", "limosalimosa")
    
    if (is_weidevogel) {
      if (huidig_spoor != "A") next
      if ("kern" %in% names(verwerkte_sporen)) next
      
      message("\n    🔹 Weidevogelmatrix-archetype uitvoeren (natte kern >= 10ha x open matrix in netwerk >= 25ha)...")
      e <- terra::ext(master_grid)
      region_str <- sprintf("%f,%f,%f,%f", e$xmin, e$xmax, e$ymin, e$ymax)
      
      # 1. Bebouwingsbuffer (GDX >= 10 ha op 100m koppeling, 100m verstoringsbuffer)
      m_bebouw_100m <- maak_bebouwingsmasker(
        min_opp_ha      = 10,
        buf_m           = 100,
        afstand_cl_m    = 100,
        ds_parquet      = ds_parquet,
        master_grid     = master_grid,
        paths           = paths,
        force_overwrite = force_overwrite
      )
      
      # 2. Groenkaart 1m -> 10m aggregatie met MAX
      f_cache_bomen_max <- file.path(paths$dir_cache, "masker_groenkaart_10m_max_bomen.tif")
      if (file.exists(f_cache_bomen_max) && !force_overwrite) {
        r_bomen_direct <- terra::rast(f_cache_bomen_max)
      } else {
        r_gr_src <- terra::rast(paths$groenkaart)
        if (round(terra::res(r_gr_src)[1]) == 1) {
          r_gr_1m_bin <- terra::ifel(!is.na(r_gr_src) & r_gr_src == 1, 1, 0)
          r_gr_10m    <- terra::aggregate(r_gr_1m_bin, fact = 10, fun = "max", na.rm = TRUE)
          r_gr_10m    <- terra::resample(r_gr_10m, master_grid, method = "near")
          r_bomen_direct <- terra::ifel(!is.na(r_gr_10m) & r_gr_10m == 1, 1L, NA_integer_)
          rm(r_gr_1m_bin, r_gr_10m)
        } else {
          f_gr_al <- hulp_align_gdal(paths$groenkaart, master_grid)
          r_gr_raw <- terra::rast(f_gr_al)
          r_bomen_direct <- terra::ifel(!is.na(r_gr_raw) & r_gr_raw == 1, 1L, NA_integer_)
          unlink(f_gr_al)
          rm(r_gr_raw)
        }
        rm(r_gr_src)
        terra::writeRaster(r_bomen_direct, f_cache_bomen_max, overwrite = TRUE, datatype = "INT2S", NAflag = -9999, gdal = "COMPRESS=LZW")
      }
      
      # 3. Boskernen (>= 12 bomen in 5x5 venster) + 50m verstoringsbuffer
      f_cache_groen50 <- file.path(paths$dir_cache, "masker_groen_50m_cluster_max.tif")
      if (file.exists(f_cache_groen50) && !force_overwrite) {
        m_groen_50m <- terra::rast(f_cache_groen50)
      } else {
        r_gr_bin <- terra::ifel(!is.na(r_bomen_direct) & r_bomen_direct == 1, 1, 0)
        w_mat5x5 <- matrix(1, nrow = 5, ncol = 5)
        f_fc_g   <- tempfile(fileext = ".tif")
        
        gr_som <- terra::focal(r_gr_bin, w = w_mat5x5, fun = "sum", na.rm = TRUE, filename = f_fc_g, overwrite = TRUE)
        rm(r_gr_bin)
        
        r_gr_cluster <- terra::ifel(!is.na(gr_som) & gr_som >= 12 & !is.na(r_bomen_direct) & r_bomen_direct == 1, 1, NA)
        rm(gr_som)
        unlink(f_fc_g)
        
        if (!is_leeg(r_gr_cluster)) {
          f_g_tmp <- tempfile(fileext = ".tif")
          terra::writeRaster(r_gr_cluster, f_g_tmp, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
          rm(r_gr_cluster)
          f_prox_g <- hulp_proximity_gdal(f_g_tmp, 50)
          unlink(f_g_tmp)
          p_g <- terra::rast(f_prox_g)
          m_groen_50m <- terra::ifel(!is.na(p_g) & p_g >= 0 & p_g <= 50, 1L, NA_integer_)
          unlink(f_prox_g)
          rm(p_g)
        } else {
          m_groen_50m <- master_grid * NA_real_
        }
        terra::writeRaster(m_groen_50m, f_cache_groen50, overwrite = TRUE, datatype = "INT2S", NAflag = -9999, gdal = "COMPRESS=LZW")
      }
      
      # Combineer storings- en uitsluitingsmaskers veilig met INT2S + NAflag = -9999
      f_stoor_tmp <- tempfile(fileext = ".tif")
      m_stoor_totaal <- terra::writeRaster(
        terra::ifel((!is.na(m_bebouw_100m) & m_bebouw_100m == 1) |
                      (!is.na(m_groen_50m) & m_groen_50m == 1) |
                      (!is.na(r_bomen_direct) & r_bomen_direct == 1), 1L, NA_integer_),
        f_stoor_tmp, overwrite = TRUE, datatype = "INT2S", NAflag = -9999
      )
      rm(m_bebouw_100m, m_groen_50m, r_bomen_direct)
      gc()
      message(sprintf("       * Totale verstoringsuitsluiting (bebouwing/bomen): %.2f ha", calc_ha_exact(m_stoor_totaal)))
      
      # 4. Natte kernen opbouwen (kern of kern1 + kern2 -> drainage -> verstoring uitsluiten -> >= 10 ha op 50m)
      message("       -> Stap 4/6: Natte kernen (>= 10 ha) filteren...")
      if ("kern" %in% names(lagen_ruw)) {
        r_kern_bruto <- lagen_ruw[["kern"]]
      } else {
        r_k1 <- if ("kern1" %in% names(lagen_ruw)) lagen_ruw[["kern1"]] else (master_grid * 0)
        r_k2 <- if ("kern2" %in% names(lagen_ruw)) lagen_ruw[["kern2"]] else (master_grid * 0)
        r_kern_bruto <- voeg_toe_geclamped(r_k1, r_k2)
        rm(r_k1, r_k2)
      }
      message(sprintf("          * Bruto kern: %.2f ha", calc_ha_exact(r_kern_bruto)))
      
      if (!is.null(paths$drainage) && file.exists(paths$drainage)) {
        drain_letters <- c("d", "e", "f", "h", "i", "g", "e-f", "h-i", "e-i")
        m_drain <- haal_categorisch_masker(paths$drainage, drain_letters, master_grid, paths$dir_cache, "drainage")
        r_kern_bruto <- terra::mask(r_kern_bruto, m_drain)
        rm(m_drain)
        message(sprintf("          * Kern na drainagefilter: %.2f ha", calc_ha_exact(r_kern_bruto)))
      }
      
      r_kern_clean <- terra::mask(r_kern_bruto, m_stoor_totaal, inverse = TRUE)
      rm(r_kern_bruto)
      message(sprintf("          * Kern na verstoringsuitsluiting: %.2f ha", calc_ha_exact(r_kern_clean)))
      
      res_kern10 <- cluster_filter_universeel(
        masker = r_kern_clean, opp_laag = r_kern_clean, min_opp_ha = 10.0,
        min_totale_opp_ha = NA, afstand_m = 50, dispersie_m = 0, werkelijk = TRUE
      )
      f_kern10_tmp <- tempfile(fileext = ".tif")
      r_kern_10ha <- terra::writeRaster(res_kern10$raster, f_kern10_tmp, overwrite = TRUE, datatype = "FLT4S")
      rm(res_kern10, r_kern_clean)
      gc()
      message(sprintf("          * Kernclusters >= 10 ha: %.2f ha", calc_ha_exact(r_kern_10ha)))
      
      # 5. Open matrix opbouwen
      message("       -> Stap 5/6: Open graslandmatrix samenstellen...")
      r_cultuur <- if ("cultuur0a" %in% names(lagen_ruw)) lagen_ruw[["cultuur0a"]] else (master_grid * NA_real_)
      r_akker   <- if ("akker" %in% names(lagen_ruw)) lagen_ruw[["akker"]] else (master_grid * NA_real_)
      r_hooi    <- if ("hooiland" %in% names(lagen_ruw)) lagen_ruw[["hooiland"]] else (master_grid * NA_real_)
      r_c1b     <- if ("cultuurhooi1b" %in% names(lagen_ruw)) lagen_ruw[["cultuurhooi1b"]] else (master_grid * NA_real_)
      r_w_bwk   <- if ("kleinwater_bwk" %in% names(lagen_ruw)) lagen_ruw[["kleinwater_bwk"]] else (master_grid * NA_real_)
      
      f_cache_grb_kw <- file.path(paths$dir_cache, "masker_grb_kleinwater_50px.tif")
      if (file.exists(f_cache_grb_kw) && !force_overwrite) {
        r_w_grb <- terra::rast(f_cache_grb_kw)
      } else if (!is.null(paths$grb_water) && file.exists(paths$grb_water)) {
        f_grb_al <- hulp_align_gdal(paths$grb_water, master_grid)
        r_grb    <- terra::rast(f_grb_al)
        r_grb_b  <- terra::ifel(!is.na(r_grb) & r_grb > 0, 1, NA)
        unlink(f_grb_al)
        rm(r_grb)
        
        cl_w <- clump_fast(r_grb_b, diag = FALSE, titel = "cl_water_klein")
        dt_w <- data.table::as.data.table(terra::freq(cl_w))
        data.table::setnames(dt_w, c("layer", "ID", "count"))
        klein_w_ids <- dt_w[!is.na(ID) & ID > 0 & count <= 50, ID]
        r_w_grb <- if (length(klein_w_ids) > 0) terra::ifel(cl_w %in% as.integer(klein_w_ids), 0.01, NA_real_) else (master_grid * NA_real_)
        unlink(attr(cl_w, "tmp_file"))
        rm(cl_w, dt_w, r_grb_b)
        
        terra::writeRaster(r_w_grb, f_cache_grb_kw, overwrite = TRUE, datatype = "FLT4S", gdal = "COMPRESS=LZW")
      } else {
        r_w_grb <- master_grid * NA_real_
      }
      
      # Combineer matrix-lagen stapsgewijs
      r_mat_som <- r_cultuur
      r_mat_som <- voeg_toe_geclamped(r_mat_som, r_hooi)
      r_mat_som <- voeg_toe_geclamped(r_mat_som, r_w_bwk)
      r_mat_som <- voeg_toe_geclamped(r_mat_som, r_w_grb)
      r_mat_som <- voeg_toe_geclamped(r_mat_som, r_c1b)
      rm(r_cultuur, r_hooi, r_w_bwk, r_w_grb, r_c1b)
      
      r_mat_som <- terra::mask(r_mat_som, m_stoor_totaal, inverse = TRUE)
      if (!is_leeg(r_akker)) {
        r_mat_som <- terra::mask(r_mat_som, terra::ifel(!is.na(r_akker) & r_akker > 0, 1L, NA_integer_), inverse = TRUE)
      }
      f_mat_tmp <- tempfile(fileext = ".tif")
      r_matrix_clean <- terra::writeRaster(r_mat_som, f_mat_tmp, overwrite = TRUE, datatype = "FLT4S")
      rm(r_mat_som, m_stoor_totaal, r_akker)
      unlink(f_stoor_tmp)
      gc()
      message(sprintf("          * Open graslandmatrix schoon: %.2f ha", calc_ha_exact(r_matrix_clean)))
      
      # 6. Integrale netwerkkoppeling (100% streaming C++ zonal)
      message("       -> Stap 6/6: Netwerkkoppeling (200m overbrugging, totaal >= 25 ha & kern >= 10 ha)...")
      if (!is_leeg(r_kern_10ha)) {
        r_comb_val <- voeg_toe_geclamped(r_kern_10ha, r_matrix_clean)
        
        f_cb_bin <- tempfile(fileext = ".tif")
        terra::writeRaster(terra::ifel(!is.na(r_comb_val) & r_comb_val > 0, 1, NA), f_cb_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
        f_prox_200 <- hulp_proximity_gdal(f_cb_bin, 100)
        unlink(f_cb_bin)
        
        p_200 <- terra::rast(f_prox_200)
        r_buf_200 <- terra::ifel(!is.na(p_200) & p_200 >= 0 & p_200 <= 100, 1, NA)
        unlink(f_prox_200)
        rm(p_200)
        
        cl_net <- clump_fast(r_buf_200, diag = FALSE, titel = "cl_grutto_net")
        rm(r_buf_200)
        
        df_tot <- terra::zonal(r_comb_val, cl_net, fun = "sum", na.rm = TRUE)
        colnames(df_tot) <- c("netwerk_id", "som_tot")
        
        df_kern <- terra::zonal(r_kern_10ha, cl_net, fun = "sum", na.rm = TRUE)
        colnames(df_kern) <- c("netwerk_id", "som_kern")
        
        dt_stats <- merge(data.table::as.data.table(df_tot), data.table::as.data.table(df_kern), by = "netwerk_id", all.x = TRUE)
        dt_stats[is.na(som_kern), som_kern := 0]
        dt_stats[, totaal_ha := som_tot * 0.01]
        dt_stats[, kern_ha   := som_kern * 0.01]
        
        valide_netwerken <- dt_stats[!is.na(netwerk_id) & netwerk_id > 0 & totaal_ha >= 25 & kern_ha >= 10, netwerk_id]
        rm(df_tot, df_kern, dt_stats)
        
        if (length(valide_netwerken) > 0) {
          m_net_ok  <- terra::ifel(cl_net %in% as.integer(valide_netwerken), 1, NA)
          final_opp <- terra::mask(r_comb_val, m_net_ok)
          kern_out  <- terra::mask(r_kern_10ha, m_net_ok)
          mat_out   <- terra::mask(r_matrix_clean, m_net_ok)
          rm(m_net_ok)
        } else {
          final_opp <- master_grid * NA_real_
          kern_out  <- master_grid * NA_real_
          mat_out   <- master_grid * NA_real_
        }
        unlink(attr(cl_net, "tmp_file"))
        rm(cl_net, r_comb_val)
      } else {
        final_opp <- master_grid * NA_real_
        kern_out  <- master_grid * NA_real_
        mat_out   <- master_grid * NA_real_
      }
      rm(r_kern_10ha, r_matrix_clean)
      unlink(c(f_kern10_tmp, f_mat_tmp))
      gc()
      
      verwerkte_sporen[["kern"]]   <- kern_out
      verwerkte_sporen[["matrix"]] <- mat_out
      
      message(sprintf("       📊 Na weidevogel-matrixanalyse: %.2f ha (Kern: %.2f ha, Matrix: %.2f ha)", 
                      calc_ha_exact(final_opp), calc_ha_exact(kern_out), calc_ha_exact(mat_out)))
      next
    }
    
    message(sprintf("\n    🔹 Verwerken van laag/spoor: %s (Spoor %s)...", tp_naam, huidig_spoor))
    r_huidig <- if (tp_naam %in% names(lagen_ruw)) lagen_ruw[[tp_naam]] else (master_grid * NA_real_)
    
    # 2-fusie: Types binnen hetzelfde spoor sommatief fuseren
    zelfde_spoor_types <- params_gesorteerd %>% 
      dplyr::filter(Spoor == huidig_spoor & 
                      !is.na(MinPatch_ha) & 
                      !grepl("^not_|_notbwk$", tolower(Type)) & 
                      !tolower(Type) %in% c("bosranden", "struweel1", "landbiotoop_bos") & 
                      !grepl("_bos$", tolower(Type))) %>% 
      dplyr::pull(Type)
    
    if (length(zelfde_spoor_types) > 1 && tp_naam == zelfde_spoor_types[1]) {
      message(sprintf("       -> Deellagen voor Spoor %s sommatief samenvoegen (%s)...", 
                      huidig_spoor, paste(zelfde_spoor_types, collapse = " + ")))
      r_som_spoor <- NULL
      for (st in zelfde_spoor_types) {
        if (st %in% names(lagen_ruw)) {
          r_som_spoor <- voeg_toe_geclamped(r_som_spoor, lagen_ruw[[st]])
        }
      }
      r_huidig <- if (is.null(r_som_spoor)) (master_grid * NA_real_) else r_som_spoor
      rm(r_som_spoor)
    } else if (length(zelfde_spoor_types) > 1 && tp_naam != zelfde_spoor_types[1]) {
      next
    }
    
    # 2-poelamfibie-kleinwater: kleine watervlakken (<= 0.40 ha op 5m) voor poelamfibieën + Boomkikker 500m koppeling
    pad_watervlakken <- get_path("watervlakken_raster", "watervlakken")
    
    if (tolower(bestandsnaam) %in% c("boomkikker", "hylaarborea",
                                     "heikikker", "ranaarvalis",
                                     "kamsalamander", "trituruscristatus",
                                     "knoflookpad", "pelobatesfuscus",
                                     "poelkikker", "pelophylaxlessonae",
                                     "rugstreeppad", "epidaleacalamita", "bufocalamita",
                                     "vroedmeesterpad", "alytesobstetricans") &&
        huidig_spoor == "A" && !is.null(pad_watervlakken)) {
      message("       -> Amfibie-specifiek: kleine watervlakken (<= 0.40 ha op 5m) toevoegen...")
      f_wv_al <- hulp_align_gdal(pad_watervlakken, master_grid)
      r_wv    <- terra::rast(f_wv_al)
      m_wv    <- terra::ifel(!is.na(r_wv) & r_wv > 0, r_wv, NA_real_)
      unlink(f_wv_al)
      rm(r_wv)
      
      if (!is_leeg(m_wv)) {
        res_wv_all <- cluster_filter_universeel(m_wv, m_wv, min_opp_ha = 0, afstand_m = 5, werkelijk = TRUE)
        cl_wv <- res_wv_all$clusters
        df_wv <- terra::zonal(m_wv, cl_wv, fun = "sum", na.rm = TRUE)
        ok_wv_ids <- df_wv[df_wv[, 2] * 0.01 <= 0.40, 1]
        if (length(ok_wv_ids) > 0) {
          m_wv_klein <- terra::mask(m_wv, cl_wv %in% as.integer(ok_wv_ids))
          r_huidig <- voeg_toe_geclamped(r_huidig, m_wv_klein)
          rm(m_wv_klein)
        }
        rm(res_wv_all, cl_wv, df_wv)
      }
      rm(m_wv)
    }
    
    if (tolower(bestandsnaam) %in% c("boomkikker", "hylaarborea") &&
        huidig_spoor == "B" && "landbiotoop_hms" %in% names(lagen_ruw) && "landbiotoop_k" %in% names(lagen_ruw)) {
      message("       -> Boomkikker-specifiek: struweel (k) en kruidenlaag (hms) wederzijds koppelen binnen 500 m...")
      res_bk500 <- hulp_wederzijdse_koppeling(
        laag_a      = lagen_ruw[["landbiotoop_hms"]],
        laag_b      = lagen_ruw[["landbiotoop_k"]],
        dist_a_m    = 500,
        dist_b_m    = 500,
        master_grid = master_grid
      )
      r_huidig <- voeg_toe_geclamped(res_bk500$laag_a_gekoppeld, res_bk500$laag_b_gekoppeld)
      rm(res_bk500)
    }
    
    # 2-lgp-basis
    lgp_als_basis <- FALSE
    if (is_leeg(r_huidig) && !is.na(get_param(rij, "LGP_Gewascodes", NA))) {
      lgp_str <- sub("^\\+\\s*", "", trimws(as.character(get_param(rij, "LGP_Gewascodes", NA))))
      codes_v <- strsplit(lgp_str, "[,; ]+")[[1]]
      message(sprintf("       -> Basislaag initialiseren vanuit LGP-teelten (%s)...", paste(codes_v, collapse = ", ")))
      r_huidig <- haal_categorisch_masker(paths$lgp_landbouw, codes_v, master_grid, paths$dir_cache, "lgp")
      lgp_als_basis <- TRUE
    }
    
    # 2-water: GRB waternabijheid
    if (isTRUE(as.logical(get_param(rij, "Gebruik_GRB_Water", FALSE)))) {
      water_dist_m <- as.numeric(get_param(rij, "Water_Afstand_m", NA))
      if (!is.na(water_dist_m) && water_dist_m > 0 && !is.null(paths$grb_water) && file.exists(paths$grb_water)) {
        message(sprintf("       -> Afstandsfilter tot open water (<= %dm)...", as.integer(water_dist_m)))
        f_grb_al  <- hulp_align_gdal(paths$grb_water, master_grid)
        r_grb     <- terra::rast(f_grb_al)
        r_grb_bin <- terra::ifel(!is.na(r_grb) & r_grb > 0, 1, NA)
        unlink(f_grb_al)
        rm(r_grb)
        
        r_wat_tot <- if ("water_bwk" %in% names(lagen_ruw)) voeg_toe_geclamped(r_grb_bin, lagen_ruw[["water_bwk"]]) else r_grb_bin
        rm(r_grb_bin)
        
        f_wb_bin  <- tempfile(fileext = ".tif")
        terra::writeRaster(terra::ifel(!is.na(r_wat_tot), 1, NA), f_wb_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
        f_prox_wb <- hulp_proximity_gdal(f_wb_bin, water_dist_m)
        unlink(f_wb_bin)
        rm(r_wat_tot)
        
        p_wb <- terra::rast(f_prox_wb)
        m_water_buf <- terra::ifel(!is.na(p_wb) & p_wb >= 0 & p_wb <= water_dist_m, 1, NA)
        unlink(f_prox_wb)
        rm(p_wb)
        
        r_huidig <- terra::mask(r_huidig, m_water_buf)
        rm(m_water_buf)
        message(sprintf("       📊 Na waternabijheidsfilter (<= %dm): %.2f ha", as.integer(water_dist_m), calc_ha_exact(r_huidig)))
      }
    }
    
    # Kwaliteitsfilter fractie
    if (grepl("waterbiotoop_bwk1", tolower(tp_naam))) {
      r_huidig <- terra::ifel(!is.na(r_huidig) & r_huidig >= 0.20, r_huidig, NA_real_)
    }
    
    # 2-akker: Akker-exclusie
    if (isTRUE(as.logical(get_param(rij, "Filter_Akker", FALSE)))) {
      message("       -> Akker-exclusie toepassen (^b codes uitsluiten)...")
      r_akker_bin <- vul_raster_vanuit_parquet(
        ds_parquet  = ds_parquet,
        codes_bevat = "b",
        master_grid = master_grid,
        waarde_vast = 1
      )
      if (!is_leeg(r_akker_bin)) {
        r_huidig <- terra::mask(r_huidig, r_akker_bin, inverse = TRUE)
      }
      rm(r_akker_bin)
      message(sprintf("       📊 Na akker-exclusie: %.2f ha", calc_ha_exact(r_huidig)))
    }
    
    # 2A. Aftrek negatieve types
    types_not <- params_soort %>%
      dplyr::filter(grepl("^not_", tolower(trimws(Type))) | 
                      grepl("_notbwk$", tolower(trimws(Type))) | 
                      (tolower(trimws(Type)) %in% c("openwater", "terril") & is.na(MinPatch_ha))) %>%
      dplyr::pull(Type)
    
    not_matches <- types_not[
      tolower(types_not) == paste0("not_", tolower(tp_naam)) |
        tolower(types_not) == gsub("_bwk1$", "_notbwk", tolower(tp_naam)) |
        (tolower(tp_naam) == "voortplanting0" & tolower(types_not) == "not_eikenbos") |
        (tolower(tp_naam) == "voortplanting" & tolower(types_not) == "not_eikenbos") |
        (params_soort$Spoor[params_soort$Type %in% types_not] == huidig_spoor)
    ]
    
    if (length(not_matches) > 0) {
      for (neg_tp in unique(not_matches)) {
        if (neg_tp %in% names(lagen_ruw)) {
          message(sprintf("       -> Aftrek negatief type '%s' van '%s'...", neg_tp, tp_naam))
          r_neg <- terra::ifel(is.na(lagen_ruw[[neg_tp]]), 0, lagen_ruw[[neg_tp]])
          r_huidig <- terra::clamp(r_huidig - r_neg, lower = 0, upper = 1.0)
          r_huidig <- terra::ifel(r_huidig > 0, r_huidig, NA_real_)
          rm(r_neg)
        }
      }
    }
    
    # 2B. Externe watervlakken
    pad_watervlakken <- get_path("watervlakken_raster", "watervlakken")
    
    if (isTRUE(as.logical(get_param(rij, "Gebruik_Watervlakken_Shp", FALSE))) && !is.null(pad_watervlakken)) {
      message("       -> Externe watervlakken toevoegen...")
      f_wv_al  <- hulp_align_gdal(pad_watervlakken, master_grid)
      r_wv     <- terra::rast(f_wv_al)
      r_wv_bin <- terra::ifel(!is.na(r_wv) & r_wv > 0, 1.0, NA_real_)
      unlink(f_wv_al)
      rm(r_wv)
      
      if (length(not_matches) > 0) {
        for (neg_tp in unique(not_matches)) {
          if (neg_tp %in% names(lagen_ruw)) {
            r_neg_wv <- terra::ifel(is.na(lagen_ruw[[neg_tp]]), 0, lagen_ruw[[neg_tp]])
            r_wv_bin <- terra::clamp(r_wv_bin - r_neg_wv, lower = 0, upper = 1.0)
            r_wv_bin <- terra::ifel(r_wv_bin > 0, r_wv_bin, NA_real_)
          }
        }
      }
      r_huidig <- voeg_toe_geclamped(r_huidig, r_wv_bin)
      rm(r_wv_bin)
    }
    
    # 2C. Huet-waterlopen
    if (isTRUE(as.logical(get_param(rij, "Gebruik_Huetzones", FALSE)))) {
      huet_lbls <- get_param(rij, "Huet_Labels", NA)
      if (!is.na(huet_lbls) && !is.null(paths$huetzon) && file.exists(paths$huetzon)) {
        huet_str <- tolower(trimws(as.character(huet_lbls)))
        
        if (huet_str %in% c("alle", "all", ">=1", ">= 1")) {
          message("       -> Alle Huet-waterloopvakken toevoegen (waarde >= 1)...")
          f_huet_al <- hulp_align_gdal(paths$huetzon, master_grid)
          r_huet    <- terra::rast(f_huet_al)
          m_huet    <- terra::ifel(!is.na(r_huet) & r_huet >= 1, 1.0, NA_real_)
          unlink(f_huet_al)
          rm(r_huet)
          r_huidig  <- voeg_toe_geclamped(r_huidig, m_huet)
          rm(m_huet)
        } else {
          codes_v <- strsplit(as.character(huet_lbls), "[,; ]+")[[1]]
          message(sprintf("       -> Huet-waterlopen toevoegen (%s)...", paste(codes_v, collapse = ", ")))
          m_huet  <- haal_categorisch_masker(paths$huetzon, codes_v, master_grid, paths$dir_cache, "huet")
          r_huidig <- voeg_toe_geclamped(r_huidig, m_huet)
          rm(m_huet)
        }
      }
    }
    
    # 2C-bis. Maas toevoegen
    if (isTRUE(as.logical(get_param(rij, "Gebruik_Maas", FALSE))) && 
        !is.null(paths$maas) && file.exists(paths$maas)) {
      message("       -> Maas toevoegen aan watercomplex...")
      is_vector <- grepl("\\.(shp|gpkg|geojson)$", paths$maas, ignore.case = TRUE)
      
      if (is_vector) {
        f_maas_rast <- tempfile(fileext = ".tif")
        e <- terra::ext(master_grid)
        qgisprocess::qgis_run_algorithm(
          "gdal:rasterize", INPUT = paths$maas, BURN = 1, UNITS = 1, WIDTH = 10, HEIGHT = 10,
          EXTENT = sprintf("%f,%f,%f,%f", e$xmin, e$xmax, e$ymin, e$ymax),
          NODATA = 255, DATA_TYPE = 0, INIT = 255, OUTPUT = f_maas_rast, .quiet = TRUE
        )
        r_maas_raw <- terra::rast(f_maas_rast)
        r_maas_ok  <- terra::ifel(!is.na(r_maas_raw) & r_maas_raw == 1, 1.0, NA_real_)
        unlink(f_maas_rast)
        rm(r_maas_raw)
      } else {
        f_maas_al  <- hulp_align_gdal(paths$maas, master_grid)
        r_maas_raw <- terra::rast(f_maas_al)
        r_maas_ok  <- terra::ifel(!is.na(r_maas_raw) & r_maas_raw >= 1, 1.0, NA_real_)
        unlink(f_maas_al)
        rm(r_maas_raw)
      }
      
      if (length(not_matches) > 0) {
        for (neg_tp in unique(not_matches)) {
          if (neg_tp %in% names(lagen_ruw)) {
            r_neg_m   <- terra::ifel(is.na(lagen_ruw[[neg_tp]]), 0, lagen_ruw[[neg_tp]])
            r_maas_ok <- terra::clamp(r_maas_ok - r_neg_m, lower = 0, upper = 1.0)
            r_maas_ok <- terra::ifel(r_maas_ok > 0, r_maas_ok, NA_real_)
          }
        }
      }
      r_huidig <- voeg_toe_geclamped(r_huidig, r_maas_ok)
      rm(r_maas_ok)
    }
    
    # 2C-ter. VHA-waterlopen (cat. 0, 1, 2) toevoegen via Excel-kolom 'Gebruik_VHA'
    if (isTRUE(as.logical(get_param(rij, "Gebruik_VHA", FALSE)))) {
      pad_vha <- get_path("vha_waterlopen")
      if (!is.null(pad_vha)) {
        message("       -> VHA-waterlopen (cat. 0, 1, 2) toevoegen...")
        f_vha_al <- hulp_align_gdal(pad_vha, master_grid)
        r_vha    <- terra::rast(f_vha_al)
        m_vha    <- terra::ifel(!is.na(r_vha) & r_vha > 0, r_vha, NA_real_)
        unlink(f_vha_al)
        rm(r_vha)
        
        r_huidig <- voeg_toe_geclamped(r_huidig, m_vha)
        
        # Bever-specifiek: watervlakken > 1000 m² binnen 500 m van VHA
        if (tolower(bestandsnaam) %in% c("europesebever", "castorfiber") &&
            !is.null(pad_watervlakken) && file.exists(pad_watervlakken) && !is_leeg(m_vha)) {
          message("       -> Bever-specifiek: watervlakken > 1000 m² binnen 500 m van VHA zoeken...")
          f_wv_al <- hulp_align_gdal(pad_watervlakken, master_grid)
          r_wv    <- terra::rast(f_wv_al)
          m_wv    <- terra::ifel(!is.na(r_wv) & r_wv > 0, r_wv, NA_real_)
          unlink(f_wv_al)
          rm(r_wv)
          
          if (!is_leeg(m_wv)) {
            # Buffer van 500 m rond VHA-waterlopen (via snelle GDAL proximity)
            f_vha_bin <- tempfile(fileext = ".tif")
            f_vha_dst <- tempfile(fileext = ".tif")
            terra::writeRaster(terra::ifel(!is.na(m_vha), 1L, 0L), f_vha_bin, overwrite = TRUE, datatype = "INT1U")
            qgisprocess::qgis_run_algorithm(
              "gdal:proximity",
              INPUT     = f_vha_bin,
              VALUES    = "1",
              UNITS     = 0, # georeferenced units (meters)
              MAX_DIST  = 500,
              NODATA    = -9999,
              DATA_TYPE = 6,
              OUTPUT    = f_vha_dst,
              .quiet    = TRUE
            )
            r_dst <- terra::rast(f_vha_dst)
            m_buf500 <- terra::ifel(!is.na(r_dst) & r_dst <= 500, 1L, NA_integer_)
            unlink(c(f_vha_bin, f_vha_dst))
            rm(r_dst)
            
            # Clump watervlakken en filter op opp > 1000 m² (0.10 ha) én overlap met 500m-buffer
            r_clumps <- terra::patches(terra::ifel(!is.na(m_wv), 1L, NA_integer_), directions = 8)
            opp_tbl  <- terra::zonal(m_wv, r_clumps, fun = "sum", na.rm = TRUE)
            # 1 cel = 100 m² -> som van fracties * 100 = oppervlakte in m²
            geldige_ids <- opp_tbl[opp_tbl[, 2] * 100 > 1000, 1]
            
            if (length(geldige_ids) > 0) {
              r_wv_groot <- terra::ifel(r_clumps %in% geldige_ids, r_clumps, NA_integer_)
              # Welke van die grote plassen raken de 500m-buffer rond VHA?
              ids_binnen_500 <- unique(terra::values(terra::mask(r_wv_groot, m_buf500), na.rm = TRUE))
              if (length(ids_binnen_500) > 0) {
                m_wv_ok <- terra::ifel(r_clumps %in% ids_binnen_500, m_wv, NA_real_)
                r_huidig <- voeg_toe_geclamped(r_huidig, m_wv_ok)
                rm(m_wv_ok)
              }
              rm(r_wv_groot)
            }
            rm(r_clumps, m_buf500)
          }
          rm(m_wv)
        }
        rm(m_vha)
        message(sprintf("       📊 Na toevoegen VHA-hydrografie: %.2f ha", calc_ha_exact(r_huidig)))
      }
    }
    
    # 2C-oever: Zelfstandige oeverlijnanalyse
    min_oever_l <- as.numeric(dplyr::coalesce(
      as.numeric(get_param(rij, "Min_Rand_Lengte_m", NA)),
      as.numeric(get_param(rij, "Min_Oever_Lengte_m", NA)),
      0
    ))
    oever_zoom_m <- as.numeric(dplyr::coalesce(
      as.numeric(get_param(rij, "Randzone_Diepte_m", NA)),
      as.numeric(get_param(rij, "Grens_Afstand_m", NA)),
      20
    ))
    
    if (!is.na(min_oever_l) && min_oever_l > 0 && is.na(get_param(rij, "Grens_Met_Type", NA)) && 
        !is.null(paths$watervlakken_raster) && file.exists(paths$watervlakken_raster)) {
      
      message(sprintf("       -> Oeveranalyse: Watergrenzen isoleren en toetsen op >= %dm lengte...", as.integer(min_oever_l)))
      f_cache_oevers <- file.path(paths$dir_cache, sprintf("oevers_gevalideerd_%dm.tif", as.integer(min_oever_l)))
      
      if (file.exists(f_cache_oevers) && !force_overwrite) {
        message(sprintf("          [CACHE] Goedgekeurde oevers (>= %dm) ingeladen...", as.integer(min_oever_l)))
        r_oever_ok <- terra::rast(f_cache_oevers)
      } else {
        f_wv_al <- hulp_align_gdal(paths$watervlakken_raster, master_grid)
        r_wv    <- terra::rast(f_wv_al)
        
        r_w <- terra::ifel(!is.na(r_wv) & r_wv > 0, 1, 0)
        r_l <- terra::ifel(is.na(r_wv) | r_wv == 0, 1, 0)
        unlink(f_wv_al)
        rm(r_wv)
        
        f_b1 <- tempfile(fileext = ".tif")
        f_b2 <- tempfile(fileext = ".tif")
        bnd_w <- terra::boundaries(r_w, classes = FALSE, inner = TRUE, filename = f_b1, overwrite = TRUE)
        bnd_l <- terra::boundaries(r_l, classes = FALSE, inner = TRUE, filename = f_b2, overwrite = TRUE)
        
        r_oev_raw <- terra::ifel((bnd_w == 1 | bnd_l == 1), 1, NA)
        unlink(c(f_b1, f_b2))
        rm(r_w, r_l, bnd_w, bnd_l)
        
        f_oev_bin <- tempfile(fileext = ".tif")
        terra::writeRaster(r_oev_raw, f_oev_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 0)
        f_prox_oev <- hulp_proximity_gdal(f_oev_bin, 25)
        unlink(f_oev_bin)
        
        p_oev <- terra::rast(f_prox_oev)
        mask_oev25 <- terra::ifel(!is.na(p_oev) & p_oev >= 0 & p_oev <= 25, 1, NA)
        unlink(f_prox_oev)
        rm(p_oev)
        
        cl_oevers   <- clump_fast(mask_oev25, diag = FALSE, titel = "cl_oevers")
        rm(mask_oev25)
        cl_oev_only <- terra::mask(cl_oevers, r_oev_raw)
        unlink(attr(cl_oevers, "tmp_file"))
        rm(r_oev_raw, cl_oevers)
        
        min_pixels_oever <- ceiling(min_oever_l / 10)
        dt_oev <- data.table::as.data.table(terra::freq(cl_oev_only))
        setnames(dt_oev, c("layer", "ID", "count"))
        valide_oev_ids <- dt_oev[!is.na(ID) & ID > 0 & count >= min_pixels_oever, ID]
        rm(dt_oev)
        
        r_oever_ok <- if (length(valide_oev_ids) > 0) terra::ifel(cl_oev_only %in% as.integer(valide_oev_ids), 1, NA) else (master_grid * NA_real_)
        rm(cl_oev_only)
        terra::writeRaster(r_oever_ok, f_cache_oevers, overwrite = TRUE, datatype = "INT1U", gdal = "COMPRESS=LZW")
      }
      
      if (!is_leeg(r_oever_ok)) {
        f_oev_bin <- tempfile(fileext = ".tif")
        terra::writeRaster(r_oever_ok, f_oev_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 0)
        f_prox_20m <- hulp_proximity_gdal(f_oev_bin, oever_zoom_m)
        unlink(f_oev_bin)
        
        p_20m <- terra::rast(f_prox_20m)
        mask_zone20 <- terra::ifel(!is.na(p_20m) & p_20m >= 0 & p_20m <= oever_zoom_m, 1, NA)
        unlink(f_prox_20m)
        rm(p_20m)
        
        r_huidig <- terra::mask(r_huidig, mask_zone20)
        rm(mask_zone20)
      } else {
        r_huidig <- master_grid * NA_real_
      }
      rm(r_oever_ok)
      message(sprintf("       📊 Na oeverlijntoets (>= %dm) en %dm oeverstrook: %.2f ha", 
                      as.integer(min_oever_l), as.integer(oever_zoom_m), calc_ha_exact(r_huidig)))
    }
    
    # 2D. LGP Gewassen
    lgp_codes <- get_param(rij, "LGP_Gewascodes", NA)
    if (!is.na(lgp_codes) && !lgp_als_basis) {
      lgp_str <- trimws(as.character(lgp_codes))
      is_additief <- grepl("^\\+", lgp_str)
      if (is_additief) lgp_str <- sub("^\\+\\s*", "", lgp_str)
      
      codes_v <- strsplit(lgp_str, "[,; ]+")[[1]]
      message(sprintf("       -> LGP-landbouwgewassen %s (%s)...", 
                      if (is_additief) "additief toevoegen" else "filteren als doorsnede", 
                      paste(codes_v, collapse = ", ")))
      
      m_lgp <- haal_categorisch_masker(paths$lgp_landbouw, codes_v, master_grid, paths$dir_cache, "lgp")
      
      if (is_additief) {
        r_huidig <- voeg_toe_geclamped(r_huidig, m_lgp)
      } else {
        r_huidig <- terra::mask(r_huidig, m_lgp)
      }
      message(sprintf("       📊 Na LGP-bewerking: %.2f ha", calc_ha_exact(r_huidig)))
    }
    
    # 2E-G. Bodem- en overstromingscriteria
    heeft_drain   <- isTRUE(as.logical(get_param(rij, "Gebruik_Drainage", FALSE))) && !is.na(get_param(rij, "Drainage_Codes", NA))
    heeft_ovstr   <- isTRUE(as.logical(get_param(rij, "Gebruik_Overstroming", FALSE))) && !is.na(get_param(rij, "Overstroming_Codes", NA))
    heeft_prof    <- !is.na(get_param(rij, "Profiel_Codes", NA))
    
    ovstr_val     <- trimws(as.character(get_param(rij, "Overstroming_Codes", "")))
    is_droog_arch <- (ovstr_val == "0")
    is_wetland    <- grepl("wetland|actwet", tolower(ovstr_val))
    
    if (is_droog_arch) {
      if (heeft_ovstr && !is.null(paths$overstroming) && file.exists(paths$overstroming)) {
        message("       -> Inundatiegevoelige gronden uitsluiten (enkel niet-overstroombaar klasse 0)...")
        f_ov_al  <- hulp_align_gdal(paths$overstroming, master_grid)
        r_ov     <- terra::rast(f_ov_al)
        m_droog  <- terra::ifel(is.na(r_ov) | r_ov == 0, 1, NA)
        unlink(f_ov_al)
        r_huidig <- terra::mask(r_huidig, m_droog)
        rm(r_ov, m_droog)
      }
      
      if (heeft_prof) {
        codes_p   <- strsplit(as.character(get_param(rij, "Profiel_Codes", NA)), "[,; ]+")[[1]]
        message(sprintf("       -> Filteren op bodemtextuur/leem (klasse %s)...", paste(codes_p, collapse = ", ")))
        pad_bodem <- if (!is.null(paths$textuur) && file.exists(paths$textuur)) paths$textuur else paths$profiel
        m_leem    <- haal_categorisch_masker(pad_bodem, codes_p, master_grid, paths$dir_cache, "textuur")
        r_huidig  <- terra::mask(r_huidig, m_leem)
        rm(m_leem)
      }
      message(sprintf("       📊 Na bodem- en overstromingsfilter: %.2f ha", calc_ha_exact(r_huidig)))
      
    } else if (is_wetland || heeft_drain || heeft_ovstr || heeft_prof) {
      m_vallei <- master_grid * NA_real_
      
      if (is_wetland) {
        heeft_wetland_file <- !is.null(paths$wetland) && file.exists(paths$wetland)
        
        if (heeft_wetland_file) {
          message("       -> Alluviaal wetlandmasker (actwet.tif) toepassen...")
          f_wet_al  <- hulp_align_gdal(paths$wetland, master_grid)
          r_wet_raw <- terra::rast(f_wet_al)
          m_wet     <- terra::ifel(!is.na(r_wet_raw) & (r_wet_raw %in% c(50:81, 120) | r_wet_raw > 0), 1, NA)
          unlink(f_wet_al)
          rm(r_wet_raw)
          m_vallei  <- terra::cover(m_vallei, m_wet)
          rm(m_wet)
        } else {
          message("       -> Geen 'actwet.tif' gevonden; fallback via 100m waterloopbuffer en natte drainage...")
          if (!is.null(paths$huetzon) && file.exists(paths$huetzon)) {
            f_hz_al  <- hulp_align_gdal(paths$huetzon, master_grid)
            r_hz     <- terra::rast(f_hz_al)
            r_hz_bin <- terra::ifel(!is.na(r_hz) & r_hz > 0, 1, NA)
            unlink(f_hz_al)
            rm(r_hz)
            
            if (!is_leeg(r_hz_bin)) {
              f_hz_tmp  <- tempfile(fileext = ".tif")
              terra::writeRaster(r_hz_bin, f_hz_tmp, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
              f_prox_hz <- hulp_proximity_gdal(f_hz_tmp, 100)
              unlink(f_hz_tmp)
              p_hz      <- terra::rast(f_prox_hz)
              m_w100    <- terra::ifel(!is.na(p_hz) & p_hz >= 0 & p_hz <= 100, 1, NA)
              unlink(f_prox_hz)
              rm(p_hz)
              m_vallei  <- terra::cover(m_vallei, m_w100)
              rm(m_w100)
            }
            rm(r_hz_bin)
          }
          
          if (!is.null(paths$drainage) && file.exists(paths$drainage)) {
            drain_codes_fallback <- c("e", "f", "g", "h", "i", "e-f", "h-i", "e-i")
            m_d_fb   <- haal_categorisch_masker(paths$drainage, drain_codes_fallback, master_grid, paths$dir_cache, "drainage")
            m_vallei <- terra::cover(m_vallei, m_d_fb)
            rm(m_d_fb)
          }
        }
      }
      
      if (heeft_drain && !is_wetland) {
        drain_codes <- get_param(rij, "Drainage_Codes", NA)
        drain_str   <- trimws(as.character(drain_codes))
        is_neg      <- grepl("^!", drain_str)
        if (is_neg) drain_str <- sub("^!\\s*", "", drain_str)
        codes_d     <- strsplit(drain_str, "[,; ]+")[[1]]
        m_d         <- haal_categorisch_masker(paths$drainage, codes_d, master_grid, paths$dir_cache, "drainage")
        if (is_neg) m_d <- terra::ifel(is.na(m_d), 1, NA)
        m_vallei    <- terra::cover(m_vallei, m_d)
        rm(m_d)
      }
      
      if (heeft_ovstr && !is_wetland) {
        codes_o  <- strsplit(as.character(get_param(rij, "Overstroming_Codes", NA)), "[,; ]+")[[1]]
        m_o      <- haal_categorisch_masker(paths$overstroming, codes_o, master_grid, paths$dir_cache, "overstroming")
        m_vallei <- terra::cover(m_vallei, m_o)
        rm(m_o)
      }
      
      if (heeft_prof) {
        codes_p  <- strsplit(as.character(get_param(rij, "Profiel_Codes", NA)), "[,; ]+")[[1]]
        m_p      <- haal_categorisch_masker(paths$profiel, codes_p, master_grid, paths$dir_cache, "profiel")
        m_vallei <- terra::cover(m_vallei, m_p)
        rm(m_p)
      }
      
      message("       -> Gecombineerde alluviale/valleitest toepassen (Wetland OR Drainage OR Overstroming OR Profiel)...")
      r_huidig <- terra::mask(r_huidig, m_vallei)
      rm(m_vallei)
      message(sprintf("       📊 Na valleitest: %.2f ha", calc_ha_exact(r_huidig)))
    }
    
    # 2H. Morfologische opening
    min_straal <- as.numeric(get_param(rij, "Min_Breedte_Straal_m", 0))
    if (!is.na(min_straal) && min_straal > 0) {
      breedte_m <- min_straal * 2
      message(sprintf("       -> Morfologische opening (min. breedte %dm)...", as.integer(breedte_m)))
      r_huidig <- hulp_morfologische_opening(r_huidig, diameter_m = breedte_m, master_grid = master_grid)
    }
    
    # 2I. Groenkaart integratie & Verbossingsfilter
    bos_drempel <- NA_real_
    bos_buf_m   <- 0
    
    if (isTRUE(as.logical(get_param(rij, "Gebruik_Groenkaart", FALSE)))) {
      bos_drempel <- as.numeric(get_param(rij, "Groenkaart_Drempel_m2", NA))
      bos_buf_m   <- as.numeric(dplyr::coalesce(
        as.numeric(get_param(rij, "Randzone_Diepte_m", NA)),
        as.numeric(get_param(rij, "Grens_Afstand_m", NA)),
        0
      ))
      
      # TAK A: Verbossingsfilter op aaneengesloten clusters
      if (!is.na(bos_drempel) && bos_drempel > 0) {
        message(sprintf("       -> Verbossingsfilter: Boscomplexen >= %d m² uitsluiten (buffer: %dm)...", 
                        as.integer(bos_drempel), as.integer(bos_buf_m)))
        f_cache_bos <- file.path(paths$dir_cache, sprintf("masker_bos_%dm2_%dm.tif", as.integer(bos_drempel), as.integer(bos_buf_m)))
        
        if (file.exists(f_cache_bos) && !force_overwrite) {
          m_bos_stoor <- terra::rast(f_cache_bos)
        } else {
          f_gr_al <- hulp_align_gdal(paths$groenkaart, master_grid)
          r_gr    <- terra::rast(f_gr_al)
          r_b1    <- terra::ifel(!is.na(r_gr) & r_gr == 1, 1, NA)
          unlink(f_gr_al)
          rm(r_gr)
          
          if ("bwkbos" %in% names(lagen_ruw)) {
            r_bwkbos <- lagen_ruw[["bwkbos"]]
            r_bwkbos_buiten <- terra::ifel(!is.na(r_huidig), NA, r_bwkbos)
            r_bwkbos_bin    <- terra::ifel(!is.na(r_bwkbos_buiten) & r_bwkbos_buiten > 0, 1, NA)
            r_b1 <- terra::cover(r_bwkbos_bin, r_b1)
            rm(r_bwkbos, r_bwkbos_buiten, r_bwkbos_bin)
          }
          
          cl_b <- clump_fast(r_b1, diag = FALSE, titel = "cl_bos_groot")
          rm(r_b1)
          
          min_cellen_bos <- ceiling(bos_drempel / 100)
          dt_b <- data.table::as.data.table(terra::freq(cl_b))
          setnames(dt_b, c("layer", "ID", "count"))
          grote_ids <- dt_b[!is.na(ID) & ID > 0 & count >= min_cellen_bos, ID]
          rm(dt_b)
          
          if (length(grote_ids) > 0) {
            r_bos_groot <- terra::ifel(cl_b %in% as.integer(grote_ids), 1, NA)
            if (bos_buf_m > 0) {
              f_bg_bin <- tempfile(fileext = ".tif")
              terra::writeRaster(r_bos_groot, f_bg_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
              f_prox_bg <- hulp_proximity_gdal(f_bg_bin, bos_buf_m)
              unlink(f_bg_bin)
              p_bg <- terra::rast(f_prox_bg)
              m_bos_stoor <- terra::ifel(!is.na(p_bg) & p_bg >= 0 & p_bg <= bos_buf_m, 1, NA)
              unlink(f_prox_bg)
              rm(p_bg, r_bos_groot)
            } else {
              m_bos_stoor <- r_bos_groot
            }
          } else {
            m_bos_stoor <- master_grid * NA_real_
          }
          unlink(attr(cl_b, "tmp_file"))
          rm(cl_b)
          terra::writeRaster(m_bos_stoor, f_cache_bos, overwrite = TRUE, datatype = "INT1U", gdal = "COMPRESS=LZW")
        }
        
        r_huidig <- terra::mask(r_huidig, m_bos_stoor, inverse = TRUE)
        message(sprintf("       📊 Na verbossingsfilter (>= %d m² + %dm buffer): %.2f ha", 
                        as.integer(bos_drempel), as.integer(bos_buf_m), calc_ha_exact(r_huidig)))
        
        # TAK B: Directe hooggroenuitsluiting
      } else if (!is.na(bos_drempel) && bos_drempel == 0) {
        message("       -> Directe hooggroenuitsluiting (alle bomen klasse 1 uitsluiten)...")
        f_gr_al <- hulp_align_gdal(paths$groenkaart, master_grid)
        r_gr    <- terra::rast(f_gr_al)
        m_bomen <- terra::ifel(!is.na(r_gr) & r_gr == 1, 1, NA)
        unlink(f_gr_al)
        
        if (bos_buf_m > 0 && !is_leeg(m_bomen)) {
          f_bm_bin  <- tempfile(fileext = ".tif")
          terra::writeRaster(m_bomen, f_bm_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
          f_prox_bm <- hulp_proximity_gdal(f_bm_bin, bos_buf_m)
          unlink(f_bm_bin)
          p_bm <- terra::rast(f_prox_bm)
          m_bomen <- terra::ifel(!is.na(p_bm) & p_bm >= 0 & p_bm <= bos_buf_m, 1, NA)
          unlink(f_prox_bm)
          rm(p_bm)
        }
        
        r_huidig <- terra::mask(r_huidig, m_bomen, inverse = TRUE)
        rm(r_gr, m_bomen)
        message(sprintf("       📊 Na directe hooggroenuitsluiting: %.2f ha", calc_ha_exact(r_huidig)))
        
        # Gladde slang-specifiek: aangrenzend hooggroen binnen 50 m van open biotoop toevoegen als bosrand
        if (tolower(bestandsnaam) %in% c("gladdeslang", "coronellaaustriaca") &&
            !is_leeg(r_huidig)) {
          message("       -> Gladde slang-specifiek: aangrenzend hooggroen binnen 50 m als bosrand toevoegen...")
          f_gr_al2  <- hulp_align_gdal(paths$groenkaart, master_grid)
          r_gr2     <- terra::rast(f_gr_al2)
          r_bomen50 <- terra::ifel(!is.na(r_gr2) & r_gr2 == 1, 1.0, NA_real_)
          unlink(f_gr_al2)
          rm(r_gr2)
          
          f_open_bin <- tempfile(fileext = ".tif")
          terra::writeRaster(terra::ifel(!is.na(r_huidig) & r_huidig > 0, 1, NA), f_open_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
          f_prox_50m <- hulp_proximity_gdal(f_open_bin, 50)
          unlink(f_open_bin)
          
          p_50m <- terra::rast(f_prox_50m)
          m_50m <- terra::ifel(!is.na(p_50m) & p_50m >= 0 & p_50m <= 50, 1, NA)
          unlink(f_prox_50m)
          rm(p_50m)
          
          r_bosrand50 <- terra::mask(r_bomen50, m_50m)
          r_huidig    <- voeg_toe_geclamped(r_huidig, r_bosrand50)
          rm(r_bomen50, m_50m, r_bosrand50)
          message(sprintf("       📊 Na toevoegen 50m Groenkaart-bosrand: %.2f ha", calc_ha_exact(r_huidig)))
        }
        
        # TAK C: Additief toevoegen (-1 = alle bomen, -2 = bomenrijen/KLE)
      } else if (!is.na(bos_drempel) && bos_drempel < 0) {
        f_gr_al <- hulp_align_gdal(paths$groenkaart, master_grid)
        r_gr    <- terra::rast(f_gr_al)
        
        if (bos_drempel == -2) {
          message("       -> Groenkaart bomenrijen/KLE (klasse 1 binnen 10m van open groen) additief toevoegen...")
          f_cache_kle <- file.path(paths$dir_cache, "masker_groenkaart_kle_bomen.tif")
          
          if (file.exists(f_cache_kle) && !force_overwrite) {
            m_bomen_add <- terra::rast(f_cache_kle)
          } else {
            r_b1 <- terra::ifel(!is.na(r_gr) & r_gr == 1, 1, NA)
            r_op <- terra::ifel(!is.na(r_gr) & (r_gr %in% c(2, 3, 4)), 1, NA)
            
            if (!is_leeg(r_b1)) {
              f_op_bin  <- tempfile(fileext = ".tif")
              terra::writeRaster(r_op, f_op_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
              f_prox_op <- hulp_proximity_gdal(f_op_bin, 10)
              unlink(f_op_bin)
              p_op <- terra::rast(f_prox_op)
              m_op10 <- terra::ifel(!is.na(p_op) & p_op >= 0 & p_op <= 10, 1, NA)
              unlink(f_prox_op)
              rm(p_op)
              
              m_bomen_add <- terra::mask(r_b1, m_op10)
              rm(m_op10, r_b1, r_op)
            } else {
              m_bomen_add <- master_grid * NA_real_
              rm(r_b1, r_op)
            }
            terra::writeRaster(m_bomen_add, f_cache_kle, overwrite = TRUE, datatype = "INT1U", gdal = "COMPRESS=LZW")
          }
        } else {
          message("       -> Groenkaart klasse 1 (alle bomen/ooibos) additief toevoegen aan substraat...")
          m_bomen_add <- terra::ifel(!is.na(r_gr) & r_gr == 1, 1.0, NA_real_)
        }
        
        unlink(f_gr_al)
        rm(r_gr)
        r_huidig <- voeg_toe_geclamped(r_huidig, m_bomen_add)
        rm(m_bomen_add)
        message(sprintf("       📊 Na additieve Groenkaart-toevoeging: %.2f ha", calc_ha_exact(r_huidig)))
        
        # TAK D: Positieve verificatie
      } else {
        message("       -> Positieve verificatie met Groenkaart (strikt klasse 1: opgaand bos)...")
        f_gr_al <- hulp_align_gdal(paths$groenkaart, master_grid)
        r_gr    <- terra::rast(f_gr_al)
        m_bos1  <- terra::ifel(!is.na(r_gr) & r_gr == 1, 1.0, NA)
        unlink(f_gr_al)
        r_huidig <- terra::mask(r_huidig, m_bos1)
        rm(r_gr, m_bos1)
      }
    }
    
    # 2J. Negatieve Bebouwingsbuffer
    buf_bebouw_m <- as.numeric(get_param(rij, "Bebouwing_Buffer_m", 0))
    if (isTRUE(as.logical(get_param(rij, "Gebruik_Bebouwing", FALSE))) && buf_bebouw_m > 0) {
      bebouw_min_ha  <- as.numeric(dplyr::coalesce(as.numeric(get_param(rij, "Bebouwing_MinOpp_ha", NA)), 1.5))
      m_bebouw_stoor <- maak_bebouwingsmasker(
        min_opp_ha      = bebouw_min_ha,
        buf_m           = buf_bebouw_m,
        afstand_cl_m    = 50,
        ds_parquet      = ds_parquet,
        master_grid     = master_grid,
        paths           = paths,
        force_overwrite = force_overwrite
      )
      r_huidig       <- terra::mask(r_huidig, m_bebouw_stoor, inverse = TRUE)
      rm(m_bebouw_stoor)
      message(sprintf("       📊 Na bebouwingsbuffer: %.2f ha", calc_ha_exact(r_huidig)))
    }
    
    # 2K. Negatieve Wegenbuffer
    buf_wegen_m <- as.numeric(get_param(rij, "Wegen_Buffer_m", 0))
    if (isTRUE(as.logical(get_param(rij, "Gebruik_Wegen", FALSE))) && buf_wegen_m > 0 && 
        !is.null(paths$grb_wegen) && file.exists(paths$grb_wegen)) {
      f_cache_wegen <- file.path(paths$dir_cache, sprintf("masker_wegen_stoor_%dm.tif", as.integer(buf_wegen_m)))
      
      is_valid_cache_w <- FALSE
      if (file.exists(f_cache_wegen) && !force_overwrite) {
        test_rw <- suppressWarnings(terra::rast(f_cache_wegen))
        if (terra::compareGeom(test_rw, master_grid, stopOnError = FALSE)) {
          is_valid_cache_w <- TRUE
          m_weg_stoor <- test_rw
          message(sprintf("       -> [CACHE] Wegenbuffer (%dm) ingeladen...", as.integer(buf_wegen_m)))
        }
      }
      
      if (!is_valid_cache_w) {
        message(sprintf("       -> Wegenbuffer genereren (snelwegen / gescheiden rijbanen, %dm)...", as.integer(buf_wegen_m)))
        f_w_al <- hulp_align_gdal(paths$grb_wegen, master_grid)
        r_w    <- terra::rast(f_w_al)
        w_cats <- terra::cats(r_w)[[1]]
        ids_w  <- integer(0)
        if (!is.null(w_cats)) {
          kol_l <- colnames(w_cats)[tolower(colnames(w_cats)) %in% c("label", "omschrijving")][1]
          kol_v <- colnames(w_cats)[tolower(colnames(w_cats)) %in% c("value", "id")][1]
          if (is.na(kol_v)) kol_v <- 1
          if (!is.na(kol_l)) {
            ids_w <- w_cats[[kol_v]][grepl("autosnelweg|gescheiden rijbanen", tolower(w_cats[[kol_l]]))]
          }
        }
        r_w_bin <- terra::ifel(!is.na(r_w) & r_w %in% ids_w, 1, NA)
        unlink(f_w_al)
        f_wb_bin <- tempfile(fileext = ".tif")
        terra::writeRaster(r_w_bin, f_wb_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
        f_prox_w <- hulp_proximity_gdal(f_wb_bin, buf_wegen_m)
        unlink(f_wb_bin)
        p_w <- terra::rast(f_prox_w)
        m_weg_stoor <- terra::ifel(!is.na(p_w) & p_w >= 0 & p_w <= buf_wegen_m, 1, NA)
        unlink(f_prox_w)
        terra::writeRaster(m_weg_stoor, f_cache_wegen, overwrite = TRUE, datatype = "INT1U", gdal = "COMPRESS=LZW")
      }
      r_huidig <- terra::mask(r_huidig, m_weg_stoor, inverse = TRUE)
      message(sprintf("       📊 Na wegenbuffer: %.2f ha", calc_ha_exact(r_huidig)))
    }
    
    # 2-vliegendhert: Saproxyle hellingbossen
    is_vliegendhert <- tolower(bestandsnaam) %in% c("vliegendhert", "lucanuscervus") || 
      (!is.null(paths$helling_4gr) && file.exists(paths$helling_4gr) && !is.na(get_param(rij, "Randzone_Diepte_m", NA)))
    
    if (is_vliegendhert && !is_leeg(r_huidig)) {
      message("       -> Saproxyl hellingbos-archetype uitvoeren (Oud bos 1km x Helling 1km x 15m rand / 40m diepte)...")
      
      f_gr_al <- hulp_align_gdal(paths$groenkaart, master_grid)
      r_gr    <- terra::rast(f_gr_al)
      r_gr_hoog <- terra::ifel(!is.na(r_gr) & r_gr == 1, 1, NA)
      r_gr_open <- terra::ifel(!is.na(r_gr) & (r_gr %in% c(2, 3, 4)), 1, NA)
      unlink(f_gr_al)
      rm(r_gr)
      
      r_naald_raw <- if ("naaldbos" %in% names(lagen_ruw)) lagen_ruw[["naaldbos"]] else (master_grid * 0)
      r_naald_cl  <- terra::ifel(is.na(r_naald_raw), 0, r_naald_raw)
      r_loof_frac <- terra::clamp(1.0 - r_naald_cl, lower = 0.0, upper = 1.0)
      
      r_loof_dom  <- terra::ifel(!is.na(r_gr_hoog) & r_naald_cl <= 0.40, r_loof_frac, NA)
      rm(r_naald_cl, r_gr_hoog, r_naald_raw)
      
      if (!is.null(paths$bosleeftijd) && file.exists(paths$bosleeftijd)) {
        f_cache_ob <- file.path(paths$dir_cache, "masker_oud_bos_1000m.tif")
        if (file.exists(f_cache_ob) && !force_overwrite) {
          m_ob1km <- terra::rast(f_cache_ob)
        } else {
          f_ob_al  <- hulp_align_gdal(paths$bosleeftijd, master_grid)
          r_ob_raw <- terra::rast(f_ob_al)
          r_ob_bin <- terra::ifel(!is.na(r_ob_raw) & (r_ob_raw %in% c(1, 2)), 1, NA)
          unlink(f_ob_al)
          rm(r_ob_raw)
          
          if (!is_leeg(r_ob_bin)) {
            f_ob_tmp <- tempfile(fileext = ".tif")
            terra::writeRaster(r_ob_bin, f_ob_tmp, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
            f_prox_ob <- hulp_proximity_gdal(f_ob_tmp, 1000)
            unlink(f_ob_tmp)
            p_ob <- terra::rast(f_prox_ob)
            m_ob1km <- terra::ifel(!is.na(p_ob) & p_ob >= 0 & p_ob <= 1000, 1, NA)
            unlink(f_prox_ob)
            rm(p_ob)
          } else {
            m_ob1km <- master_grid * 1
          }
          rm(r_ob_bin)
          terra::writeRaster(m_ob1km, f_cache_ob, overwrite = TRUE, datatype = "INT1U", gdal = "COMPRESS=LZW")
        }
        r_loof_dom <- terra::mask(r_loof_dom, m_ob1km)
        rm(m_ob1km)
      }
      
      if (!is.null(paths$helling_4gr) && file.exists(paths$helling_4gr)) {
        f_cache_hel <- file.path(paths$dir_cache, "masker_helling_4gr_1000m.tif")
        if (file.exists(f_cache_hel) && !force_overwrite) {
          m_hel1km <- terra::rast(f_cache_hel)
        } else {
          f_hel_al  <- hulp_align_gdal(paths$helling_4gr, master_grid)
          r_hel_raw <- terra::rast(f_hel_al)
          r_hel_bin <- terra::ifel(!is.na(r_hel_raw) & r_hel_raw == 1, 1, NA)
          unlink(f_hel_al)
          rm(r_hel_raw)
          
          if (!is_leeg(r_hel_bin)) {
            f_hel_tmp <- tempfile(fileext = ".tif")
            terra::writeRaster(r_hel_bin, f_hel_tmp, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
            f_prox_hel <- hulp_proximity_gdal(f_hel_tmp, 1000)
            unlink(f_hel_tmp)
            p_hel <- terra::rast(f_prox_hel)
            m_hel1km <- terra::ifel(!is.na(p_hel) & p_hel >= 0 & p_hel <= 1000, 1, NA)
            unlink(f_prox_hel)
            rm(p_hel)
          } else {
            m_hel1km <- master_grid * 1
          }
          rm(r_hel_bin)
          terra::writeRaster(m_hel1km, f_cache_hel, overwrite = TRUE, datatype = "INT1U", gdal = "COMPRESS=LZW")
        }
        r_loof_dom <- terra::mask(r_loof_dom, m_hel1km)
        rm(m_hel1km)
      }
      
      min_bos_c_ha <- as.numeric(dplyr::coalesce(as.numeric(get_param(rij, "Kern_MinOpp_ha", NA)), 10.0))
      res_b10 <- cluster_filter_universeel(
        masker = r_loof_dom, opp_laag = r_loof_dom, min_opp_ha = min_bos_c_ha,
        min_totale_opp_ha = NA, afstand_m = 50, dispersie_m = 0, werkelijk = FALSE
      )
      vh_bos_10ha <- res_b10$raster
      rm(res_b10, r_loof_dom)
      
      res_op5 <- cluster_filter_universeel(
        masker = r_gr_open, opp_laag = r_gr_open, min_opp_ha = 5.0,
        min_totale_opp_ha = NA, afstand_m = 50, dispersie_m = 0, werkelijk = FALSE
      )
      vh_open_5ha <- res_op5$raster
      rm(res_op5, r_gr_open)
      
      if (!is_leeg(vh_bos_10ha) && 
          !is_leeg(vh_open_5ha)) {
        f_op_bin <- tempfile(fileext = ".tif")
        terra::writeRaster(terra::ifel(!is.na(vh_open_5ha) & vh_open_5ha > 0, 1, NA), f_op_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
        f_prox_op <- hulp_proximity_gdal(f_op_bin, 15)
        unlink(f_op_bin)
        p_op <- terra::rast(f_prox_op)
        m_op15 <- terra::ifel(!is.na(p_op) & p_op >= 0 & p_op <= 15, 1, NA)
        unlink(f_prox_op)
        rm(p_op)
        
        r_rand_raw <- terra::mask(vh_bos_10ha, m_op15)
        rm(m_op15)
        
        res_r025 <- cluster_filter_universeel(
          masker = r_rand_raw, opp_laag = r_rand_raw, min_opp_ha = 0.25,
          min_totale_opp_ha = NA, afstand_m = 10, dispersie_m = 0, werkelijk = FALSE
        )
        vh_rand_ok <- res_r025$raster
        rm(res_r025, r_rand_raw)
        
        if (!is_leeg(vh_rand_ok)) {
          f_rd_bin <- tempfile(fileext = ".tif")
          terra::writeRaster(terra::ifel(!is.na(vh_rand_ok) & vh_rand_ok > 0, 1, NA), f_rd_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
          f_prox_rd <- hulp_proximity_gdal(f_rd_bin, 40)
          unlink(f_rd_bin)
          p_rd <- terra::rast(f_prox_rd)
          m_rd40 <- terra::ifel(!is.na(p_rd) & p_rd >= 0 & p_rd <= 40, 1, NA)
          unlink(f_prox_rd)
          rm(p_rd)
          
          vh_bos_diep <- terra::mask(vh_bos_10ha, m_rd40)
          rm(m_rd40)
          
          vh_bruto <- terra::cover(vh_rand_ok, vh_bos_diep)
          r_huidig <- terra::mask(r_loof_frac, vh_bruto)
          rm(vh_bos_diep, vh_bruto)
        } else {
          r_huidig <- master_grid * NA_real_
        }
        rm(vh_rand_ok)
      } else {
        r_huidig <- master_grid * NA_real_
      }
      rm(vh_bos_10ha, vh_open_5ha, r_loof_frac)
      message(sprintf("       📊 Na saproxyle hellingbos-analyse: %.2f ha", calc_ha_exact(r_huidig)))
    }
    
    # 2L-pre. Richtinggevende voorfilter
    grens_moment_val <- tolower(trimws(get_param(rij, "Grens_Moment", "")))
    if (grens_moment_val == "voor_clustering") {
      ref_tp_raw <- get_param(rij, "Grens_Met_Type", NA)
      grens_m_raw <- get_param(rij, "Grens_Afstand_m", 0)
      
      if (!is.na(ref_tp_raw)) {
        ref_types <- trimws(strsplit(as.character(ref_tp_raw), "[,;]+")[[1]])
        afstanden <- as.numeric(trimws(strsplit(as.character(grens_m_raw), "[,;]+")[[1]]))
        if (length(afstanden) < length(ref_types)) afstanden <- rep(afstanden[1], length(ref_types))
        
        for (k in seq_along(ref_types)) {
          ref_tp  <- ref_types[k]
          grens_m <- afstanden[k]
          r_ref_bron <- if (ref_tp %in% names(verwerkte_sporen)) verwerkte_sporen[[ref_tp]] else (if (ref_tp %in% names(lagen_ruw)) lagen_ruw[[ref_tp]] else NULL)
          
          if (!is.null(r_ref_bron) && grens_m > 0) {
            message(sprintf("       -> Voorfilter: '%s' beperken tot binnen %dm van '%s'...", tp_naam, as.integer(grens_m), ref_tp))
            r_ref_bin <- terra::ifel(!is.na(r_ref_bron) & r_ref_bron > 0, 1, NA)
            if (!is_leeg(r_ref_bin)) {
              f_ref_bin <- tempfile(fileext = ".tif")
              terra::writeRaster(r_ref_bin, f_ref_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
              f_prox_ref <- hulp_proximity_gdal(f_ref_bin, grens_m)
              unlink(f_ref_bin)
              
              p_ref <- terra::rast(f_prox_ref)
              m_ref <- terra::ifel(!is.na(p_ref) & p_ref >= 0 & p_ref <= grens_m, 1, NA)
              unlink(f_prox_ref)
              r_huidig <- terra::mask(r_huidig, m_ref)
              rm(p_ref, m_ref)
            } else {
              r_huidig <- master_grid * NA_real_
            }
            message(sprintf("       📊 Na voorfilter (%s binnen %dm %s): %.2f ha", tp_naam, as.integer(grens_m), ref_tp, calc_ha_exact(r_huidig)))
          }
        }
      }
    }
    
    # 2-rietmoeras: 3-schillen ecotoontoets (Roerdomp)
    rij_openw <- params_gesorteerd %>% dplyr::filter(tolower(Type) == "openwater") %>% dplyr::slice(1)
    rij_riet  <- params_gesorteerd %>% dplyr::filter(tolower(Type) == "rietland") %>% dplyr::slice(1)
    
    is_rietmoeras_complex <- nrow(rij_openw) > 0 && nrow(rij_riet) > 0 && 
      !is.na(get_param(rij_openw, "Kern_MinOpp_ha", NA)) && 
      !is.na(get_param(rij_riet, "Kern_MinOpp_ha", NA)) && 
      !is.na(get_param(rij, "Water_Afstand_m", NA))
    
    if (is_rietmoeras_complex && "openwater" %in% names(lagen_ruw) && "rietland" %in% names(lagen_ruw)) {
      message("       -> 3-schillen rietmoerasanalyse uitvoeren (water >= 10ha, rietkern >= 0.5ha, ecotoon >= 0.5ha)...")
      
      drempel_w_ha  <- as.numeric(get_param(rij_openw, "Kern_MinOpp_ha", 10.0))
      drempel_rk_ha <- as.numeric(get_param(rij_riet, "Kern_MinOpp_ha", 0.5))
      schil_buf_m   <- as.numeric(get_param(rij, "Water_Afstand_m", 100))
      
      r_w_raw <- lagen_ruw[["openwater"]]
      res_w10 <- cluster_filter_universeel(
        masker = r_w_raw, opp_laag = r_w_raw, min_opp_ha = drempel_w_ha,
        min_totale_opp_ha = NA, afstand_m = 50, dispersie_m = 0, werkelijk = TRUE
      )
      r_w10 <- res_w10$raster
      rm(res_w10)
      
      if (!is_leeg(r_w10)) {
        f_w_bin <- tempfile(fileext = ".tif")
        terra::writeRaster(terra::ifel(!is.na(r_w10) & r_w10 > 0, 1, NA), f_w_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
        f_prox_w <- hulp_proximity_gdal(f_w_bin, schil_buf_m)
        unlink(f_w_bin)
        p_w <- terra::rast(f_prox_w)
        m_w100 <- terra::ifel(!is.na(p_w) & p_w >= 0 & p_w <= schil_buf_m, 1, NA)
        unlink(f_prox_w)
        rm(p_w)
        r_huidig <- terra::mask(r_huidig, m_w100)
        rm(m_w100)
      } else {
        r_huidig <- master_grid * NA_real_
      }
      rm(r_w10)
      
      r_r_raw <- lagen_ruw[["rietland"]]
      r_r_bin <- terra::ifel(!is.na(r_r_raw) & r_r_raw > 0, 1, NA)
      
      if (!is_leeg(r_r_bin)) {
        f_r_bin <- tempfile(fileext = ".tif")
        terra::writeRaster(r_r_bin, f_r_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
        f_prox_r10 <- hulp_proximity_gdal(f_r_bin, 10)
        unlink(f_r_bin)
        p_r10 <- terra::rast(f_prox_r10)
        m_r10 <- terra::ifel(!is.na(p_r10) & p_r10 >= 0 & p_r10 <= 10, 1, NA)
        unlink(f_prox_r10)
        rm(p_r10)
        
        r_wk <- terra::mask(r_w_raw, m_r10)
        rm(m_r10)
        
        res_wk <- cluster_filter_universeel(
          masker = r_wk, opp_laag = r_wk, min_opp_ha = drempel_rk_ha,
          min_totale_opp_ha = NA, afstand_m = 50, dispersie_m = 0, werkelijk = TRUE
        )
        r_wk_ok <- res_wk$raster
        rm(res_wk, r_wk)
        
        if (!is_leeg(r_wk_ok)) {
          f_wk_bin <- tempfile(fileext = ".tif")
          terra::writeRaster(terra::ifel(!is.na(r_wk_ok) & r_wk_ok > 0, 1, NA), f_wk_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
          f_prox_wk <- hulp_proximity_gdal(f_wk_bin, schil_buf_m)
          unlink(f_wk_bin)
          p_wk <- terra::rast(f_prox_wk)
          m_wk100 <- terra::ifel(!is.na(p_wk) & p_wk >= 0 & p_wk <= schil_buf_m, 1, NA)
          unlink(f_prox_wk)
          rm(p_wk)
          r_huidig <- terra::mask(r_huidig, m_wk100)
          rm(m_wk100)
        } else {
          r_huidig <- master_grid * NA_real_
        }
        rm(r_wk_ok)
      } else {
        r_huidig <- master_grid * NA_real_
      }
      
      res_rk <- cluster_filter_universeel(
        masker = r_r_raw, opp_laag = r_r_raw, min_opp_ha = drempel_rk_ha,
        min_totale_opp_ha = NA, afstand_m = 50, dispersie_m = 0, werkelijk = TRUE
      )
      r_rk_ok <- res_rk$raster
      rm(res_rk)
      
      if (!is_leeg(r_rk_ok)) {
        f_rk_bin <- tempfile(fileext = ".tif")
        terra::writeRaster(terra::ifel(!is.na(r_rk_ok) & r_rk_ok > 0, 1, NA), f_rk_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
        f_prox_rk <- hulp_proximity_gdal(f_rk_bin, schil_buf_m)
        unlink(f_rk_bin)
        p_rk <- terra::rast(f_prox_rk)
        m_rk100 <- terra::ifel(!is.na(p_rk) & p_rk >= 0 & p_rk <= schil_buf_m, 1, NA)
        unlink(f_prox_rk)
        rm(p_rk)
        r_huidig <- terra::mask(r_huidig, m_rk100)
        rm(m_rk100)
      } else {
        r_huidig <- master_grid * NA_real_
      }
      rm(r_rk_ok)
      
      message(sprintf("       📊 Na 3 rietmoeras-schillen (elk <= %dm): %.2f ha", 
                      as.integer(schil_buf_m), calc_ha_exact(r_huidig)))
    }
    
    # 2-venlibellen: Boreo-alpiene venlibellen & venjuffers (water-first 500m + bospercentage-eisen)
    is_venlibel <- tolower(bestandsnaam) %in% c("speerwaterjuffer", "coenagrionhastulatum",
                                                "venglazenmaker", "aeshnajuncea",
                                                "venwitsnuitlibel", "leucorrhiniadubia",
                                                "zwarteheidelibel", "sympetrumdanae",
                                                "maanwaterjuffer", "coenagrionlunulatum",
                                                "noordsewitsnuitlibel", "leucorrhiniarubicunda",
                                                "hoogveenglanslibel", "somatochloraarctica")
    
    if (is_venlibel && grepl("^waterbiotoop", tolower(tp_naam)) && !grepl("not", tolower(tp_naam))) {
      message("       -> Boreo-alpien venlibellen-archetype uitvoeren (water-first 500m + bospercentage-eisen)...")
      is_hvg <- tolower(bestandsnaam) %in% c("hoogveenglanslibel", "somatochloraarctica")
      
      # 1. Spoor A: Water schonen en clusteren
      r_wat_raw <- lagen_ruw[[tp_naam]]
      not_wat_tp <- names(lagen_ruw)[grepl("waterbiotoop.*not", tolower(names(lagen_ruw)))][1]
      r_not_wat  <- if (!is.na(not_wat_tp)) lagen_ruw[[not_wat_tp]] else (master_grid * 0)
      r_wat_cl   <- terra::clamp(r_wat_raw - terra::ifel(is.na(r_not_wat), 0, r_not_wat), lower = 0, upper = 1.0)
      r_wat_cl   <- terra::ifel(r_wat_cl >= 0.20, r_wat_cl, NA_real_)
      
      w_drempel_ha <- as.numeric(dplyr::coalesce(as.numeric(get_param(rij, "MinPatch_ha", NA)), 0.01))
      w_zoek_m     <- as.numeric(dplyr::coalesce(as.numeric(get_param(rij, "AfstandBiotopen_m", NA)), 50))
      
      res_w <- cluster_filter_universeel(
        masker = r_wat_cl, opp_laag = r_wat_cl, min_opp_ha = w_drempel_ha,
        min_totale_opp_ha = NA, afstand_m = w_zoek_m, dispersie_m = 0, werkelijk = TRUE
      )
      r_wat_ok <- res_w$raster
      rm(res_w, r_wat_cl, r_not_wat, r_wat_raw)
      
      # 2. Spoor B: Nabij land schonen & boseis toetsen (voor hoogveenglanslibel >= 5 ha, 10m, min 60% bos)
      lnd_tp <- params_gesorteerd %>% 
        dplyr::filter(!is.na(Spoor) & Spoor == "B" & !grepl("not", tolower(Type)) & !grepl("bos", tolower(Type))) %>% 
        dplyr::slice(1)
      
      lnd_tp_naam <- trimws(lnd_tp$Type[1])
      r_lnd_raw   <- if (lnd_tp_naam %in% names(lagen_ruw)) lagen_ruw[[lnd_tp_naam]] else (master_grid * 0)
      not_lnd_tp  <- names(lagen_ruw)[grepl("landbiotoop.*not", tolower(names(lagen_ruw)))][1]
      r_not_lnd   <- if (!is.na(not_lnd_tp)) lagen_ruw[[not_lnd_tp]] else (master_grid * 0)
      r_lnd_cl    <- terra::clamp(r_lnd_raw - terra::ifel(is.na(r_not_lnd), 0, r_not_lnd), lower = 0, upper = 1.0)
      r_lnd_cl    <- terra::ifel(r_lnd_cl > 0, r_lnd_cl, NA_real_)
      
      l_drempel_ha <- as.numeric(dplyr::coalesce(as.numeric(get_param(lnd_tp, "MinPatch_ha", NA)), if (is_hvg) 5.0 else 10.0))
      l_zoek_m     <- as.numeric(dplyr::coalesce(as.numeric(get_param(lnd_tp, "AfstandBiotopen_m", NA)), if (is_hvg) 10 else 50))
      l_bos_pct    <- as.numeric(dplyr::coalesce(as.numeric(get_param(lnd_tp, "Min_Struweel_Pct", NA)), if (is_hvg) 60.0 else NA_real_))
      
      res_l <- cluster_filter_universeel(
        masker = r_lnd_cl, opp_laag = r_lnd_cl, min_opp_ha = l_drempel_ha,
        min_totale_opp_ha = NA, afstand_m = l_zoek_m, dispersie_m = 0, werkelijk = TRUE
      )
      r_lnd_10ha  <- res_l$raster
      cl_lnd_10ha <- res_l$clusters
      rm(res_l, r_lnd_cl, r_not_lnd, r_lnd_raw)
      
      if (!is_leeg(cl_lnd_10ha)) {
        bos_lnd_kand <- names(lagen_ruw)[grepl("landbiotoop.*(bos|struweel)", tolower(names(lagen_ruw))) & !grepl("ruim", tolower(names(lagen_ruw)))]
        
        if (!is.na(l_bos_pct) && l_bos_pct > 0 && length(bos_lnd_kand) > 0) {
          r_bos_lnd <- lagen_ruw[[bos_lnd_kand[1]]]
          r_bos_lnd_clean <- terra::ifel(!is.na(r_bos_lnd) & r_bos_lnd > 0, r_bos_lnd, 0)
          
          if (tolower(bestandsnaam) %in% c("speerwaterjuffer", "coenagrionhastulatum")) {
            w_50m <- terra::focalMat(master_grid, 50, type = "circle")
            w_50m[w_50m > 0] <- 1
            cl_calc_zone <- terra::focal(cl_lnd_10ha, w = w_50m, fun = "max", na.rm = TRUE)
            rm(w_50m)
          } else {
            cl_calc_zone <- cl_lnd_10ha
          }
          
          df_l_tot <- terra::zonal(r_lnd_10ha, cl_lnd_10ha, fun = "sum", na.rm = TRUE)
          colnames(df_l_tot) <- c("ID", "Totaal_Sum")
          
          bos_in_z <- terra::mask(r_bos_lnd_clean, cl_calc_zone)
          df_b_tot <- terra::zonal(bos_in_z, cl_calc_zone, fun = "sum", na.rm = TRUE)
          colnames(df_b_tot) <- c("ID", "Bos_Sum")
          rm(bos_in_z, r_bos_lnd_clean, r_bos_lnd, cl_calc_zone)
          
          dt_l <- merge(data.table::as.data.table(df_l_tot), data.table::as.data.table(df_b_tot), by = "ID", all.x = TRUE)
          dt_l[is.na(Bos_Sum), Bos_Sum := 0]
          dt_l[, Pct_Bos := (Bos_Sum / Totaal_Sum) * 100]
          valide_l_ids <- dt_l[!is.na(ID) & Pct_Bos >= l_bos_pct, ID]
          rm(df_l_tot, df_b_tot, dt_l)
          
          r_lnd_ok <- if (length(valide_l_ids) > 0) terra::mask(r_lnd_10ha, cl_lnd_10ha %in% as.integer(valide_l_ids)) else (master_grid * NA_real_)
        } else {
          r_lnd_ok <- r_lnd_10ha
        }
      } else {
        r_lnd_ok <- master_grid * NA_real_
      }
      rm(r_lnd_10ha, cl_lnd_10ha)
      
      # 3. Spoor C: Ruim land (>= 50 ha op 50 m voor hoogveenglanslibel zonder boseis; >= 30 ha + bos% voor overige)
      rm_tp <- params_gesorteerd %>% 
        dplyr::filter(!is.na(Spoor) & Spoor == "C" & !grepl("bos", tolower(Type))) %>% 
        dplyr::slice(1)
      
      r_rm_ok <- master_grid * 1
      if (nrow(rm_tp) > 0) {
        rm_tp_naam <- trimws(rm_tp$Type[1])
        r_ruim_raw <- if (rm_tp_naam %in% names(lagen_ruw)) lagen_ruw[[rm_tp_naam]] else (master_grid * 0)
        
        rm_drempel_ha <- as.numeric(dplyr::coalesce(as.numeric(get_param(rm_tp, "MinPatch_ha", NA)), if (is_hvg) 50.0 else 30.0))
        rm_zoek_m     <- as.numeric(dplyr::coalesce(as.numeric(get_param(rm_tp, "AfstandBiotopen_m", NA)), 50))
        rm_bos_pct    <- as.numeric(dplyr::coalesce(as.numeric(get_param(rm_tp, "Min_Struweel_Pct", NA)), if (is_hvg) NA_real_ else 40.0))
        
        res_rm <- cluster_filter_universeel(
          masker = r_ruim_raw, opp_laag = r_ruim_raw, min_opp_ha = rm_drempel_ha,
          min_totale_opp_ha = NA, afstand_m = rm_zoek_m, dispersie_m = 0, werkelijk = TRUE
        )
        r_rm_ha  <- res_rm$raster
        cl_rm_ha <- res_rm$clusters
        rm(res_rm, r_ruim_raw)
        
        if (!is_leeg(cl_rm_ha)) {
          if (!is.na(rm_bos_pct) && rm_bos_pct > 0) {
            bos_rm_kand <- names(lagen_ruw)[grepl("ruim.*bos|bos.*ruim", tolower(names(lagen_ruw)))]
            r_bos_rm <- if (length(bos_rm_kand) > 0) lagen_ruw[[bos_rm_kand[1]]] else (master_grid * 0)
            r_bos_rm_clean <- terra::ifel(!is.na(r_bos_rm) & r_bos_rm > 0, r_bos_rm, 0)
            
            df_rm_tot <- terra::zonal(r_rm_ha, cl_rm_ha, fun = "sum", na.rm = TRUE)
            colnames(df_rm_tot) <- c("ID", "Totaal_Sum")
            
            bos_in_rm <- terra::mask(r_bos_rm_clean, cl_rm_ha)
            df_b_rm   <- terra::zonal(bos_in_rm, cl_rm_ha, fun = "sum", na.rm = TRUE)
            colnames(df_b_rm) <- c("ID", "Bos_Sum")
            rm(bos_in_rm, r_bos_rm_clean, r_bos_rm)
            
            dt_rm <- merge(data.table::as.data.table(df_rm_tot), data.table::as.data.table(df_b_rm), by = "ID", all.x = TRUE)
            dt_rm[is.na(Bos_Sum), Bos_Sum := 0]
            dt_rm[, Pct_Bos := (Bos_Sum / Totaal_Sum) * 100]
            valide_rm_ids <- dt_rm[!is.na(ID) & Pct_Bos >= rm_bos_pct, ID]
            rm(df_rm_tot, df_b_rm, dt_rm)
            
            r_rm_ok <- if (length(valide_rm_ids) > 0) terra::mask(r_rm_ha, cl_rm_ha %in% as.integer(valide_rm_ids)) else (master_grid * NA_real_)
          } else {
            # Geen bospercentage-eis op ruim land (bv. Hoogveenglanslibel: 50 ha open heide/veen/boscomplex)
            r_rm_ok <- r_rm_ha
          }
        } else {
          r_rm_ok <- master_grid * NA_real_
        }
        rm(r_rm_ha, cl_rm_ha)
      }
      
      # 4. Interacties:
      heeft_w <- !is_leeg(r_wat_ok)
      heeft_l <- !is_leeg(r_lnd_ok)
      heeft_r <- !is_leeg(r_rm_ok)
      
      if (is_hvg) {
        # HOOGVEENGLANSLIBEL:
        # 1) Water (0.01 ha, 50m) moet binnen <= 500 m van ruim land (50 ha, 50m) liggen
        # 2) Aanpalend bosrijk land (5 ha, 10m, >= 60% bos) wordt afgesneden op <= 500 m van het water
        if (heeft_w && heeft_r) {
          f_rm_tmp <- tempfile(fileext = ".tif")
          terra::writeRaster(terra::ifel(!is.na(r_rm_ok) & r_rm_ok > 0, 1, NA), f_rm_tmp, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
          f_prox_rm500 <- hulp_proximity_gdal(f_rm_tmp, 500)
          unlink(f_rm_tmp)
          p_rm500 <- terra::rast(f_prox_rm500)
          m_rm500 <- terra::ifel(!is.na(p_rm500) & p_rm500 >= 0 & p_rm500 <= 500, 1, NA)
          unlink(f_prox_rm500)
          rm(p_rm500)
          
          water_finaal <- terra::mask(r_wat_ok, m_rm500)
          rm(m_rm500)
        } else {
          water_finaal <- master_grid * NA_real_
        }
        
        if (heeft_l && !is_leeg(water_finaal)) {
          f_w_tmp <- tempfile(fileext = ".tif")
          terra::writeRaster(terra::ifel(!is.na(water_finaal) & water_finaal > 0, 1, NA), f_w_tmp, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
          f_prox_w500 <- hulp_proximity_gdal(f_w_tmp, 500)
          unlink(f_w_tmp)
          p_w500 <- terra::rast(f_prox_w500)
          m_w500 <- terra::ifel(!is.na(p_w500) & p_w500 >= 0 & p_w500 <= 500, 1, NA)
          unlink(f_prox_w500)
          rm(p_w500)
          
          land_finaal <- terra::mask(r_lnd_ok, m_w500)
          rm(m_w500)
        } else {
          land_finaal <- master_grid * NA_real_
        }
        
      } else {
        # OVERIGE VENLIBELLEN: oevercontact (10m of 20m) & water-first afkapping (500m)
        oever_dist_m <- as.numeric(dplyr::coalesce(as.numeric(get_param(lnd_tp, "Grens_Afstand_m", NA)), 20))
        
        if (heeft_w && heeft_l && heeft_r) {
          f_w_tmp <- tempfile(fileext = ".tif")
          terra::writeRaster(terra::ifel(!is.na(r_wat_ok) & r_wat_ok > 0, 1, NA), f_w_tmp, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
          f_prox_w <- hulp_proximity_gdal(f_w_tmp, oever_dist_m)
          p_w <- terra::rast(f_prox_w)
          m_w <- terra::ifel(!is.na(p_w) & p_w >= 0 & p_w <= oever_dist_m, 1, NA)
          unlink(f_prox_w)
          rm(p_w)
          
          f_l_tmp <- tempfile(fileext = ".tif")
          terra::writeRaster(terra::ifel(!is.na(r_lnd_ok) & r_lnd_ok > 0, 1, NA), f_l_tmp, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
          f_prox_l <- hulp_proximity_gdal(f_l_tmp, oever_dist_m)
          p_l <- terra::rast(f_prox_l)
          m_l <- terra::ifel(!is.na(p_l) & p_l >= 0 & p_l <= oever_dist_m, 1, NA)
          unlink(f_prox_l)
          rm(p_l)
          
          r_wat_gekoppeld <- terra::mask(r_wat_ok, m_l)
          r_lnd_gekoppeld <- terra::mask(r_lnd_ok, m_w)
          rm(m_w, m_l)
          
          if (nrow(rm_tp) > 0 && tolower(bestandsnaam) %in% c("speerwaterjuffer", "coenagrionhastulatum", 
                                                              "venwitsnuitlibel", "leucorrhiniadubia",
                                                              "zwarteheidelibel", "sympetrumdanae")) {
            f_rm_tmp <- tempfile(fileext = ".tif")
            terra::writeRaster(terra::ifel(!is.na(r_rm_ok) & r_rm_ok > 0, 1, NA), f_rm_tmp, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
            f_prox_rm500 <- hulp_proximity_gdal(f_rm_tmp, 500)
            unlink(f_rm_tmp)
            p_rm500 <- terra::rast(f_prox_rm500)
            m_rm500 <- terra::ifel(!is.na(p_rm500) & p_rm500 >= 0 & p_rm500 <= 500, 1, NA)
            unlink(f_prox_rm500)
            rm(p_rm500)
            
            r_wat_gekoppeld <- terra::mask(r_wat_gekoppeld, m_rm500)
            r_lnd_gekoppeld <- terra::mask(r_lnd_gekoppeld, m_rm500)
            rm(m_rm500)
          }
          
          if (!is_leeg(r_wat_gekoppeld)) {
            f_wi_tmp <- tempfile(fileext = ".tif")
            terra::writeRaster(terra::ifel(!is.na(r_wat_gekoppeld) & r_wat_gekoppeld > 0, 1, NA), f_wi_tmp, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
            f_prox_wi500 <- hulp_proximity_gdal(f_wi_tmp, 500)
            unlink(f_wi_tmp)
            p_wi500 <- terra::rast(f_prox_wi500)
            m_wi500 <- terra::ifel(!is.na(p_wi500) & p_wi500 >= 0 & p_wi500 <= 500, 1, NA)
            unlink(f_prox_wi500)
            rm(p_wi500)
            
            water_finaal <- r_wat_gekoppeld
            land_finaal  <- terra::mask(r_lnd_gekoppeld, m_wi500)
            rm(m_wi500)
          } else {
            water_finaal <- master_grid * NA_real_
            land_finaal  <- master_grid * NA_real_
          }
          unlink(c(f_w_tmp, f_l_tmp))
          rm(r_wat_gekoppeld, r_lnd_gekoppeld)
        } else {
          water_finaal <- master_grid * NA_real_
          land_finaal  <- master_grid * NA_real_
        }
      }
      rm(r_wat_ok, r_lnd_ok, r_rm_ok)
      
      r_huidig <- voeg_toe_geclamped(water_finaal, land_finaal)
      verwerkte_sporen[[tp_naam]]     <- water_finaal
      verwerkte_sporen[[lnd_tp_naam]] <- land_finaal
      
      message(sprintf("       📊 Na boreo-alpiene venlibellen-analyse: %.2f ha (Water: %.2f ha, Land: %.2f ha)", 
                      calc_ha_exact(r_huidig), calc_ha_exact(water_finaal), calc_ha_exact(land_finaal)))
      next
    }
    
    # 2-waterspitsmuis: Semi-aquatische oeverzoogdieren (Geoptimaliseerd voor Vlaanderen)
    is_waterspitsmuis <- tolower(bestandsnaam) %in% c("waterspitsmuis", "neomysfodiens")
    
    if (is_waterspitsmuis && tp_naam == zelfde_spoor_types[1] && huidig_spoor == "A") {
      message("       -> Symmetrische oeveranalyse uitvoeren (30% contact water/land + 50m oeverstrook)...")
      
      r_wat_tot <- r_huidig
      r_lnd_raw <- if ("bwk_land" %in% names(lagen_ruw)) lagen_ruw[["bwk_land"]] else (master_grid * NA_real_)
      
      if (!is_leeg(r_wat_tot) && 
          !is_leeg(r_lnd_raw)) {
        
        e <- terra::ext(master_grid)
        reg_str <- sprintf("%f,%f,%f,%f", e$xmin, e$xmax, e$ymin, e$ymax)
        
        # Hulpfunctie voor razendsnelle 8-richtingen r.clump via GRASS
        fast_clump8 <- function(r_in, titel) {
          f_in  <- tempfile(fileext = ".tif")
          f_out <- tempfile(fileext = ".tif")
          terra::writeRaster(terra::ifel(!is.na(r_in) & r_in > 0, 1L, NA_integer_), f_in, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
          qgisprocess::qgis_run_algorithm(
            "grass:r.clump", input = f_in, output = f_out, title = titel,
            "-d" = TRUE, GRASS_REGION_PARAMETER = reg_str,
            GRASS_REGION_CELLSIZE_PARAMETER = 10, GRASS_RASTER_FORMAT_OPT = "COMPRESS=LZW", .quiet = TRUE
          )
          unlink(f_in)
          cl <- terra::rast(f_out)
          try(terra::set.cats(cl, layer = 1, value = NULL), silent = TRUE)
          try(terra::coltab(cl) <- NULL, silent = TRUE)
          return(list(rast = cl, file = f_out))
        }
        
        r_bin_w <- terra::ifel(!is.na(r_wat_tot) & r_wat_tot > 0, 1L, 0L)
        r_bin_l <- terra::ifel(!is.na(r_lnd_raw) & r_lnd_raw > 0, 1L, 0L)
        m_buren <- matrix(c(1, 1, 1, 1, 0, 1, 1, 1, 1), nrow = 3, ncol = 3)
        
        # 1. Waterclusters testen op >= 30% geschikt land langs oever
        f_land_w <- terra::focal(r_bin_l, w = m_buren, fun = "sum", na.rm = TRUE)
        f_oev_w  <- terra::focal(1L - r_bin_w, w = m_buren, fun = "sum", na.rm = TRUE)
        
        cl_w_obj <- fast_clump8(r_wat_tot, "cl_wsm_w")
        dt_w <- data.table::data.table(
          ID   = as.integer(terra::values(cl_w_obj$rast, mat = FALSE)),
          Tot  = as.numeric(terra::values(f_oev_w, mat = FALSE)),
          Land = as.numeric(terra::values(f_land_w, mat = FALSE))
        )[!is.na(ID) & ID > 0]
        rm(f_land_w, f_oev_w)
        
        dt_w_agg <- dt_w[, .(Tot_Oever = sum(Tot, na.rm = TRUE), Land_Contact = sum(Land, na.rm = TRUE)), by = ID]
        ok_w_ids <- dt_w_agg[Tot_Oever > 0 & (Land_Contact / Tot_Oever) >= 0.30, ID]
        rm(dt_w, dt_w_agg)
        
        water_ok <- if (length(ok_w_ids) > 0) terra::mask(r_wat_tot, cl_w_obj$rast %in% ok_w_ids) else (master_grid * NA_real_)
        unlink(cl_w_obj$file)
        rm(cl_w_obj)
        
        # 2. Landclusters testen op >= 30% water langs rand
        f_wat_l  <- terra::focal(r_bin_w, w = m_buren, fun = "sum", na.rm = TRUE)
        f_rand_l <- terra::focal(1L - r_bin_l, w = m_buren, fun = "sum", na.rm = TRUE)
        rm(r_bin_w, r_bin_l)
        
        cl_l_obj <- fast_clump8(r_lnd_raw, "cl_wsm_l")
        dt_l <- data.table::data.table(
          ID  = as.integer(terra::values(cl_l_obj$rast, mat = FALSE)),
          Tot = as.numeric(terra::values(f_rand_l, mat = FALSE)),
          Wat = as.numeric(terra::values(f_wat_l, mat = FALSE))
        )[!is.na(ID) & ID > 0]
        rm(f_wat_l, f_rand_l)
        
        dt_l_agg <- dt_l[, .(Tot_Rand = sum(Tot, na.rm = TRUE), Wat_Contact = sum(Wat, na.rm = TRUE)), by = ID]
        ok_l_ids <- dt_l_agg[Tot_Rand > 0 & (Wat_Contact / Tot_Rand) >= 0.30, ID]
        rm(dt_l, dt_l_agg)
        
        land_ok_direct <- if (length(ok_l_ids) > 0) terra::mask(r_lnd_raw, cl_l_obj$rast %in% ok_l_ids) else (master_grid * NA_real_)
        unlink(cl_l_obj$file)
        rm(cl_l_obj)
        
        # 3. Aaneengesloten 50m oeverstrook rond goedgekeurd water
        if (!is_leeg(water_ok)) {
          f_w_bin   <- tempfile(fileext = ".tif")
          terra::writeRaster(terra::ifel(!is.na(water_ok) & water_ok > 0, 1L, NA_integer_), f_w_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
          f_prox_50 <- hulp_proximity_gdal(f_w_bin, 50)
          f_prox_10 <- hulp_proximity_gdal(f_w_bin, 10)
          unlink(f_w_bin)
          
          p_w50 <- terra::rast(f_prox_50)
          m_w50 <- terra::ifel(!is.na(p_w50) & p_w50 >= 0 & p_w50 <= 50, 1L, NA_integer_)
          unlink(f_prox_50)
          rm(p_w50)
          
          land_oever_raw <- terra::mask(r_lnd_raw, m_w50)
          rm(m_w50)
          
          p_w10 <- terra::rast(f_prox_10)
          m_w10 <- terra::ifel(!is.na(p_w10) & p_w10 >= 0 & p_w10 <= 10, 1L, NA_integer_)
          unlink(f_prox_10)
          rm(p_w10)
          
          cl_oev_obj <- fast_clump8(land_oever_raw, "cl_wsm_oev")
          cl_hit     <- terra::mask(cl_oev_obj$rast, m_w10)
          ids_touch  <- unique(na.omit(as.integer(terra::values(cl_hit, mat = FALSE))))
          rm(m_w10, cl_hit)
          
          land_oever_50m <- if (length(ids_touch) > 0) terra::mask(land_oever_raw, cl_oev_obj$rast %in% ids_touch) else (master_grid * NA_real_)
          unlink(cl_oev_obj$file)
          rm(cl_oev_obj, land_oever_raw)
        } else {
          land_oever_50m <- master_grid * NA_real_
        }
        
        oever_totaal <- voeg_toe_geclamped(land_oever_50m, land_ok_direct)
        r_huidig     <- voeg_toe_geclamped(water_ok, oever_totaal)
        
        verwerkte_sporen[["water_bwk"]] <- water_ok
        verwerkte_sporen[["bwk_land"]]  <- oever_totaal
        
        rm(land_oever_50m, land_ok_direct, oever_totaal, water_ok)
        gc()
      } else {
        r_huidig <- master_grid * NA_real_
      }
      message(sprintf("       📊 Na symmetrische oeveranalyse: %.2f ha", calc_ha_exact(r_huidig)))
      next
    }
    
    if (is_waterspitsmuis && huidig_spoor == "B") {
      next
    }
    
    # 2-kemphaan: Plasdras-weidevogels
    is_kemphaan <- tolower(bestandsnaam) %in% c("kemphaan", "kemphaan_wv", "calidrispugnax", "philomachuspugnax")
    
    if (is_kemphaan && tp_naam == zelfde_spoor_types[1] && huidig_spoor == "A") {
      message("       -> Kemphaan plasdras-archetype uitvoeren (grasland <= 100m van water/moeras + water/moeras)...")
      
      if (!is.null(paths$huetzon) && file.exists(paths$huetzon)) {
        f_hz_al <- hulp_align_gdal(paths$huetzon, master_grid)
        r_hz    <- terra::rast(f_hz_al)
        r_hz_ok <- terra::ifel(!is.na(r_hz) & r_hz > 0, 1.0, NA_real_)
        unlink(f_hz_al)
        rm(r_hz)
      } else {
        r_hz_ok <- master_grid * NA_real_
      }
      
      r_bwk_wm <- vul_raster_vanuit_parquet(
        ds_parquet  = ds_parquet,
        codes_bevat = c("31", "2190", "mr", "v", "kn", "a", "da"),
        master_grid = master_grid,
        waarde_vast = 1
      )
      
      r_water_moeras <- voeg_toe_geclamped(r_bwk_wm, r_hz_ok)
      rm(r_bwk_wm, r_hz_ok)
      
      if (!is_leeg(r_water_moeras)) {
        f_wm_bin <- tempfile(fileext = ".tif")
        terra::writeRaster(terra::ifel(!is.na(r_water_moeras) & r_water_moeras > 0, 1, NA), f_wm_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
        f_prox_wm <- hulp_proximity_gdal(f_wm_bin, 100)
        unlink(f_wm_bin)
        p_wm <- terra::rast(f_prox_wm)
        m_wm100 <- terra::ifel(!is.na(p_wm) & p_wm >= 0 & p_wm <= 100, 1, NA)
        unlink(f_prox_wm)
        rm(p_wm)
        
        r_gras_gekoppeld <- terra::mask(r_huidig, m_wm100)
        rm(m_wm100)
        r_huidig <- voeg_toe_geclamped(r_water_moeras, r_gras_gekoppeld)
        rm(r_gras_gekoppeld)
      }
      rm(r_water_moeras)
      message(sprintf("       📊 Na plasdras-koppeling (<= 100m water/moeras): %.2f ha", calc_ha_exact(r_huidig)))
    }
    
    # 2-smient: Watervogels met grasland-water-asymmetrie
    is_smient <- tolower(bestandsnaam) %in% c("smient", "smient_wv", "marecapenelope", "anaspenelope")
    
    if (is_smient && tp_naam == zelfde_spoor_types[1] && huidig_spoor == "A") {
      message("       -> Smient-archetype uitvoeren (water + Huet/watervlakken <-> grasland <= 100m, losstaand water >= 1ha)...")
      
      r_wat_tot <- r_huidig
      
      if (!is.null(paths$huetzon) && file.exists(paths$huetzon)) {
        m_huet <- haal_categorisch_masker(paths$huetzon, c("Mlz", "O1brak", "O1o", "O2zout", "Rg"), master_grid, paths$dir_cache, "huet")
        r_wat_tot <- voeg_toe_geclamped(r_wat_tot, m_huet)
        rm(m_huet)
      }
      
      pad_wv <- get_path("watervlakken_raster", "watervlakken")
      if (!is.null(pad_wv) && file.exists(pad_wv)) {
        f_wv_al <- hulp_align_gdal(pad_wv, master_grid)
        r_wv    <- terra::rast(f_wv_al)
        r_wv_b  <- terra::ifel(!is.na(r_wv) & r_wv > 0, 1.0, NA_real_)
        unlink(f_wv_al)
        rm(r_wv)
        r_wat_tot <- voeg_toe_geclamped(r_wat_tot, r_wv_b)
        rm(r_wv_b)
      }
      
      r_gras_raw <- if ("foerageer_grasland" %in% names(lagen_ruw)) lagen_ruw[["foerageer_grasland"]] else (master_grid * NA_real_)
      
      if (!is_leeg(r_wat_tot) && !is_leeg(r_gras_raw)) {
        f_w_bin <- tempfile(fileext = ".tif")
        terra::writeRaster(terra::ifel(!is.na(r_wat_tot) & r_wat_tot > 0, 1, NA), f_w_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
        f_prox_w100 <- hulp_proximity_gdal(f_w_bin, 100)
        unlink(f_w_bin)
        p_w100 <- terra::rast(f_prox_w100)
        m_w100 <- terra::ifel(!is.na(p_w100) & p_w100 >= 0 & p_w100 <= 100, 1, NA)
        unlink(f_prox_w100)
        rm(p_w100)
        
        r_gras_ok <- terra::mask(r_gras_raw, m_w100)
        rm(m_w100)
      } else {
        r_gras_ok <- master_grid * NA_real_
      }
      
      if (!is_leeg(r_gras_ok)) {
        f_g_bin <- tempfile(fileext = ".tif")
        terra::writeRaster(terra::ifel(!is.na(r_gras_ok) & r_gras_ok > 0, 1, NA), f_g_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
        f_prox_g20 <- hulp_proximity_gdal(f_g_bin, 20)
        unlink(f_g_bin)
        p_g20 <- terra::rast(f_prox_g20)
        m_g20 <- terra::ifel(!is.na(p_g20) & p_g20 >= 0 & p_g20 <= 20, 1, NA)
        unlink(f_prox_g20)
        rm(p_g20)
        
        r_wat_aan_gras <- terra::mask(r_wat_tot, m_g20)
        r_wat_los      <- terra::mask(r_wat_tot, m_g20, inverse = TRUE)
        rm(m_g20)
      } else {
        r_wat_aan_gras <- master_grid * NA_real_
        r_wat_los      <- r_wat_tot
      }
      
      drempel_w_ha <- as.numeric(dplyr::coalesce(as.numeric(get_param(rij, "MinPatch_ha", NA)), 1.0))
      res_wlos <- cluster_filter_universeel(
        masker = r_wat_los, opp_laag = r_wat_los, min_opp_ha = drempel_w_ha,
        min_totale_opp_ha = NA, afstand_m = 50, dispersie_m = 0, werkelijk = TRUE
      )
      r_wat_los_ok <- res_wlos$raster
      rm(res_wlos, r_wat_los, r_wat_tot)
      
      r_wat_finaal <- voeg_toe_geclamped(r_wat_aan_gras, r_wat_los_ok)
      r_huidig     <- voeg_toe_geclamped(r_wat_finaal, r_gras_ok)
      
      verwerkte_sporen[["foerageer_water"]]    <- r_wat_finaal
      verwerkte_sporen[["foerageer_grasland"]] <- r_gras_ok
      
      rm(r_wat_aan_gras, r_wat_los_ok, r_wat_finaal, r_gras_ok, r_gras_raw)
      message(sprintf("       📊 Na Smient-analyse: %.2f ha", calc_ha_exact(r_huidig)))
      next
    }
    
    if (is_smient && huidig_spoor == "B") {
      next
    }
    
    # 2-kleinerietgans: Foerageercomplex >= 50 ha op 50 m WAARVAN minstens 80% grasland
    is_kleinerietgans <- tolower(bestandsnaam) %in% c("kleinerietgans", "kleinerietgans_wv", "anserbrachyrhynchus")
    
    if (is_kleinerietgans && tp_naam == zelfde_spoor_types[1] && huidig_spoor == "A") {
      message("       -> Kleine rietgans-archetype uitvoeren (>= 50 ha op 50 m waarvan >= 80% grasland)...")
      
      r_gras_sub <- if ("grasland_bwk" %in% names(lagen_ruw)) lagen_ruw[["grasland_bwk"]] else (master_grid * NA_real_)
      
      if (!is_leeg(r_huidig)) {
        f_bin <- tempfile(fileext = ".tif")
        terra::writeRaster(terra::ifel(!is.na(r_huidig) & r_huidig > 0, 1L, NA_integer_), f_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
        f_buf <- hulp_proximity_gdal(f_bin, 25)
        unlink(f_bin)
        p_buf <- terra::rast(f_buf)
        m_buf <- terra::ifel(!is.na(p_buf) & p_buf >= 0 & p_buf <= 25, 1L, NA_integer_)
        unlink(f_buf)
        rm(p_buf)
        
        f_buf_bin <- tempfile(fileext = ".tif")
        f_cl_net  <- tempfile(fileext = ".tif")
        terra::writeRaster(m_buf, f_buf_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
        rm(m_buf)
        
        e <- terra::ext(master_grid)
        qgisprocess::qgis_run_algorithm(
          "grass:r.clump", input = f_buf_bin, output = f_cl_net, title = "cl_kleinerietgans",
          "-d" = FALSE, GRASS_REGION_PARAMETER = sprintf("%f,%f,%f,%f", e$xmin, e$xmax, e$ymin, e$ymax),
          GRASS_REGION_CELLSIZE_PARAMETER = 10, GRASS_RASTER_FORMAT_OPT = "COMPRESS=LZW", .quiet = TRUE
        )
        unlink(f_buf_bin)
        
        cl_net <- terra::rast(f_cl_net)
        try(terra::set.cats(cl_net, layer = 1, value = NULL), silent = TRUE)
        try(terra::coltab(cl_net) <- NULL, silent = TRUE)
        cl_hab <- terra::mask(cl_net, r_huidig)
        rm(cl_net)
        
        df_tot <- terra::zonal(r_huidig, cl_hab, fun = "sum", na.rm = TRUE)
        colnames(df_tot) <- c("ID", "Som_Tot")
        dt_tot <- data.table::as.data.table(df_tot)[!is.na(ID) & ID > 0]
        dt_tot[, Area_Tot_ha := Som_Tot * 0.01]
        
        if (!is_leeg(r_gras_sub)) {
          df_gras <- terra::zonal(r_gras_sub, cl_hab, fun = "sum", na.rm = TRUE)
          colnames(df_gras) <- c("ID", "Som_Gras")
          dt_gras <- data.table::as.data.table(df_gras)[!is.na(ID) & ID > 0]
        } else {
          dt_gras <- data.table::data.table(ID = dt_tot$ID, Som_Gras = 0)
        }
        
        dt_stats <- merge(dt_tot, dt_gras, by = "ID", all.x = TRUE)
        dt_stats[is.na(Som_Gras), Som_Gras := 0]
        dt_stats[, Pct_Gras := ifelse(Som_Tot > 0, Som_Gras / Som_Tot, 0)]
        
        drempel_ha  <- as.numeric(dplyr::coalesce(as.numeric(get_param(rij, "MinPatch_ha", NA)), 50.0))
        geldige_ids <- dt_stats[Area_Tot_ha >= drempel_ha & Pct_Gras >= 0.80, ID]
        
        if (length(geldige_ids) > 0) {
          m_ok <- terra::ifel(cl_hab %in% as.integer(geldige_ids), 1L, NA_integer_)
          r_huidig <- terra::mask(r_huidig, m_ok)
          rm(m_ok)
        } else {
          r_huidig <- master_grid * NA_real_
        }
        unlink(f_cl_net)
        rm(cl_hab, df_tot, dt_tot, dt_gras, dt_stats)
      }
      
      verwerkte_sporen[["foerageer_bwk"]] <- r_huidig
      rm(r_gras_sub)
      gc()
      message(sprintf("       📊 Na Kleine rietgans-filter (>= 50 ha & >= 80%% grasland): %.2f ha", calc_ha_exact(r_huidig)))
      next
    }
    
    # 2-nauwekorfslak: Duinmodel (ecoregio == 1, 50m) + Binnenlandmodel (voortplanting_leem + ab <= 100m, 100m), min 0.04 ha
    is_nauwekorfslak <- tolower(bestandsnaam) %in% c("nauwekorfslak", "vertigoangustior")
    
    if (is_nauwekorfslak && tp_naam == zelfde_spoor_types[1] && huidig_spoor == "A") {
      message("       -> Nauwe korfslak-archetype uitvoeren (Duin 50m + Binnenland leem/zegge + ab <= 100m, min 0.04 ha)...")
      
      r_duin_raw <- r_huidig
      r_leem_raw <- if ("voortplanting_leem" %in% names(lagen_ruw)) lagen_ruw[["voortplanting_leem"]] else (master_grid * NA_real_)
      r_ab_raw   <- if ("ab_biotoop" %in% names(lagen_ruw)) lagen_ruw[["ab_biotoop"]] else (master_grid * NA_real_)
      
      pad_eco <- get_path("ecoregio", fallback = here::here("data/input/Raster_Vlaanderen/vlaanderen_ecoregio_10m.tif"))
      if (!is.null(pad_eco)) {
        f_eco_al <- hulp_align_gdal(pad_eco, master_grid)
        r_eco    <- terra::rast(f_eco_al)
        r_duin_raw <- terra::mask(r_duin_raw, terra::ifel(r_eco == 1, 1, NA))
        r_leem_raw <- terra::mask(r_leem_raw, terra::ifel(r_eco > 1, 1, NA))
        r_ab_raw   <- terra::mask(r_ab_raw,   terra::ifel(r_eco > 1, 1, NA))
        unlink(f_eco_al)
        rm(r_eco)
      }
      
      if (!is_leeg(r_duin_raw)) {
        res_duin <- cluster_filter_universeel(
          masker = r_duin_raw, opp_laag = r_duin_raw, min_opp_ha = 0.04,
          min_totale_opp_ha = NA, afstand_m = 50, dispersie_m = 0, werkelijk = TRUE
        )
        r_duin_ok <- res_duin$raster
        rm(res_duin)
      } else {
        r_duin_ok <- master_grid * NA_real_
      }
      
      if (!is_leeg(r_leem_raw)) {
        if (!is_leeg(r_ab_raw)) {
          f_l_bin <- tempfile(fileext = ".tif")
          terra::writeRaster(terra::ifel(!is.na(r_leem_raw) & r_leem_raw > 0, 1, NA), f_l_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
          f_prox_l100 <- hulp_proximity_gdal(f_l_bin, 100)
          unlink(f_l_bin)
          p_l100 <- terra::rast(f_prox_l100)
          m_l100 <- terra::ifel(!is.na(p_l100) & p_l100 >= 0 & p_l100 <= 100, 1, NA)
          unlink(f_prox_l100)
          rm(p_l100)
          
          r_ab_ok <- terra::mask(r_ab_raw, m_l100)
          rm(m_l100)
          r_leem_comb <- voeg_toe_geclamped(r_leem_raw, r_ab_ok)
          rm(r_ab_ok)
        } else {
          r_leem_comb <- r_leem_raw
        }
        
        res_leem <- cluster_filter_universeel(
          masker = r_leem_comb, opp_laag = r_leem_comb, min_opp_ha = 0.04,
          min_totale_opp_ha = NA, afstand_m = 100, dispersie_m = 0, werkelijk = TRUE
        )
        r_leem_ok <- res_leem$raster
        rm(res_leem, r_leem_comb)
      } else {
        r_leem_ok <- master_grid * NA_real_
      }
      
      r_totaal_nk <- voeg_toe_geclamped(r_duin_ok, r_leem_ok)
      if (!is_leeg(r_totaal_nk)) {
        res_fin <- cluster_filter_universeel(
          masker = r_totaal_nk, opp_laag = r_totaal_nk, min_opp_ha = 0.04,
          min_totale_opp_ha = NA, afstand_m = 100, dispersie_m = 0, werkelijk = TRUE
        )
        r_huidig <- res_fin$raster
        rm(res_fin)
      } else {
        r_huidig <- master_grid * NA_real_
      }
      
      verwerkte_sporen[["voortplanting_duin"]] <- r_duin_ok
      verwerkte_sporen[["voortplanting_leem"]] <- r_leem_ok
      verwerkte_sporen[["ab_biotoop"]]         <- master_grid * NA_real_
      
      rm(r_duin_raw, r_leem_raw, r_ab_raw, r_duin_ok, r_leem_ok, r_totaal_nk)
      message(sprintf("       📊 Na Nauwe korfslak-analyse: %.2f ha", calc_ha_exact(r_huidig)))
      next
    }
    
    if (is_nauwekorfslak && huidig_spoor %in% c("B", "C")) {
      next
    }
    
    # 2-platteschijfhoren: Spoor A (water BWK) + Spoor B ((voortplanting_h + Huet) & bodemtype) -> 50m, min 1.0 ha
    is_platteschijfhoren <- tolower(bestandsnaam) %in% c("platteschijfhoren", "anisusvorticulus")
    
    if (is_platteschijfhoren && tp_naam == zelfde_spoor_types[1] && huidig_spoor == "A") {
      message("       -> Platte schijfhoren-archetype uitvoeren (Spoor A water + Spoor B grasland/Huet op bodemtype, 50m, min 1 ha)...")
      
      r_spoor_a <- r_huidig
      r_spoor_b <- if ("voortplanting_h" %in% names(lagen_ruw)) lagen_ruw[["voortplanting_h"]] else (master_grid * NA_real_)
      
      if (!is.null(paths$huetzon) && file.exists(paths$huetzon)) {
        m_huet <- haal_categorisch_masker(paths$huetzon, c("Bk", "Bkk", "Bkl", "Bkz", "P", "Rk"), master_grid, paths$dir_cache, "huet")
        r_spoor_b <- voeg_toe_geclamped(r_spoor_b, m_huet)
        rm(m_huet)
      }
      
      pad_bodem <- get_path("bodem_bodemtype", fallback = here::here("data/input/Raster_Vlaanderen/vlaanderen_bodem_bodemtype_10m.tif"))
      if (!is.null(pad_bodem) && !is_leeg(r_spoor_b)) {
        f_bod_al <- hulp_align_gdal(pad_bodem, master_grid)
        r_bod    <- terra::rast(f_bod_al)
        bodem_ids <- c(89, 1460, 2830, 2913, 3089, 3381, 3513, 3640, 3687, 4011, 4036)
        m_bod    <- terra::ifel(r_bod %in% bodem_ids, 1, NA)
        r_spoor_b <- terra::mask(r_spoor_b, m_bod)
        unlink(f_bod_al)
        rm(r_bod, m_bod)
      } else {
        warning("⚠️ bodem_bodemtype raster niet gevonden voor Platte schijfhoren; Spoor B wordt zonder bodemfilter meegenomen.")
      }
      
      r_comb <- voeg_toe_geclamped(r_spoor_a, r_spoor_b)
      drempel_ps_ha <- as.numeric(dplyr::coalesce(as.numeric(get_param(rij, "MinPatch_ha", NA)), 1.0))
      
      if (!is_leeg(r_comb)) {
        res_ps <- cluster_filter_universeel(
          masker = r_comb, opp_laag = r_comb, min_opp_ha = drempel_ps_ha,
          min_totale_opp_ha = NA, afstand_m = 50, dispersie_m = 0, werkelijk = TRUE
        )
        r_huidig <- res_ps$raster
        rm(res_ps)
      } else {
        r_huidig <- master_grid * NA_real_
      }
      
      verwerkte_sporen[["voortplanting_a"]] <- r_huidig
      verwerkte_sporen[["voortplanting_h"]] <- master_grid * NA_real_
      
      rm(r_spoor_a, r_spoor_b, r_comb)
      message(sprintf("       📊 Na Platte schijfhoren-analyse: %.2f ha", calc_ha_exact(r_huidig)))
      next
    }
    
    if (is_platteschijfhoren && huidig_spoor == "B") {
      next
    }
    
    # 2-kempenseheidelibel
    is_kempenseheidelibel <- tolower(bestandsnaam) %in% c("kempenseheidelibel", "sympetrumdepressiusculum")
    
    if (is_kempenseheidelibel && tp_naam == zelfde_spoor_types[1] && huidig_spoor == "A") {
      message("       -> Kempense heidelibel-archetype uitvoeren (water >= 1 ha met >= 0.8 ha moeras/ruigte <= 10m <-> land >= 2.5 ha <= 100m)...")
      
      r_w1_raw <- r_huidig
      r_w2_raw <- if ("waterbiotoop_bwk2" %in% names(lagen_ruw)) lagen_ruw[["waterbiotoop_bwk2"]] else (master_grid * NA_real_)
      r_lb_raw <- if ("landbiotoop_bwk" %in% names(lagen_ruw))   lagen_ruw[["landbiotoop_bwk"]]   else (master_grid * NA_real_)
      
      if (!is_leeg(r_w1_raw) &&
          !is_leeg(r_w2_raw) &&
          !is_leeg(r_lb_raw)) {
        
        res_w1 <- cluster_filter_universeel(
          masker = r_w1_raw, opp_laag = r_w1_raw, min_opp_ha = 1.0,
          min_totale_opp_ha = NA, afstand_m = 20, dispersie_m = 0, werkelijk = TRUE
        )
        r_w1_ok  <- res_w1$raster
        cl_w1_ok <- res_w1$clusters
        
        res_w2 <- cluster_filter_universeel(
          masker = r_w2_raw, opp_laag = r_w2_raw, min_opp_ha = 0.8,
          min_totale_opp_ha = NA, afstand_m = 20, dispersie_m = 0, werkelijk = TRUE
        )
        r_w2_ok <- res_w2$raster
        rm(res_w1, res_w2)
        
        res_lb <- cluster_filter_universeel(
          masker = r_lb_raw, opp_laag = r_lb_raw, min_opp_ha = 2.5,
          min_totale_opp_ha = NA, afstand_m = 50, dispersie_m = 0, werkelijk = TRUE
        )
        r_lb_ok <- res_lb$raster
        rm(res_lb)
        
        if (!is_leeg(r_w1_ok) &&
            !is_leeg(r_w2_ok) &&
            !is_leeg(r_lb_ok)) {
          
          # Stap A: hele waterclusters (>= 1 ha) behouden die binnen <= 10m van moeras/ruigte (>= 0.8 ha) liggen
          f_w2_bin <- tempfile(fileext = ".tif")
          terra::writeRaster(terra::ifel(!is.na(r_w2_ok) & r_w2_ok > 0, 1, NA), f_w2_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
          f_prox_w2 <- hulp_proximity_gdal(f_w2_bin, 10)
          unlink(f_w2_bin)
          p_w2 <- terra::rast(f_prox_w2)
          m_w2 <- terra::ifel(!is.na(p_w2) & p_w2 >= 0 & p_w2 <= 10, 1, NA)
          unlink(f_prox_w2)
          rm(p_w2)
          
          cl_w1_hit <- terra::mask(cl_w1_ok, m_w2)
          rm(m_w2)
          geldige_w1_ids <- unique(na.omit(terra::values(cl_w1_hit, mat = FALSE)))
          rm(cl_w1_hit)
          
          if (length(geldige_w1_ids) > 0) {
            r_w3 <- terra::mask(r_w1_ok, terra::ifel(cl_w1_ok %in% geldige_w1_ids, 1, NA))
            
            # Stap B: wederzijdse koppeling <= 100m tussen r_w3 (water) en r_lb_ok (land >= 2.5 ha)
            f_lb_bin <- tempfile(fileext = ".tif")
            terra::writeRaster(terra::ifel(!is.na(r_lb_ok) & r_lb_ok > 0, 1, NA), f_lb_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
            f_prox_lb <- hulp_proximity_gdal(f_lb_bin, 100)
            unlink(f_lb_bin)
            p_lb <- terra::rast(f_prox_lb)
            m_lb <- terra::ifel(!is.na(p_lb) & p_lb >= 0 & p_lb <= 100, 1, NA)
            unlink(f_prox_lb)
            rm(p_lb)
            
            water_finaal <- terra::mask(r_w3, m_lb)
            rm(m_lb)
            
            if (!is_leeg(water_finaal)) {
              f_w3_bin <- tempfile(fileext = ".tif")
              terra::writeRaster(terra::ifel(!is.na(water_finaal) & water_finaal > 0, 1, NA), f_w3_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
              f_prox_w3 <- hulp_proximity_gdal(f_w3_bin, 100)
              unlink(f_w3_bin)
              p_w3 <- terra::rast(f_prox_w3)
              m_w3 <- terra::ifel(!is.na(p_w3) & p_w3 >= 0 & p_w3 <= 100, 1, NA)
              unlink(f_prox_w3)
              rm(p_w3)
              
              land_finaal <- terra::mask(r_lb_ok, m_w3)
              rm(m_w3)
            } else {
              land_finaal <- master_grid * NA_real_
            }
            rm(r_w3)
          } else {
            water_finaal <- master_grid * NA_real_
            land_finaal  <- master_grid * NA_real_
          }
        } else {
          water_finaal <- master_grid * NA_real_
          land_finaal  <- master_grid * NA_real_
        }
        rm(r_w1_ok, cl_w1_ok, r_w2_ok, r_lb_ok)
      } else {
        water_finaal <- master_grid * NA_real_
        land_finaal  <- master_grid * NA_real_
      }
      
      r_huidig <- voeg_toe_geclamped(water_finaal, land_finaal)
      
      verwerkte_sporen[["waterbiotoop_bwk1"]] <- water_finaal
      verwerkte_sporen[["landbiotoop_bwk"]]   <- land_finaal
      
      rm(r_w1_raw, r_w2_raw, r_lb_raw, water_finaal, land_finaal)
      message(sprintf("       📊 Na Kempense heidelibel-analyse (Water + Land): %.2f ha (Water: %.2f ha | Land: %.2f ha)", 
                      calc_ha_exact(r_huidig),
                      calc_ha_exact(verwerkte_sporen[["waterbiotoop_bwk1"]]),
                      calc_ha_exact(verwerkte_sporen[["landbiotoop_bwk"]])))
      next
    }
    
    if (is_kempenseheidelibel && huidig_spoor %in% c("B", "C")) {
      next
    }
    
    # SCENARIO 2h: Ecoregio-filter (hard afsnijden op ecoregiogrens vóór MinPatch clustering)
    eco_val <- trimws(as.character(get_param(rij, "Ecoregio", get_param(rij, "Ecoregio_Codes", ""))))
    
    if (eco_val != "" && !is.na(eco_val) && eco_val != "NA") {
      pad_eco <- get_path("ecoregio", fallback = here::here("data/input/Raster_Vlaanderen/vlaanderen_ecoregio_10m.tif"))
      if (!is.null(pad_eco)) {
        eco_targets <- as.integer(trimws(unlist(strsplit(eco_val, ","))))
        eco_targets <- eco_targets[!is.na(eco_targets)]
        
        if (length(eco_targets) > 0 && !is_leeg(r_huidig)) {
          message(sprintf("       -> Ecoregio-filter toepassen (toegestane NR-codes: %s)...", paste(eco_targets, collapse = ", ")))
          
          f_eco_al <- hulp_align_gdal(pad_eco, master_grid)
          r_eco    <- terra::rast(f_eco_al)
          try(terra::set.cats(r_eco, layer = 1, value = NULL), silent = TRUE)
          
          m_eco_ok <- terra::ifel(r_eco %in% eco_targets, 1L, NA_integer_)
          r_huidig <- terra::mask(r_huidig, m_eco_ok)
          
          unlink(f_eco_al)
          rm(r_eco, m_eco_ok)
          gc()
          message(sprintf("       📊 Na ecoregio-filter: %.2f ha", calc_ha_exact(r_huidig)))
        }
      } else {
        warning(sprintf("       ⚠️ Ecoregio-bestand niet gevonden (%s); filter overgeslagen.", as.character(pad_eco)))
      }
    }
    
    # 2L. Clustering en oppervlakte-drempel per spoor / deelhabitat
    min_opp_ha  <- as.numeric(dplyr::coalesce(
      as.numeric(get_param(rij, "MinPatch_ha", NA)),
      as.numeric(get_param(rij, "Min_Patch_ha", NA)),
      as.numeric(get_param(rij, "MinOpp_ha", NA)),
      0
    ))
    kern_opp_ha   <- as.numeric(get_param(rij, "Kern_MinOpp_ha", NA))
    min_buren_val <- as.numeric(get_param(rij, "Min_Aantal_Buren", NA))
    zoekdist_m    <- as.numeric(get_param(rij, "AfstandBiotopen_m", 0))
    actieradius_m <- as.numeric(get_param(rij, "Actieradius_m", NA))
    
    rand_diepte_val  <- as.numeric(get_param(rij, "Randzone_Diepte_m", NA))
    is_inwaarts_rand <- !is.na(rand_diepte_val) && rand_diepte_val > 0 && grens_moment_val == "na_clustering"
    
    if (!is.na(kern_opp_ha) && kern_opp_ha > min_opp_ha && !is_inwaarts_rand) {
      
      res_cl <- cluster_filter_universeel(
        masker            = r_huidig,
        opp_laag          = r_huidig,
        min_opp_ha        = 0,
        min_totale_opp_ha = NA,
        afstand_m         = zoekdist_m,
        dispersie_m       = 0,
        werkelijk         = TRUE
      )
      cl_raw <- res_cl$clusters
      
      if (!is_leeg(cl_raw)) {
        df_opp <- terra::zonal(r_huidig, cl_raw, fun = "sum", na.rm = TRUE)
        colnames(df_opp) <- c("ID", "Opp_Sum")
        dt_agg <- data.table::as.data.table(df_opp)
        dt_agg[, Opp_ha := Opp_Sum * 0.01]
        
        ids_groot <- dt_agg[!is.na(ID) & Opp_ha >= kern_opp_ha, ID]
        ids_klein <- dt_agg[!is.na(ID) & Opp_ha >= min_opp_ha & Opp_ha < kern_opp_ha, ID]
        
        r_groot <- if (length(ids_groot) > 0) terra::mask(r_huidig, cl_raw %in% as.integer(ids_groot)) else (master_grid * NA_real_)
        r_klein <- if (length(ids_klein) > 0) terra::mask(r_huidig, cl_raw %in% as.integer(ids_klein)) else (master_grid * NA_real_)
        
        heeft_groot <- !is_leeg(r_groot)
        heeft_klein <- !is_leeg(r_klein)
        
        if (!is.na(min_buren_val) && min_buren_val > 0) {
          dispersie_m   <- as.numeric(get_param(rij, "Dispersiecap_m", NA))
          buur_straal_m <- dplyr::coalesce(actieradius_m, dispersie_m, 500)
          message(sprintf("       -> Buurtoets: Patches (%.1f-%.1f ha) toetsen op >= %d patches (>= %.1f ha) binnen %dm...",
                          min_opp_ha, kern_opp_ha, as.integer(min_buren_val), min_opp_ha, as.integer(buur_straal_m)))
          
          if (heeft_klein) {
            ids_alle_min <- c(ids_groot, ids_klein)
            cl_alle_min  <- terra::mask(cl_raw, cl_raw %in% as.integer(ids_alle_min))
            p_poly_all   <- terra::as.polygons(cl_alle_min, dissolve = TRUE)
            
            if (!is.null(p_poly_all) && nrow(p_poly_all) > 0) {
              sf_poly_all   <- sf::st_as_sf(p_poly_all)
              sf_buf_all    <- sf::st_buffer(sf_poly_all, dist = buur_straal_m)
              int_list      <- sf::st_intersects(sf_buf_all, sf_poly_all)
              buren_telling <- lengths(int_list)
              
              meta_kern_idx <- which(buren_telling >= min_buren_val)
              
              if (length(meta_kern_idx) > 0) {
                # Alle patches >= min_opp_ha binnen buur_straal_m van een metapopulatiekern (min_buren_val) zijn geldig
                bereikbare_idx   <- unique(unlist(int_list[meta_kern_idx]))
                alle_ids_vec     <- as.vector(terra::values(p_poly_all)[[1]])
                valide_klein_ids <- intersect(alle_ids_vec[bereikbare_idx], ids_klein)
                
                r_klein_ok <- if (length(valide_klein_ids) > 0) {
                  terra::mask(r_klein, cl_raw %in% as.integer(valide_klein_ids))
                } else {
                  master_grid * NA_real_
                }
              } else {
                r_klein_ok <- master_grid * NA_real_
              }
              rm(sf_poly_all, sf_buf_all, int_list, buren_telling, p_poly_all)
            } else {
              r_klein_ok <- master_grid * NA_real_
            }
            rm(cl_alle_min)
          } else {
            r_klein_ok <- master_grid * NA_real_
          }
          
        } else if (!is.na(actieradius_m) && actieradius_m > 0) {
          message(sprintf("       -> Kern/Satelliet netwerk (Kernen >= %.0f ha, Satellieten >= %.0f ha, actieradius %dm)...",
                          kern_opp_ha, min_opp_ha, as.integer(actieradius_m)))
          if (heeft_klein) {
            sat_bij_groot <- master_grid * NA_real_
            if (heeft_groot) {
              f_g_bin  <- tempfile(fileext = ".tif")
              terra::writeRaster(terra::ifel(!is.na(r_groot) & r_groot > 0, 1, NA), f_g_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
              f_prox_g <- hulp_proximity_gdal(f_g_bin, actieradius_m)
              unlink(f_g_bin)
              
              p_g <- terra::rast(f_prox_g)
              m_g4k <- terra::ifel(!is.na(p_g) & p_g >= 0 & p_g <= actieradius_m, 1, NA)
              unlink(f_prox_g)
              sat_bij_groot <- terra::mask(r_klein, m_g4k)
              rm(p_g, m_g4k)
            }
            
            res_sat_cl <- cluster_filter_universeel(
              masker            = r_klein,
              opp_laag          = r_klein,
              min_opp_ha        = (min_opp_ha * 2),
              min_totale_opp_ha = NA,
              afstand_m         = actieradius_m,
              dispersie_m       = 0,
              werkelijk         = TRUE
            )
            sat_bij_sat <- res_sat_cl$raster
            rm(res_sat_cl)
            
            r_klein_ok <- voeg_toe_geclamped(sat_bij_groot, sat_bij_sat)
            rm(sat_bij_groot, sat_bij_sat)
          } else {
            r_klein_ok <- master_grid * NA_real_
          }
        } else {
          r_klein_ok <- master_grid * NA_real_
        }
        
        r_huidig <- voeg_toe_geclamped(r_groot, r_klein_ok)
        rm(r_groot, r_klein, r_klein_ok, cl_raw, dt_agg)
        
        res_fin_cl <- cluster_filter_universeel(
          masker            = r_huidig,
          opp_laag          = r_huidig,
          min_opp_ha        = min_opp_ha,
          min_totale_opp_ha = NA,
          afstand_m         = zoekdist_m,
          dispersie_m       = 0,
          werkelijk         = TRUE
        )
        r_huidig  <- res_fin_cl$raster
        cl_huidig <- res_fin_cl$clusters
        rm(res_fin_cl)
      } else {
        r_huidig  <- master_grid * NA_real_
        cl_huidig <- master_grid * NA_real_
      }
      message(sprintf("       📊 Na Kern/Netwerk selectie (>= %.0f ha): %.2f ha", min_opp_ha, calc_ha_exact(r_huidig)))
      
    } else if (!is.na(min_opp_ha) && min_opp_ha > 0) {
      message(sprintf("       -> Deelhabitat clusteren (MinPatch): >= %.2f ha (zoekafstand %dm)...", 
                      min_opp_ha, as.integer(zoekdist_m)))
      res_cl <- cluster_filter_universeel(
        masker            = r_huidig,
        opp_laag          = r_huidig,
        min_opp_ha        = min_opp_ha,
        min_totale_opp_ha = NA,
        afstand_m         = zoekdist_m,
        dispersie_m       = 0,
        werkelijk         = TRUE
      )
      r_huidig  <- res_cl$raster
      cl_huidig <- res_cl$clusters
      message(sprintf("       📊 Na MinPatch clustering (>= %.2f ha): %.2f ha", min_opp_ha, calc_ha_exact(r_huidig)))
      
      min_raaklijn_m <- as.numeric(get_param(rij, "Min_Rand_Lengte_m", NA))
      bos_zoom_tp    <- trimws(get_param(rij, "Grens_Met_Type", NA))
      bos_zoom_dist  <- as.numeric(get_param(rij, "Grens_Afstand_m", 20))
      
      if (!is.na(bos_zoom_tp) && bos_zoom_tp %in% names(lagen_ruw)) {
        r_bos_opp <- lagen_ruw[[bos_zoom_tp]]
        
        if (isTRUE(as.logical(get_param(rij, "Gebruik_Groenkaart", FALSE))) && 
            !is.null(paths$groenkaart) && file.exists(paths$groenkaart)) {
          f_gr_al <- hulp_align_gdal(paths$groenkaart, master_grid)
          r_gr    <- terra::rast(f_gr_al)
          r_b1    <- terra::ifel(!is.na(r_gr) & r_gr == 1, 1.0, NA_real_)
          unlink(f_gr_al)
          rm(r_gr)
          r_bos_opp <- voeg_toe_geclamped(r_bos_opp, r_b1)
          rm(r_b1)
        }
        
        if (!is.na(min_raaklijn_m) && min_raaklijn_m > 0 && grens_moment_val == "na_clustering" &&
            !is_leeg(cl_huidig)) {
          message(sprintf("       -> Ecotonale bosrandtoets (>= %dm raaklijn met '%s') & %dm boszoom...", 
                          as.integer(min_raaklijn_m), bos_zoom_tp, as.integer(bos_zoom_dist)))
          
          f_hab_bin <- tempfile(fileext = ".tif")
          terra::writeRaster(terra::ifel(!is.na(r_huidig) & r_huidig > 0, 1, NA), f_hab_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
          f_prox_hab10 <- hulp_proximity_gdal(f_hab_bin, 10)
          unlink(f_hab_bin)
          
          p_hab10 <- terra::rast(f_prox_hab10)
          r_bos_bin <- terra::ifel(!is.na(r_bos_opp) & r_bos_opp > 0, 1, NA)
          r_contact_pixels <- terra::ifel(!is.na(p_hab10) & p_hab10 >= 0 & p_hab10 <= 10 & !is.na(r_bos_bin), 1, NA)
          unlink(f_prox_hab10)
          rm(p_hab10, r_bos_bin)
          
          if (!is_leeg(r_contact_pixels)) {
            f_con_bin <- tempfile(fileext = ".tif")
            f_cl_con  <- tempfile(fileext = ".tif")
            terra::writeRaster(r_contact_pixels, f_con_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
            rm(r_contact_pixels)
            
            e <- terra::ext(master_grid)
            qgisprocess::qgis_run_algorithm(
              "grass:r.clump", input = f_con_bin, output = f_cl_con, title = "cl_contact_lijn",
              "-d" = FALSE, GRASS_REGION_PARAMETER = sprintf("%f,%f,%f,%f", e$xmin, e$xmax, e$ymin, e$ymax),
              GRASS_REGION_CELLSIZE_PARAMETER = 10, .quiet = TRUE
            )
            unlink(f_con_bin)
            
            cl_con <- terra::rast(f_cl_con)
            try(terra::set.cats(cl_con, layer = 1, value = NULL), silent = TRUE)
            try(terra::coltab(cl_con) <- NULL, silent = TRUE)
            
            dt_con <- data.table::as.data.table(terra::freq(cl_con))
            setnames(dt_con, c("layer", "ID", "count"))
            min_cellen_lijn <- ceiling(min_raaklijn_m / 10)
            ok_con_ids <- dt_con[!is.na(ID) & ID > 0 & count >= min_cellen_lijn, ID]
            rm(dt_con)
            
            if (length(ok_con_ids) > 0) {
              r_raaklijn_ok <- terra::ifel(cl_con %in% as.integer(ok_con_ids), 1, NA)
              f_rl_bin <- tempfile(fileext = ".tif")
              terra::writeRaster(r_raaklijn_ok, f_rl_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
              unlink(f_cl_con)
              rm(r_raaklijn_ok, cl_con)
              
              f_prox_rl15 <- hulp_proximity_gdal(f_rl_bin, 15)
              unlink(f_rl_bin)
              p_rl15 <- terra::rast(f_prox_rl15)
              r_touch <- terra::ifel(!is.na(p_rl15) & p_rl15 >= 0 & p_rl15 <= 15, 1, NA)
              cells_touch <- terra::cells(r_touch)
              unlink(f_prox_rl15)
              rm(p_rl15, r_touch)
              
              voldoet_cl_ids <- unique(na.omit(as.integer(cl_huidig[cells_touch][[1]])))
              rm(cells_touch)
              
              if (length(voldoet_cl_ids) > 0) {
                cl_huidig <- terra::ifel(cl_huidig %in% as.integer(voldoet_cl_ids), cl_huidig, NA)
                r_huidig  <- terra::mask(r_huidig, cl_huidig)
                
                f_hab_ok_bin <- tempfile(fileext = ".tif")
                terra::writeRaster(terra::ifel(!is.na(r_huidig), 1, NA), f_hab_ok_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
                f_prox_bz <- hulp_proximity_gdal(f_hab_ok_bin, bos_zoom_dist)
                unlink(f_hab_ok_bin)
                p_bz <- terra::rast(f_prox_bz)
                r_bos_zoom <- terra::mask(r_bos_opp, terra::ifel(!is.na(p_bz) & p_bz >= 0 & p_bz <= bos_zoom_dist, 1, NA))
                unlink(f_prox_bz)
                rm(p_bz)
                
                r_huidig <- voeg_toe_geclamped(r_huidig, r_bos_zoom)
                rm(r_bos_zoom)
              } else {
                cl_huidig <- master_grid * NA_real_
                r_huidig  <- master_grid * NA_real_
              }
            } else {
              unlink(f_cl_con)
              cl_huidig <- master_grid * NA_real_
              r_huidig  <- master_grid * NA_real_
            }
          } else {
            cl_huidig <- master_grid * NA_real_
            r_huidig  <- master_grid * NA_real_
          }
          message(sprintf("       📊 Na bosrandraaklijntoets (>= %dm) + boszoom: %.2f ha", as.integer(min_raaklijn_m), calc_ha_exact(r_huidig)))
          
        } else if (!is.na(bos_zoom_dist) && bos_zoom_dist > 0 && grepl("bos|struweel", tolower(bos_zoom_tp))) {
          message(sprintf("       -> Directe bosrandzoom (%dm) meenemen vanuit '%s'...", as.integer(bos_zoom_dist), bos_zoom_tp))
          f_hab_bin <- tempfile(fileext = ".tif")
          terra::writeRaster(terra::ifel(!is.na(r_huidig) & r_huidig > 0, 1, NA), f_hab_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
          f_prox_bz <- hulp_proximity_gdal(f_hab_bin, bos_zoom_dist)
          unlink(f_hab_bin)
          
          p_bz <- terra::rast(f_prox_bz)
          mask_bz <- terra::ifel(!is.na(p_bz) & p_bz >= 0 & p_bz <= bos_zoom_dist, 1, NA)
          unlink(f_prox_bz)
          rm(p_bz)
          
          r_bos_zoom <- terra::mask(r_bos_opp, mask_bz)
          rm(mask_bz)
          
          r_huidig <- voeg_toe_geclamped(r_huidig, r_bos_zoom)
          rm(r_bos_zoom)
          message(sprintf("       📊 Na toevoegen van %dm bosrandzoom: %.2f ha", as.integer(bos_zoom_dist), calc_ha_exact(r_huidig)))
        }
      }
      
      max_opp_ha_val <- as.numeric(get_param(rij, "MaxOpp_ha", NA))
      if (!is.na(max_opp_ha_val) && max_opp_ha_val > 0 && !is_leeg(cl_huidig)) {
        df_opp_chk <- terra::zonal(r_huidig, cl_huidig, fun = "sum", na.rm = TRUE)
        colnames(df_opp_chk) <- c("ID", "Opp_Sum")
        dt_opp_chk <- data.table::as.data.table(df_opp_chk)
        dt_opp_chk[, Opp_ha := Opp_Sum * 0.01]
        
        ok_max_ids <- dt_opp_chk[!is.na(ID) & Opp_ha <= max_opp_ha_val, ID]
        if (length(ok_max_ids) > 0) {
          cl_huidig <- terra::ifel(cl_huidig %in% as.integer(ok_max_ids), cl_huidig, NA)
          r_huidig  <- terra::mask(r_huidig, cl_huidig)
        } else {
          cl_huidig <- master_grid * NA_real_
          r_huidig  <- master_grid * NA_real_
        }
        rm(df_opp_chk, dt_opp_chk)
      }
      
      min_buren_val <- as.numeric(get_param(rij, "Min_Aantal_Buren", NA))
      if (!is.na(min_buren_val) && min_buren_val > 0 && is.na(kern_opp_ha) && !is_leeg(cl_huidig)) {
        buur_straal_m <- if (!is.na(actieradius_m) && actieradius_m > 0) actieradius_m else 200
        message(sprintf("       -> Dichtheidstoets: Clusters toetsen op >= %d buren binnen %dm...", 
                        as.integer(min_buren_val), as.integer(buur_straal_m)))
        
        p_poly <- terra::as.polygons(cl_huidig, dissolve = TRUE)
        if (!is.null(p_poly) && nrow(p_poly) > 0) {
          sf_poly  <- sf::st_as_sf(p_poly)
          sf_buf   <- sf::st_buffer(sf_poly, dist = buur_straal_m)
          int_list <- sf::st_intersects(sf_buf, sf_poly)
          buren_telling <- lengths(int_list)
          
          ok_idx <- which(buren_telling >= min_buren_val)
          
          if (length(ok_idx) > 0) {
            valide_ids <- as.vector(terra::values(p_poly)[[1]])[ok_idx]
            cl_huidig  <- terra::ifel(cl_huidig %in% as.integer(valide_ids), cl_huidig, NA)
            r_huidig   <- terra::mask(r_huidig, cl_huidig)
          } else {
            cl_huidig <- master_grid * NA_real_
            r_huidig  <- master_grid * NA_real_
          }
          rm(sf_poly, sf_buf, int_list, buren_telling, p_poly)
        } else {
          cl_huidig <- master_grid * NA_real_
          r_huidig  <- master_grid * NA_real_
        }
        message(sprintf("       📊 Na dichtheidstoets (>= %d poelen): %.2f ha", as.integer(min_buren_val), calc_ha_exact(r_huidig)))
      }
      
      if (is_inwaarts_rand && !is_leeg(r_huidig)) {
        message(sprintf("       -> Inwaartse randzone van %dm uitsnijden...", as.integer(rand_diepte_val)))
        r_niet_b_bin <- terra::ifel(is.na(r_huidig), 1, NA)
        f_nb_bin <- tempfile(fileext = ".tif")
        terra::writeRaster(r_niet_b_bin, f_nb_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 0)
        f_prox_nb <- hulp_proximity_gdal(f_nb_bin, rand_diepte_val)
        unlink(f_nb_bin)
        rm(r_niet_b_bin)
        
        p_nb <- terra::rast(f_prox_nb)
        m_inwaarts <- terra::ifel(!is.na(p_nb) & p_nb >= 0 & p_nb <= rand_diepte_val, 1, NA)
        unlink(f_prox_nb)
        rm(p_nb)
        
        r_huidig <- terra::mask(r_huidig, m_inwaarts)
        rm(m_inwaarts)
        message(sprintf("       📊 Na inwaartse randzone-uitsnijding (%dm): %.2f ha", as.integer(rand_diepte_val), calc_ha_exact(r_huidig)))
      }
      
      min_bos_pct <- as.numeric(get_param(rij, "Min_Struweel_Pct", NA))
      if (!is.na(min_bos_pct) && min_bos_pct > 0 && !is_leeg(cl_huidig)) {
        bos_kandidaten <- c(
          if (grepl("ruim", tolower(tp_naam))) "landbiotoop_bwk_ruim_bos" else "landbiotoop_bos",
          "landbiotoop_bwk_ruim_bos", "landbiotoop_bos", "landbiotoop_nabij_bos", 
          "landbiotoop_nabij_bos_bwk", "struweel1"
        )
        bos_zoek_naam <- intersect(bos_kandidaten, names(lagen_ruw))[1]
        r_bos_bron <- if (!is.na(bos_zoek_naam)) lagen_ruw[[bos_zoek_naam]] else NULL
        
        if (!is.null(r_bos_bron)) {
          message(sprintf("       -> Structuurtoets: Cluster toetsen op >= %.0f%% bos/struweel ('%s')...", min_bos_pct, bos_zoek_naam))
          r_bos_clean <- terra::ifel(!is.na(r_bos_bron) & r_bos_bron > 0, 1.0, 0)
          
          df_tot_cl <- terra::zonal(r_huidig, cl_huidig, fun = "sum", na.rm = TRUE)
          colnames(df_tot_cl) <- c("ID", "Tot_Val")
          
          df_bos_cl <- terra::zonal(r_bos_clean, cl_huidig, fun = "sum", na.rm = TRUE)
          colnames(df_bos_cl) <- c("ID", "Bos_Val")
          
          dt_str <- merge(data.table::as.data.table(df_tot_cl), data.table::as.data.table(df_bos_cl), by = "ID", all.x = TRUE)
          dt_str[is.na(Bos_Val), Bos_Val := 0]
          dt_str[, Pct_Bos := (Bos_Val / Tot_Val) * 100]
          
          ok_cl_ids <- if (min_bos_pct == 50) dt_str[!is.na(ID) & Pct_Bos > min_bos_pct, ID] else dt_str[!is.na(ID) & Pct_Bos >= min_bos_pct, ID]
          
          if (length(ok_cl_ids) > 0) {
            cl_huidig <- terra::ifel(cl_huidig %in% as.integer(ok_cl_ids), cl_huidig, NA)
            r_huidig  <- terra::mask(r_huidig, cl_huidig)
          } else {
            cl_huidig <- master_grid * NA_real_
            r_huidig  <- master_grid * NA_real_
          }
          rm(r_bos_clean, df_tot_cl, df_bos_cl, dt_str)
          message(sprintf("       📊 Na structuurtoets (>= %.0f%% bos): %.2f ha", min_bos_pct, calc_ha_exact(r_huidig)))
        }
      }
      
      min_water_pct_val <- as.numeric(get_param(rij, "Min_Water_Pct", NA))
      if (!is.na(min_water_pct_val) && min_water_pct_val > 0 && !is_leeg(cl_huidig)) {
        r_water_bron <- if ("openwater" %in% names(lagen_ruw)) lagen_ruw[["openwater"]] else NULL
        
        if (!is.null(r_water_bron)) {
          message(sprintf("       -> Wateraandeeltoets: Cluster strikt toetsen op >= %.0f%% open water...", min_water_pct_val))
          r_wat_clean <- terra::ifel(is.na(r_water_bron), 0, r_water_bron)
          
          df_foer_cl <- terra::zonal(r_huidig, cl_huidig, fun = "sum", na.rm = TRUE)
          colnames(df_foer_cl) <- c("ID", "Foerageer_ha")
          df_foer_cl$Foerageer_ha <- df_foer_cl$Foerageer_ha * 0.01
          
          df_water_cl <- terra::zonal(r_wat_clean, cl_huidig, fun = "sum", na.rm = TRUE)
          colnames(df_water_cl) <- c("ID", "Water_ha")
          df_water_cl$Water_ha <- df_water_cl$Water_ha * 0.01
          
          dt_wat <- merge(data.table::as.data.table(df_foer_cl), data.table::as.data.table(df_water_cl), by = "ID", all.x = TRUE)
          dt_wat[is.na(Water_ha), Water_ha := 0]
          dt_wat[, Percentage_Water := (Water_ha / Foerageer_ha) * 100]
          
          drempel_pct <- min_water_pct_val
          ok_ids <- dt_wat[!is.na(ID) & Percentage_Water >= drempel_pct, ID]
          
          message(sprintf("       📊 Clusters met >= %.0f%% open water: %d van de %d", 
                          drempel_pct, length(ok_ids), nrow(dt_wat)))
          
          if (length(ok_ids) > 0) {
            cl_huidig <- terra::ifel(cl_huidig %in% as.integer(ok_ids), cl_huidig, NA)
            r_huidig  <- terra::mask(r_huidig, cl_huidig)
          } else {
            cl_huidig <- master_grid * NA_real_
            r_huidig  <- master_grid * NA_real_
          }
          rm(r_wat_clean, df_foer_cl, df_water_cl, dt_wat)
          message(sprintf("       📊 Na zuivere waterfractietoets (>= %.0f%%): %.2f ha", drempel_pct, calc_ha_exact(r_huidig)))
        }
      }
      
      # SCENARIO 2d-grasland: Zuivere procentuele graslandaandeel-toets (Kwartelkoning)
      min_gras_pct_val <- as.numeric(get_param(rij, "Min_Grasland_Pct", NA))
      if (!is.na(min_gras_pct_val) && min_gras_pct_val > 0 && !is_leeg(cl_huidig)) {
        r_gras_bron <- if ("grasland" %in% names(lagen_ruw)) lagen_ruw[["grasland"]] else NULL
        
        if (!is.null(r_gras_bron)) {
          message(sprintf("       -> Graslandaandeeltoets: Cluster strikt toetsen op >= %.0f%% grasland...", min_gras_pct_val))
          r_gras_clean <- terra::ifel(is.na(r_gras_bron), 0, r_gras_bron)
          
          df_bio_cl <- terra::zonal(r_huidig, cl_huidig, fun = "sum", na.rm = TRUE)
          colnames(df_bio_cl) <- c("ID", "Biotoop_ha")
          df_bio_cl$Biotoop_ha <- df_bio_cl$Biotoop_ha * 0.01
          
          df_gras_cl <- terra::zonal(r_gras_clean, cl_huidig, fun = "sum", na.rm = TRUE)
          colnames(df_gras_cl) <- c("ID", "Gras_ha")
          df_gras_cl$Gras_ha <- df_gras_cl$Gras_ha * 0.01
          
          dt_gras <- merge(data.table::as.data.table(df_bio_cl), data.table::as.data.table(df_gras_cl), by = "ID", all.x = TRUE)
          dt_gras[is.na(Gras_ha), Gras_ha := 0]
          dt_gras[, Percentage_Gras := (Gras_ha / Biotoop_ha) * 100]
          
          drempel_gras_pct <- min_gras_pct_val
          ok_gras_ids <- dt_gras[!is.na(ID) & Percentage_Gras >= drempel_gras_pct, ID]
          
          message(sprintf("       📊 Clusters met >= %.0f%% grasland: %d van de %d", 
                          drempel_gras_pct, length(ok_gras_ids), nrow(dt_gras)))
          
          if (length(ok_gras_ids) > 0) {
            cl_huidig <- terra::ifel(cl_huidig %in% as.integer(ok_gras_ids), cl_huidig, NA)
            r_huidig  <- terra::mask(r_huidig, cl_huidig)
          } else {
            cl_huidig <- master_grid * NA_real_
            r_huidig  <- master_grid * NA_real_
          }
          rm(r_gras_clean, df_bio_cl, df_gras_cl, dt_gras)
          message(sprintf("       📊 Na zuivere graslandfractietoets (>= %.0f%%): %.2f ha", drempel_gras_pct, calc_ha_exact(r_huidig)))
        }
      }
      
      # SCENARIO 2e: Inbedding in overkoepelend type
      inbedding_tp <- trimws(get_param(rij, "Inbedding_In_Type", NA))
      if (!is.na(inbedding_tp) && inbedding_tp %in% names(verwerkte_sporen) && !is_leeg(r_huidig)) {
        message(sprintf("       -> Inbeddingstoets: '%s' controleren op inbedding in '%s'...", tp_naam, inbedding_tp))
        r_inbed_mask <- terra::ifel(!is.na(verwerkte_sporen[[inbedding_tp]]) & verwerkte_sporen[[inbedding_tp]] > 0, 1, NA)
        
        inbed_buf_m <- as.numeric(dplyr::coalesce(as.numeric(get_param(rij, "Grens_Afstand_m", NA)), 0))
        if (inbed_buf_m > 0 && !is_leeg(r_inbed_mask)) {
          f_ib_bin  <- tempfile(fileext = ".tif")
          terra::writeRaster(r_inbed_mask, f_ib_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
          f_prox_ib <- hulp_proximity_gdal(f_ib_bin, inbed_buf_m)
          unlink(f_ib_bin)
          p_ib <- terra::rast(f_prox_ib)
          r_inbed_mask <- terra::ifel(!is.na(p_ib) & p_ib >= 0 & p_ib <= inbed_buf_m, 1, NA)
          unlink(f_prox_ib)
          rm(p_ib)
        }
        
        if (!is_leeg(cl_huidig)) {
          df_over <- terra::zonal(r_inbed_mask, cl_huidig, fun = "max", na.rm = TRUE)
          colnames(df_over) <- c("ID", "In_Bed")
          ok_ids <- df_over$ID[!is.na(df_over$In_Bed) & df_over$In_Bed > 0]
          if (length(ok_ids) > 0) {
            cl_huidig <- terra::ifel(cl_huidig %in% as.integer(ok_ids), cl_huidig, NA)
            r_huidig  <- terra::mask(r_huidig, cl_huidig)
          } else {
            cl_huidig <- master_grid * NA_real_
            r_huidig  <- master_grid * NA_real_
          }
          rm(df_over)
        } else {
          r_huidig <- terra::mask(r_huidig, r_inbed_mask)
        }
        rm(r_inbed_mask)
        message(sprintf("       📊 Na inbeddingstoets: %.2f ha", calc_ha_exact(r_huidig)))
      }
      
      gebruik_grb_water <- isTRUE(as.logical(get_param(rij, "Gebruik_GRB_Water", FALSE)))
      kern_drempel_ha   <- as.numeric(get_param(rij, "Kern_MinOpp_ha", NA))
      
      if (gebruik_grb_water && !is.na(kern_drempel_ha) && kern_drempel_ha > 0 && 
          !is_leeg(cl_huidig)) {
        
        message("       -> Landschapskwaliteit toetsen (< 75% bos, < 90% water) & terrestrische kerncheck...")
        
        if (!is.null(paths$grb_water) && file.exists(paths$grb_water)) {
          f_grb_al  <- hulp_align_gdal(paths$grb_water, master_grid)
          r_grb     <- terra::rast(f_grb_al)
          r_wat_bin <- terra::ifel(!is.na(r_grb) & r_grb > 0, 1.0, NA)
          unlink(f_grb_al)
          rm(r_grb)
        } else {
          r_wat_bin <- master_grid * NA_real_
        }
        
        f_gr_al   <- hulp_align_gdal(paths$groenkaart, master_grid)
        r_gr      <- terra::rast(f_gr_al)
        r_bos_bin <- terra::ifel(!is.na(r_gr) & r_gr == 1, 1.0, NA)
        unlink(f_gr_al)
        rm(r_gr)
        
        df_tot <- terra::zonal(r_huidig, cl_huidig, fun = "sum", na.rm = TRUE)
        colnames(df_tot) <- c("Cluster_ID", "Biotoop_Som")
        dt_tot <- data.table::as.data.table(df_tot)
        dt_tot[, Totaal_ha := Biotoop_Som * 0.01]
        
        bos_in_zone <- terra::mask(r_bos_bin, cl_huidig)
        df_bos <- terra::zonal(bos_in_zone, cl_huidig, fun = "sum", na.rm = TRUE)
        colnames(df_bos) <- c("Cluster_ID", "Bos_Som")
        dt_bos <- data.table::as.data.table(df_bos)
        dt_bos[, Bos_ha := Bos_Som * 0.01]
        rm(r_bos_bin, bos_in_zone, df_bos)
        
        wat_in_zone <- terra::mask(r_wat_bin, cl_huidig)
        df_wat <- terra::zonal(wat_in_zone, cl_huidig, fun = "sum", na.rm = TRUE)
        colnames(df_wat) <- c("Cluster_ID", "Water_Som")
        dt_wat <- data.table::as.data.table(df_wat)
        dt_wat[, Water_ha := Water_Som * 0.01]
        rm(wat_in_zone, df_wat)
        
        dt_merged <- merge(dt_tot, dt_bos, by = "Cluster_ID", all.x = TRUE)
        dt_merged <- merge(dt_merged, dt_wat, by = "Cluster_ID", all.x = TRUE)
        dt_merged[is.na(Bos_ha), Bos_ha := 0]
        dt_merged[is.na(Water_ha), Water_ha := 0]
        dt_merged[, Pct_Bos   := (Bos_ha / Totaal_ha) * 100]
        dt_merged[, Pct_Water := (Water_ha / Totaal_ha) * 100]
        
        ids_kwaliteit_ok <- dt_merged[!is.na(Cluster_ID) & Pct_Bos < 75 & Pct_Water < 90, Cluster_ID]
        message(sprintf("       📊 Clusters goedgekeurd op openheid (<75%% bos, <90%% water): %d van de %d", 
                        length(ids_kwaliteit_ok), nrow(dt_merged)))
        rm(dt_tot, dt_bos, dt_wat, dt_merged)
        
        if (length(ids_kwaliteit_ok) > 0) {
          r_terr_biotoop <- terra::mask(r_huidig, r_wat_bin, inverse = TRUE)
          rm(r_wat_bin)
          
          res_terr <- cluster_filter_universeel(
            masker            = r_terr_biotoop,
            opp_laag          = r_terr_biotoop,
            min_opp_ha        = kern_drempel_ha,
            min_totale_opp_ha = NA,
            afstand_m         = 0,
            dispersie_m       = 0,
            werkelijk         = TRUE
          )
          cl_terr_10ha <- res_terr$clusters
          rm(res_terr, r_terr_biotoop)
          
          if (!is_leeg(cl_terr_10ha)) {
            ct <- terra::crosstab(c(cl_huidig, cl_terr_10ha), long = TRUE)
            colnames(ct) <- c("Cluster_ID", "Terr_ID", "n")
            dt_ct <- data.table::as.data.table(ct)[!is.na(Cluster_ID) & !is.na(Terr_ID) & n > 0]
            valide_cluster_ids <- unique(dt_ct[Cluster_ID %in% ids_kwaliteit_ok, Cluster_ID])
            rm(ct, dt_ct, cl_terr_10ha)
          } else {
            valide_cluster_ids <- integer(0)
            rm(cl_terr_10ha)
          }
          
          if (length(valide_cluster_ids) > 0) {
            cl_huidig <- terra::ifel(cl_huidig %in% as.integer(valide_cluster_ids), cl_huidig, NA)
            r_huidig  <- terra::mask(r_huidig, cl_huidig)
          } else {
            cl_huidig <- master_grid * NA_real_
            r_huidig  <- master_grid * NA_real_
          }
        } else {
          rm(r_wat_bin)
          cl_huidig <- master_grid * NA_real_
          r_huidig  <- master_grid * NA_real_
        }
        message(sprintf("       📊 Na terrestrische kerncheck (>= %.0f ha land): %.2f ha", kern_drempel_ha, calc_ha_exact(r_huidig)))
      }
      
    } else {
      cl_huidig <- master_grid * NA_real_
    }
    
    # SCENARIO 2g: Dynamische kustweidevogel-toets
    heeft_slik <- "slik_schor" %in% names(lagen_ruw)
    heeft_ven  <- "vennen_water" %in% names(lagen_ruw)
    
    if ((heeft_slik || heeft_ven) && !is_leeg(cl_huidig)) {
      message("       -> Dynamische oppervlaktetoets (5-10 ha) en vennenfilter (< 95%) berekenen...")
      
      r_slik_clean <- if (heeft_slik) terra::ifel(is.na(lagen_ruw[["slik_schor"]]), 0, lagen_ruw[["slik_schor"]]) else (master_grid * 0)
      r_ven_clean  <- if (heeft_ven)  terra::ifel(is.na(lagen_ruw[["vennen_water"]]), 0, lagen_ruw[["vennen_water"]]) else (master_grid * 0)
      
      df_tot_cl <- terra::zonal(r_huidig, cl_huidig, fun = "sum", na.rm = TRUE)
      colnames(df_tot_cl) <- c("ID", "Totaal_Val")
      df_tot_cl$Totaal_ha <- df_tot_cl$Totaal_Val * 0.01
      
      df_slik_cl <- terra::zonal(r_slik_clean, cl_huidig, fun = "sum", na.rm = TRUE)
      colnames(df_slik_cl) <- c("ID", "Slik_Val")
      df_slik_cl$Slik_ha <- df_slik_cl$Slik_Val * 0.01
      
      df_ven_cl <- terra::zonal(r_ven_clean, cl_huidig, fun = "sum", na.rm = TRUE)
      colnames(df_ven_cl) <- c("ID", "Ven_Val")
      df_ven_cl$Ven_ha <- df_ven_cl$Ven_Val * 0.01
      
      dt_dyn <- merge(data.table::as.data.table(df_tot_cl), data.table::as.data.table(df_slik_cl), by = "ID", all.x = TRUE)
      dt_dyn <- merge(dt_dyn, data.table::as.data.table(df_ven_cl), by = "ID", all.x = TRUE)
      dt_dyn[is.na(Slik_ha), Slik_ha := 0]
      dt_dyn[is.na(Ven_ha), Ven_ha := 0]
      
      dt_dyn[, Slik_Fractie := Slik_ha / Totaal_ha]
      dt_dyn[, MinOpp_Req   := 10.0 - (5.0 * pmin(1.0, Slik_Fractie))]
      dt_dyn[, Ven_Pct      := (Ven_ha / Totaal_ha) * 100]
      
      ok_dyn_ids <- dt_dyn[!is.na(ID) & (Totaal_ha >= MinOpp_Req) & (Ven_Pct < 95), ID]
      message(sprintf("       📊 Clusters goedgekeurd na dynamische opp/vennentoets: %d van de %d", 
                      length(ok_dyn_ids), nrow(dt_dyn)))
      
      if (length(ok_dyn_ids) > 0) {
        cl_huidig <- terra::ifel(cl_huidig %in% as.integer(ok_dyn_ids), cl_huidig, NA)
        r_huidig  <- terra::mask(r_huidig, cl_huidig)
      } else {
        cl_huidig <- master_grid * NA_real_
        r_huidig  <- master_grid * NA_real_
      }
      rm(r_slik_clean, r_ven_clean, df_tot_cl, df_slik_cl, df_ven_cl, dt_dyn)
    }
    
    # 2M. Lineaire KLE-toets
    min_kle_m <- as.numeric(get_param(rij, "Min_KLE_Lengte_m", 0))
    if (!is.na(min_kle_m) && min_kle_m > 0 && !is_leeg(cl_huidig)) {
      min_kle_cellen <- round(min_kle_m / 10)
      message(sprintf("       -> KLE-toets (min. %dm / %d cellen houtkanten/struwelen)...", as.integer(min_kle_m), as.integer(min_kle_cellen)))
      f_gr_al <- hulp_align_gdal(paths$groenkaart, master_grid)
      r_gr    <- terra::rast(f_gr_al)
      r_kle_b <- terra::ifel(!is.na(r_gr) & (r_gr %in% c(1, 2)), 1, 0)
      unlink(f_gr_al)
      
      df_kle <- terra::zonal(r_kle_b, cl_huidig, fun = "sum", na.rm = TRUE)
      colnames(df_kle) <- c("ID", "KLE_pixels")
      dt_kle <- data.table::as.data.table(df_kle)
      ok_ids <- dt_kle[!is.na(ID) & ID > 0 & KLE_pixels >= min_kle_cellen, ID]
      
      if (length(ok_ids) > 0) {
        cl_huidig <- terra::ifel(cl_huidig %in% as.integer(ok_ids), cl_huidig, NA)
        r_huidig  <- terra::mask(r_huidig, cl_huidig)
      } else {
        cl_huidig <- master_grid * NA_real_
        r_huidig  <- master_grid * NA_real_
      }
      message(sprintf("       📊 Na KLE-toets (>= %dm): %.2f ha", as.integer(min_kle_m), calc_ha_exact(r_huidig)))
    }
    
    verwerkte_sporen[[tp_naam]] <- r_huidig
    verwerkte_sporen[[paste0(tp_naam, ":", huidig_spoor)]] <- r_huidig
    cl_sporen[[tp_naam]]        <- cl_huidig
  }
  
  # ----------------------------------------------------------------------------
  # STAP 3: KOPPELINGEN TUSSEN TYPES
  # ----------------------------------------------------------------------------
  kol_rand_l <- intersect(c("Min_Rand_Lengte_m", "Min_Oever_Lengte_m"), colnames(params_soort))[1]
  waarden_rand_l <- if (!is.na(kol_rand_l)) suppressWarnings(as.numeric(params_soort[[kol_rand_l]])) else numeric(0)
  min_rand_m <- if (any(!is.na(waarden_rand_l))) max(waarden_rand_l, na.rm = TRUE) else NA_real_
  heeft_rand_lengte_toets <- isTRUE(!is.na(min_rand_m) && min_rand_m > 0)
  
  waarden_diepte <- if ("Randzone_Diepte_m" %in% colnames(params_soort)) suppressWarnings(as.numeric(params_soort$Randzone_Diepte_m)) else numeric(0)
  rand_diepte_m <- if (any(!is.na(waarden_diepte))) max(waarden_diepte, na.rm = TRUE) else NA_real_
  
  grens_momenten_v <- if ("Grens_Moment" %in% colnames(params_soort)) tolower(trimws(as.character(params_soort$Grens_Moment))) else character(0)
  heeft_ecotoon_archetype <- isTRUE(!is.na(rand_diepte_m) && rand_diepte_m > 0 && !("na_clustering" %in% grens_momenten_v))
  
  wed_vals <- if ("Gebruik_Wederzijds" %in% colnames(params_soort)) as.logical(params_soort$Gebruik_Wederzijds) else logical(0)
  is_wederzijds <- isTRUE(any(wed_vals == TRUE, na.rm = TRUE))
  
  if (heeft_rand_lengte_toets && length(verwerkte_sporen) >= 2) {
    rij_rl     <- params_soort %>% 
      dplyr::filter(!is.na(.data[[kol_rand_l]]) & suppressWarnings(as.numeric(.data[[kol_rand_l]])) > 0) %>% 
      dplyr::slice(1)
    open_tp    <- trimws(rij_rl$Type[1])
    bos_tp     <- trimws(get_param(rij_rl, "Grens_Met_Type", NA))
    groep_buf  <- as.numeric(get_param(rij_rl, "Grens_Afstand_m", 25))
    min_pixels <- round(min_rand_m / 10)
    
    if (!is.na(bos_tp) && open_tp %in% names(verwerkte_sporen) && bos_tp %in% names(verwerkte_sporen)) {
      message(sprintf("\n    📏 Bosrandlengte-toets: '%s' toetsen op >= %dm (%d pixels) bosrand met '%s'...", 
                      open_tp, as.integer(min_rand_m), as.integer(min_pixels), bos_tp))
      
      r_open <- verwerkte_sporen[[open_tp]]
      r_bos  <- verwerkte_sporen[[bos_tp]]
      
      if (!is_leeg(r_open) && 
          !is_leeg(r_bos)) {
        r_bos_bin <- terra::ifel(!is.na(r_bos) & r_bos > 0, 1L, 0L)
        w_3x3     <- matrix(1, nrow = 3, ncol = 3)
        r_b_focal <- terra::focal(r_bos_bin, w = w_3x3, fun = "max", na.rm = TRUE)
        rm(r_bos_bin)
        
        r_open_bin <- terra::ifel(!is.na(r_open) & r_open > 0, 1L, 0L)
        r_rand_bin <- terra::ifel(r_open_bin == 1L & r_b_focal == 1L, 1L, 0L)
        rm(r_b_focal)
        
        r_open_alleen <- terra::ifel(r_open_bin == 1L, 1L, NA_integer_)
        rm(r_open_bin)
        
        f_op_bin  <- tempfile(fileext = ".tif")
        terra::writeRaster(r_open_alleen, f_op_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
        f_prox_op <- hulp_proximity_gdal(f_op_bin, groep_buf)
        unlink(f_op_bin)
        
        p_op <- terra::rast(f_prox_op)
        mask_op_buf <- terra::ifel(!is.na(p_op) & p_op >= 0 & p_op <= groep_buf, 1L, NA_integer_)
        unlink(f_prox_op)
        rm(p_op)
        
        f_buf_bin <- tempfile(fileext = ".tif")
        f_cl_c    <- tempfile(fileext = ".tif")
        terra::writeRaster(mask_op_buf, f_buf_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
        rm(mask_op_buf)
        
        e <- terra::ext(master_grid)
        qgisprocess::qgis_run_algorithm(
          "grass:r.clump", input = f_buf_bin, output = f_cl_c, title = paste0("cl_", open_tp),
          "-d" = FALSE, GRASS_REGION_PARAMETER = sprintf("%f,%f,%f,%f", e$xmin, e$xmax, e$ymin, e$ymax),
          GRASS_REGION_CELLSIZE_PARAMETER = 10, .quiet = TRUE
        )
        unlink(f_buf_bin)
        
        cl_comp <- terra::rast(f_cl_c)
        try(terra::set.cats(cl_comp, layer = 1, value = NULL), silent = TRUE)
        try(terra::coltab(cl_comp) <- NULL, silent = TRUE)
        cl_comp_open <- terra::mask(cl_comp, r_open_alleen)
        
        df_rand_som <- terra::zonal(r_rand_bin, cl_comp_open, fun = "sum", na.rm = TRUE)
        colnames(df_rand_som) <- c("ID", "Totaal_Rand")
        dt_stats <- data.table::as.data.table(df_rand_som)
        ok_ids   <- dt_stats[!is.na(ID) & ID > 0 & Totaal_Rand >= min_pixels, ID]
        rm(df_rand_som, dt_stats, r_rand_bin, r_open_alleen, cl_comp)
        
        message(sprintf("       📊 Complexen met >= %dm bosrand: %d", as.integer(min_rand_m), length(ok_ids)))
        if (length(ok_ids) > 0) {
          r_open_ok <- terra::mask(r_open, cl_comp_open %in% as.integer(ok_ids))
          
          # Boszoom (groep_buf, bv. 25m) vanuit het gekoppelde bos meenemen
          f_hab_ok <- tempfile(fileext = ".tif")
          terra::writeRaster(terra::ifel(!is.na(r_open_ok), 1L, NA_integer_), f_hab_ok, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
          f_prox_bz <- hulp_proximity_gdal(f_hab_ok, groep_buf)
          unlink(f_hab_ok)
          p_bz <- terra::rast(f_prox_bz)
          r_bos_zoom <- terra::mask(r_bos, terra::ifel(!is.na(p_bz) & p_bz >= 0 & p_bz <= groep_buf, 1L, NA_integer_))
          unlink(f_prox_bz)
          rm(p_bz)
          
          verwerkte_sporen[[open_tp]] <- voeg_toe_geclamped(r_open_ok, r_bos_zoom)
          rm(r_open_ok, r_bos_zoom)
        } else {
          verwerkte_sporen[[open_tp]] <- master_grid * NA_real_
        }
        unlink(f_cl_c)
        rm(cl_comp_open)
      } else {
        verwerkte_sporen[[open_tp]] <- master_grid * NA_real_
      }
      
      verwerkte_sporen[[bos_tp]] <- master_grid * NA_real_
      message(sprintf("       📊 Open habitat na bosrandlengtetoets: %.2f ha", calc_ha_exact(verwerkte_sporen[[open_tp]])))
    }
    
  } else if (heeft_ecotoon_archetype && length(verwerkte_sporen) >= 2) {
    message(sprintf("\n    🌿 Ecotoon-analyse: Scherpe contactzone bepalen (<= %dm) en bufferen...", as.integer(rand_diepte_m)))
    
    tp_namen <- names(verwerkte_sporen)
    tp_a <- tp_namen[1]
    tp_b <- tp_namen[2]
    
    laag_a <- verwerkte_sporen[[tp_a]]
    laag_b <- verwerkte_sporen[[tp_b]]
    
    heeft_a <- !is_leeg(laag_a)
    heeft_b <- !is_leeg(laag_b)
    
    if (heeft_a && heeft_b) {
      f_b_bin  <- tempfile(fileext = ".tif")
      terra::writeRaster(terra::ifel(!is.na(laag_b) & laag_b > 0, 1, NA), f_b_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
      f_prox_b <- hulp_proximity_gdal(f_b_bin, rand_diepte_m)
      unlink(f_b_bin)
      
      p_b <- terra::rast(f_prox_b)
      mask_a_rand <- terra::ifel(!is.na(p_b) & p_b >= 0 & p_b <= rand_diepte_m, 1, NA)
      unlink(f_prox_b)
      r_rand_a <- terra::mask(laag_a, mask_a_rand)
      rm(p_b, mask_a_rand)
      
      f_a_bin  <- tempfile(fileext = ".tif")
      terra::writeRaster(terra::ifel(!is.na(laag_a) & laag_a > 0, 1, NA), f_a_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
      f_prox_a <- hulp_proximity_gdal(f_a_bin, rand_diepte_m)
      unlink(f_a_bin)
      
      p_a <- terra::rast(f_prox_a)
      mask_b_rand <- terra::ifel(!is.na(p_a) & p_a >= 0 & p_a <= rand_diepte_m, 1, NA)
      unlink(f_prox_a)
      r_rand_b <- terra::mask(laag_b, mask_b_rand)
      rm(p_a, mask_b_rand)
      
      r_rand_tot <- voeg_toe_geclamped(r_rand_a, r_rand_b)
      rm(r_rand_a, r_rand_b)
      
      drempel_rand_ha <- suppressWarnings(min(as.numeric(params_soort$Kern_MinOpp_ha), na.rm = TRUE))
      if (is.infinite(drempel_rand_ha) || is.na(drempel_rand_ha) || drempel_rand_ha <= 0) drempel_rand_ha <- 0.25
      
      message(sprintf("       -> Contactzone filteren op aaneengesloten lengte (>= %.2f ha, 0m netwerk)...", drempel_rand_ha))
      res_rand <- cluster_filter_universeel(
        masker            = r_rand_tot,
        opp_laag          = r_rand_tot,
        min_opp_ha        = drempel_rand_ha,
        min_totale_opp_ha = NA,
        afstand_m         = 0,
        dispersie_m       = 0,
        werkelijk         = TRUE
      )
      r_rand_val <- res_rand$raster
      rm(res_rand, r_rand_tot)
      
      ha_rand_val <- calc_ha_exact(r_rand_val)
      message(sprintf("       📊 Goedgekeurde ecotoon-randen: %.2f ha", ha_rand_val))
      
      buf_leef_m <- suppressWarnings(max(as.numeric(params_soort$Grens_Afstand_m), na.rm = TRUE))
      if (is.infinite(buf_leef_m) || is.na(buf_leef_m) || buf_leef_m <= 0) buf_leef_m <- 100
      
      if (ha_rand_val > 0) {
        f_r_bin <- tempfile(fileext = ".tif")
        terra::writeRaster(terra::ifel(!is.na(r_rand_val) & r_rand_val > 0, 1, NA), f_r_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
        f_prox_r <- hulp_proximity_gdal(f_r_bin, buf_leef_m)
        unlink(f_r_bin)
        
        p_r <- terra::rast(f_prox_r)
        mask_leef_buf <- terra::ifel(!is.na(p_r) & p_r >= 0 & p_r <= buf_leef_m, 1, NA)
        unlink(f_prox_r)
        rm(p_r)
        
        for (tp in tp_namen) {
          verwerkte_sporen[[tp]] <- terra::mask(verwerkte_sporen[[tp]], mask_leef_buf)
        }
        verwerkte_sporen[["bosrand_ecotoon"]] <- r_rand_val
        rm(mask_leef_buf, r_rand_val)
      } else {
        for (tp in tp_namen) {
          verwerkte_sporen[[tp]] <- master_grid * NA_real_
        }
      }
    }
  
} else if (!isTRUE(tolower(bestandsnaam) %in% c("kempenseheidelibel", "sympetrumdepressiusculum"))) {
  if (is_wederzijds && length(verwerkte_sporen) >= 2) {
      rij_a <- params_soort %>% dplyr::filter(!is.na(Gebruik_Wederzijds) & as.logical(Gebruik_Wederzijds) & Spoor == "A") %>% dplyr::slice(1)
      rij_b <- params_soort %>% dplyr::filter(!is.na(Gebruik_Wederzijds) & as.logical(Gebruik_Wederzijds) & Spoor == "B") %>% dplyr::slice(1)
      
      tp_a <- if (nrow(rij_a) > 0) trimws(rij_a$Type[1]) else names(verwerkte_sporen)[1]
      tp_b <- if (nrow(rij_b) > 0) trimws(rij_b$Type[1]) else setdiff(names(verwerkte_sporen), tp_a)[1]
      
      dist_a <- if (nrow(rij_a) > 0) as.numeric(dplyr::coalesce(as.numeric(get_param(rij_a, "Grens_Afstand_m", NA)),
                                                                as.numeric(get_param(rij_a, "Wederzijdse_Afstand_m", NA)), 50)) else 50
      dist_b <- if (nrow(rij_b) > 0) as.numeric(dplyr::coalesce(as.numeric(get_param(rij_b, "Grens_Afstand_m", NA)),
                                                                as.numeric(get_param(rij_b, "Wederzijdse_Afstand_m", NA)), dist_a)) else dist_a
      
      if (length(tp_a) == 1 && length(tp_b) == 1 && !is.na(tp_a) && !is.na(tp_b) &&
          tp_a %in% names(verwerkte_sporen) && tp_b %in% names(verwerkte_sporen)) {
        message(sprintf("\n    🔗 Wederzijdse koppeling: '%s' (<= %dm) <-> '%s' (<= %dm)...", 
                        tp_a, as.integer(dist_a), tp_b, as.integer(dist_b)))
        
        res_wed <- hulp_wederzijdse_koppeling(
          laag_a      = verwerkte_sporen[[tp_a]], 
          laag_b      = verwerkte_sporen[[tp_b]], 
          dist_a_m    = dist_a, 
          dist_b_m    = dist_b, 
          master_grid = master_grid
        )
        verwerkte_sporen[[tp_a]] <- res_wed$laag_a_gekoppeld
        verwerkte_sporen[[tp_b]] <- res_wed$laag_b_gekoppeld
        message(sprintf("       📊 Na wederzijdse koppeling: '%s' = %.2f ha | '%s' = %.2f ha", 
                        tp_a, calc_ha_exact(verwerkte_sporen[[tp_a]]), 
                        tp_b, calc_ha_exact(verwerkte_sporen[[tp_b]])))
      }
    }
    
    for (i in 1:nrow(params_soort)) {
      rij          <- params_soort[i, ]
      doel_tp      <- trimws(rij$Type)
      huidig_sp    <- trimws(as.character(get_param(rij, "Spoor", "")))
      doel_key     <- if (huidig_sp != "" && paste0(doel_tp, ":", huidig_sp) %in% names(verwerkte_sporen)) {
        paste0(doel_tp, ":", huidig_sp)
      } else {
        doel_tp
      }
      
      kern_tp_raw  <- get_param(rij, "Grens_Met_Type", NA)
      grens_m_raw  <- get_param(rij, "Grens_Afstand_m", NA)
      grens_mom    <- tolower(trimws(get_param(rij, "Grens_Moment", "")))
      
      if (!is.na(kern_tp_raw) && !is.na(grens_m_raw) && grens_mom != "" && grens_mom != "voor_clustering" && doel_key %in% names(verwerkte_sporen)) {
        kern_types <- trimws(strsplit(as.character(kern_tp_raw), "[,;]+")[[1]])
        afstanden  <- as.numeric(trimws(strsplit(as.character(grens_m_raw), "[,;]+")[[1]]))
        if (length(afstanden) < length(kern_types)) afstanden <- rep(afstanden[1], length(kern_types))
        
        for (k in seq_along(kern_types)) {
          kern_tp <- kern_types[k]
          grens_m <- afstanden[k]
          
          if (kern_tp %in% names(verwerkte_sporen) && grens_m > 0) {
            r_kern_bin <- terra::ifel(!is.na(verwerkte_sporen[[kern_tp]]) & verwerkte_sporen[[kern_tp]] > 0, 1, NA)
            
            if (!is_leeg(r_kern_bin)) {
              f_k_bin  <- tempfile(fileext = ".tif")
              terra::writeRaster(r_kern_bin, f_k_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
              f_prox_k <- hulp_proximity_gdal(f_k_bin, grens_m)
              unlink(f_k_bin)
              
              p_k <- terra::rast(f_prox_k)
              m_grens <- terra::ifel(!is.na(p_k) & p_k >= 0 & p_k <= grens_m, 1, NA)
              unlink(f_prox_k)
              
              if (grens_mom == "cluster_select" && doel_tp %in% names(cl_sporen) && 
                  !is_leeg(cl_sporen[[doel_tp]])) {
                message(sprintf("\n    🎯 Clusterselectie: Hele clusters van '%s' behouden die binnen %dm van '%s' liggen...", 
                                doel_key, as.integer(grens_m), kern_tp))
                cl_d <- cl_sporen[[doel_tp]]
                cl_contact <- terra::mask(cl_d, m_grens)
                contact_ids <- unique(na.omit(as.integer(terra::values(cl_contact, mat = FALSE))))
                rm(cl_contact)
                
                if (length(contact_ids) > 0) {
                  cl_sporen[[doel_tp]] <- terra::ifel(cl_d %in% as.integer(contact_ids), cl_d, NA)
                  verwerkte_sporen[[doel_key]] <- terra::mask(verwerkte_sporen[[doel_key]], cl_sporen[[doel_tp]])
                } else {
                  verwerkte_sporen[[doel_key]] <- master_grid * NA_real_
                  cl_sporen[[doel_tp]]         <- master_grid * NA_real_
                }
                rm(cl_d)
                
              } else {
                message(sprintf("\n    🔗 Koppeling na clustering: '%s' beperken tot binnen %dm van '%s'...", 
                                doel_key, as.integer(grens_m), kern_tp))
                verwerkte_sporen[[doel_key]] <- terra::mask(verwerkte_sporen[[doel_key]], m_grens)
              }
              verwerkte_sporen[[doel_tp]] <- verwerkte_sporen[[doel_key]]
              rm(p_k, m_grens)
            } else {
              message(sprintf("       ⚠️ Referentielaag '%s' is leeg! Doel '%s' wordt op NA gezet.", kern_tp, doel_key))
              verwerkte_sporen[[doel_key]] <- master_grid * NA_real_
              verwerkte_sporen[[doel_tp]]  <- master_grid * NA_real_
            }
            message(sprintf("       📊 Gekoppeld habitat '%s': %.2f ha", doel_key, calc_ha_exact(verwerkte_sporen[[doel_key]])))
          }
        }
      }
    }
  }
  
  # ----------------------------------------------------------------------------
  # STAP 4: BRUTO COMPLEX SAMENVOEGEN, METACLUSTER & ACTIERADIUSTOETS
  # ----------------------------------------------------------------------------
  message("\n    🧩 Stap 4: Samenstellen van bruto complex (optelling van fracties)...")
  
  alle_keys <- names(verwerkte_sporen)
  basis_types <- alle_keys[!grepl(":[A-Z0-9]+$", alle_keys)]
  types_te_fuseren <- unlist(lapply(basis_types, function(bt) {
    spoor_variants <- alle_keys[startsWith(alle_keys, paste0(bt, ":"))]
    if (length(spoor_variants) > 1) spoor_variants else bt
  }))
  context_types <- c("landbiotoop_bwk_ruim", "landbiotoop_bwk_ruim_bos", "land_ruim", 
                     "landbiotoop_rest_bwk", "landbiotoop_rest", "land_rest", "openwater")
  types_te_fuseren <- setdiff(types_te_fuseren, context_types)
  
  final_opp <- NULL
  for (nm in types_te_fuseren) {
    r_sp <- verwerkte_sporen[[nm]]
    if (!is_leeg(r_sp)) {
      final_opp <- voeg_toe_geclamped(final_opp, r_sp)
    }
  }
  
  if (is.null(final_opp)) {
    final_opp <- master_grid * NA_real_
  }
  
  ha_bruto_complex <- calc_ha_exact(final_opp)
  message(sprintf("       📊 Bruto samengesteld complex vóór metacluster: %.2f ha", ha_bruto_complex))
  
  # 4B. Toets op MinComplex_ha
  kol_c <- intersect(c("MinComplex_ha", "Min_Complex_ha", "MinTotaleOpp_ha"), colnames(params_soort))[1]
  waarden_complex <- if (!is.na(kol_c)) params_soort[[kol_c]] else NA
  min_tot_ha <- suppressWarnings(max(as.numeric(waarden_complex), na.rm = TRUE))
  
  zoek_tot_m <- suppressWarnings(max(as.numeric(params_soort$AfstandBiotopen_m), na.rm = TRUE))
  if (is.infinite(zoek_tot_m) || is.na(zoek_tot_m) || zoek_tot_m <= 0) zoek_tot_m <- 50
  
  cl_opp <- NULL
  
  if (!is.infinite(min_tot_ha) && !is.na(min_tot_ha) && min_tot_ha > 0) {
    message(sprintf("    🏝️ Metacluster-toets (MinComplex): Complex >= %.2f ha (zoekafstand %dm)...", 
                    min_tot_ha, as.integer(zoek_tot_m)))
    
    res_tot <- cluster_filter_universeel(
      masker            = final_opp, 
      opp_laag          = final_opp, 
      min_opp_ha        = min_tot_ha, 
      min_totale_opp_ha = NA, 
      afstand_m         = zoek_tot_m, 
      dispersie_m       = 0, 
      werkelijk         = TRUE
    )
    final_opp <- res_tot$raster
    cl_opp    <- res_tot$clusters
    ha_na_complex <- calc_ha_exact(final_opp)
    message(sprintf("       📊 Na MinComplex clustering (>= %.1f ha): %.2f ha", min_tot_ha, ha_na_complex))
    
    # 4C. Actieradiuskoppeling
    actieradius_comp <- suppressWarnings(max(as.numeric(params_soort$Actieradius_m), na.rm = TRUE))
    
    if (!is.infinite(actieradius_comp) && !is.na(actieradius_comp) && actieradius_comp > 0 && 
        "voortplanting_bwk" %in% names(verwerkte_sporen) && ha_na_complex > 0) {
      
      message(sprintf("    🎯 Actieradius-toets: Broedbos valideren binnen %d m van complex...", as.integer(actieradius_comp)))
      
      r_comp_bin <- terra::ifel(!is.na(final_opp) & final_opp > 0, 1, NA)
      f_c_bin    <- tempfile(fileext = ".tif")
      terra::writeRaster(r_comp_bin, f_c_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
      f_prox_c   <- hulp_proximity_gdal(f_c_bin, actieradius_comp)
      unlink(f_c_bin)
      rm(r_comp_bin)
      
      p_c <- terra::rast(f_prox_c)
      m_disp <- terra::ifel(!is.na(p_c) & p_c >= 0 & p_c <= actieradius_comp, 1, NA)
      unlink(f_prox_c)
      rm(p_c)
      
      verwerkte_sporen[["voortplanting_bwk"]] <- terra::mask(verwerkte_sporen[["voortplanting_bwk"]], m_disp)
      rm(m_disp)
      
      ha_bos_actieradius <- calc_ha_exact(verwerkte_sporen[["voortplanting_bwk"]])
      message(sprintf("       📊 Broedbos binnen %d m actieradius: %.2f ha", as.integer(actieradius_comp), ha_bos_actieradius))
      
      if ("foerageer_open" %in% names(verwerkte_sporen)) {
        verwerkte_sporen[["foerageer_open"]] <- terra::mask(verwerkte_sporen[["foerageer_open"]], final_opp)
      }
      
      final_opp <- NULL
      for (nm in types_te_fuseren) {
        r_sp <- verwerkte_sporen[[nm]]
        if (!is_leeg(r_sp)) {
          final_opp <- voeg_toe_geclamped(final_opp, r_sp)
        }
      }
      if (is.null(final_opp)) final_opp <- master_grid * NA_real_
      
      # Reset cl_opp zodat de geactualiseerde final_opp hieronder in 4D opnieuw geclusterd wordt
      cl_opp <- NULL
    }
  }
  
  # 4D. Finale cluster-ID's opbouwen (indien nog niet gezet door MinComplex)
  if (is.null(cl_opp)) {
    if (length(types_te_fuseren) == 1 && types_te_fuseren[1] %in% names(cl_sporen) &&
        !is_leeg(cl_sporen[[types_te_fuseren[1]]])) {
      cl_opp <- cl_sporen[[types_te_fuseren[1]]]
    } else if (!is_leeg(final_opp)) {
      r_fin_bin <- terra::ifel(!is.na(final_opp) & final_opp > 0, 1, NA)
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
      rm(r_fin_bin)
      gc()
    } else {
      cl_opp <- master_grid * NA_real_
    }
  }
  
  # ----------------------------------------------------------------------------
  # STAP 5: MATRIXFILTER OP INTENSIEF GRASLAND (^hp < 90%)
  # ----------------------------------------------------------------------------
  hp_waarden <- suppressWarnings(as.numeric(params_soort$Max_Intensief_Grasland_Pct))
  max_hp_pct <- if (any(!is.na(hp_waarden))) max(hp_waarden, na.rm = TRUE) else NA_real_
  
  if (!is.na(max_hp_pct) && max_hp_pct > 0 && !is_leeg(cl_opp)) {
    message(sprintf("    📊 Matrixfilter: Aandeel intensief cultuurgrasland (hp) toetsen (< %.0f%%)...", max_hp_pct))
    
    r_hp <- vul_raster_vanuit_parquet(
      ds_parquet  = ds_parquet,
      codes_bevat = "hp",
      master_grid = master_grid,
      waarde_vast = NULL
    )
    r_hp_clean <- terra::ifel(is.na(r_hp), 0, r_hp)
    rm(r_hp)
    
    df_tot <- terra::zonal(final_opp, cl_opp, fun = "sum", na.rm = TRUE)
    colnames(df_tot) <- c("ID", "Tot_Val")
    
    df_hp <- terra::zonal(r_hp_clean, cl_opp, fun = "sum", na.rm = TRUE)
    colnames(df_hp) <- c("ID", "HP_Val")
    rm(r_hp_clean)
    
    dt_hp <- merge(data.table::as.data.table(df_tot), data.table::as.data.table(df_hp), by = "ID", all.x = TRUE)
    dt_hp[is.na(HP_Val), HP_Val := 0]
    dt_hp[, Fractie_HP := HP_Val / Tot_Val]
    
    ok_complex_ids <- dt_hp[!is.na(ID) & Fractie_HP < (max_hp_pct / 100), ID]
    
    if (length(ok_complex_ids) > 0) {
      cl_opp    <- terra::ifel(cl_opp %in% as.integer(ok_complex_ids), cl_opp, NA)
      final_opp <- terra::mask(final_opp, cl_opp)
    } else {
      final_opp <- master_grid * NA_real_
      cl_opp    <- master_grid * NA_real_
    }
    message(sprintf("       📊 Na matrixfilter (< %.0f%% hp): %.2f ha", max_hp_pct, calc_ha_exact(final_opp)))
  }
  
  ha_fin <- calc_ha_exact(final_opp)
  message(sprintf("    ✅ [OK] Definitief berekend leefgebied voor %s: %.2f ha", soort_naam, ha_fin))
  
  res_uitvoer <- list(final_opp = final_opp, cl_opp = cl_opp)
  for (nm in types_te_fuseren) {
    res_uitvoer[[nm]] <- terra::mask(verwerkte_sporen[[nm]], final_opp)
  }
  
  if ("bwk_strikt" %in% names(verwerkte_sporen) && "bwk_inclusief" %in% names(verwerkte_sporen)) {
    message("    🎨 Gradatieraster samenstellen (Klasse 2 = Strikt droog, Klasse 1 = Inclusief/suboptimaal)...")
    r_bs <- terra::ifel(!is.na(verwerkte_sporen[["bwk_strikt"]]) & verwerkte_sporen[["bwk_strikt"]] > 0, 1, 0)
    r_bi <- terra::ifel(!is.na(verwerkte_sporen[["bwk_inclusief"]]) & verwerkte_sporen[["bwk_inclusief"]] > 0, 1, 0)
    res_uitvoer$gradatie <- terra::ifel(r_bs == 1, 2, terra::ifel(r_bi == 1, 1, NA_integer_))
    rm(r_bs, r_bi)
  }
  
  return(res_uitvoer)
}
