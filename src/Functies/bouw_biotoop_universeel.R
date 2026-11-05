# ==============================================================================
# bouw_biotoop_universeel.R
# Universele data-gedreven motor voor alle ecologische archetypes (Productieversie)
# Ondersteunt: Paapje, Modderkruiper, Bever, Klauwier, Wespendief, Zwarte specht, 
#              Boompieper, Nachtzwaluw, Maanwaterjuffer, Kwak, Kleine zilverreiger,
#              Watersnip, Zadelsprinkhaan
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
  # 0. INTERNE HELPERS
  # ----------------------------------------------------------------------------
  get_param <- function(df_rij, col_naam, default_val = NA) {
    if (!(col_naam %in% colnames(df_rij)) || nrow(df_rij) == 0) return(default_val)
    val <- df_rij[[col_naam]][1]
    if (is.null(val) || is.na(val) || val == "" || val == "NA") return(default_val)
    return(val)
  }
  
  # Veilige additieve fusie van fracties (voorkomt NA-besmetting en clamp op 1.0)
  voeg_toe_geclamped <- function(r_basis, r_nieuw) {
    heeft_b <- !is.null(r_basis) && !all(is.na(suppressWarnings(terra::minmax(r_basis))))
    heeft_n <- !is.null(r_nieuw) && !all(is.na(suppressWarnings(terra::minmax(r_nieuw))))
    
    if (!heeft_n) return(r_basis)
    if (!heeft_b) return(r_nieuw)
    
    b_nul <- terra::ifel(is.na(r_basis), 0, r_basis)
    n_nul <- terra::ifel(is.na(r_nieuw), 0, r_nieuw)
    som   <- terra::clamp(b_nul + n_nul, upper = 1.0)
    return(terra::ifel(som > 0, som, NA_real_))
  }
  
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
      r_tp <- extract_bwk_laag(soort_naam, tp, ds_parquet, df_bwk_codes, master_grid)
      if (is.null(r_tp)) r_tp <- master_grid * NA_real_
      terra::writeRaster(r_tp, f_tif, overwrite = TRUE, datatype = "FLT4S", gdal = "COMPRESS=LZW")
      lagen_ruw[[tp]] <- terra::rast(f_tif)
    }
    message(sprintf("       📊 Ruwe oppervlakte '%s': %.2f ha", tp, calc_ha_exact(lagen_ruw[[tp]])))
  }
  
  verwerkte_sporen <- list()
  cl_sporen        <- list()
  
  # ----------------------------------------------------------------------------
  # STAP 2: VERWERKING PER TYPE (RIJ IN EXCEL)
  # ----------------------------------------------------------------------------
  rij_volgorde <- order(params_soort$Spoor == "A", params_soort$Spoor == "B", decreasing = FALSE)
  params_gesorteerd <- params_soort[rij_volgorde, ]
  
  for (i in 1:nrow(params_gesorteerd)) {
    rij     <- params_gesorteerd[i, ]
    tp_naam <- trimws(rij$Type)
    
    # Exclusielagen en pure referentielagen zonder oppervlakte overslaan
    min_patch_val <- get_param(rij, "MinPatch_ha", NA)
    if (grepl("^not_", tolower(tp_naam)) || tolower(tp_naam) == "bwkbos" || 
        grepl("_bos$", tolower(tp_naam)) || tolower(tp_naam) == "openwater" || is.na(min_patch_val)) {
      next
    }
    
    message(sprintf("\n    🔹 Verwerken van laag/spoor: %s (Spoor %s)...", tp_naam, get_param(rij, "Spoor", "A")))
    r_huidig <- if (tp_naam %in% names(lagen_ruw)) lagen_ruw[[tp_naam]] else (master_grid * NA_real_)
    
    # 2-lgp-basis: Indien laag primair gedefinieerd wordt door teelten (bv. Europese hamster)
    if (all(is.na(suppressWarnings(terra::minmax(r_huidig)))) && !is.na(get_param(rij, "LGP_Gewascodes", NA))) {
      lgp_codes <- get_param(rij, "LGP_Gewascodes", NA)
      codes_v <- strsplit(as.character(lgp_codes), "[,; ]+")[[1]]
      message(sprintf("       -> Basislaag initialiseren vanuit LGP-teelten (%s)...", paste(codes_v, collapse = ", ")))
      r_huidig <- haal_categorisch_masker(paths$lgp_landbouw, codes_v, master_grid, paths$dir_cache, "lgp")
    }
    
    # 2-fusie: Types binnen hetzelfde spoor sommatief fuseren (bv. eiafzet_bwk + eiafzet_lineair)
    zelfde_spoor_types <- params_gesorteerd %>% 
      dplyr::filter(Spoor == huidig_spoor & !grepl("^not_", tolower(Type)) & 
                      !tolower(Type) %in% c("bosranden", "openwater", "water_bwk") &
                      !grepl("_bos$", tolower(Type))) %>% 
      dplyr::pull(Type)
    
    if (length(zelfde_spoor_types) > 1 && tp_naam == zelfde_spoor_types[1]) {
      message(sprintf("       -> Deellagen voor Spoor %s sommatief samenvoegen (%s)...", 
                      huidig_spoor, paste(zelfde_spoor_types, collapse = " + ")))
      r_som_spoor <- master_grid * 0
      for (st in zelfde_spoor_types) {
        if (st %in% names(lagen_ruw)) {
          r_som_spoor <- voeg_toe_geclamped(r_som_spoor, lagen_ruw[[st]])
        }
      }
      r_huidig <- r_som_spoor
      rm(r_som_spoor)
    } else if (length(zelfde_spoor_types) > 1 && tp_naam != zelfde_spoor_types[1]) {
      # Reeds samengevoegd onder het eerste type van dit spoor
      next
    }
    
    # 2-water: Filter op afstand tot open water (bv. Bont dikkopje <= 200m van water)
    if (isTRUE(as.logical(get_param(rij, "Gebruik_GRB_Water", FALSE)))) {
      water_dist_m <- as.numeric(get_param(rij, "Water_Afstand_m", NA))
      if (!is.na(water_dist_m) && water_dist_m > 0 && !is.null(paths$grb_water) && file.exists(paths$grb_water)) {
        message(sprintf("       -> Afstandsfilter tot open water (<= %dm)...", as.integer(water_dist_m)))
        f_grb_al  <- hulp_align_gdal(paths$grb_water, master_grid)
        r_grb     <- terra::rast(f_grb_al)
        r_grb_bin <- terra::ifel(!is.na(r_grb) & r_grb > 0, 1, NA)
        unlink(f_grb_al)
        rm(r_grb)
        
        # Optioneel aanvullen met specifiek BWK-water indien geladen
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
    
    # Kwaliteitsfilter op initiële fractie (bv. >= 0.20 bij vennen)
    if (grepl("waterbiotoop_bwk1", tolower(tp_naam))) {
      r_huidig <- terra::ifel(!is.na(r_huidig) & r_huidig >= 0.20, r_huidig, NA_real_)
    }
    
    # 2-akker: Negatieve akker-exclusie (Filter_Akker = TRUE: ^b codes uitsluiten)
    if (isTRUE(as.logical(get_param(rij, "Filter_Akker", FALSE)))) {
      message("       -> Akker-exclusie toepassen (^b codes uitsluiten)...")
      r_akker_bin <- vul_raster_vanuit_parquet(
        ds_parquet  = ds_parquet,
        codes_bevat = "b",
        master_grid = master_grid,
        waarde_vast = 1
      )
      if (!all(is.na(suppressWarnings(terra::minmax(r_akker_bin))))) {
        r_huidig <- terra::mask(r_huidig, r_akker_bin, inverse = TRUE)
      }
      rm(r_akker_bin)
      message(sprintf("       📊 Na akker-exclusie: %.2f ha", calc_ha_exact(r_huidig)))
    }
    
    # 2A. Automatische aftrek van types die beginnen met 'not_'
    huidig_spoor <- get_param(rij, "Spoor", "A")
    types_not <- params_soort %>%
      dplyr::filter(grepl("^not_", tolower(trimws(Type))) | grepl("_notbwk$", tolower(trimws(Type)))) %>%
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
    
    # 2B. Externe watervlakken toevoegen
    if (isTRUE(as.logical(get_param(rij, "Gebruik_Watervlakken_Shp", FALSE))) && 
        !is.null(paths$watervlakken_raster) && file.exists(paths$watervlakken_raster)) {
      message("       -> Externe watervlakken toevoegen...")
      f_wv_al  <- hulp_align_gdal(paths$watervlakken_raster, master_grid)
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
    
    # 2C. Huet-waterlopen toevoegen
    if (isTRUE(as.logical(get_param(rij, "Gebruik_Huetzones", FALSE)))) {
      huet_lbls <- get_param(rij, "Huet_Labels", NA)
      if (!is.na(huet_lbls)) {
        codes_v <- strsplit(as.character(huet_lbls), "[,; ]+")[[1]]
        message(sprintf("       -> Huet-waterlopen toevoegen (%s)...", paste(codes_v, collapse = ", ")))
        m_huet  <- haal_categorisch_masker(paths$huetzon, codes_v, master_grid, paths$dir_cache, "huet")
        r_huidig <- voeg_toe_geclamped(r_huidig, m_huet)
      }
    }
    
    # 2C-bis. Maas toevoegen
    if (isTRUE(as.logical(get_param(rij, "Gebruik_Maas", FALSE))) && 
        !is.null(paths$maas) && file.exists(paths$maas)) {
      message("       -> Maas toevoegen aan watercomplex...")
      f_maas_al <- hulp_align_gdal(paths$maas, master_grid)
      r_maas    <- terra::rast(f_maas_al)
      r_maas_ok <- terra::ifel(!is.na(r_maas) & r_maas >= 1, 1.0, NA_real_)
      unlink(f_maas_al)
      
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
      rm(r_maas, r_maas_ok)
    }
    
    # 2D. Landbouwgewassen (LGP): toevoegen (bij '+' prefix) of filteren als doorsnede
    lgp_codes <- get_param(rij, "LGP_Gewascodes", NA)
    if (!is.na(lgp_codes)) {
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
    heeft_drain <- isTRUE(as.logical(get_param(rij, "Gebruik_Drainage", FALSE))) && !is.na(get_param(rij, "Drainage_Codes", NA))
    heeft_ovstr <- isTRUE(as.logical(get_param(rij, "Gebruik_Overstroming", FALSE))) && !is.na(get_param(rij, "Overstroming_Codes", NA))
    heeft_prof  <- !is.na(get_param(rij, "Profiel_Codes", NA))
    
    # Check of het gaat om niet-overstroombaar (klasse 0 = strikte exclusie van inundatie, bv. Hamster)
    ovstr_val <- trimws(as.character(get_param(rij, "Overstroming_Codes", "")))
    is_droog_archetype <- (ovstr_val == "0")
    
    if (is_droog_archetype) {
      # Strikte AND-doorsnede: Niet-overstroombaar en leembodem (klasse 6)
      if (heeft_ovstr && !is.null(paths$overstroming) && file.exists(paths$overstroming)) {
        message("       -> Inundatiegevoelige gronden uitsluiten (enkel niet-overstroombaar klasse 0)...")
        f_ov_al <- hulp_align_gdal(paths$overstroming, master_grid)
        r_ov    <- terra::rast(f_ov_al)
        m_droog <- terra::ifel(is.na(r_ov) | r_ov == 0, 1, NA)
        unlink(f_ov_al)
        r_huidig <- terra::mask(r_huidig, m_droog)
        rm(r_ov, m_droog)
      }
      
      if (heeft_prof) {
        codes_p <- strsplit(as.character(get_param(rij, "Profiel_Codes", NA)), "[,; ]+")[[1]]
        message(sprintf("       -> Filteren op bodemtextuur/leem (klasse %s)...", paste(codes_p, collapse = ", ")))
        pad_bodem <- if (!is.null(paths$textuur) && file.exists(paths$textuur)) paths$textuur else paths$profiel
        m_leem <- haal_categorisch_masker(pad_bodem, codes_p, master_grid, paths$dir_cache, "textuur")
        r_huidig <- terra::mask(r_huidig, m_leem)
        rm(m_leem)
      }
      message(sprintf("       📊 Na bodem- en overstromingsfilter: %.2f ha", calc_ha_exact(r_huidig)))
      
    } else if (heeft_drain || heeft_ovstr || heeft_prof) {
      # Valleicontour (OR): Drainage OR Overstroming OR Bodemprofiel (bv. Wulp, Modderkruiper)
      m_vallei <- master_grid * NA_real_
      
      if (heeft_drain) {
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
      
      if (heeft_ovstr) {
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
      
      message("       -> Gecombineerde valleitest toepassen (Drainage OR Overstroming OR Profiel)...")
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
    
    # 2I. Groenkaart integratie
    if (isTRUE(as.logical(get_param(rij, "Gebruik_Groenkaart", FALSE)))) {
      bos_drempel <- as.numeric(get_param(rij, "Groenkaart_Drempel_m2", NA))
      
      if (!is.na(bos_drempel) && bos_drempel > 0) {
        message(sprintf("       -> Verbossingsfilter: Boscomplexen >= %d m² uitsluiten...", as.integer(bos_drempel)))
        f_cache_bos <- file.path(paths$dir_cache, sprintf("masker_bos_%dm2.tif", as.integer(bos_drempel)))
        
        if (file.exists(f_cache_bos) && !force_overwrite) {
          m_bos_groot <- terra::rast(f_cache_bos)
        } else {
          f_gr_al <- hulp_align_gdal(paths$groenkaart, master_grid)
          r_gr    <- terra::rast(f_gr_al)
          r_b1    <- terra::ifel(!is.na(r_gr) & r_gr == 1, 1, NA)
          unlink(f_gr_al)
          rm(r_gr)
          
          f_b_bin <- tempfile(fileext = ".tif")
          f_cl_b  <- tempfile(fileext = ".tif")
          terra::writeRaster(r_b1, f_b_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
          rm(r_b1)
          
          e <- terra::ext(master_grid)
          qgisprocess::qgis_run_algorithm(
            "grass:r.clump", input = f_b_bin, output = f_cl_b, title = "cl_bos_groot",
            "-d" = FALSE, GRASS_REGION_PARAMETER = sprintf("%f,%f,%f,%f", e$xmin, e$xmax, e$ymin, e$ymax),
            GRASS_REGION_CELLSIZE_PARAMETER = 10, .quiet = TRUE
          )
          unlink(f_b_bin)
          
          cl_b <- terra::rast(f_cl_b)
          try(terra::set.cats(cl_b, layer = 1, value = NULL), silent = TRUE)
          try(terra::coltab(cl_b) <- NULL, silent = TRUE)
          
          min_cellen_bos <- ceiling(bos_drempel / 100) # 2000m² = 20 cellen
          dt_b <- data.table::as.data.table(terra::freq(cl_b))
          setnames(dt_b, c("layer", "ID", "count"))
          grote_ids <- dt_b[!is.na(ID) & ID > 0 & count >= min_cellen_bos, ID]
          rm(dt_b)
          
          m_bos_groot <- if (length(grote_ids) > 0) terra::ifel(cl_b %in% as.integer(grote_ids), 1, NA) else (master_grid * NA_real_)
          unlink(f_cl_b)
          rm(cl_b)
          terra::writeRaster(m_bos_groot, f_cache_bos, overwrite = TRUE, datatype = "INT1U", gdal = "COMPRESS=LZW")
        }
        
        r_huidig <- terra::mask(r_huidig, m_bos_groot, inverse = TRUE)
        message(sprintf("       📊 Na verbossingsfilter (>= %d m²): %.2f ha", as.integer(bos_drempel), calc_ha_exact(r_huidig)))
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
      f_cache_bebouw <- file.path(paths$dir_cache, sprintf("masker_bebouwing_groot_%dm.tif", as.integer(buf_bebouw_m)))
      
      is_valid_cache <- FALSE
      if (file.exists(f_cache_bebouw) && !force_overwrite) {
        test_r <- suppressWarnings(terra::rast(f_cache_bebouw))
        if (terra::compareGeom(test_r, master_grid, stopOnError = FALSE)) {
          is_valid_cache <- TRUE
          m_bebouw_stoor <- test_r
          message(sprintf("       -> [CACHE] Bebouwingsbuffer (%dm) ingeladen...", as.integer(buf_bebouw_m)))
        }
      }
      
      if (!is_valid_cache) {
        message(sprintf("       -> Bebouwingsbuffer genereren (GDX >= 10 ha, %dm buffer)...", as.integer(buf_bebouw_m)))
        gdx_codes <- c("u", "ua", "ua-", "uc", "ud", "ui", "ui-", "un", "un-", "ur", "ur-", "ur+", "uv", "uv-", "uv+")
        
        r_bebouw_bin <- vul_raster_vanuit_parquet(
          ds_parquet  = ds_parquet,
          codes_exact = gdx_codes,
          master_grid = master_grid,
          waarde_vast = 1
        )
        
        res_b <- cluster_filter_universeel(
          masker            = r_bebouw_bin, 
          opp_laag          = r_bebouw_bin, 
          min_opp_ha        = 10, 
          min_totale_opp_ha = NA, 
          afstand_m         = 25, 
          dispersie_m       = 0, 
          werkelijk         = FALSE
        )
        
        if (!all(is.na(suppressWarnings(terra::minmax(res_b$raster))))) {
          f_b_bin <- tempfile(fileext = ".tif")
          terra::writeRaster(res_b$raster, f_b_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
          f_prox_b <- hulp_proximity_gdal(f_b_bin, buf_bebouw_m)
          unlink(f_b_bin)
          
          p_b <- terra::rast(f_prox_b)
          m_bebouw_stoor <- terra::ifel(!is.na(p_b) & p_b >= 0 & p_b <= buf_bebouw_m, 1, NA)
          unlink(f_prox_b)
          rm(p_b)
        } else {
          m_bebouw_stoor <- master_grid * NA_real_
        }
        rm(res_b, r_bebouw_bin)
        terra::writeRaster(m_bebouw_stoor, f_cache_bebouw, overwrite = TRUE, datatype = "INT1U", gdal = "COMPRESS=LZW")
        gc()
      }
      r_huidig <- terra::mask(r_huidig, m_bebouw_stoor, inverse = TRUE)
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
    
    # 2L-pre. Richtinggevende voorfilter
    grens_moment_val <- tolower(trimws(get_param(rij, "Grens_Moment", "")))
    if (grens_moment_val == "voor_clustering") {
      ref_tp  <- trimws(get_param(rij, "Grens_Met_Type", NA))
      grens_m <- as.numeric(get_param(rij, "Grens_Afstand_m", 0))
      
      r_ref_bron <- if (ref_tp %in% names(verwerkte_sporen)) verwerkte_sporen[[ref_tp]] else (if (ref_tp %in% names(lagen_ruw)) lagen_ruw[[ref_tp]] else NULL)
      
      if (!is.null(r_ref_bron) && grens_m > 0) {
        message(sprintf("       -> Voorfilter: '%s' beperken tot binnen %dm van '%s'...", tp_naam, as.integer(grens_m), ref_tp))
        r_ref_bin <- terra::ifel(!is.na(r_ref_bron) & r_ref_bron > 0, 1, NA)
        if (!all(is.na(suppressWarnings(terra::minmax(r_ref_bin))))) {
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
    
    # SCENARIO 1: Kernen vs. Satellieten / Buur-netwerken
    rand_diepte_val <- as.numeric(get_param(rij, "Randzone_Diepte_m", NA))
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
      
      if (!all(is.na(suppressWarnings(terra::minmax(cl_raw))))) {
        df_opp <- terra::zonal(r_huidig, cl_raw, fun = "sum", na.rm = TRUE)
        colnames(df_opp) <- c("ID", "Opp_Sum")
        dt_agg <- data.table::as.data.table(df_opp)
        dt_agg[, Opp_ha := Opp_Sum * 0.01]
        
        ids_groot <- dt_agg[!is.na(ID) & Opp_ha >= kern_opp_ha, ID]
        ids_klein <- dt_agg[!is.na(ID) & Opp_ha >= min_opp_ha & Opp_ha < kern_opp_ha, ID]
        
        r_groot <- if (length(ids_groot) > 0) terra::mask(r_huidig, cl_raw %in% as.integer(ids_groot)) else (master_grid * NA_real_)
        r_klein <- if (length(ids_klein) > 0) terra::mask(r_huidig, cl_raw %in% as.integer(ids_klein)) else (master_grid * NA_real_)
        
        heeft_groot <- !all(is.na(suppressWarnings(terra::minmax(r_groot))))
        heeft_klein <- !all(is.na(suppressWarnings(terra::minmax(r_klein))))
        
        # TAK A: Buur-netwerktoets (bv. Zadelsprinkhaan: 1-5 ha met >= 6 buren binnen 500m)
        if (!is.na(min_buren_val) && min_buren_val > 0) {
          buur_straal_m <- if (!is.na(actieradius_m) && actieradius_m > 0) actieradius_m else 500
          message(sprintf("       -> Buurtoets: Patches (%.1f-%.1f ha) toetsen op >= %d buren binnen %dm...",
                          min_opp_ha, kern_opp_ha, as.integer(min_buren_val), as.integer(buur_straal_m)))
          
          if (heeft_klein) {
            cl_klein <- terra::mask(cl_raw, cl_raw %in% as.integer(ids_klein))
            p_poly   <- terra::as.polygons(cl_klein, dissolve = TRUE)
            
            if (!is.null(p_poly) && nrow(p_poly) > 0) {
              pts_cent <- terra::centroids(p_poly, inside = TRUE)
              buf_buur <- terra::buffer(pts_cent, width = buur_straal_m)
              rel_mat  <- terra::is.related(buf_buur, pts_cent, "intersects")
              
              pts_cent$buren_count <- rowSums(rel_mat)
              ok_pts <- pts_cent[pts_cent$buren_count >= min_buren_val, ]
              
              if (nrow(ok_pts) > 0) {
                id_col <- names(p_poly)[1]
                valide_klein_ids <- ok_pts[[id_col]]
                r_klein_ok <- terra::mask(r_klein, cl_raw %in% as.integer(valide_klein_ids))
              } else {
                r_klein_ok <- master_grid * NA_real_
              }
              rm(pts_cent, buf_buur, rel_mat, p_poly)
            } else {
              r_klein_ok <- master_grid * NA_real_
            }
            rm(cl_klein)
          } else {
            r_klein_ok <- master_grid * NA_real_
          }
          
          # TAK B: Satelliet/kern interactie op actieradius (bv. Zwarte specht binnen 4000m)
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
      
      # SCENARIO 2: Standaard clustering op enkele drempel
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
      
      # SCENARIO 2-ecotoon: Bosrand-raaklijntoets (>= 50m) en 20m boszoom (bv. Dagvlinders)
      min_raaklijn_m <- as.numeric(get_param(rij, "Min_Rand_Lengte_m", NA))
      bos_zoom_tp    <- trimws(get_param(rij, "Grens_Met_Type", NA))
      
      if (!is.na(min_raaklijn_m) && min_raaklijn_m > 0 && !is.na(bos_zoom_tp) && 
          bos_zoom_tp %in% names(lagen_ruw) && !all(is.na(suppressWarnings(terra::minmax(cl_huidig))))) {
        
        message(sprintf("       -> Ecotonale bosrandtoets (>= %dm raaklijn met '%s') & 20m boszoom...", 
                        as.integer(min_raaklijn_m), bos_zoom_tp))
        
        r_bos_opp <- lagen_ruw[[bos_zoom_tp]]
        f_hab_bin <- tempfile(fileext = ".tif")
        terra::writeRaster(terra::ifel(!is.na(r_huidig) & r_huidig > 0, 1, NA), f_hab_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
        f_prox_hab10 <- hulp_proximity_gdal(f_hab_bin, 10)
        unlink(f_hab_bin)
        
        p_hab10 <- terra::rast(f_prox_hab10)
        r_bos_bin <- terra::ifel(!is.na(r_bos_opp) & r_bos_opp > 0, 1, NA)
        r_contact_pixels <- terra::ifel(!is.na(p_hab10) & p_hab10 >= 0 & p_hab10 <= 10 & !is.na(r_bos_bin), 1, NA)
        unlink(f_prox_hab10)
        rm(p_hab10, r_bos_bin)
        
        if (!all(is.na(suppressWarnings(terra::minmax(r_contact_pixels))))) {
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
              
              # 20m boszoom toevoegen aan het habitat
              f_hab_ok_bin <- tempfile(fileext = ".tif")
              terra::writeRaster(terra::ifel(!is.na(r_huidig), 1, NA), f_hab_ok_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
              f_prox_bz20 <- hulp_proximity_gdal(f_hab_ok_bin, 20)
              unlink(f_hab_ok_bin)
              p_bz20 <- terra::rast(f_prox_bz20)
              r_bos_zoom <- terra::mask(r_bos_opp, terra::ifel(!is.na(p_bz20) & p_bz20 >= 0 & p_bz20 <= 20, 1, NA))
              unlink(f_prox_bz20)
              rm(p_bz20)
              
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
      }
      
      # SCENARIO 2-max: Maximale oppervlaktefilter (bv. poelen > 1.0 ha uitsluiten)
      max_opp_ha_val <- as.numeric(get_param(rij, "MaxOpp_ha", NA))
      if (!is.na(max_opp_ha_val) && max_opp_ha_val > 0 && !all(is.na(suppressWarnings(terra::minmax(cl_huidig))))) {
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
      
      # SCENARIO 2-buren: Dichtheidseis / Poelencluster (bv. Gevlekte witsnuitlibel: >= 4 poelen binnen 200m)
      min_buren_val <- as.numeric(get_param(rij, "Min_Aantal_Buren", NA))
      if (!is.na(min_buren_val) && min_buren_val > 0 && is.na(kern_opp_ha) && !all(is.na(suppressWarnings(terra::minmax(cl_huidig))))) {
        buur_straal_m <- if (!is.na(actieradius_m) && actieradius_m > 0) actieradius_m else 200
        message(sprintf("       -> Dichtheidstoets: Clusters toetsen op >= %d buren binnen %dm...", 
                        as.integer(min_buren_val), as.integer(buur_straal_m)))
        
        p_poly <- terra::as.polygons(cl_huidig, dissolve = TRUE)
        if (!is.null(p_poly) && nrow(p_poly) > 0) {
          # Centroïden omzetten naar sf voor veilige buur-telling
          pts_cent <- terra::centroids(p_poly, inside = TRUE)
          sf_pts   <- sf::st_as_sf(pts_cent)
          
          # Veilige intersectie-telling via sf (crasht nooit op matrix-dimensies)
          sf_buf   <- sf::st_buffer(sf_pts, dist = buur_straal_m)
          int_list <- sf::st_intersects(sf_buf, sf_pts)
          buren_telling <- lengths(int_list)
          
          ok_idx <- which(buren_telling >= min_buren_val)
          
          if (length(ok_idx) > 0) {
            id_col <- names(p_poly)[1]
            valide_ids <- p_poly[[id_col]][ok_idx]
            cl_huidig  <- terra::ifel(cl_huidig %in% as.integer(valide_ids), cl_huidig, NA)
            r_huidig   <- terra::mask(r_huidig, cl_huidig)
          } else {
            cl_huidig <- master_grid * NA_real_
            r_huidig  <- master_grid * NA_real_
          }
          rm(pts_cent, sf_pts, sf_buf, int_list, buren_telling, p_poly)
        } else {
          cl_huidig <- master_grid * NA_real_
          r_huidig  <- master_grid * NA_real_
        }
        message(sprintf("       📊 Na dichtheidstoets (>= %d poelen): %.2f ha", as.integer(min_buren_val), calc_ha_exact(r_huidig)))
      }
      
      # SCENARIO 2b: Inwaartse randzone uitsnijden
      if (is_inwaarts_rand && !all(is.na(suppressWarnings(terra::minmax(r_huidig))))) {
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
      
      # SCENARIO 2c: Bospercentage / Structuureis
      min_bos_pct <- as.numeric(get_param(rij, "Min_Struweel_Pct", NA))
      if (!is.na(min_bos_pct) && min_bos_pct > 0 && !all(is.na(suppressWarnings(terra::minmax(cl_huidig))))) {
        # Flexibel identificeren van de bijbehorende bos-/struweelreferentielaag
        bos_kandidaten <- c(
          if (grepl("ruim", tolower(tp_naam))) "landbiotoop_bwk_ruim_bos" else "landbiotoop_bos",
          "landbiotoop_bwk_ruim_bos", "landbiotoop_bos", "landbiotoop_nabij_bos", "struweel1"
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
      
      # SCENARIO 2d: Zuivere procentuele wateraandeel-toets (bv. Kwak & Kleine zilverreiger >= 20% water)
      min_water_pct_val <- as.numeric(get_param(rij, "Min_Water_Pct", NA))
      if (!is.na(min_water_pct_val) && min_water_pct_val > 0 && !all(is.na(suppressWarnings(terra::minmax(cl_huidig))))) {
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
      
      # SCENARIO 2e: Inbedding in overkoepelend type
      inbedding_tp <- trimws(get_param(rij, "Inbedding_In_Type", NA))
      if (!is.na(inbedding_tp) && inbedding_tp %in% names(verwerkte_sporen) && !all(is.na(suppressWarnings(terra::minmax(r_huidig))))) {
        message(sprintf("       -> Inbeddingstoets: '%s' beperken tot binnen goedgekeurd '%s'...", tp_naam, inbedding_tp))
        r_inbed_mask <- terra::ifel(!is.na(verwerkte_sporen[[inbedding_tp]]) & verwerkte_sporen[[inbedding_tp]] > 0, 1, NA)
        r_huidig <- terra::mask(r_huidig, r_inbed_mask)
        cl_huidig <- terra::mask(cl_huidig, r_inbed_mask)
        rm(r_inbed_mask)
        message(sprintf("       📊 Na inbeddingstoets: %.2f ha", calc_ha_exact(r_huidig)))
      }
      
      # SCENARIO 2f: Steltloper-kwaliteitstoets & Terrestrische kerncheck (bv. Watersnip)
      gebruik_grb_water <- isTRUE(as.logical(get_param(rij, "Gebruik_GRB_Water", FALSE)))
      kern_drempel_ha   <- as.numeric(get_param(rij, "Kern_MinOpp_ha", NA))
      
      if (gebruik_grb_water && !is.na(kern_drempel_ha) && kern_drempel_ha > 0 && 
          !all(is.na(suppressWarnings(terra::minmax(cl_huidig))))) {
        
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
          
          if (!all(is.na(suppressWarnings(terra::minmax(cl_terr_10ha))))) {
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
    
    # SCENARIO 2g: Dynamische kustweidevogel-toets & Vennenfilter (bv. Tureluur)
    heeft_slik <- "slik_schor" %in% names(lagen_ruw)
    heeft_ven  <- "vennen_water" %in% names(lagen_ruw)
    
    if ((heeft_slik || heeft_ven) && !all(is.na(suppressWarnings(terra::minmax(cl_huidig))))) {
      message("       -> Dynamische oppervlaktetoets (5-10 ha) en vennenfilter (< 95%) berekenen...")
      
      r_slik_clean <- if (heeft_slik) terra::ifel(is.na(lagen_ruw[["slik_schor"]]), 0, lagen_ruw[["slik_schor"]]) else (master_grid * 0)
      r_ven_clean  <- if (heeft_ven)  terra::ifel(is.na(lagen_ruw[["vennen_water"]]), 0, lagen_ruw[["vennen_water"]]) else (master_grid * 0)
      
      # 1. Sommaties per cluster
      df_tot_cl <- terra::zonal(r_huidig, cl_huidig, fun = "sum", na.rm = TRUE)
      colnames(df_tot_cl) <- c("ID", "Totaal_Val")
      df_tot_cl$Totaal_ha <- df_tot_cl$Totaal_Val * 0.01
      
      df_slik_cl <- terra::zonal(r_slik_clean, cl_huidig, fun = "sum", na.rm = TRUE)
      colnames(df_slik_cl) <- c("ID", "Slik_Val")
      df_slik_cl$Slik_ha <- df_slik_cl$Slik_Val * 0.01
      
      df_ven_cl <- terra::zonal(r_ven_clean, cl_huidig, fun = "sum", na.rm = TRUE)
      colnames(df_ven_cl) <- c("ID", "Ven_Val")
      df_ven_cl$Ven_ha <- df_ven_cl$Ven_Val * 0.01
      
      # 2. Tabel samenvoegen en dynamische formule berekenen
      dt_dyn <- merge(data.table::as.data.table(df_tot_cl), data.table::as.data.table(df_slik_cl), by = "ID", all.x = TRUE)
      dt_dyn <- merge(dt_dyn, data.table::as.data.table(df_ven_cl), by = "ID", all.x = TRUE)
      dt_dyn[is.na(Slik_ha), Slik_ha := 0]
      dt_dyn[is.na(Ven_ha), Ven_ha := 0]
      
      # Formule: 10 - 5 * (Slik_ha / Totaal_ha)
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
    
    # SCENARIO 2h: Ecoregio-filter (bv. Tureluur: Polders, Schelde, IJzervallei)
    ecoregio_val <- get_param(rij, "Ecoregio", NA)
    if (!is.na(ecoregio_val) && !is.null(paths$ecoregio) && file.exists(paths$ecoregio) && 
        !all(is.na(suppressWarnings(terra::minmax(cl_huidig))))) {
      message("       -> Ecoregio-filter toepassen...")
      f_eco_al <- hulp_align_gdal(paths$ecoregio, master_grid)
      r_eco    <- terra::rast(f_eco_al)
      m_eco    <- terra::ifel(!is.na(r_eco) & r_eco >= 1, 1, NA)
      unlink(f_eco_al)
      
      # Cluster moet overlappen met goedgekeurde ecoregio
      ct_eco <- terra::crosstab(c(cl_huidig, m_eco), long = TRUE)
      colnames(ct_eco) <- c("Cluster_ID", "Eco_Val", "n")
      dt_eco <- data.table::as.data.table(ct_eco)[!is.na(Cluster_ID) & !is.na(Eco_Val) & n > 0]
      ok_eco_ids <- unique(dt_eco$Cluster_ID)
      
      if (length(ok_eco_ids) > 0) {
        cl_huidig <- terra::ifel(cl_huidig %in% as.integer(ok_eco_ids), cl_huidig, NA)
        r_huidig  <- terra::mask(r_huidig, cl_huidig)
      } else {
        cl_huidig <- master_grid * NA_real_
        r_huidig  <- master_grid * NA_real_
      }
      rm(r_eco, m_eco, ct_eco, dt_eco)
      message(sprintf("       📊 Na Ecoregio-filter: %.2f ha", calc_ha_exact(r_huidig)))
    }
    
    # 2M. Lineaire KLE-toets
    min_kle_m <- as.numeric(get_param(rij, "Min_KLE_Lengte_m", 0))
    if (!is.na(min_kle_m) && min_kle_m > 0 && !all(is.na(suppressWarnings(terra::minmax(cl_huidig))))) {
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
    cl_sporen[[tp_naam]]        <- cl_huidig
  }
  
  # ----------------------------------------------------------------------------
  # STAP 3: KOPPELINGEN TUSSEN TYPES (Randlengte, Ecotoon, Wederzijds of Richtinggevend)
  # ----------------------------------------------------------------------------
  
  # 3A. Bosrandlengtetoets
  kol_rand_l <- intersect(c("Min_Rand_Lengte_m", "Min_Oever_Lengte_m"), colnames(params_soort))[1]
  waarden_rand_l <- if (!is.na(kol_rand_l)) params_soort[[kol_rand_l]] else NA
  min_rand_m <- suppressWarnings(max(as.numeric(waarden_rand_l), na.rm = TRUE))
  heeft_rand_lengte_toets <- !is.infinite(min_rand_m) && !is.na(min_rand_m) && min_rand_m > 0
  
  # 3B. Ecotoon / Bosrand-contactzone
  rand_diepte_m <- suppressWarnings(max(as.numeric(params_soort$Randzone_Diepte_m), na.rm = TRUE))
  heeft_ecotoon_archetype <- !is.infinite(rand_diepte_m) && !is.na(rand_diepte_m) && rand_diepte_m > 0 && 
    !("na_clustering" %in% tolower(trimws(params_soort$Grens_Moment)))
  
  # 3C. Wederzijdse koppeling
  is_wederzijds <- any(as.logical(params_soort$Gebruik_Wederzijds) == TRUE, na.rm = TRUE)
  
  if (heeft_rand_lengte_toets && length(verwerkte_sporen) >= 2) {
    rij_rl     <- params_soort %>% 
      dplyr::filter(!is.na(.data[[kol_rand_l]]) & as.numeric(.data[[kol_rand_l]]) > 0) %>% 
      dplyr::slice(1)
    open_tp    <- trimws(rij_rl$Type)
    bos_tp     <- trimws(get_param(rij_rl, "Grens_Met_Type", NA))
    groep_buf  <- as.numeric(get_param(rij_rl, "Grens_Afstand_m", 25))
    min_pixels <- round(min_rand_m / 10)
    
    if (!is.na(bos_tp) && open_tp %in% names(verwerkte_sporen) && bos_tp %in% names(verwerkte_sporen)) {
      message(sprintf("\n    📏 Bosrandlengte-toets: '%s' toetsen op >= %dm (%d pixels) bosrand met '%s'...",
                      open_tp, as.integer(min_rand_m), as.integer(min_pixels), bos_tp))
      
      r_open <- verwerkte_sporen[[open_tp]]
      r_bos  <- verwerkte_sporen[[bos_tp]]
      
      if (!all(is.na(suppressWarnings(terra::minmax(r_open)))) && !all(is.na(suppressWarnings(terra::minmax(r_bos))))) {
        r_bos_bin <- terra::ifel(!is.na(r_bos) & r_bos > 0, 1, 0)
        w_3x3     <- matrix(1, nrow = 3, ncol = 3)
        f_fc_b    <- tempfile(fileext = ".tif")
        r_b_focal <- terra::focal(r_bos_bin, w = w_3x3, fun = "max", na.rm = TRUE, filename = f_fc_b, overwrite = TRUE)
        unlink(f_fc_b)
        rm(r_bos_bin)
        
        r_open_bin <- terra::ifel(!is.na(r_open) & r_open > 0, 1, 0)
        r_rand_bin <- terra::ifel(r_open_bin == 1 & r_b_focal == 1, 1, 0)
        rm(r_b_focal)
        
        r_open_alleen <- terra::ifel(r_open_bin == 1, 1, NA)
        rm(r_open_bin)
        
        f_op_bin   <- tempfile(fileext = ".tif")
        terra::writeRaster(r_open_alleen, f_op_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 0)
        f_prox_op  <- hulp_proximity_gdal(f_op_bin, groep_buf)
        unlink(f_op_bin)
        
        p_op <- terra::rast(f_prox_op)
        mask_op_buf <- terra::ifel(!is.na(p_op) & p_op >= 0 & p_op <= groep_buf, 1, NA)
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
        unlink(f_cl_c)
        rm(r_open_alleen)
        
        cells_c <- terra::cells(!is.na(cl_comp_open) & cl_comp_open > 0)
        if (length(cells_c) > 0) {
          dt_rand <- data.table::data.table(
            ID         = as.integer(cl_comp_open[cells_c][[1]]),
            Rand_Pixel = as.integer(r_rand_bin[cells_c][[1]])
          )
          dt_stats <- dt_rand[, .(Totaal_Rand = sum(Rand_Pixel, na.rm = TRUE)), by = ID]
          ok_ids   <- dt_stats[Totaal_Rand >= min_pixels, ID]
          rm(dt_rand, dt_stats)
          
          message(sprintf("       📊 Complexen met >= %dm bosrand: %d", as.integer(min_rand_m), length(ok_ids)))
          if (length(ok_ids) > 0) {
            verwerkte_sporen[[open_tp]] <- terra::mask(r_open, cl_comp_open %in% as.integer(ok_ids))
          } else {
            verwerkte_sporen[[open_tp]] <- master_grid * NA_real_
          }
        } else {
          verwerkte_sporen[[open_tp]] <- master_grid * NA_real_
        }
        rm(cells_c, r_rand_bin, cl_comp, cl_comp_open)
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
    
    heeft_a <- !all(is.na(suppressWarnings(terra::minmax(laag_a))))
    heeft_b <- !all(is.na(suppressWarnings(terra::minmax(laag_b))))
    
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
    
    # 3C. Wederzijdse koppeling (symmetrisch en asymmetrisch via hulp_wederzijdse_koppeling)
  } else if (is_wederzijds && length(verwerkte_sporen) >= 2) {
    rij_a <- params_soort %>% dplyr::filter(as.logical(Gebruik_Wederzijds) == TRUE & Spoor == "A") %>% dplyr::slice(1)
    rij_b <- params_soort %>% dplyr::filter(as.logical(Gebruik_Wederzijds) == TRUE & Spoor == "B") %>% dplyr::slice(1)
    
    tp_a <- trimws(rij_a$Type)
    tp_b <- trimws(get_param(rij_a, "Grens_Met_Type", rij_b$Type))
    
    dist_a <- as.numeric(dplyr::coalesce(as.numeric(get_param(rij_a, "Grens_Afstand_m", NA)),
                                         as.numeric(get_param(rij_a, "Wederzijdse_Afstand_m", NA)), 50))
    dist_b <- as.numeric(dplyr::coalesce(as.numeric(get_param(rij_b, "Grens_Afstand_m", NA)),
                                         as.numeric(get_param(rij_b, "Wederzijdse_Afstand_m", NA)), dist_a))
    
    if (!is.na(tp_b) && tp_a %in% names(verwerkte_sporen) && tp_b %in% names(verwerkte_sporen)) {
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
    
  } else if (!is_wederzijds) {
    for (i in 1:nrow(params_soort)) {
      rij     <- params_soort[i, ]
      doel_tp <- trimws(rij$Type)
      kern_tp_raw <- get_param(rij, "Grens_Met_Type", NA)
      grens_m_raw <- get_param(rij, "Grens_Afstand_m", NA)
      grens_mom   <- tolower(trimws(get_param(rij, "Grens_Moment", "")))
      
      if (!is.na(kern_tp_raw) && !is.na(grens_m_raw) && grens_mom != "voor_clustering" && doel_tp %in% names(verwerkte_sporen)) {
        kern_types <- trimws(strsplit(as.character(kern_tp_raw), "[,;]+")[[1]])
        afstanden  <- as.numeric(trimws(strsplit(as.character(grens_m_raw), "[,;]+")[[1]]))
        if (length(afstanden) < length(kern_types)) afstanden <- rep(afstanden[1], length(kern_types))
        
        for (k in seq_along(kern_types)) {
          kern_tp <- kern_types[k]
          grens_m <- afstanden[k]
          
          if (kern_tp %in% names(verwerkte_sporen) && grens_m > 0) {
            message(sprintf("\n    🔗 Koppeling na clustering: '%s' beperken tot binnen %dm van '%s'...", doel_tp, as.integer(grens_m), kern_tp))
            r_kern_bin <- terra::ifel(!is.na(verwerkte_sporen[[kern_tp]]) & verwerkte_sporen[[kern_tp]] > 0, 1, NA)
            
            if (!all(is.na(suppressWarnings(terra::minmax(r_kern_bin))))) {
              f_k_bin  <- tempfile(fileext = ".tif")
              terra::writeRaster(r_kern_bin, f_k_bin, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
              f_prox_k <- hulp_proximity_gdal(f_k_bin, grens_m)
              unlink(f_k_bin)
              
              p_k <- terra::rast(f_prox_k)
              m_grens <- terra::ifel(!is.na(p_k) & p_k >= 0 & p_k <= grens_m, 1, NA)
              unlink(f_prox_k)
              
              verwerkte_sporen[[doel_tp]] <- terra::mask(verwerkte_sporen[[doel_tp]], m_grens)
              rm(p_k, m_grens)
            } else {
              message(sprintf("       ⚠️ Referentielaag '%s' is leeg! Doel '%s' wordt op NA gezet.", kern_tp, doel_tp))
              verwerkte_sporen[[doel_tp]] <- master_grid * NA_real_
            }
            message(sprintf("       📊 Gekoppeld habitat '%s': %.2f ha", doel_tp, calc_ha_exact(verwerkte_sporen[[doel_tp]])))
          }
        }
      }
    }
  }
  
  # ----------------------------------------------------------------------------
  # STAP 4: BRUTO COMPLEX SAMENVOEGEN, METACLUSTER & ACTIERADIUSTOETS
  # ----------------------------------------------------------------------------
  message("\n    🧩 Stap 4: Samenstellen van bruto complex (optelling van fracties)...")
  
  # Sluit types uit die enkel als ruime contextschil dienden
  types_te_fuseren <- names(verwerkte_sporen)
  context_types <- c("landbiotoop_bwk_ruim", "landbiotoop_bwk_ruim_bos", "land_ruim")
  types_te_fuseren <- setdiff(types_te_fuseren, context_types)
  
  final_opp <- NULL
  for (nm in types_te_fuseren) {
    r_sp <- verwerkte_sporen[[nm]]
    if (!is.null(r_sp) && !all(is.na(suppressWarnings(terra::minmax(r_sp))))) {
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
        if (!is.null(r_sp) && !all(is.na(suppressWarnings(terra::minmax(r_sp))))) {
          final_opp <- voeg_toe_geclamped(final_opp, r_sp)
        }
      }
      if (is.null(final_opp)) final_opp <- master_grid * NA_real_
      
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
        rm(r_fin_bin)
      } else {
        cl_opp <- master_grid * NA_real_
      }
    }
  } else {
    # 4D. Geen MinComplex vereist
    if (length(types_te_fuseren) == 1 && types_te_fuseren[1] %in% names(cl_sporen)) {
      cl_opp <- cl_sporen[[types_te_fuseren[1]]]
    } else {
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
        rm(r_fin_bin)
      } else {
        cl_opp <- master_grid * NA_real_
      }
      gc()
    }
  }
  
  # ----------------------------------------------------------------------------
  # STAP 5: MATRIXFILTER OP INTENSIEF GRASLAND (^hp < 90%)
  # ----------------------------------------------------------------------------
  max_hp_pct <- max(as.numeric(params_soort$Max_Intensief_Grasland_Pct), na.rm = TRUE)
  if (!is.infinite(max_hp_pct) && !is.na(max_hp_pct) && max_hp_pct > 0 && !all(is.na(suppressWarnings(terra::minmax(cl_opp))))) {
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
  for (nm in names(verwerkte_sporen)) {
    res_uitvoer[[nm]] <- terra::mask(verwerkte_sporen[[nm]], final_opp)
  }
  
  # Indien zowel een strikt als een inclusief spoor aanwezig is (bv. Veldkrekel): gradatie toevoegen
  if ("bwk_strikt" %in% names(verwerkte_sporen) && "bwk_inclusief" %in% names(verwerkte_sporen)) {
    message("    🎨 Gradatieraster samenstellen (Klasse 2 = Strikt droog, Klasse 1 = Inclusief/suboptimaal)...")
    r_bs <- terra::ifel(!is.na(verwerkte_sporen[["bwk_strikt"]]) & verwerkte_sporen[["bwk_strikt"]] > 0, 1, 0)
    r_bi <- terra::ifel(!is.na(verwerkte_sporen[["bwk_inclusief"]]) & verwerkte_sporen[["bwk_inclusief"]] > 0, 1, 0)
    res_uitvoer$gradatie <- terra::ifel(r_bs == 1, 2, terra::ifel(r_bi == 1, 1, NA_integer_))
    rm(r_bs, r_bi)
  }
  
  return(res_uitvoer)
}
