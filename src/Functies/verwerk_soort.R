# ==============================================================================
# verwerk_soort.R - DATA-GEDREVEN CONTROLLER VOOR ALLE SOORTEN
# ==============================================================================

verwerk_soort <- function(soort_naam, params_soort, master_grid, ds_parquet, df_bwk_codes, paths, force_overwrite = FALSE) {
  
  bestandsnaam <- tolower(gsub(" ", "", soort_naam))
  f_opp  <- file.path(paths$out_opp,  paste0("Potentieel_Leefgebied_", bestandsnaam, ".tif"))
  f_id   <- file.path(paths$out_id,   paste0("ID_Netwerken_", bestandsnaam, ".tif"))
  f_arpl <- file.path(paths$out_arpl, paste0("Habitat_ARPL_", bestandsnaam, ".tif"))
  
  message(paste0("\n>>> Start verwerking voor: ", toupper(soort_naam)))
  
  # ----------------------------------------------------------------------------
  # DEEL 1: CLUSTERING & POTENTIEEL LEEFGEBIED (00 & 01)
  # ----------------------------------------------------------------------------
  heeft_basis_rasters <- file.exists(f_opp) && file.exists(f_id)
  
  if (!force_overwrite && heeft_basis_rasters) {
    message("    ⏩ [OVERSLAAN BASIS] ID- en oppervlakte-rasters bestaan al.")
    cl_opp    <- terra::rast(f_id)
    final_opp <- terra::rast(f_opp)
  } else {
    
    # Bepaal parameters voor Spoor A
    spoor_a <- params_soort %>% dplyr::filter(Spoor == "A")
    if (nrow(spoor_a) == 0) spoor_a <- params_soort
    
    # Zoek de hoofddrempels
    hoofd_rij <- spoor_a %>% 
      dplyr::filter(!is.na(MinOpp_ha) & suppressWarnings(as.numeric(MinOpp_ha)) >= 0.01) %>% 
      dplyr::slice(1)
    if (nrow(hoofd_rij) == 0) hoofd_rij <- spoor_a %>% dplyr::slice(1)
    
    min_opp_ha     <- suppressWarnings(as.numeric(hoofd_rij$MinOpp_ha[1]))
    if (is.na(min_opp_ha)) min_opp_ha <- 0.1
    
    min_totale_opp <- suppressWarnings(as.numeric(hoofd_rij$MinTotaleOpp_ha[1]))
    afstand_m      <- suppressWarnings(as.numeric(hoofd_rij$AfstandBiotopen_m[1]))
    if (is.na(afstand_m)) afstand_m <- 0
    dispersie_m    <- suppressWarnings(as.numeric(hoofd_rij$Dispersiecap_m[1]))
    if (is.na(dispersie_m)) dispersie_m <- 0
    
    message(paste0("    Parameters: MinOpp=", min_opp_ha, " ha",
                   " | MinTotaleOpp=", ifelse(is.na(min_totale_opp), "Geen", paste0(min_totale_opp, " ha")),
                   " | Afstand=", afstand_m, " m",
                   " | Dispersie=", dispersie_m, " m"))
    
    dir_sub <- file.path(paths$out_opp, "Deellagen")
    dir.create(dir_sub, showWarnings = FALSE, recursive = TRUE)
    
    res_model <- NULL
    
    # ==========================================================================
    # ROUTE A: SPECIFIEKE ECOLOGISCHE ARCHETYPES
    # ==========================================================================
    if (bestandsnaam == "kommavlinder") {
      message("    🦋 Specifieke workflow voor kommavlinder...")
      res_model <- bouw_biotoop_kern_met_overgang(params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (bestandsnaam == "heivlinder") {
      message("    🦋 Specifieke workflow voor heivlinder...")
      res_model <- bouw_biotoop_mozaiek_met_substraat(params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (bestandsnaam == "boompieper") {
      message(paste0("    🐦 Bosrand-gradiënt workflow voor: ", soort_naam))
      res_model <- bouw_biotoop_bosrand_gradient(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (bestandsnaam == "bruinekiekendief") {
      message(paste0("    🦅 Actieradius & kernbiotoop workflow voor: ", soort_naam))
      res_model <- bouw_biotoop_actieradius_kern(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (grepl("beekrombout", tolower(soort_naam))) {
      source(here::here("src", "Functies", "bouw_biotoop_stromend_water_libel.R"))
      res <- bouw_biotoop_stromend_water_libel(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (grepl("gevlektewitsnuitlibel", tolower(soort_naam))) {
      source(here::here("src", "Functies", "bouw_biotoop_vennencomplex_libel.R"))
      res <- bouw_biotoop_vennencomplex_libel(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (grepl("hoogveenglanslibel", tolower(soort_naam))) {
      source(here::here("src", "Functies", "bouw_biotoop_hoogveen_libel.R"))
      res <- bouw_biotoop_hoogveen_libel(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (grepl("zwarteheidelibel", tolower(soort_naam))) {
      source(here::here("src", "Functies", "bouw_biotoop_hoogveen_libel.R"))
      res <- bouw_biotoop_zwarteheidelibel(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (grepl("maanwaterjuffer", tolower(soort_naam))) {
      source(here::here("src", "Functies", "bouw_biotoop_bosven_waterjuffer.R"))
      res <- bouw_biotoop_bosven_waterjuffer(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (grepl("noordsewitsnuitlibel", tolower(soort_naam))) {
      source(here::here("src", "Functies", "bouw_biotoop_dystroof_ven_libel.R"))
      res <- bouw_biotoop_dystroof_ven_libel(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (grepl("sierlijkewitsnuitlibel", tolower(soort_naam))) {
      source(here::here("src", "Functies", "bouw_biotoop_verlandingsplas_libel.R"))
      res <- bouw_biotoop_verlandingsplas_libel(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (grepl("speerwaterjuffer", tolower(soort_naam))) {
      source(here::here("src", "Functies", "bouw_biotoop_speerwaterjuffer.R"))
      res <- bouw_biotoop_speerwaterjuffer(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (grepl("venglazenmaker", tolower(soort_naam))) {
      source(here::here("src", "Functies", "bouw_biotoop_venglazenmaker.R"))
      res <- bouw_biotoop_venglazenmaker(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (grepl("venwitsnuitlibel", tolower(soort_naam))) {
      source(here::here("src", "Functies", "bouw_biotoop_venwitsnuitlibel.R"))
      res <- bouw_biotoop_venwitsnuitlibel(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (bestandsnaam == "grauweklauwier") {
      message(paste0("    🐦 KLE-mozaïek & storingsbuffer workflow voor: ", soort_naam))
      res_model <- bouw_biotoop_kle_mozaiek_met_storing(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (bestandsnaam == "grutto") {
      message(paste0("    🪶 Weidevogel-matrix workflow voor: ", soort_naam))
      res_model <- bouw_biotoop_weidevogel_matrix(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (bestandsnaam %in% c("eikelmuis", "hazelmuis")) {
      message(paste0("    🐭 Kleinschalig struweel- en loofbos workflow voor: ", soort_naam))
      res_model <- bouw_biotoop_kleinschalig_struweel_loofbos(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (bestandsnaam == "ijsvogel") {
      message(paste0("    🐟 Oever-gradiënt workflow voor: ", soort_naam))
      res_model <- bouw_biotoop_oever_gradient(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (bestandsnaam == "waterspitsmuis") {
      message(paste0("    🌊 Semi-aquatisch oeverzoogdier workflow voor: ", soort_naam))
      res_model <- bouw_biotoop_semi_aquatisch_oeverzoogdier(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (bestandsnaam == "europeseotter") {
      message(paste0("    🦦 Rivierlandschap & matrix workflow voor: ", soort_naam))
      res_model <- bouw_biotoop_otter_rivierlandschap(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (bestandsnaam %in% c("kluut", "steltkluut", "strandplevier")) {
      message(paste0("    🪶 Kuststeltloper & slik workflow voor: ", soort_naam))
      res_model <- bouw_biotoop_kuststeltloper_slik(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (bestandsnaam %in% c("kwak", "kleinezilverreiger", "lepelaar", "purperreiger")) {
      message(paste0("    🪶 Commuter water-moeras workflow voor: ", soort_naam))
      res_model <- bouw_biotoop_commuter_watermoeras(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (bestandsnaam == "tureluur") {
      message(paste0("    🪶 Dynamisch kustgrasland workflow voor: ", soort_naam))
      res_model <- bouw_biotoop_dynamisch_kustgrasland(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (bestandsnaam %in% c("europesebever", "bever")) {
      message(paste0("    🦫 Beekdal- en ooibos workflow voor: ", soort_naam))
      res_model <- bouw_biotoop_bever_beekdal(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (bestandsnaam %in% c("europesehamster", "hamster")) {
      message(paste0("    🐹 Leemakker- en graancomplex workflow voor: ", soort_naam))
      res_model <- bouw_biotoop_hamster_leemakker(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (bestandsnaam == "kwartelkoning") {
      message(paste0("    🌾 Vallei-grasland workflow voor: ", soort_naam))
      res_model <- bouw_biotoop_vallei_grasland_mozaiek(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (bestandsnaam %in% c("paapje", "graspieper")) {
      message(paste0("    🌾 Open grasland (bosrandvrij) workflow voor: ", soort_naam))
      res_model <- bouw_biotoop_open_grasland_bosrandvrij(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (bestandsnaam == "watersnip") {
      message(paste0("    🪶 Drasland-steltloper workflow voor: ", soort_naam))
      res_model <- bouw_biotoop_drasland_steltloper(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (bestandsnaam == "nachtzwaluw") {
      message(paste0("    🦉 Bosrand-lengte workflow voor: ", soort_naam))
      res_model <- bouw_biotoop_bosrand_lengte(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (bestandsnaam == "porseleinhoen") {
      message(paste0("    🌾 Open inundatiemoeras workflow voor: ", soort_naam))
      res_model <- bouw_biotoop_open_inundatiemoeras(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (bestandsnaam == "roerdomp") {
      message(paste0("    🪶 Rietmoeras-ecotoon complex workflow voor: ", soort_naam))
      res_model <- bouw_biotoop_rietmoeras_ecotoon_complex(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (bestandsnaam == "wespendief") {
      message(paste0("    🦅 Bos-landschap actieradius workflow voor: ", soort_naam))
      res_model <- bouw_biotoop_bos_landschap_actieradius(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (bestandsnaam == "woudaap") {
      message(paste0("    🪶 Rietmoeras met waterkern workflow voor: ", soort_naam))
      res_model <- bouw_biotoop_rietmoeras_waterkern(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (bestandsnaam == "wulp") {
      message(paste0("    🪶 Grootschalige weidevogel workflow voor: ", soort_naam))
      res_model <- bouw_biotoop_grootschalige_weidevogel(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (bestandsnaam == "zomertortel") {
      message(paste0("    🕊️ Kleinschalig cultuurland driehoek workflow voor: ", soort_naam))
      res_model <- bouw_biotoop_kleinschalig_cultuurland_driehoek(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (bestandsnaam == "zwartespecht") {
      message(paste0("    🪵 Bos-metapopulatie workflow voor: ", soort_naam))
      res_model <- bouw_biotoop_bos_metapopulatie_netwerk(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (bestandsnaam == "zwartkopmeeuw") {
      message(paste0("    🪶 Kolonievogel kust & polder workflow voor: ", soort_naam))
      res_model <- bouw_biotoop_kolonievogel_kust_polder(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (grepl("drijvendewaterweegbree|waterplant|hydrofiet", tolower(soort_naam))) {
      source(here::here("src", "Functies", "bouw_biotoop_waterplant_hydrograaf.R"))
      res <- bouw_biotoop_waterplant_hydrograaf(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (grepl("grotemodderkruiper|modderkruiper|vis", tolower(soort_naam))) {
      source(here::here("src", "Functies", "bouw_biotoop_alluviale_vis.R"))
      res <- bouw_biotoop_alluviale_vis(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
    } else if (grepl("spaansevlag|nachtvlinder", tolower(soort_naam))) {
      source(here::here("src", "Functies", "bouw_biotoop_bosrand_ecotoon.R"))
      res <- bouw_biotoop_bosrand_ecotoon(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (bestandsnaam == "duinpieper") {
      message(paste0("    🌾 Open zand & stuifduin workflow voor: ", soort_naam))
      res_model <- bouw_biotoop_open_zand_stuifduin(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (bestandsnaam %in% c("dwergstern", "grotestern", "visdief")) {
      message(paste0("    🪶 Kust-stern workflow voor: ", soort_naam))
      res_model <- bouw_biotoop_kust_stern(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (bestandsnaam == "grauwekiekendief") {
      message(paste0("    🦅 Grootschalig akkercomplex workflow voor: ", soort_naam))
      res_model <- bouw_biotoop_grootschalig_akkercomplex(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (grepl("veldkrekel", tolower(soort_naam))) {
      source(here::here("src", "Functies", "bouw_biotoop_droge_metapopulatie_sprinkhaan.R"))
      res <- bouw_biotoop_droge_metapopulatie_sprinkhaan(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (grepl("zadelsprinkhaan", tolower(soort_naam))) {
      source(here::here("src", "Functies", "bouw_biotoop_droge_heide_sprinkhaan.R"))
      res <- bouw_biotoop_droge_heide_sprinkhaan(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (grepl("zompsprinkhaan|moerassprinkhaan", tolower(soort_naam))) {
      source(here::here("src", "Functies", "bouw_biotoop_natte_metapopulatie_sprinkhaan.R"))
      res <- bouw_biotoop_natte_metapopulatie_sprinkhaan(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (grepl("vliegendhert", tolower(soort_naam))) {
      source(here::here("src", "Functies", "bouw_biotoop_saproxyl_hellingbos.R"))
      res <- bouw_biotoop_saproxyl_hellingbos(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (grepl("juchtleerkever|osmoderma", tolower(soort_naam))) {
      source(here::here("src", "Functies", "bouw_biotoop_saproxyl_kle_park.R"))
      res <- bouw_biotoop_saproxyl_kle_park(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (grepl("vermiljoenkever|cucujus", tolower(soort_naam))) {
      source(here::here("src", "Functies", "bouw_biotoop_saproxyl_alluviaal_bos.R"))
      res <- bouw_biotoop_saproxyl_alluviaal_bos(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (bestandsnaam %in% c("bontdikkopje", "bruineeikenpage") || "bosranden" %in% tolower(params_soort$Type)) {
      message(paste0("    🦋 Ecotonale bosrand-workflow voor: ", soort_naam))
      res_model <- bouw_biotoop_ecotoon(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths)
      
    } else if (nrow(spoor_a) > 1) {
      # ========================================================================
      # ROUTE A-MEERVOUDIG: HISTORISCHE AMPHIBIEEN/MEERVOUDIGE BIOTOPEN
      # ========================================================================
      message("    🧩 Meervoudig biotoop gedetecteerd. Complex biotoop opbouwen...")
      comp_res <- bouw_biotoop_complex(soort_naam, spoor_a, ds_parquet, df_bwk_codes, master_grid, paths)
      
      r_opp   <- comp_res$r_gecombineerd
      ha_voor <- calc_ha_exact(r_opp)
      message(paste0("    Oppervlakte vóór clustering: ", round(ha_voor, 2), " ha"))
      
      drempel_netwerk <- if (!is.na(min_totale_opp) && min_totale_opp > 0) min_totale_opp else min_opp_ha
      
      res <- cluster_filter_universeel(
        masker            = r_opp,
        opp_laag          = r_opp,
        min_opp_ha        = drempel_netwerk,
        min_totale_opp_ha = NA,
        afstand_m         = dispersie_m,
        dispersie_m       = 0,
        werkelijk         = TRUE
      )
      final_opp <- res$raster
      cl_opp    <- res$clusters
      
      # Optionele poelentoets
      heeft_plassentoets <- any(as.logical(spoor_a$Gebruik_Watervlakken_Shp), na.rm = TRUE) &&
        any(!is.na(spoor_a$Min_Aantal_Poelen) & suppressWarnings(as.numeric(spoor_a$Min_Aantal_Poelen)) > 0)
      
      if (heeft_plassentoets && !is.null(comp_res$kleine_plassen_v) && !all(is.na(suppressWarnings(terra::minmax(cl_opp))))) {
        min_poelen  <- max(suppressWarnings(as.numeric(spoor_a$Min_Aantal_Poelen)), na.rm = TRUE)
        plassen_pts <- terra::centroids(comp_res$kleine_plassen_v)
        plassen_ext <- terra::extract(cl_opp, plassen_pts, buffer = 50, fun = "max", na.rm = TRUE)
        plassen_pts$cluster_ID <- as.integer(plassen_ext[[2]])
        
        dt_pl     <- data.table::as.data.table(as.data.frame(plassen_pts))
        valid_ids <- dt_pl[!is.na(cluster_ID) & cluster_ID > 0, .N, by = cluster_ID][N >= min_poelen, cluster_ID]
        
        if (length(valid_ids) > 0) {
          try(terra::set.cats(cl_opp, layer = 1, value = NULL), silent = TRUE)
          if (terra::has.colors(cl_opp)) try(terra::coltab(cl_opp) <- NULL, silent = TRUE)
          cl_opp    <- terra::ifel(cl_opp %in% as.integer(valid_ids), cl_opp, NA)
          final_opp <- terra::mask(final_opp, cl_opp)
        } else {
          cl_opp    <- master_grid * NA
          final_opp <- master_grid * NA
        }
      }
      
      res_model <- list(
        final_opp     = final_opp,
        cl_opp        = cl_opp,
        landbiotoop   = terra::mask(comp_res$r_land_deel, cl_opp),
        waterbiotoop  = terra::mask(comp_res$r_water_deel, cl_opp)
      )
    }
    
    # ==========================================================================
    # ROUTE B: STANDAARD ROUTE (Middelste bonte specht, Blauwborst, etc.)
    # ==========================================================================
    if (is.null(res_model)) {
      message(paste0("    🌲 Standaard Route B workflow voor: ", soort_naam))
      type_naam <- hoofd_rij$Type
      r_opp     <- extract_bwk_laag(soort_naam, type_naam, ds_parquet, df_bwk_codes, master_grid)
      
      if (is.null(r_opp) || all(is.na(suppressWarnings(terra::minmax(r_opp))))) {
        message("    ⚠️ Geen BWK data gevonden. Leeg raster wordt weggeschreven.")
        final_opp <- master_grid * NA
        cl_opp    <- master_grid * NA
      } else {
        ha_voor <- calc_ha_exact(r_opp)
        message(paste0("    Oppervlakte vóór clustering: ", round(ha_voor, 2), " ha"))
        
        res <- cluster_filter_universeel(
          masker            = r_opp,
          opp_laag          = r_opp,
          min_opp_ha        = min_opp_ha,
          min_totale_opp_ha = min_totale_opp,
          afstand_m         = afstand_m,
          dispersie_m       = dispersie_m,
          werkelijk         = TRUE
        )
        final_opp <- res$raster
        cl_opp    <- res$clusters
      }
    } else {
      # Verwerk resultaten van Route A
      final_opp <- res_model$final_opp
      cl_opp    <- res_model$cl_opp
      
      # Automatische, universele export van deellagen indien aanwezig
      deellaag_namen <- setdiff(names(res_model), c("final_opp", "cl_opp"))
      for (dl in deellaag_namen) {
        if (!is.null(res_model[[dl]])) {
          f_sub <- file.path(dir_sub, paste0(bestandsnaam, "_", dl, ".tif"))
          terra::writeRaster(res_model[[dl]], f_sub, overwrite = TRUE, datatype = "FLT4S", gdal = "COMPRESS=LZW")
        }
      }
      if (length(deellaag_namen) > 0) {
        message(paste0("    ✅ [OK] Deellagen opgeslagen in Deellagen/: ", paste(deellaag_namen, collapse = ", ")))
      }
    }
    
    ha_na <- calc_ha_exact(final_opp)
    message(paste0("    Oppervlakte ná clustering: ", round(ha_na, 2), " ha"))
    
    # Exporteren basisrasters
    terra::writeRaster(final_opp, f_opp, overwrite = TRUE, gdal = c("COMPRESS=LZW", "BIGTIFF=YES"), datatype = "FLT4S", NAflag = -9999)
    terra::writeRaster(cl_opp,    f_id,  overwrite = TRUE, gdal = c("COMPRESS=LZW", "BIGTIFF=YES"), datatype = "INT4U", NAflag = 0)
    message(paste0("    ✅ [OK] Opgeslagen:\n      - ", basename(f_opp), "\n      - ", basename(f_id)))
    
    # Automatische opkuis van persistente cache na succesvolle export
    cache_soort <- file.path(paths$dir_cache, bestandsnaam)
    if (dir.exists(cache_soort)) {
      unlink(cache_soort, recursive = TRUE)
      message("    🧹 [CACHE] Tussenstappen succesvol opgeruimd.")
    }
  }
  
  # ----------------------------------------------------------------------------
  # DEEL 2: ARPL BEREKENEN (02_ARPL)
  # ----------------------------------------------------------------------------
  if (exists("bereken_en_exporteer_arpl", mode = "function")) {
    if (!force_overwrite && file.exists(f_arpl)) {
      message(paste0("    ⏩ [OVERSLAAN ARPL] ", basename(f_arpl), " bestaat reeds."))
    } else {
      rij         <- params_soort %>% dplyr::slice(1)
      dispersie_m <- suppressWarnings(as.numeric(rij$Dispersiecap_m[1]))
      is_wv       <- grepl("_wv$", tolower(soort_naam))
      
      bereken_en_exporteer_arpl(
        soort_naam      = soort_naam,
        r_id            = cl_opp,
        r_opp           = final_opp,
        dispersie_m     = dispersie_m,
        paths           = paths,
        is_wintervogel  = is_wv
      )
    }
  }
  
  rm(cl_opp, final_opp)
  gc()
  return(invisible(TRUE))
}
