# ==============================================================================
# bouw_biotoop_hamster_leemakker.R
# Archetype voor Europese hamster: Graangewassen x Bodemtextuur Leem (6) x Niet-overstroombaar (>= 125 ha op 50m)
# Conform officiële INBO-standaard ANB2019_EuropeseHamster2018.asc
# ==============================================================================

bouw_biotoop_hamster_leemakker <- function(soort_naam, params_soort, ds_parquet, df_bwk_codes, master_grid, paths) {
  
  bestandsnaam <- tolower(gsub(" ", "", soort_naam))
  soort_cache  <- file.path(paths$dir_cache, bestandsnaam)
  dir.create(soort_cache, showWarnings = FALSE, recursive = TRUE)
  
  f_cache_rds <- file.path(soort_cache, "01_raw_parquet_layers.rds")
  
  # Parameters ophalen
  rij_params    <- params_soort %>% dplyr::slice(1)
  min_opp_ha    <- as.numeric(ifelse(is.na(rij_params$MinOpp_ha), 125, rij_params$MinOpp_ha))
  zoekafstand_m <- as.numeric(ifelse(is.na(rij_params$AfstandBiotopen_m), 50, rij_params$AfstandBiotopen_m))
  
  # ----------------------------------------------------------------------------
  # 1. GEWASSENLAAG (LGP / EPR GRAANGEWASSEN)
  # ----------------------------------------------------------------------------
  message("    🌾 Landbouwgewassen filteren op wintergranen en vlinderbloemigen...")
  
  codes_graan <- c(6, 19, 21, 36, 38, 45, 49, 63, 74, 100, 174, 219, 225, 241)
  if (!is.na(rij_params$LGP_Gewascodes)) {
    codes_graan <- as.numeric(trimws(unlist(strsplit(as.character(rij_params$LGP_Gewascodes), ","))))
  }
  
  f_cache_graan <- file.path(paths$dir_cache, "masker_graangewassen_hamster.tif")
  
  if (file.exists(f_cache_graan)) {
    r_graan_bin <- terra::rast(f_cache_graan)
  } else if (!is.null(paths$lgp) && file.exists(paths$lgp)) {
    f_lgp_al    <- hulp_align_gdal(paths$lgp, master_grid)
    r_lgp       <- terra::rast(f_lgp_al)
    r_graan_bin <- terra::ifel(!is.na(r_lgp) & (r_lgp %in% codes_graan), 1, NA)
    unlink(f_lgp_al)
    rm(r_lgp)
    terra::writeRaster(r_graan_bin, f_cache_graan, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
  } else {
    # Fallback: BWK akkerfractie
    message("       ℹ️ Geen LGP-raster gevonden op paths$lgp; fallback naar akkers via BWK/Parquet...")
    r_akker <- extract_bwk_laag(soort_naam, "akker", ds_parquet, df_bwk_codes, master_grid)
    r_graan_bin <- if (!is.null(r_akker)) terra::ifel(!is.na(r_akker) & r_akker > 0, 1, NA) else (master_grid * NA)
  }
  
  # ----------------------------------------------------------------------------
  # 2. BODEMTEXTUUR (LEEM = KLASSE 6)
  # ----------------------------------------------------------------------------
  message("    🌱 Bodemtextuur filteren op Leem (klasse 6)...")
  
  f_cache_leem <- file.path(paths$dir_cache, "masker_bodem_leem_kl6.tif")
  
  if (file.exists(f_cache_leem)) {
    r_leem_bin <- terra::rast(f_cache_leem)
  } else if (!is.null(paths$textuur) && file.exists(paths$textuur)) {
    f_text_al  <- hulp_align_gdal(paths$textuur, master_grid)
    r_text     <- terra::rast(f_text_al)
    
    # Klasse 6 is Leem conform Vlaamse bodemkaart
    r_leem_bin <- terra::ifel(!is.na(r_text) & r_text == 6, 1, NA)
    unlink(f_text_al)
    rm(r_text)
    terra::writeRaster(r_leem_bin, f_cache_leem, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
  } else if (!is.null(paths$bodem) && file.exists(paths$bodem)) {
    # Alternatief pad naar bodemkaart
    f_bod_al   <- hulp_align_gdal(paths$bodem, master_grid)
    r_bod      <- terra::rast(f_bod_al)
    r_leem_bin <- terra::ifel(!is.na(r_bod) & r_bod == 6, 1, NA)
    unlink(f_bod_al)
    rm(r_bod)
    terra::writeRaster(r_leem_bin, f_cache_leem, overwrite = TRUE, datatype = "INT1U", NAflag = 255)
  } else {
    message("       ⚠️ Geen bodemtextuurkaart gevonden; leemconditie niet toegepast.")
    r_leem_bin <- master_grid * 1
  }
  
  # ----------------------------------------------------------------------------
  # 3. OVERSTROMINGSFILTER (UITSLUITEND NIET-OVERSTROOMBAAR = KLASSE 0)
  # ----------------------------------------------------------------------------
  message("    🚫 Inundatiegevoelige gronden uitsluiten (enkel niet-overstroombaar)...")
  
  if (!is.null(paths$overstroming) && file.exists(paths$overstroming)) {
    f_ovstr_al <- hulp_align_gdal(paths$overstroming, master_grid)
    r_ovstr    <- terra::rast(f_ovstr_al)
    # Alleen klasse 0 of NoData/0 overstromingskans is toegestaan
    mask_niet_overstroombaar <- terra::ifel(is.na(r_ovstr) | r_ovstr == 0, 1, NA)
    unlink(f_ovstr_al)
    rm(r_ovstr)
  } else {
    mask_niet_overstroombaar <- master_grid * 1
  }
  
  # ----------------------------------------------------------------------------
  # 4. SAMENVOEGEN BASISBIOTOOP (GRAAN & LEEM & NIET-OVERSTROOMBAAR)
  # ----------------------------------------------------------------------------
  r_hamster_basis <- terra::ifel(!is.na(r_graan_bin) & !is.na(r_leem_bin) & !is.na(mask_niet_overstroombaar), 1.0, NA)
  rm(r_graan_bin, r_leem_bin, mask_niet_overstroombaar)
  gc()
  
  ha_basis <- calc_ha_exact(r_hamster_basis)
  message(sprintf("       -> Potentieel hamsterareaal (graan op leem, droog): %.2f ha", ha_basis))
  
  if (ha_basis == 0) {
    message("       ⚠️ Geen cellen gevonden die voldoen aan de combinatie graan x leem.")
    return(list(final_opp = master_grid * NA, cl_opp = master_grid * NA))
  }
  
  # ----------------------------------------------------------------------------
  # 5. CLUSTERING OP 50M EN MINIMUMOPPERVLAKTE (>= 125 HA)
  # ----------------------------------------------------------------------------
  message(sprintf("    🏝️ Akkercomplexen clusteren op %dm en filteren op >= %.1f ha...", as.integer(zoekafstand_m), min_opp_ha))
  
  res_cluster <- cluster_filter_universeel(
    masker            = r_hamster_basis,
    opp_laag          = r_hamster_basis,
    min_opp_ha        = min_opp_ha,
    min_totale_opp_ha = NA,
    afstand_m         = zoekafstand_m,
    dispersie_m       = 0,
    werkelijk         = TRUE
  )
  
  final_opp <- res_cluster$raster
  cl_opp    <- res_cluster$clusters
  rm(res_cluster, r_hamster_basis)
  gc()
  
  ha_finaal <- calc_ha_exact(final_opp)
  message(sprintf("       -> Definitief leefgebied Europese hamster: %.2f ha", ha_finaal))
  
  return(list(
    final_opp = final_opp,
    cl_opp    = cl_opp
  ))
}
