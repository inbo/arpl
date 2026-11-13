# ==============================================================================
# extract_bwk_laag.R
# Geheugenveilige extractie van BWK-Parquet gegevens
# Ondersteunt zowel heel Vlaanderen als willekeurige subgebieden (crop / AOI)
# ==============================================================================

extract_bwk_laag <- function(soort_naam, tp, ds_parquet, df_bwk_codes, master_grid, paths = NULL, force_overwrite = FALSE) {
  
  bestandsnaam <- tolower(gsub(" ", "", soort_naam))
  tp_clean     <- tolower(trimws(tp))
  
  # ----------------------------------------------------------------------------
  # 1. BESTANDSPADEN EN SCHIJFCACHE
  # ----------------------------------------------------------------------------
  if (!is.null(paths$dir_cache)) {
    soort_cache_dir <- file.path(paths$dir_cache, bestandsnaam)
    dir.create(soort_cache_dir, showWarnings = FALSE, recursive = TRUE)
    f_out_tif <- file.path(soort_cache_dir, sprintf("raw_bwk_%s.tif", tp_clean))
  } else {
    f_out_tif <- tempfile(fileext = ".tif")
  }
  
  # Indien reeds gecached en compatibel met huidig master_grid, direct inlezen
  if (file.exists(f_out_tif) && !force_overwrite) {
    r_cached <- suppressWarnings(terra::rast(f_out_tif))
    if (terra::compareGeom(r_cached, master_grid, stopOnError = FALSE)) {
      message(sprintf("       📦 [CACHE] Ruw Parquet-raster ingeladen van schijf: %s", basename(f_out_tif)))
      return(r_cached)
    }
  }
  
  message(sprintf("       ⏳ Parquet-laag '%s' extraheren voor %s...", tp, soort_naam))
  
  # ----------------------------------------------------------------------------
  # 2. BWK-CODES OPZOEKEN IN DE REFERENTIETABEL
  # ----------------------------------------------------------------------------
  codes_voor_type <- df_bwk_codes %>%
    dplyr::mutate(Soort_clean = tolower(trimws(Soort)),
                  Type_clean  = tolower(trimws(Type))) %>%
    dplyr::filter(Soort_clean == bestandsnaam & Type_clean == tp_clean)
  
  if (nrow(codes_voor_type) == 0) {
    message(sprintf("       ⚠️ Geen BWK-codes gevonden voor type '%s' bij soort '%s'. Leeg raster geretourneerd.", tp, soort_naam))
    return(master_grid * NA_real_)
  }
  
  exact_codes <- tolower(trimws(codes_voor_type$Code[codes_voor_type$Match == "exact"]))
  bevat_codes <- tolower(trimws(codes_voor_type$Code[codes_voor_type$Match == "bevat"]))
  
  # ----------------------------------------------------------------------------
  # 3. FILTEREN UIT DE ARROW / PARQUET DATASET
  # ----------------------------------------------------------------------------
  if (length(bevat_codes) > 0) {
    bevat_codes_escaped  <- gsub("([\\.\\^\\$\\*\\+\\?\\(\\)\\[\\{\\\\\\|])", "\\\\\\1", bevat_codes)
    bevat_codes_anchored <- paste0("^", bevat_codes_escaped)
    regex_term           <- paste0(bevat_codes_anchored, collapse = "|")
    
    tabel_sub <- ds_parquet %>%
      dplyr::filter(tolower(CODE) %in% exact_codes | stringr::str_detect(tolower(CODE), regex_term)) %>%
      dplyr::collect() %>%
      data.table::as.data.table()
  } else {
    tabel_sub <- ds_parquet %>%
      dplyr::filter(tolower(CODE) %in% exact_codes) %>%
      dplyr::collect() %>%
      data.table::as.data.table()
  }
  
  if (nrow(tabel_sub) == 0) {
    message(sprintf("       ℹ️ Geen records aangetroffen in Parquet voor '%s'.", tp))
    return(master_grid * NA_real_)
  }
  
  # Aggregeer oppervlakte per cel_id (maximaal 1.0 aandeel per rastercel)
  tabel_sub[, CODE := tolower(trimws(CODE))]
  tabel_unique  <- unique(tabel_sub, by = c("cel_id", "CODE"))
  tabel_cel_som <- tabel_unique[, .(Oppervlakte = pmin(sum(BWK_FRAC, na.rm = TRUE), 1.0)), by = .(cel_id)]
  rm(tabel_unique, tabel_sub)
  gc()
  
  # ----------------------------------------------------------------------------
  # 4. ROBUUSTE CEL-INDEXERING: VLAANDEREN VS SUBGEBIED
  # ----------------------------------------------------------------------------
  # Maak een schone raster-header aan zónder schijfbelasting
  r_opp_type <- terra::rast(master_grid)
  terra::values(r_opp_type) <- NA_real_
  
  # Referentiekader van heel Vlaanderen ophalen
  pad_vl_template <- if (!is.null(paths$master_grid_vlaanderen) && file.exists(paths$master_grid_vlaanderen)) {
    paths$master_grid_vlaanderen
  } else if (!is.null(paths$master_grid) && file.exists(paths$master_grid)) {
    paths$master_grid
  } else {
    NULL
  }
  
  r_vl_ref <- if (!is.null(pad_vl_template)) {
    terra::rast(pad_vl_template)[[1]]
  } else {
    terra::rast(
      xmin = 20000, xmax = 260000, ymin = 150000, ymax = 250000,
      resolution = terra::res(master_grid),
      crs = terra::crs(master_grid)
    )
  }
  
  is_subgebied <- !terra::compareGeom(master_grid, r_vl_ref, stopOnError = FALSE)
  
  if (is_subgebied) {
    message("       🔍 Subgebied gedetecteerd: cel_id's transformeren naar lokale uitsnede...")
    coords_xy <- terra::xyFromCell(r_vl_ref, tabel_cel_som$cel_id)
    lokale_cel_ids <- terra::cellFromXY(master_grid, coords_xy)
    rm(coords_xy)
    
    geldige_idx <- which(!is.na(lokale_cel_ids))
    if (length(geldige_idx) > 0) {
      r_opp_type[lokale_cel_ids[geldige_idx]] <- tabel_cel_som$Oppervlakte[geldige_idx]
    }
    rm(lokale_cel_ids, geldige_idx)
  } else {
    max_cells <- terra::ncell(master_grid)
    tabel_cel_som <- tabel_cel_som[!is.na(cel_id) & cel_id >= 1 & cel_id <= max_cells]
    if (nrow(tabel_cel_som) > 0) {
      r_opp_type[tabel_cel_som$cel_id] <- tabel_cel_som$Oppervlakte
    }
  }
  
  rm(tabel_cel_som)
  gc()
  
  # ----------------------------------------------------------------------------
  # 5. WEGSCHRIJVEN NAAR PERSISTENTE SCHIJFCACHE (LZW GECOMPRIMEERD)
  # ----------------------------------------------------------------------------
  dir.create(dirname(f_out_tif), showWarnings = FALSE, recursive = TRUE)
  r_stabiel <- terra::writeRaster(
    r_opp_type,
    filename  = f_out_tif,
    overwrite = TRUE,
    datatype  = "FLT4S",
    gdal      = c("COMPRESS=LZW", "TILED=YES"),
    NAflag    = -9999
  )
  rm(r_opp_type)
  gc()
  
  ha_totaal <- calc_ha_exact(r_stabiel)
  message(sprintf("       📊 Ruwe oppervlakte '%s': %.2f ha (veilig bewaard op schijf)", tp, ha_totaal))
  
  return(r_stabiel)
}
