# ==============================================================================
# bereken_arpl.R - ACTUEEL EN RECENT POTENTIEEL LEEFGEBIED (ARPL) MAKER
# ==============================================================================

# Vaste lijst met broedvogels waarvoor territoria primeren
BROEDVOGELS_LIJST <- c(
  "blauwborst", "boomleeuwerik", "boompieper", "bruinekiekendief", "fluiter",
  "grauweklauwier", "grutto", "ijsvogel", "kwak", "kwartelkoning", "matkop",
  "middelstebontespecht", "nachtegaal", "nachtzwaluw", "paapje", "porseleinhoen", 
  "roerdomp", "tapuit", "watersnip", "wespendief", "wielewaal", "woudaap", 
  "wulp", "zomertortel", "zwartespecht", "zwartkopmeeuw"
)

#' Laadt waarnemingen of territoria voor een soort als SpatVector punten
laad_waarnemingen_soort <- function(soort_naam, dir_waarnemingen) {
  if (!dir.exists(dir_waarnemingen)) {
    message("    ⚠️ Map met waarnemingen niet gevonden: ", dir_waarnemingen)
    return(NULL)
  }
  
  s_clean <- tolower(gsub("[_ ]", "", trimws(soort_naam)))
  is_broedvogel <- s_clean %in% BROEDVOGELS_LIJST
  
  alle_bestanden <- list.files(dir_waarnemingen, pattern = "\\.csv$", full.names = TRUE)
  bestanden_clean <- tolower(gsub("[_ ]", "", basename(alle_bestanden)))
  
  # Filter alle bestanden die bij deze soort horen
  matches <- alle_bestanden[grepl(s_clean, bestanden_clean)]
  
  if (length(matches) == 0) {
    message("    ℹ️ Geen waarnemingen- of territoria-CSV gevonden voor: ", soort_naam)
    return(NULL)
  }
  
  # Selecteer bestand: voor broedvogels eerst zoeken naar 'territoria'
  gekozen_bestand <- NULL
  if (is_broedvogel) {
    match_terr <- matches[grepl("territori", tolower(basename(matches)))]
    if (length(match_terr) > 0) gekozen_bestand <- match_terr[1]
  }
  
  # Indien geen territoria gevonden (of geen broedvogel), neem de reguliere waarnemingen
  if (is.null(gekozen_bestand)) {
    gekozen_bestand <- matches[1]
  }
  
  # Inlezen via data.table (flexibel met delimiters)
  dt <- tryCatch({
    data.table::fread(gekozen_bestand, fill = TRUE, encoding = "UTF-8")
  }, error = function(e) {
    message("    ❌ Fout bij inlezen CSV: ", basename(gekozen_bestand), " (", e$message, ")")
    return(NULL)
  })
  
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  
  data.table::setnames(dt, tolower(names(dt)))
  
  # Coördinatenkolommen detecteren
  x_col <- intersect(c("x", "x_coord", "lon", "easting"), names(dt))[1]
  y_col <- intersect(c("y", "y_coord", "lat", "northing"), names(dt))[1]
  
  if (is.na(x_col) || is.na(y_col)) {
    message("    ⚠️ Geen x/y coördinaten gevonden in: ", basename(gekozen_bestand))
    return(NULL)
  }
  
  dt_clean <- dt[!is.na(get(x_col)) & !is.na(get(y_col))]
  if (nrow(dt_clean) == 0) return(NULL)
  
  pts <- terra::vect(
    as.matrix(dt_clean[, c(x_col, y_col), with = FALSE]), 
    type = "points", 
    crs  = "EPSG:31370"
  )
  
  label_type <- if (grepl("territori", tolower(basename(gekozen_bestand)))) "Territoria" else "Waarnemingen"
  message(paste0("    📍 ", label_type, " geladen: ", length(pts), " punten (", basename(gekozen_bestand), ")"))
  return(pts)
}

#' Berekent en exporteert het ARPL raster
bereken_en_exporteer_arpl <- function(soort_naam, r_id, r_opp, dispersie_m, paths, is_wintervogel = FALSE) {
  
  bestandsnaam <- tolower(gsub(" ", "", soort_naam))
  f_out_arpl   <- file.path(paths$out_arpl, paste0("Habitat_ARPL_", bestandsnaam, ".tif"))
  
  # 1. Controle op leeg netwerk-raster
  if (is.null(r_id) || all(is.na(suppressWarnings(terra::minmax(r_id))))) {
    r_leeg <- r_id * NA
    terra::writeRaster(r_leeg, f_out_arpl, overwrite = TRUE, datatype = "INT1U", NAflag = 255, gdal = "COMPRESS=LZW")
    return(f_out_arpl)
  }
  
  # 2. Wintervogel: werkelijk biotoop is direct ARPL
  if (is_wintervogel) {
    message("    ❄️ Wintervogel: werkelijk habitat wordt direct als ARPL geëxporteerd.")
    r_arpl <- terra::ifel(!is.na(r_opp) & r_opp > 0, 1, NA)
    terra::writeRaster(r_arpl, f_out_arpl, overwrite = TRUE, datatype = "INT1U", NAflag = 255, gdal = "COMPRESS=LZW")
    message(paste0("    ✅ [ARPL] Opgeslagen: ", basename(f_out_arpl)))
    return(f_out_arpl)
  }
  
  # 3. Waarnemingen ophalen
  v_pts <- laad_waarnemingen_soort(soort_naam, paths$dir_waarnemingen)
  
  if (is.null(v_pts) || length(v_pts) == 0) {
    message("    ℹ️ Geen waarnemingen gevonden: ARPL raster blijft leeg.")
    r_leeg <- r_id * NA
    terra::writeRaster(r_leeg, f_out_arpl, overwrite = TRUE, datatype = "INT1U", NAflag = 255, gdal = "COMPRESS=LZW")
    return(f_out_arpl)
  }
  
  # 4. Zeer mobiele soorten (Dispersiecap >= 250 km)
  if (!is.na(dispersie_m) && dispersie_m >= 250000) {
    message("    🦅 Zeer mobiele soort (Dispersiecap >= 250 km): alle actieve clusters worden ARPL.")
    r_arpl <- terra::ifel(!is.na(r_id) & r_id > 0, 1, NA)
    terra::writeRaster(r_arpl, f_out_arpl, overwrite = TRUE, datatype = "INT1U", NAflag = 255, gdal = "COMPRESS=LZW")
    message(paste0("    ✅ [ARPL] Opgeslagen: ", basename(f_out_arpl)))
    return(f_out_arpl)
  }
  
  # 5. Standaard soorten: Waarnemingen bufferen met dispersiecapaciteit
  buffer_afstand <- ifelse(is.na(dispersie_m) || dispersie_m <= 0, 500, dispersie_m)
  message(paste0("    -> Buffering waarnemingen met dispersiecapaciteit: ", buffer_afstand, " m..."))
  
  v_blobs <- tryCatch({
    terra::aggregate(terra::buffer(v_pts, width = buffer_afstand))
  }, error = function(e) NULL)
  
  if (is.null(v_blobs)) {
    r_leeg <- r_id * NA
    terra::writeRaster(r_leeg, f_out_arpl, overwrite = TRUE, datatype = "INT1U", NAflag = 255, gdal = "COMPRESS=LZW")
    return(f_out_arpl)
  }
  
  # 6. Intersectie met netwerk ID's
  message("    -> Netwerk-ID's koppelen aan gebufferde waarnemingen...")
  ext_blobs   <- terra::extract(r_id, v_blobs, ID = FALSE)
  id_col      <- names(ext_blobs)[1]
  ids_arpl    <- unique(na.omit(ext_blobs[[id_col]]))
  
  if (length(ids_arpl) > 0) {
    message(paste0("    -> ", length(ids_arpl), " netwerkcluster(s) geselecteerd voor ARPL."))
    r_arpl <- terra::ifel(r_id %in% ids_arpl, 1, NA)
  } else {
    message("    ℹ️ Geen netwerken geraakt door de gebufferde waarnemingen.")
    r_arpl <- r_id * NA
  }
  
  # 7. Exporteren als 1-byte LZW raster
  terra::writeRaster(r_arpl, f_out_arpl, overwrite = TRUE, datatype = "INT1U", NAflag = 255, gdal = "COMPRESS=LZW")
  message(paste0("    ✅ [ARPL] Succesvol opgeslagen: ", basename(f_out_arpl)))
  
  rm(v_pts, v_blobs, r_arpl)
  gc()
  return(f_out_arpl)
}
