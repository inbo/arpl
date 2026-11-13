# ==============================================================================
# haal_categorisch_masker.R
# Universele extractor voor categorische hulpkaarten (drainage, profiel, huetzon, etc.)
# ==============================================================================

haal_categorisch_masker <- function(pad_raster, 
                                    codes_selectie, 
                                    master_grid, 
                                    dir_cache, 
                                    naam_prefix = "masker") {
  
  if (is.null(pad_raster) || !file.exists(pad_raster)) {
    return(master_grid * NA_real_)
  }
  
  # Genereer een unieke cache-naam op basis van gevraagde codes
  codes_clean <- gsub("[^a-zA-Z0-9]", "_", tolower(trimws(as.character(codes_selectie))))
  hash_codes  <- paste(sort(codes_clean), collapse = "_")
  # Beperk bestandsnaamlengte indien er veel codes zijn
  if (nchar(hash_codes) > 50) {
    hash_codes <- digest::digest(hash_codes, algo = "crc32")
  }
  
  f_cache <- file.path(dir_cache, sprintf("%s_%s.tif", naam_prefix, hash_codes))
  
  if (file.exists(f_cache)) {
    return(terra::rast(f_cache))
  }
  
  message(sprintf("    🛠️ Genereren van %s voor codes: %s...", naam_prefix, paste(codes_selectie, collapse = ", ")))
  
  # 1. Inspecteer de categorieëntabel van het bronbestand vóór warping
  r_orig  <- terra::rast(pad_raster)
  leg     <- terra::cats(r_orig)[[1]]
  ok_ids  <- integer(0)
  
  codes_zoek_chr <- tolower(trimws(as.character(codes_selectie)))
  
  if (!is.null(leg)) {
    kol_val   <- colnames(leg)[tolower(colnames(leg)) %in% c("value", "id")][1]
    kol_label <- colnames(leg)[tolower(colnames(leg)) %in% c("label", "omschrijving", "naam", "profiel", "drainage", "huet")][1]
    if (is.na(kol_val)) kol_val <- 1
    
    # Check zowel match op tekstlabel als numerieke value
    if (!is.na(kol_label)) {
      labels_raw <- tolower(trimws(as.character(leg[[kol_label]])))
      ok_ids <- c(ok_ids, leg[[kol_val]][labels_raw %in% codes_zoek_chr])
    }
    values_raw <- tolower(trimws(as.character(leg[[kol_val]])))
    ok_ids <- c(ok_ids, leg[[kol_val]][values_raw %in% codes_zoek_chr])
    ok_ids <- unique(na.omit(ok_ids))
  }
  
  # Fallback: als er geen VAT-tabel is, filter direct op de getallen
  if (length(ok_ids) == 0) {
    ok_ids <- suppressWarnings(as.integer(codes_zoek_chr[!is.na(as.numeric(codes_zoek_chr))]))
  }
  rm(r_orig)
  
  if (length(ok_ids) == 0) {
    message(sprintf("    ⚠️ Geen overeenkomende IDs gevonden in %s!", basename(pad_raster)))
    return(master_grid * NA_real_)
  }
  
  # 2. Warp via GDAL
  f_aligned <- hulp_align_gdal(pad_raster, master_grid)
  r_geom    <- terra::rast(f_aligned)
  
  mask_res  <- terra::ifel(!is.na(r_geom) & r_geom %in% ok_ids, 1.0, NA_real_)
  unlink(f_aligned)
  rm(r_geom)
  
  # 3. Persistent bewaren
  terra::writeRaster(mask_res, f_cache, overwrite = TRUE, datatype = "INT1U", NAflag = 255, gdal = "COMPRESS=LZW")
  
  return(mask_res)
}