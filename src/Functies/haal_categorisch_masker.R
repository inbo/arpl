haal_categorisch_masker <- function(pad_raster, codes, master_grid, dir_cache = NULL, prefix = "cat") {
  codes_clean <- trimws(as.character(codes))
  codes_clean <- codes_clean[!is.na(codes_clean) & codes_clean != ""]
  if (length(codes_clean) == 0 || is.null(pad_raster) || !file.exists(pad_raster)) {
    return(master_grid * NA_real_)
  }
  
  f_al <- hulp_align_gdal(pad_raster, master_grid)
  r    <- terra::rast(f_al)
  
  ct <- terra::cats(r)[[1]]
  if (!is.null(ct) && nrow(ct) > 0) {
    id_kol  <- colnames(ct)[1]
    lbl_kol <- if (ncol(ct) >= 2) colnames(ct)[2] else colnames(ct)[1]
    
    # Match op Label (case-insensitive) of op ID
    match_lbl <- tolower(trimws(as.character(ct[[lbl_kol]]))) %in% tolower(codes_clean)
    match_id  <- trimws(as.character(ct[[id_kol]])) %in% codes_clean
    target_ids <- ct[[id_kol]][match_lbl | match_id]
  } else {
    # Puur numeriek raster zonder attribuuttabel (zoals Vlaanderen_textuur_10m.tif)
    target_ids <- suppressWarnings(as.numeric(codes_clean))
    target_ids <- target_ids[!is.na(target_ids)]
  }
  
  if (length(target_ids) > 0) {
    m_out <- terra::ifel(r %in% target_ids, 1L, NA_integer_)
  } else {
    m_out <- master_grid * NA_real_
  }
  
  unlink(f_al)
  rm(r)
  return(m_out)
}
