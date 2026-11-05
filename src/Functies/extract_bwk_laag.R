extract_bwk_laag <- function(soort_naam, type_laag, ds_parquet, df_bwk_codes, master_grid) {
  # Maak doelraster
  r_doel <- terra::rast(
    extent  = terra::ext(master_grid),
    res     = terra::res(master_grid),
    crs     = terra::crs(master_grid),
    vals    = NA_real_
  )
  
  # Filter codes voor soort en type
  df_codes_sub <- df_bwk_codes %>%
    dplyr::filter(tolower(trimws(Soort)) == tolower(trimws(soort_naam)),
                  tolower(trimws(Type))  == tolower(trimws(type_laag)))
  
  if (nrow(df_codes_sub) == 0) return(r_doel)
  
  codes_exact <- tolower(trimws(df_codes_sub$Code[tolower(df_codes_sub$Match) == "exact"]))
  codes_bevat <- tolower(trimws(df_codes_sub$Code[tolower(df_codes_sub$Match) == "bevat"]))
  
  # Parquet extractie met veilige regex escaping
  t_sub <- if (length(codes_bevat) > 0) {
    bevat_codes_escaped  <- gsub("([.\\^$*+?()\\[{\\\\|])", "\\\\\\1", codes_bevat)
    regex_str           <- paste0("^", paste(bevat_codes_escaped, collapse = "|^"))
    
    ds_parquet %>%
      dplyr::filter(tolower(CODE) %in% codes_exact | grepl(regex_str, tolower(CODE))) %>%
      dplyr::collect() %>%
      data.table::as.data.table()
  } else {
    ds_parquet %>%
      dplyr::filter(tolower(CODE) %in% codes_exact) %>%
      dplyr::collect() %>%
      data.table::as.data.table()
  }
  
  if (nrow(t_sub) == 0) return(r_doel)
  
  t_sub[, CODE := tolower(trimws(CODE))]
  t_cel <- unique(t_sub, by = c("cel_id", "CODE"))[, .(Opp = pmin(sum(BWK_FRAC, na.rm = TRUE), 1.0)), by = cel_id]
  t_cel <- t_cel[Opp >= 0.01]
  
  if (nrow(t_cel) == 0) return(r_doel)
  
  # ----------------------------------------------------------------------------
  # CEL-ID VERTALING: GLOBALE VLAANDEREN MASTERGRID ID -> LOKALE CELLEN
  # ----------------------------------------------------------------------------
  # Vlaanderen MasterGrid 10m oorsprong & resolutie:
  # xmin: 22000, ymax: 245000, resolutie: 10m
  # Aantal kolommen Vlaanderen grid = (xmax - xmin) / 10 = (260000 - 22000) / 10 = 23800
  # Indien globale cel-ID's gebruikt worden:
  
  vg_xmin <- 22000
  vg_ymax <- 245000
  vg_res  <- 10
  
  # Controleer of cel_id groter is dan het lokale raster (dan zijn het globale IDs)
  if (max(t_cel$cel_id, na.rm = TRUE) > terra::ncell(r_doel)) {
    # Bereken X en Y van globale Vlaanderen cel
    # cel_id = (row - 1) * ncol + col
    # We gebruiken het globale Vlaanderen raster als referentie (of uit paths$master_grid)
    if (exists("vlaanderen_master_grid") && !is.null(vlaanderen_master_grid)) {
      xy_pts <- terra::xyFromCell(vlaanderen_master_grid, t_cel$cel_id)
    } else {
      # Direct via wiskundige projectie van het 10m Vlaanderen grid
      # Breedte Vlaanderen raster: 238000 m / 10 = 23800 cols (of haal ncol uit paths)
      r_vlaanderen_ref <- terra::rast(here::here("data/input/Raster_Vlaanderen/Vlaanderen_MasterGrid_10m.tif"))[[1]]
      xy_pts <- terra::xyFromCell(r_vlaanderen_ref, t_cel$cel_id)
    }
    
    # Bepaal welke cellen binnen de uitsnede van r_doel vallen
    lokale_cellen <- terra::cellFromXY(r_doel, xy_pts)
    valide_idx <- !is.na(lokale_cellen)
    
    if (any(valide_idx)) {
      r_doel[lokale_cellen[valide_idx]] <- t_cel$Opp[valide_idx]
    }
  } else {
    # Cel_ids zijn al lokaal
    r_doel[t_cel$cel_id] <- t_cel$Opp
  }
  
  return(r_doel)
}
