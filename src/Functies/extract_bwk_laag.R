extract_bwk_laag <- function(soort_naam, type_laag, ds_parquet, df_bwk_codes, master_grid) {
  # Maak een maagdelijk leeg raster puur op basis van dimensies/CRS (geen disk-pointer!)
  r_doel <- terra::rast(
    extent  = terra::ext(master_grid),
    res     = terra::res(master_grid),
    crs     = terra::crs(master_grid),
    vals    = NA
  )
  
  # Filter de codes voor deze soort en type
  df_codes_sub <- df_bwk_codes %>%
    dplyr::filter(tolower(trimws(Soort)) == tolower(trimws(soort_naam)),
                  tolower(trimws(Type))  == tolower(trimws(type_laag)))
  
  if (nrow(df_codes_sub) == 0) return(r_doel)
  
  codes_exact <- tolower(trimws(df_codes_sub$Code[tolower(df_codes_sub$Match) == "exact"]))
  codes_bevat <- tolower(trimws(df_codes_sub$Code[tolower(df_codes_sub$Match) == "bevat"]))
  
  # Parquet query
  t_sub <- if (length(codes_bevat) > 0) {
    regex_str <- paste0("^", paste(codes_bevat, collapse = "|^"))
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
  
  if (nrow(t_cel) > 0) {
    # Alleen cellen binnen het actieve grid schrijven
    n_tot <- terra::ncell(r_doel)
    t_valide <- t_cel[cel_id >= 1 & cel_id <= n_tot]
    if (nrow(t_valide) > 0) {
      r_doel[t_valide$cel_id] <- t_valide$Opp
    }
  }
  
  return(r_doel)
}