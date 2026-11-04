# ==============================================================================
# extract_bwk_laag.R
# ==============================================================================

extract_bwk_laag <- function(soort_naam, type_naam, ds_parquet, df_bwk_codes, master_grid) {
  codes_subset <- df_bwk_codes %>%
    filter(tolower(trimws(Soort)) == tolower(trimws(soort_naam)),
           tolower(trimws(Type))  == tolower(trimws(type_naam)))
  
  if (nrow(codes_subset) == 0) return(NULL)
  
  exact_codes <- tolower(trimws(codes_subset$Code[codes_subset$Match == "exact"]))
  bevat_codes <- tolower(trimws(codes_subset$Code[codes_subset$Match == "bevat"]))
  
  if (length(bevat_codes) > 0) {
    bevat_escaped <- gsub("([\\.\\^\\$\\*\\+\\?\\(\\)\\[\\{\\\\\\|])", "\\\\\\1", bevat_codes)
    regex_term    <- paste0("^", bevat_escaped, collapse = "|")
    
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
  
  if (nrow(tabel_sub) == 0) return(NULL)
  
  tabel_sub[, CODE := tolower(trimws(CODE))]
  tabel_unique  <- unique(tabel_sub, by = c("cel_id", "CODE"))
  tabel_cel_som <- tabel_unique[, .(Oppervlakte = pmin(sum(BWK_FRAC, na.rm = TRUE), 1.0)), by = .(cel_id)]
  
  r_opp <- master_grid
  r_opp[tabel_cel_som$cel_id] <- tabel_cel_som$Oppervlakte
  return(r_opp)
}
