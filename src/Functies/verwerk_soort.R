# ==============================================================================
# verwerk_soort.R
# ==============================================================================

verwerk_soort <- function(soort_naam, 
                          df_excel_afstanden = NULL, 
                          params_soort       = NULL,
                          ds_parquet, 
                          df_bwk_codes, 
                          master_grid, 
                          paths, 
                          force_overwrite    = FALSE,
                          ...) {
  
  message(sprintf("\n>>> Start verwerking voor: %s", toupper(soort_naam)))
  
  # Vang op of params_soort direct binnenkomt, of gefilterd moet worden uit de excel-tabel
  if (is.null(params_soort)) {
    if (is.null(df_excel_afstanden)) {
      stop("Geen parameters gevonden: geef ofwel 'params_soort' ofwel 'df_excel_afstanden' mee.")
    }
    params_soort <- df_excel_afstanden %>%
      dplyr::filter(tolower(trimws(Soort)) == tolower(trimws(soort_naam)))
  }
  
  if (nrow(params_soort) == 0) {
    stop(sprintf("❌ Soort '%s' niet gevonden in Excel afstanden-tabel!", soort_naam))
  }
  
  source(here::here("src/Functies/bouw_biotoop_universeel.R"))
  
  res <- bouw_biotoop_universeel(
    soort_naam      = soort_naam,
    params_soort    = params_soort,
    ds_parquet      = ds_parquet,
    df_bwk_codes    = df_bwk_codes,
    master_grid     = master_grid,
    paths           = paths,
    force_overwrite = force_overwrite
  )
  
  # Exporteer resultaten (indien aanwezig in je pipeline)
  if (exists("exporteer_soort_resultaten", mode = "function")) {
    exporteer_soort_resultaten(soort_naam, res, master_grid, paths)
  }
  
  return(invisible(res))
}
