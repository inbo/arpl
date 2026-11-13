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
  
  bestandsnaam <- tolower(gsub(" ", "", soort_naam))
  message(sprintf("\n>>> Start verwerking voor: %s", toupper(soort_naam)))
  
  # 1. Parameters ophalen
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
  
  # 2. Zorg dat paths$dir_output gegarandeerd ingevuld is
  dir_uitvoer <- dplyr::coalesce(
    paths$dir_output,
    paths$output_dir,
    paths$dir_out,
    paths$output,
    paths$base_dir,
    paths$rasters_soorten
  )
  
  if (is.null(dir_uitvoer) || is.na(dir_uitvoer)) {
    dir_uitvoer <- here::here("data", "output", "Vlaanderen", "Rasters_Soorten")
  }
  paths$dir_output <- dir_uitvoer
  
  # 3. Model draaien
  res <- bouw_biotoop_universeel(
    soort_naam      = soort_naam,
    params_soort    = params_soort,
    ds_parquet      = ds_parquet,
    df_bwk_codes    = df_bwk_codes,
    master_grid     = master_grid,
    paths           = paths,
    force_overwrite = force_overwrite
  )
  
  # 4. Directe exportverificatie (vangnet indien bouw_biotoop niet schreef)
  scenario_map <- if (!is.null(paths$scenario_naam) && paths$scenario_naam != "") paths$scenario_naam else NULL
  base_dir     <- if (!is.null(scenario_map)) file.path(paths$dir_output, scenario_map) else paths$dir_output
  
  dir_id   <- file.path(base_dir, "00_ID_Rasters")
  dir_leef <- file.path(base_dir, "01_Potentieel_Leefgebied")
  dir_deel <- file.path(dir_leef, "Deellagen") # 📂 Map binnen 01_Potentieel_Leefgebied
  
  f_out_id  <- file.path(dir_id, paste0("ID_Netwerken_", bestandsnaam, ".tif"))
  f_out_tot <- file.path(dir_leef, paste0("Potentieel_Leefgebied_", bestandsnaam, ".tif"))
  
  # Alleen wegschrijven als bouw_biotoop dat nog niet gedaan heeft
  if (!file.exists(f_out_tot) && !is.null(res$final_opp) && !all(is.na(suppressWarnings(terra::minmax(res$final_opp))))) {
    dir.create(dir_id, showWarnings = FALSE, recursive = TRUE)
    dir.create(dir_leef, showWarnings = FALSE, recursive = TRUE)
    dir.create(dir_deel, showWarnings = FALSE, recursive = TRUE)
    
    # Export Netwerk ID
    if (!is.null(res$cl_opp) && !all(is.na(suppressWarnings(terra::minmax(res$cl_opp))))) {
      terra::writeRaster(
        res$cl_opp,
        filename  = f_out_id,
        overwrite = TRUE,
        gdal      = c("COMPRESS=LZW"),
        datatype  = "INT4U",
        NAflag    = 0
      )
      message(sprintf("    💾 [EXPORT] Netwerk-ID raster bewaard: %s", basename(f_out_id)))
    }
    
    # Export Totaal Leefgebied
    terra::writeRaster(
      res$final_opp,
      filename  = f_out_tot,
      overwrite = TRUE,
      gdal      = c("COMPRESS=LZW"),
      datatype  = "FLT4S",
      NAflag    = -9999
    )
    message(sprintf("    💾 [EXPORT] Totaal leefgebied bewaard: %s", basename(f_out_tot)))
    
    # Export Deellagen naar 01_Potentieel_Leefgebied/Deellagen/
    deellagen_namen <- setdiff(names(res), c("final_opp", "cl_opp", "gradatie"))
    if (length(deellagen_namen) > 0) {
      for (nm in deellagen_namen) {
        r_deel <- res[[nm]]
        if (!is.null(r_deel) && !all(is.na(suppressWarnings(terra::minmax(r_deel))))) {
          f_out_deel <- file.path(dir_deel, paste0(bestandsnaam, "_", nm, ".tif"))
          terra::writeRaster(
            r_deel,
            filename  = f_out_deel,
            overwrite = TRUE,
            gdal      = c("COMPRESS=LZW"),
            datatype  = "FLT4S",
            NAflag    = -9999
          )
          message(sprintf("    💾 [EXPORT] Deellaag bewaard: %s", basename(f_out_deel)))
        }
      }
    }
  }
  
  return(invisible(res))
}
