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
  
  # 2. Uitvoermap bepalen (NULL-veilig voor dplyr::coalesce)
  kandidaat_paden <- c(
    paths$dir_output,
    paths$output_dir,
    paths$dir_out,
    paths$output,
    paths$base_dir,
    paths$rasters_soorten
  )
  kandidaat_paden <- kandidaat_paden[!is.na(kandidaat_paden) & kandidaat_paden != ""]
  
  dir_uitvoer <- if (length(kandidaat_paden) > 0) {
    kandidaat_paden[1]
  } else {
    here::here("data", "output", "Vlaanderen", "Rasters_Soorten")
  }
  paths$dir_output <- dir_uitvoer
  
  # 3. Model draaien (enkel rekenwerk, geen export)
  res <- bouw_biotoop_universeel(
    soort_naam      = soort_naam,
    params_soort    = params_soort,
    ds_parquet      = ds_parquet,
    df_bwk_codes    = df_bwk_codes,
    master_grid     = master_grid,
    paths           = paths,
    force_overwrite = force_overwrite
  )
  
  # 4. Centrale export naar disk (1 / NA binaire kaarten)
  scenario_map <- if (!is.null(paths$scenario_naam) && paths$scenario_naam != "") paths$scenario_naam else NULL
  base_dir     <- if (!is.null(scenario_map)) file.path(paths$dir_output, scenario_map) else paths$dir_output
  
  dir_id   <- file.path(base_dir, "00_ID_Rasters")
  dir_leef <- file.path(base_dir, "01_Potentieel_Leefgebied")
  dir_deel <- file.path(dir_leef, "Deellagen")
  
  dir.create(dir_id,   showWarnings = FALSE, recursive = TRUE)
  dir.create(dir_leef, showWarnings = FALSE, recursive = TRUE)
  dir.create(dir_deel, showWarnings = FALSE, recursive = TRUE)
  
  # 4A. Netwerk-ID raster
  if (!is_leeg(res$cl_opp)) {
    f_out_id <- file.path(dir_id, paste0("ID_Netwerken_", bestandsnaam, ".tif"))
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
  
  # 4B. Totaal Potentieel Leefgebied (1 / NA)
  if (!is_leeg(res$final_opp)) {
    f_out_tot <- file.path(dir_leef, paste0("Potentieel_Leefgebied_", bestandsnaam, ".tif"))
    r_tot_bin <- terra::ifel(!is.na(res$final_opp) & res$final_opp > 0, 1L, NA_integer_)
    terra::writeRaster(
      r_tot_bin,
      filename  = f_out_tot,
      overwrite = TRUE,
      gdal      = c("COMPRESS=LZW"),
      datatype  = "INT1U",
      NAflag    = 255
    )
    rm(r_tot_bin)
    message(sprintf("    💾 [EXPORT] Totaal leefgebied bewaard (1/NA): %s", basename(f_out_tot)))
  }
  
  # 4C. Deellagen (enkel indien > 1 deellaag, elk als 1 / NA)
  deellagen_namen <- setdiff(names(res), c("final_opp", "cl_opp", "gradatie"))
  if (length(deellagen_namen) > 1) {
    for (nm in deellagen_namen) {
      r_deel <- res[[nm]]
      if (!is_leeg(r_deel)) {
        f_out_deel <- file.path(dir_deel, paste0(bestandsnaam, "_", gsub(":", "_", nm), ".tif"))
        r_deel_bin <- terra::ifel(!is.na(r_deel) & r_deel > 0, 1L, NA_integer_)
        terra::writeRaster(
          r_deel_bin,
          filename  = f_out_deel,
          overwrite = TRUE,
          gdal      = c("COMPRESS=LZW"),
          datatype  = "INT1U",
          NAflag    = 255
        )
        rm(r_deel_bin)
        message(sprintf("    💾 [EXPORT] Deellaag bewaard (1/NA): %s", basename(f_out_deel)))
      }
    }
  }
  
  return(invisible(res))
}
