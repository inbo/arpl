# ==============================================================================
# vul_raster_vanuit_parquet.R
# Universele en veilige Parquet-naar-Raster toewijzing voor Vlaanderen én subsets
# ==============================================================================

library(terra)
library(data.table)
library(dplyr)
library(stringr)

vul_raster_vanuit_parquet <- function(ds_parquet, 
                                      codes_exact = character(0), 
                                      codes_bevat = character(0), 
                                      master_grid, 
                                      pad_master_referentie = NULL,
                                      veld_waarde = "BWK_FRAC",
                                      aggregeer_som = TRUE,
                                      waarde_vast = NULL) {
  
  if (length(codes_exact) == 0 && length(codes_bevat) == 0) {
    return(master_grid * NA)
  }
  
  codes_exact_clean <- tolower(trimws(codes_exact))
  codes_bevat_clean <- tolower(trimws(codes_bevat))
  
  # 1. Bepaal of we op een subset / testgebied werken
  is_subset <- terra::ncell(master_grid) < 200000000 # Vlaanderen telt ~300M cellen
  
  dt_map <- NULL
  if (is_subset) {
    # Referentiegrid van Vlaanderen ophalen (gecached in global env voor snelheid)
    if (!exists(".vlaanderen_ref_grid", envir = .GlobalEnv) || is.null(.GlobalEnv$.vlaanderen_ref_grid)) {
      ref_pad <- if (!is.null(pad_master_referentie)) {
        pad_master_referentie
      } else if (exists("PATHS") && !is.null(PATHS$master_grid)) {
        PATHS$master_grid
      } else {
        here::here("data", "input", "Raster_Vlaanderen", "Vlaanderen_MasterGrid_10m.tif")
      }
      .GlobalEnv$.vlaanderen_ref_grid <- terra::rast(ref_pad)[[1]]
    }
    
    # Bereken de globale cel-indices binnen de uitsnede
    coords_sub <- terra::xyFromCell(master_grid, seq_len(terra::ncell(master_grid)))
    globale_ids_sub <- terra::cellFromXY(.GlobalEnv$.vlaanderen_ref_grid, coords_sub)
    
    dt_map <- data.table::data.table(
      lokale_id = seq_len(terra::ncell(master_grid)),
      cel_id    = globale_ids_sub
    )[!is.na(cel_id)]
  }
  
  # 2. Filter opstellen voor Parquet via Arrow
  if (length(codes_bevat_clean) > 0) {
    bevat_escaped <- gsub("([\\.\\^\\$\\*\\+\\?\\(\\)\\[\\{\\\\\\|])", "\\\\\\1", codes_bevat_clean)
    regex_term    <- paste0("^", bevat_escaped, collapse = "|")
    
    query <- ds_parquet %>%
      dplyr::filter(tolower(CODE) %in% codes_exact_clean | stringr::str_detect(tolower(CODE), regex_term))
  } else {
    query <- ds_parquet %>%
      dplyr::filter(tolower(CODE) %in% codes_exact_clean)
  }
  
  # Als we in subsetmodus zitten: filter direct in Parquet (enorme snelheidswinst!)
  if (is_subset) {
    query <- query %>% dplyr::filter(cel_id %in% dt_map$cel_id)
  }
  
  tabel_sub <- query %>%
    dplyr::collect() %>%
    data.table::as.data.table()
  
  r_opp <- master_grid * NA
  if (nrow(tabel_sub) == 0) return(r_opp)
  
  # 3. Aggregatie berekenen
  tabel_sub[, CODE := tolower(trimws(CODE))]
  tabel_uniek <- unique(tabel_sub, by = c("cel_id", "CODE"))
  
  if (!is.null(waarde_vast)) {
    dt_agg <- unique(tabel_uniek, by = "cel_id")[, .(Waarde = as.numeric(waarde_vast)), by = .(cel_id)]
  } else if (aggregeer_som) {
    dt_agg <- tabel_uniek[, .(Waarde = pmin(sum(get(veld_waarde), na.rm = TRUE), 1.0)), by = .(cel_id)]
  } else {
    dt_agg <- tabel_uniek[, .(Waarde = max(get(veld_waarde), na.rm = TRUE)), by = .(cel_id)]
  }
  
  # 4. Veilig toewijzen aan het raster
  if (is_subset) {
    dt_gemapt <- merge(dt_agg, dt_map, by = "cel_id", all.x = TRUE)[!is.na(lokale_id)]
    if (nrow(dt_gemapt) > 0) {
      r_opp[dt_gemapt$lokale_id] <- dt_gemapt$Waarde
    }
  } else {
    r_opp[dt_agg$cel_id] <- dt_agg$Waarde
  }
  
  return(r_opp)
}
