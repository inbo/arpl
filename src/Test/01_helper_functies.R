library(terra)
library(data.table)
library(dplyr)

# ------------------------------------------------------------------------------
# 1. RAZENDSNELLE CLUSTER FILTER VOOR HEEL VLAANDEREN
# ------------------------------------------------------------------------------
cluster_filter_vlaanderen <- function(masker, opp_laag, drempel_m2, dist_m, werkelijk = FALSE) {
  if (all(is.na(suppressWarnings(terra::minmax(masker))))) {
    return(list(raster = masker * NA, clusters = masker * NA))
  }
  
  # Netwerk opbouwen op basis van dispersiecapaciteit
  if (dist_m > 0) {
    r_binair <- terra::ifel(!is.na(masker) & masker > 0, 1, NA)
    r_buffered <- terra::buffer(r_binair, width = dist_m / 2)
    cl_network <- terra::patches(r_buffered, directions = 4, zeroAsNA = TRUE)
    cl_biotoop <- terra::mask(cl_network, masker)
  } else {
    cl_biotoop <- terra::patches(masker, directions = 8, zeroAsNA = TRUE)
  }
  
  # Zonale statistieken berekenen (m2 bepalen)
  if (werkelijk) {
    stats_df <- terra::zonal(opp_laag, cl_biotoop, fun = "sum", na.rm = TRUE)
    colnames(stats_df) <- c("ID", "Waarde")
    stats_df$Area_m2 <- stats_df$Waarde * 100 
  } else {
    f <- terra::freq(cl_biotoop)
    stats_df <- data.frame(ID = f$value, Waarde = f$count)
    stats_df$Area_m2 <- stats_df$Waarde * 100 
  }
  
  stats_df <- stats_df[!is.na(stats_df$ID), ]
  if (nrow(stats_df) == 0) return(list(raster = masker * NA, clusters = masker * NA))
  
  voldoet_ids <- stats_df$ID[stats_df$Area_m2 >= drempel_m2]
  if (length(voldoet_ids) == 0) return(list(raster = masker * NA, clusters = masker * NA))
  
  final_network_mask <- terra::ifel(cl_biotoop %in% voldoet_ids, 1, NA)
  
  r_finaal  <- terra::mask(masker, final_network_mask)
  cl_finaal <- terra::mask(cl_biotoop, r_finaal) 
  
  return(list(raster = r_finaal, clusters = cl_finaal))
}

# ------------------------------------------------------------------------------
# 2. UNIVERSELE EXPORT VOOR VLAANDEREN
# ------------------------------------------------------------------------------
exporteer_resultaten_vlaanderen <- function(final_max, final_opp, cl_opp, soort, scenario_naam, master_grid) {
  base_dir <- here::here("data/output/Vlaanderen/Rasters_Soorten", scenario_naam)
  
  folders <- list(
    potentie  = file.path(base_dir, "01_Maximale_Potentie"),
    werkelijk = file.path(base_dir, "02_Werkelijke_Oppervlaktes"),
    id_raster = file.path(base_dir, "00_ID_Rasters")
  )
  purrr::walk(folders, ~if (!dir.exists(.x)) dir.create(.x, showWarnings = FALSE, recursive = TRUE))
  
  potentie_export_rast <- if (!is.null(final_max) && !all(is.na(suppressWarnings(terra::minmax(final_max))))) {
    terra::ifel(!is.na(final_max) & final_max > 0, 1, NA)
  } else {
    master_grid * NA
  }
  
  werkelijk_export_rast <- if (!is.null(final_opp) && !all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
    terra::ifel(!is.na(final_opp) & final_opp > 0, 1, NA)
  } else {
    master_grid * NA
  }
  
  if (!is.null(cl_opp) && !all(is.na(suppressWarnings(terra::minmax(cl_opp))))) {
    id_export_rast <- cl_opp
  } else if (!is.null(final_opp) && !all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
    id_export_rast <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE)
  } else {
    id_export_rast <- master_grid * NA
  }
  
  bestands_naam <- tolower(gsub(" ", "", soort))
  
  export_config <- list(
    list(rast = potentie_export_rast,  naam = "Maximale_Potentie",        folder = folders$potentie),
    list(rast = werkelijk_export_rast, naam = "Werkelijke_Oppervlaktes", folder = folders$werkelijk)
  )
  
  for(item in export_config) {
    file_path <- file.path(item$folder, paste0("Habitat_", item$naam, "_", bestands_naam, ".tif"))
    terra::writeRaster(item$rast, filename = file_path, overwrite = TRUE, gdal = c("COMPRESS=LZW"), datatype = "INT1U", NAflag = 255)
    message(paste("    [OK] Geëxporteerd:", basename(file_path)))
  }
  
  file_path_id <- file.path(folders$id_raster, paste0("ID_Netwerken_", bestands_naam, ".tif"))
  terra::writeRaster(id_export_rast, filename = file_path_id, overwrite = TRUE, gdal = c("COMPRESS=LZW"), datatype = "INT4U", NAflag = 0)
  message(paste("    [OK] Analytisch ID-raster geëxporteerd:", basename(file_path_id)))
}

# ------------------------------------------------------------------------------
# CLUSTER FILTER MET TWEE DREMPELS (VLAANDEREN LEVEL)
# ------------------------------------------------------------------------------
cluster_filter_twee_drempels_vlaanderen <- function(masker, opp_laag, min_opp_ha, min_totale_opp_ha, dist_m, werkelijk = FALSE) {
  if (all(is.na(suppressWarnings(terra::minmax(masker))))) {
    return(list(raster = masker * NA, clusters = masker * NA))
  }
  
  min_opp_m2    <- min_opp_ha * 10000
  min_totale_m2 <- min_totale_opp_ha * 10000
  
  # 1. Individuele clusters identificeren en filteren op min_opp_ha
  cl_individueel <- terra::patches(masker, directions = 8, zeroAsNA = TRUE)
  
  if (werkelijk) {
    stats_indiv <- terra::zonal(opp_laag, cl_individueel, fun = "sum", na.rm = TRUE)
    colnames(stats_indiv) <- c("ID", "Waarde")
    stats_indiv$Area_m2 <- stats_indiv$Waarde * 100
  } else {
    f_indiv <- terra::freq(cl_individueel)
    stats_indiv <- data.frame(ID = f_indiv$value, Waarde = f_indiv$count)
    stats_indiv$Area_m2 <- stats_indiv$Waarde * 100
  }
  
  voldoet_indiv_ids <- stats_indiv$ID[!is.na(stats_indiv$ID) & stats_indiv$Area_m2 >= min_opp_m2]
  if (length(voldoet_indiv_ids) == 0) return(list(raster = masker * NA, clusters = masker * NA))
  
  masker_gefilterd <- cl_individueel %in% voldoet_indiv_ids
  masker_gefilterd <- terra::ifel(masker_gefilterd == 1, masker, NA)
  
  # 2. Netwerkvorming via dispersiebuffer
  if (dist_m > 0) {
    r_binair <- terra::ifel(!is.na(masker_gefilterd) & masker_gefilterd > 0, 1, NA)
    r_buffered <- terra::buffer(r_binair, width = dist_m / 2)
    cl_network <- terra::patches(r_buffered, directions = 4, zeroAsNA = TRUE)
    cl_biotoop_only <- terra::mask(cl_network, masker_gefilterd)
  } else {
    cl_network <- terra::patches(masker_gefilterd, directions = 8, zeroAsNA = TRUE)
    cl_biotoop_only <- cl_network
  }
  
  # 3. Oppervlakte per netwerk filteren op min_totale_opp_ha
  if (werkelijk) {
    stats_net <- terra::zonal(opp_laag, cl_biotoop_only, fun = "sum", na.rm = TRUE)
    colnames(stats_net) <- c("ID", "Waarde")
    stats_net$Area_m2 <- stats_net$Waarde * 100
  } else {
    f_net <- terra::freq(cl_biotoop_only)
    stats_net <- data.frame(ID = f_net$value, Waarde = f_net$count)
    stats_net$Area_m2 <- stats_net$Waarde * 100
  }
  
  stats_net <- stats_net[!is.na(stats_net$ID), ]
  if (nrow(stats_net) == 0) return(list(raster = masker * NA, clusters = masker * NA))
  
  voldoet_net_ids <- stats_net$ID[stats_net$Area_m2 >= min_totale_m2]
  if (length(voldoet_net_ids) == 0) return(list(raster = masker * NA, clusters = masker * NA))
  
  final_network_mask <- cl_biotoop_only %in% voldoet_net_ids
  final_network_mask <- terra::ifel(final_network_mask == 1, 1, NA)
  
  r_finaal  <- terra::mask(masker_gefilterd, final_network_mask)
  cl_finaal <- terra::mask(cl_biotoop_only, r_finaal)
  
  return(list(raster = r_finaal, clusters = cl_finaal))
}
