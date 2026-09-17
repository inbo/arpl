library(knitr)
library(here)
library(tidyverse)
library(sf)
library(terra)
library(readxl)
library(tidyterra)
library(leaflet)
library(kableExtra)
library(data.table)

conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::select)
conflicted::conflicts_prefer(terra::intersect)

# --- FUNCTIES ---

calc_ha_exact <- function(r) {
  if(is.null(r)) return(0)
  if(all(is.na(terra::values(r, mat=FALSE)))) return(0)
  area_raster <- r * terra::cellSize(r, unit = "ha")
  val <- terra::global(area_raster, "sum", na.rm = TRUE)[[1]]
  return(as.numeric(val))
}

cluster_filter_compleet <- function(masker, opp_laag, drempel_m2, dist_m, werkelijk = FALSE) {
  if (terra::global(is.na(masker), "sum")[[1]] == terra::ncell(masker)) {
    return(list(raster = masker * NA, clusters = masker * NA))
  }
  
  if (dist_m > 0) {
    r_binair <- terra::ifel(!is.na(masker) & masker > 0, 1, NA)
    r_buffered <- terra::buffer(r_binair, width = dist_m / 2)
    cl_network <- terra::patches(r_buffered, directions = 4, zeroAsNA = TRUE)
    cl_biotoop_only <- terra::mask(cl_network, masker)
  } else {
    cl_network <- terra::patches(masker, directions = 8, zeroAsNA = TRUE)
    cl_biotoop_only <- cl_network
  }
  
  if(werkelijk) {
    stats_df <- terra::zonal(opp_laag, cl_biotoop_only, fun = "sum", na.rm = TRUE)
    colnames(stats_df) <- c("ID", "Waarde")
    stats_df$Area_m2 <- stats_df$Waarde * 100 
  } else {
    f <- terra::freq(cl_biotoop_only)
    stats_df <- data.frame(ID = f$value, Waarde = f$count)
    stats_df$Area_m2 <- stats_df$Waarde * 100 
  }
  
  stats_df <- stats_df[!is.na(stats_df$ID), ]
  if(nrow(stats_df) == 0) return(list(raster = masker * NA, clusters = masker * NA))
  
  voldoet_ids <- stats_df$ID[stats_df$Area_m2 >= drempel_m2]
  if(length(voldoet_ids) == 0) return(list(raster = masker * NA, clusters = masker * NA))
  
  masker_binair <- cl_biotoop_only %in% voldoet_ids
  final_network_mask <- terra::ifel(masker_binair == 1, 1, NA)
  
  r_finaal  <- terra::mask(masker, final_network_mask)
  cl_finaal <- terra::mask(cl_biotoop_only, r_finaal) 
  
  return(list(raster = r_finaal, clusters = cl_finaal))
}

# Handmatige vastlegging voor Levendbarende hagedis
soort          <- "levendbarendehagedis"

# Scenario pad en naam bepalen

# --- DYNAMISCHE SCENARIO PARAMETER CHECK ---
if (!exists("params") || is.null(params$scenario_rds_path)) {
  scenario_rds_path <- "data/input/Scenario_rds/MH_Scenario_BWK_2025.rds"
} else {
  scenario_rds_path <- params$scenario_rds_path
}

p_raw <- gsub("^([.][.]/)+", "", scenario_rds_path)
scenario_path <- here::here(p_raw)


if (!file.exists(scenario_path)) {
  stop(paste("❌ FOUT: Scenario RDS bestand NIET gevonden op:", scenario_path))
}

scen_volledig <- basename(scenario_path)
scenario_naam <- gsub("^MH_Scenario_|^Scenario_|.rds$", "", scen_volledig)

message(paste("Verwerken van soort:", soort, "binnen scenario:", scenario_naam))

oppervlakte_ha <- 5    # Minimale vlek-oppervlakte (MinOpp)
min_totale_ha  <- 10   # Minimale netwerk-oppervlakte (Min totale Opp)
afstand_m      <- 0    # Afstand tussen biotopen
buffer_m       <- 100  # Dispersiecapaciteit

cat("Soort:", soort, "\n")
cat("MinOpp vlek:", oppervlakte_ha, "ha | MinOpp netwerk:", min_totale_ha, "ha\n")
cat("Dispersieafstand:", buffer_m, "m\n")

area_shape   <- vect(here("data/input/Mechelse_Heide.shp"))
master_grid  <- rast(here("data/input/Raster_Vlaanderen/Vlaanderen_MasterGrid_10m.tif"))[[1]]

area_shape_proj <- project(area_shape, crs(master_grid))
area_buffer_fix <- buffer(area_shape_proj, width = buffer_m)

message("-> Vertaalraster voor globale/lokale cellen opbouwen...")

id_raster_MH <- crop(master_grid, area_buffer_fix, snap = "near")

globale_id_raster <- master_grid
values(globale_id_raster) <- 1:ncell(globale_id_raster)
id_raster_MH_globale_values <- crop(globale_id_raster, area_buffer_fix, snap = "near")

extractie_ids <- terra::extract(id_raster_MH_globale_values, area_buffer_fix, cells = TRUE)

vertaal_df <- as.data.table(extractie_ids)
setnames(vertaal_df, c("cell", names(id_raster_MH_globale_values)), c("lokale_id", "globale_id"))
vertaal_df <- vertaal_df[!is.na(globale_id)]

studiegebied_globale_ids <- unique(vertaal_df$globale_id)

values(id_raster_MH) <- NA
# Reparatie in de Gebied-chunk:
template_MH <- terra::rasterize(area_buffer_fix, id_raster_MH, field = 1, background = 0)

rm(globale_id_raster, id_raster_MH_globale_values)
gc()

# Biotoopfiltering (BWK-codes)
# 1. LAAD DE BRON-CROSSWALK IN
df_nieuw <- read_csv(here("data/input/Excel_files/Resultaten_Totaal_Samengevoegd.csv"), show_col_types = FALSE)

resultaten_gegroepeerd <- df_nieuw %>%
  mutate(Soort_clean = tolower(trimws(Soort))) %>%
  filter(Soort_clean == soort) %>%
  group_by(Type) %>%
  nest(Data = c(Code, Match))

# Lees de scenario-RDS in
tabel_vlaanderen <- readRDS(scenario_path)
setDT(tabel_vlaanderen)
tabel_vlaanderen[, CODE := tolower(trimws(CODE))]

lijst_matches      <- list()
lijst_oppervlaktes <- list()

for(i in 1:nrow(resultaten_gegroepeerd)) {
  h_type <- resultaten_gegroepeerd$Type[i]
  h_data <- resultaten_gegroepeerd$Data[[i]]
  
  exact_codes <- tolower(trimws(h_data$Code[h_data$Match == "exact"]))
  bevat_codes <- tolower(trimws(h_data$Code[h_data$Match == "bevat"]))
  
  if(length(bevat_codes) > 0) {
    bevat_codes_escaped <- gsub("([\\.\\^\\$\\*\\+\\?\\(\\)\\[\\{\\\\\\|])", "\\\\\\1", bevat_codes)
    bevat_codes_anchored <- paste0("^", bevat_codes_escaped)
    regex_term <- paste0(bevat_codes_anchored, collapse = "|")
    
    tabel_gefilterd <- tabel_vlaanderen[CODE %in% exact_codes | grepl(regex_term, CODE)]
  } else {
    tabel_gefilterd <- tabel_vlaanderen[CODE %in% exact_codes]
  }
  
  tabel_MH <- tabel_gefilterd[cel_id %in% studiegebied_globale_ids]
  tabel_MH_unique <- unique(tabel_MH, by = c("cel_id", "CODE"))

  tabel_cel_som <- tabel_MH_unique[, .(Oppervlakte = pmin(sum(BWK_FRAC, na.rm = TRUE), 1.0)), by = .(cel_id)]
  
  r_match_type <- id_raster_MH * NA
  r_opp_type   <- id_raster_MH * NA
  
  if(nrow(tabel_cel_som) > 0) {
    tabel_cel_som[, Oppervlakte := pmin(Oppervlakte, 1)]
    tabel_cel_som[, Match := ifelse(Oppervlakte >= 0.01, 1, 0)]
    
    setnames(tabel_cel_som, "cel_id", "globale_id")
    tabel_finaal_mapping <- merge(tabel_cel_som, vertaal_df, by = "globale_id", all.x = TRUE)
    tabel_finaal_mapping <- tabel_finaal_mapping[!is.na(lokale_id)]
    
    if(nrow(tabel_finaal_mapping) > 0) {
      tabel_matches_clean <- tabel_finaal_mapping[Match == 1]
      if(nrow(tabel_matches_clean) > 0) {
        r_match_type[tabel_matches_clean$lokale_id] <- tabel_matches_clean$Match
      }
      r_opp_type[tabel_finaal_mapping$lokale_id] <- tabel_finaal_mapping$Oppervlakte
      
      lijst_matches[[h_type]]      <- r_match_type
      lijst_oppervlaktes[[h_type]] <- r_opp_type
    }
  }
}

lijst_matches      <- purrr::compact(lijst_matches)
lijst_oppervlaktes <- purrr::compact(lijst_oppervlaktes)

if (length(lijst_matches) == 0 || length(lijst_oppervlaktes) == 0) {
  message("⚠️ KRIJT: Geen enkele doel-BWK code voor '", soort, "' aangetroffen in de database. Failsafe ingeschakeld.")
  masker_bron <- template_MH * NA
  opp_bron    <- template_MH * NA
} else {
  masker_bron <- lijst_matches[[1]]
  opp_bron    <- lijst_oppervlaktes[[1]]
}

names(masker_bron) <- "Match"
names(opp_bron)    <- "Oppervlakte"

suppressWarnings(rm(tabel_vlaanderen, vertaal_df, lijst_matches, lijst_oppervlaktes))
gc()

drempel_vlek_m2    <- oppervlakte_ha * 10000  # 5 ha = 50.000 m²
drempel_netwerk_m2 <- min_totale_ha * 10000   # 10 ha = 100.000 m²

if(!all(is.na(terra::values(masker_bron, mat = FALSE)))) {
  masker_bron[masker_bron == 0] <- NA
}

# ------------------------------------------------------------------------------
# STAP 1: VOORFILTER - Verwijder losse biotopen < 5 ha (MinOpp vlek)
# ------------------------------------------------------------------------------
if (!all(is.na(terra::values(masker_bron, mat = FALSE)))) {
  vlek_patches <- terra::patches(masker_bron, directions = 8, zeroAsNA = TRUE)
  vlek_freq    <- as.data.frame(terra::freq(vlek_patches))
  
  gevalideerde_vlek_ids <- vlek_freq$value[(vlek_freq$count * 100) >= drempel_vlek_m2]
  
  if (length(gevalideerde_vlek_ids) > 0) {
    masker_gefilterd <- terra::mask(masker_bron, vlek_patches %in% gevalideerde_vlek_ids)
    opp_gefilterd    <- terra::mask(opp_bron, vlek_patches %in% gevalideerde_vlek_ids)
  } else {
    masker_gefilterd <- masker_bron * NA
    opp_gefilterd    <- opp_bron * NA
  }
} else {
  masker_gefilterd <- masker_bron * NA
  opp_gefilterd    <- opp_bron * NA
}

# ------------------------------------------------------------------------------
# STAP 2: NETWERKCLUSTERING - Verbind vlekken binnen 100 m & check netwerk >= 10 ha
# ------------------------------------------------------------------------------
res_max <- cluster_filter_compleet(
  masker     = masker_gefilterd, 
  opp_laag   = opp_gefilterd, 
  drempel_m2 = drempel_netwerk_m2, 
  dist_m     = buffer_m, 
  werkelijk  = FALSE
)

final_max <- res_max$raster
cl_max    <- res_max$clusters

# Evalueer werkelijke oppervlakte binnen de netwerken
if (!all(is.na(terra::values(cl_max, mat = FALSE)))) {
  zonal_werkelijk <- terra::zonal(opp_gefilterd, cl_max, fun = "sum", na.rm = TRUE)
  colnames(zonal_werkelijk) <- c("ID", "Werk_m2")
  zonal_werkelijk$Werk_m2 <- zonal_werkelijk$Werk_m2 * 100
  
  ids_werkelijk_voldoet <- zonal_werkelijk$ID[zonal_werkelijk$Werk_m2 >= drempel_netwerk_m2]
  
  if(length(ids_werkelijk_voldoet) > 0) {
    final_opp <- terra::mask(final_max, cl_max %in% ids_werkelijk_voldoet)
    cl_opp    <- terra::mask(cl_max, cl_max %in% ids_werkelijk_voldoet)
  } else {
    final_opp <- final_max * NA
    cl_opp    <- cl_max * NA
  }
} else {
  final_opp <- final_max * NA
  cl_opp    <- cl_max * NA
}

grens_merged  <- terra::aggregate(area_shape_proj)
buffer_lijn   <- terra::buffer(grens_merged, width = buffer_m)
binnen_masker <- terra::rasterize(grens_merged, template_MH, field = 1)

# Statusraster opbouwen
if (!all(is.na(terra::values(final_max, mat = FALSE)))) {
  is_loss   <- !is.na(final_max) & is.na(final_opp)
  is_kept   <- !is.na(final_opp)
  is_inside <- !is.na(binnen_masker)

  r_status <- terra::ifel(is_loss & is_inside, 1, NA)
  r_status <- terra::cover(r_status, terra::ifel(is_loss & !is_inside, 2, NA))
  r_status <- terra::cover(r_status, terra::ifel(is_kept & is_inside, 3, NA))
  r_status <- terra::cover(r_status, terra::ifel(is_kept & !is_inside, 4, NA))
  
  r_status <- terra::as.factor(r_status)
  status_labels <- data.frame(ID = c(1, 2, 3, 4), 
                              Label = c("Maximale Potentie (Binnen)", "Maximale Potentie (Buiten)", 
                                        "Werkelijk Habitat (Binnen)", "Werkelijk Habitat (Buiten)"))
  levels(r_status) <- status_labels
} else {
  r_status <- terra::rast(template_MH, vals = NA)
  message("Let op: Geen geschikte clusters gevonden voor deze soort.")
}


# ==============================================================================
# SCHONE EXPORT BIOTOOP EN ANALYTISCH ID-RASTER (VOOR SCRIPT 2 / ARPL)
# ==============================================================================
base_dir <- here::here("data/output/Mechelse_Heide/Rasters_Soorten", scenario_naam)

folders <- list(
  potentie  = file.path(base_dir, "01_Maximale_Potentie"),
  werkelijk = file.path(base_dir, "02_Werkelijke_Oppervlaktes"),
  id_raster = file.path(base_dir, "00_ID_Rasters")
)
purrr::walk(folders, ~if (!dir.exists(.x)) dir.create(.x, showWarnings = FALSE, recursive = TRUE))

# 1. Bepaal Maximale Potentie Raster
potentie_export_rast <- if (exists("grutto_kaart_A") && !is.null(grutto_kaart_A) && !all(is.na(terra::values(grutto_kaart_A, mat=FALSE)))) {
  terra::ifel(grutto_kaart_A > 0, 1, NA)
} else if (exists("final_max") && !all(is.na(terra::values(final_max, mat=FALSE)))) {
  terra::ifel(!is.na(final_max) & final_max > 0, 1, NA)
} else {
  terra::rast(template_MH, vals = NA)
}

# 2. Bepaal Werkelijke Oppervlakte Raster
werkelijk_export_rast <- if (exists("resB_strikt") && !is.null(resB_strikt) && (!all(is.na(terra::values(resB_strikt$kern, mat=FALSE))) || !all(is.na(terra::values(resB_strikt$bouw, mat=FALSE))))) {
  r_net_totaal_opp <- terra::cover(resB_strikt$kern, resB_strikt$bouw)
  terra::ifel(!is.na(r_net_totaal_opp) & r_net_totaal_opp > 0, 1, NA)
} else if (exists("final_opp") && !all(is.na(terra::values(final_opp, mat=FALSE)))) {
  terra::ifel(!is.na(final_opp) & final_opp > 0, 1, NA)
} else {
  terra::rast(template_MH, vals = NA)
}

# 3. Bepaal Analytisch Metacluster ID-raster
if (exists("resB_strikt") && !is.null(resB_strikt) && (!all(is.na(terra::values(resB_strikt$kern, mat=FALSE))) || !all(is.na(terra::values(resB_strikt$bouw, mat=FALSE))))) {
  r_net_totaal <- terra::cover(resB_strikt$kern, resB_strikt$bouw)
  r_buf        <- terra::buffer(r_net_totaal, width = 100)
  id_export_rast <- terra::mask(terra::patches(r_buf, directions = 8, zeroAsNA = TRUE), r_net_totaal)
} else if (exists("cl_max") && !all(is.na(terra::values(cl_max, mat=FALSE)))) {
  id_export_rast <- cl_max
} else if (exists("cl_opp") && !all(is.na(terra::values(cl_opp, mat=FALSE)))) {
  id_export_rast <- cl_opp
} else {
  id_export_rast <- terra::rast(template_MH, vals = NA)
}

# --- EXPORT LUS VOOR VISUELE BINAIR RASTERS ---
export_config <- list(
  list(rast = potentie_export_rast,  naam = "Maximale_Potentie",        folder = folders$potentie),
  list(rast = werkelijk_export_rast, naam = "Werkelijke_Oppervlaktes", folder = folders$werkelijk)
)

for(item in export_config) {
  suffix <- if (exists("zoek_sleutel") && grepl("_wv$", zoek_sleutel)) "_wv.tif" else ".tif"
  file_path <- file.path(item$folder, paste0("Habitat_", item$naam, "_", soort, suffix))
  
  export_rast <- terra::deepcopy(item$rast)
  export_rast <- terra::extend(export_rast, master_grid, fill = NA)
  
  terra::writeRaster(
    export_rast, 
    filename = file_path, 
    overwrite = TRUE, 
    gdal = c("COMPRESS=LZW"), 
    datatype = "INT1U", 
    NAflag = 255
  )
  message(paste("    [OK] Geëxporteerd:", basename(file_path)))
}

# --- EXPORT VOOR ANALYTISCH ID-RASTER ---
suffix <- if (exists("zoek_sleutel") && grepl("_wv$", zoek_sleutel)) "_wv.tif" else ".tif"
file_path_id <- file.path(folders$id_raster, paste0("ID_Netwerken_", soort, suffix))

export_id_rast <- terra::extend(id_export_rast, master_grid, fill = NA)

terra::writeRaster(
  export_id_rast, 
  filename = file_path_id, 
  overwrite = TRUE, 
  gdal = c("COMPRESS=LZW"), 
  datatype = "INT4U", 
  NAflag = 0
)
message(paste("    [OK] Analytisch ID-raster geëxporteerd:", basename(file_path_id)))

suppressWarnings(
  rm(potentie_export_rast, werkelijk_export_rast, id_export_rast, export_config, export_rast, export_id_rast)
)
gc()

message(paste("🏁 SCENARIO EXPORT VOLLEDIG AFGEROND VOOR:", toupper(soort)))

