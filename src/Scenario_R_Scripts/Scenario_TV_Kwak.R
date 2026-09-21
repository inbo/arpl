library(here)
library(knitr)
library(tidyverse)
library(sf)
library(terra)
library(readxl)
library(tidyterra)
library(data.table)

conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::select)
conflicted::conflicts_prefer(dplyr::first)
conflicted::conflicts_prefer(terra::intersect)
conflicted::conflicts_prefer(terra::any)

calc_ha_exact <- function(r) {
  if(is.null(r)) return(0)
  if(all(is.na(terra::values(r, mat=FALSE)))) return(0)
  area_raster <- r * terra::cellSize(r, unit = "ha")
  val <- terra::global(area_raster, "sum", na.rm = TRUE)[[1]]
  return(as.numeric(val))
}

maak_kaart_laag <- function(res_list) {
  if(is.null(res_list)) return(NULL)
  kernen <- terra::ifel(!is.na(res_list$kern), 1, NA)
  if(!is.null(res_list$bouw)) {
    bouw <- terra::ifel(!is.na(res_list$bouw), 2, NA)
    finaal <- terra::cover(kernen, bouw)
  } else {
    finaal = kernen
  }
  return(finaal)
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

terraOptions(
  memfrac = 0.8,
  tempdir = tempdir(),
  verbose = FALSE
)

soort <- "kwak"

# --- DYNAMISCHE SCENARIO PARAMETER CHECK ---
if (exists("SCENARIO_RDS_PAD") && !is.null(SCENARIO_RDS_PAD)) {
  scenario_rds_path <- SCENARIO_RDS_PAD
} else if (exists("params") && !is.null(params$scenario_rds_path)) {
  scenario_rds_path <- params$scenario_rds_path
} else {
  scenario_rds_path <- "data/input/Scenario_rds/TV_Scenario_BWK_2025.rds"
}

p_raw <- gsub("^([.][.]/)+", "", scenario_rds_path)
scenario_path <- here::here(p_raw)

if (!file.exists(scenario_path)) {
  stop(paste("❌ FOUT: Scenario RDS bestand NIET gevonden op:", scenario_path))
}

scen_volledig <- basename(scenario_path)
scenario_naam <- gsub("^TV_Scenario_|^Scenario_|.rds$", "", scen_volledig)

message(paste("Verwerken van soort:", soort, "binnen scenario:", scenario_naam))

df <- read_excel(here::here("data/input/Excel_files/Soorten_bwk_afstanden.xlsx"))
resultaat <- df %>%
  filter(tolower(trimws(Soort)) == soort) %>%
  select(Type, MinOpp_ha, AfstandBiotopen_m, Dispersiecap_m)

oppervlakte_ha <- resultaat$MinOpp_ha[1]
afstand_m      <- resultaat$AfstandBiotopen_m[1]
buffer_m       <- 10000 

rm(df, resultaat)

area_shape  <- vect(here("data/input/Turnhouts_Vennegebied.shp"))
master_grid <- rast(here("data/input/Raster_Vlaanderen/Vlaanderen_MasterGrid_10m.tif"))[[1]]

df_namen_sleutel <- read_csv(here("data/input/Excel_files/BWK_Laag_Namen_2025.csv"), show_col_types = FALSE)
gouden_namenlijst <- tolower(trimws(df_namen_sleutel$Laagnaam))

area_shape_proj <- project(area_shape, crs(master_grid))
area_buffer_fix <- buffer(area_shape_proj, width = buffer_m)

message("-> Vertaalraster voor globale/lokale cellen opbouwen via snelle MASK methode...")
id_raster_TV <- crop(master_grid, area_buffer_fix, snap = "near")

globale_id_raster <- master_grid
globale_id_raster <- terra::init(globale_id_raster, fun = "cell")

id_raster_TV_globale_values <- crop(globale_id_raster, area_buffer_fix, snap = "near")
id_raster_TV_masked <- mask(id_raster_TV_globale_values, area_buffer_fix)

message("-> Vertaaltabel bliksemsnel opbouwen via C++ dataframe extractie...")

df_extractie <- as.data.frame(id_raster_TV_masked, cells = TRUE)
vertaal_df <- as.data.table(df_extractie)
setnames(vertaal_df, c(1, 2), c("lokale_id", "globale_id"))

vertaal_df <- vertaal_df[!is.na(globale_id)]
studiegebied_globale_ids <- unique(vertaal_df$globale_id)

values(id_raster_TV) <- NA
template_TV <- terra::rasterize(area_buffer_fix, id_raster_TV, field = 1, background = 0)

grens_web <- sf::st_as_sf(terra::project(area_shape, "EPSG:4326"))

rm(globale_id_raster, id_raster_TV_globale_values, id_raster_TV_masked, df_extractie)
gc()

# 1. LAAD DE BRON-CROSSWALK/DICTIONARY IN
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
    bevat_codes_escaped  <- gsub("([\\.\\^\\$\\*\\+\\?\\(\\)\\[\\{\\\\\\|])", "\\\\\\1", bevat_codes)
    bevat_codes_anchored <- paste0("^", bevat_codes_escaped)
    regex_term           <- paste0(bevat_codes_anchored, collapse = "|")
    
    tabel_gefilterd <- tabel_vlaanderen[CODE %in% exact_codes | grepl(regex_term, CODE)]
  } else {
    tabel_gefilterd <- tabel_vlaanderen[CODE %in% exact_codes]
  }
  
  tabel_TV <- tabel_gefilterd[cel_id %in% studiegebied_globale_ids]
  tabel_cel_som <- tabel_TV[, .(Oppervlakte = sum(BWK_FRAC, na.rm = TRUE)), by = .(cel_id)]
  
  r_match_type <- id_raster_TV * NA
  r_opp_type   <- id_raster_TV * NA
  
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

openwater_max     <- lijst_matches[["openwater"]]
openwater_opp     <- lijst_oppervlaktes[["openwater"]]
voortplanting0_max <- lijst_matches[["voortplanting0"]]
voortplanting0_opp <- lijst_oppervlaktes[["voortplanting0"]]
not_eikenbos_max  <- lijst_matches[["not_eikenbos"]]
not_eikenbos_opp  <- lijst_oppervlaktes[["not_eikenbos"]]
moeras_max        <- lijst_matches[["moeras"]]
moeras_opp        <- lijst_oppervlaktes[["moeras"]]

rm(vertaal_df, lijst_matches, lijst_oppervlaktes, df_nieuw, resultaten_gegroepeerd)
gc()

# ==============================================================================
# KWAK STAP 1: VOORTPLANTINGSHABITAT
# ==============================================================================
message("-> Spoor 1: Broedplaatsen filteren en clusteren...")

r_kwak_voortplanting3_max <- template_TV * NA
r_kwak_voortplanting3_opp <- template_TV * NA

if (exists("voortplanting0_max") && !all(is.na(suppressWarnings(terra::minmax(voortplanting0_max))))) {
  
  # --- SPOOR A: MAX ---
  cl_voortp_max <- cluster_filter_compleet(
    masker     = voortplanting0_max,
    opp_laag   = voortplanting0_opp,
    drempel_m2 = 5 * 10000, 
    dist_m     = 100,       
    werkelijk  = FALSE
  )
  r_kwak_voortplanting2_max <- cl_voortp_max$raster
  
  if (!all(is.na(suppressWarnings(terra::minmax(r_kwak_voortplanting2_max))))) {
    if (exists("not_eikenbos_max") && !all(is.na(suppressWarnings(terra::minmax(not_eikenbos_max))))) {
      r_kwak_voortplanting3_max <- terra::mask(r_kwak_voortplanting2_max, not_eikenbos_max, inverse = TRUE)
    } else {
      r_kwak_voortplanting3_max <- r_kwak_voortplanting2_max
    }
  }
  
  # --- SPOOR B: OPP ---
  r_binair_voortp_opp <- terra::ifel(!is.na(voortplanting0_opp) & voortplanting0_opp > 0, 1, NA)
  cl_voortp_opp <- cluster_filter_compleet(
    masker     = r_binair_voortp_opp,
    opp_laag   = voortplanting0_opp,
    drempel_m2 = 5 * 10000, 
    dist_m     = 100,       
    werkelijk  = TRUE
  )
  r_kwak_voortplanting2_opp <- cl_voortp_opp$raster
  
  if (!all(is.na(suppressWarnings(terra::minmax(r_kwak_voortplanting2_opp))))) {
    if (exists("not_eikenbos_max") && !all(is.na(suppressWarnings(terra::minmax(not_eikenbos_max))))) {
      r_kwak_voortplanting3_opp <- terra::mask(r_kwak_voortplanting2_opp, not_eikenbos_max, inverse = TRUE)
    } else {
      r_kwak_voortplanting3_opp <- r_kwak_voortplanting2_opp
    }
  }
}

suppressWarnings(rm(cl_voortp_max, cl_voortp_opp, r_kwak_voortplanting2_max, r_kwak_voortplanting2_opp, r_binair_voortp_opp))
gc()

# ==============================================================================
# KWAK STAP 2: FOERAGEERHABITAT
# ==============================================================================
message("-> Spoor 2: Foerageergebieden (Procentuele moerastoets) verwerken...")

r_kwak_foerageer1_max <- template_TV * NA
r_kwak_foerageer1_opp <- template_TV * NA

if (exists("moeras_max") && !all(is.na(suppressWarnings(terra::minmax(moeras_max))))) {
  
  # --- SPOOR A: MAX ---
  cl_moeras_max <- cluster_filter_compleet(
    masker     = moeras_max,
    opp_laag   = moeras_opp,
    drempel_m2 = 10 * 10000, 
    dist_m     = 50,          
    werkelijk  = FALSE
  )
  r_moeras1_ids_max <- cl_moeras_max$clusters 
  
  if (!all(is.na(suppressWarnings(terra::minmax(r_moeras1_ids_max)))) && exists("openwater_opp")) {
    clean_moeras_opp <- terra::deepcopy(moeras_opp)
    clean_moeras_opp[is.na(clean_moeras_opp)] <- 0
    clean_water_opp  <- terra::deepcopy(openwater_opp)
    clean_water_opp[is.na(clean_water_opp)] <- 0
    
    som_moeras_max <- terra::zonal(clean_moeras_opp, r_moeras1_ids_max, fun = "sum", na.rm = TRUE)
    som_water_max  <- terra::zonal(clean_water_opp, r_moeras1_ids_max, fun = "sum", na.rm = TRUE)
    
    dt_toets_max <- data.table(ClusterID = som_moeras_max[[1]], Moeras_ha = som_moeras_max[[2]] * 0.01, Water_ha = som_water_max[[2]] * 0.01)
    dt_toets_max[, Water_Aandeel := Water_ha / Moeras_ha]
    goedgekeurde_foer_ids_max <- dt_toets_max[Water_Aandeel >= 0.20, ClusterID]
    
    if (length(goedgekeurde_foer_ids_max) > 0) {
      r_kwak_foerageer1_max <- terra::ifel(r_moeras1_ids_max %in% goedgekeurde_foer_ids_max, 1, NA)
    }
  }
  
  # --- SPOOR B: OPP ---
  r_binair_moeras_opp <- terra::ifel(!is.na(moeras_opp) & moeras_opp > 0, 1, NA)
  cl_moeras_opp <- cluster_filter_compleet(
    masker     = r_binair_moeras_opp,
    opp_laag   = moeras_opp,
    drempel_m2 = 10 * 10000, 
    dist_m     = 50,          
    werkelijk  = TRUE
  )
  r_moeras1_ids_opp <- cl_moeras_opp$clusters
  
  if (!all(is.na(suppressWarnings(terra::minmax(r_moeras1_ids_opp)))) && exists("openwater_opp")) {
    if (!exists("clean_moeras_opp")) {
      clean_moeras_opp <- terra::deepcopy(moeras_opp)
      clean_moeras_opp[is.na(clean_moeras_opp)] <- 0
      clean_water_opp  <- terra::deepcopy(openwater_opp)
      clean_water_opp[is.na(clean_water_opp)] <- 0
    }
    
    som_moeras_opp <- terra::zonal(clean_moeras_opp, r_moeras1_ids_opp, fun = "sum", na.rm = TRUE)
    som_water_opp  <- terra::zonal(clean_water_opp, r_moeras1_ids_opp, fun = "sum", na.rm = TRUE)
    
    dt_toets_opp <- data.table(ClusterID = som_moeras_opp[[1]], Moeras_ha = som_moeras_opp[[2]] * 0.01, Water_ha = som_water_opp[[2]] * 0.01)
    dt_toets_opp[, Water_Aandeel := Water_ha / Moeras_ha]
    goedgekeurde_foer_ids_opp <- dt_toets_opp[Water_Aandeel >= 0.20, ClusterID]
    
    if (length(goedgekeurde_foer_ids_opp) > 0) {
      r_kwak_foerageer1_opp <- terra::mask(moeras_opp, r_moeras1_ids_opp %in% goedgekeurde_foer_ids_opp)
    }
  }
}

suppressWarnings(rm(clean_moeras_opp, clean_water_opp, cl_moeras_max, cl_moeras_opp, r_moeras1_ids_max, r_moeras1_ids_opp, 
                    som_moeras_max, som_water_max, dt_toets_max, som_moeras_opp, som_water_opp, dt_toets_opp, r_binair_moeras_opp))
gc()

# ==============================================================================
# KWAK STAP 3: KOPPELING (5 KM) EN FINALE CLUSTERING
# ==============================================================================
message("-> Interactie-analyse: Wederzijdse 5 km commuter-toets toepassen...")

kwak_finaal_max <- template_TV * NA
kwak_finaal_opp <- template_TV * NA

r_tot_opp <- terra::cover(moeras_opp, voortplanting0_opp)

# --- SPOOR A: MAX ---
if (!all(is.na(suppressWarnings(terra::minmax(r_kwak_voortplanting3_max)))) && !all(is.na(suppressWarnings(terra::minmax(r_kwak_foerageer1_max))))) {
  dist_tot_foer_max   <- terra::distance(terra::ifel(r_kwak_foerageer1_max == 1, 1, NA))
  dist_tot_voortp_max <- terra::distance(terra::ifel(r_kwak_voortplanting3_max == 1, 1, NA))
  
  r_voortplanting1_max <- terra::ifel(dist_tot_foer_max <= 5000 & r_kwak_voortplanting3_max == 1, 1, NA)
  r_foerageer1_max     <- terra::ifel(dist_tot_voortp_max <= 5000 & r_kwak_foerageer1_max == 1, 1, NA)
  
  r_kwak_leefgebied1_max <- terra::cover(r_voortplanting1_max, r_foerageer1_max)
  
  cl_finaal_max <- cluster_filter_compleet(
    masker     = r_kwak_leefgebied1_max,
    opp_laag   = r_tot_opp,
    drempel_m2 = 10 * 10000, 
    dist_m     = 50,          
    werkelijk  = FALSE
  )
  kwak_finaal_max <- cl_finaal_max$raster %>% terra::crop(template_TV)
  rm(dist_tot_foer_max, dist_tot_voortp_max, r_voortplanting1_max, r_foerageer1_max, r_kwak_leefgebied1_max, cl_finaal_max)
}

# --- SPOOR B: OPP ---
if (!all(is.na(suppressWarnings(terra::minmax(r_kwak_voortplanting3_opp)))) && !all(is.na(suppressWarnings(terra::minmax(r_kwak_foerageer1_opp))))) {
  r_binair_vp3_opp  <- terra::ifel(!is.na(r_kwak_voortplanting3_opp) & r_kwak_voortplanting3_opp > 0, 1, NA)
  r_binair_foer_opp <- terra::ifel(!is.na(r_kwak_foerageer1_opp) & r_kwak_foerageer1_opp > 0, 1, NA)
  
  dist_tot_foer_opp   <- terra::distance(r_binair_foer_opp)
  dist_tot_voortp_opp <- terra::distance(r_binair_vp3_opp)
  
  r_voortplanting1_opp <- terra::mask(r_kwak_voortplanting3_opp, terra::ifel(dist_tot_foer_opp <= 5000, 1, NA))
  r_foerageer1_opp     <- terra::mask(r_kwak_foerageer1_opp, terra::ifel(dist_tot_voortp_opp <= 5000, 1, NA))
  
  r_kwak_leefgebied1_opp <- terra::cover(r_voortplanting1_opp, r_foerageer1_opp)
  r_binair_leef1_opp     <- terra::ifel(!is.na(r_kwak_leefgebied1_opp) & r_kwak_leefgebied1_opp > 0, 1, NA)
  
  cl_finaal_opp <- cluster_filter_compleet(
    masker     = r_binair_leef1_opp,
    opp_laag   = r_kwak_leefgebied1_opp,
    drempel_m2 = 10 * 10000, 
    dist_m     = 50,          
    werkelijk  = TRUE
  )
  kwak_finaal_opp <- cl_finaal_opp$raster %>% terra::crop(template_TV)
  rm(r_binair_vp3_opp, r_binair_foer_opp, dist_tot_foer_opp, dist_tot_voortp_opp, r_voortplanting1_opp, r_foerageer1_opp, r_kwak_leefgebied1_opp, r_binair_leef1_opp, cl_finaal_opp)
}

rm(r_kwak_voortplanting3_max, r_kwak_voortplanting3_opp, r_kwak_foerageer1_max, r_kwak_foerageer1_opp, r_tot_opp)
gc()

# DEFINITIEVE MODELUITGANGEN
final_max <- kwak_finaal_max
final_opp <- kwak_finaal_opp

if (!all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  cl_opp <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE)
} else {
  cl_opp <- template_TV * NA
}

if (!all(is.na(suppressWarnings(terra::minmax(final_max))))) {
  cl_max <- terra::patches(final_max, directions = 8, zeroAsNA = TRUE)
} else {
  cl_max <- template_TV * NA
}

# ==============================================================================
# SCHONE EXPORT BIOTOOP EN ANALYTISCH ID-RASTER (VOOR SCRIPT 2 / ARPL)
# ==============================================================================
base_dir <- here::here("data/output/Turnhouts_Vennegebied/Rasters_Soorten", scenario_naam)

folders <- list(
  potentie  = file.path(base_dir, "01_Maximale_Potentie"),
  werkelijk = file.path(base_dir, "02_Werkelijke_Oppervlaktes"),
  id_raster = file.path(base_dir, "00_ID_Rasters")
)
purrr::walk(folders, ~if (!dir.exists(.x)) dir.create(.x, showWarnings = FALSE, recursive = TRUE))

# 1. Bepaal Maximale Potentie Raster
potentie_export_rast <- if (exists("final_max") && !is.null(final_max) && !all(is.na(suppressWarnings(terra::minmax(final_max))))) {
  terra::ifel(!is.na(final_max) & final_max > 0, 1, NA)
} else {
  terra::rast(template_TV, vals = NA)
}

# 2. Bepaal Werkelijke Oppervlakte Raster
werkelijk_export_rast <- if (exists("final_opp") && !is.null(final_opp) && !all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  terra::ifel(!is.na(final_opp) & final_opp > 0, 1, NA)
} else {
  terra::rast(template_TV, vals = NA)
}

# 3. Bepaal Analytisch Metacluster ID-raster (EXCLUSIEF OP BASIS VAN WERKELIJKE OPPERVLAKTES)
if (exists("cl_opp") && !is.null(cl_opp) && !all(is.na(suppressWarnings(terra::minmax(cl_opp))))) {
  id_export_rast <- cl_opp
} else if (exists("final_opp") && !is.null(final_opp) && !all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  id_export_rast <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE)
} else {
  id_export_rast <- terra::rast(template_TV, vals = NA)
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
