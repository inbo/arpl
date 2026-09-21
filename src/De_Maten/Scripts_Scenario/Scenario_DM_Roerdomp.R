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

soort <- "roerdomp"

# --- DYNAMISCHE SCENARIO PARAMETER CHECK ---
if (exists("SCENARIO_RDS_PAD") && !is.null(SCENARIO_RDS_PAD)) {
  scenario_rds_path <- SCENARIO_RDS_PAD
} else if (exists("params") && !is.null(params$scenario_rds_path)) {
  scenario_rds_path <- params$scenario_rds_path
} else {
  scenario_rds_path <- "data/input/Scenario_rds/DM_Scenario_BWK_2025.rds"
}

p_raw <- gsub("^([.][.]/)+", "", scenario_rds_path)
scenario_path <- here::here(p_raw)

if (!file.exists(scenario_path)) {
  stop(paste("❌ FOUT: Scenario RDS bestand NIET gevonden op:", scenario_path))
}

scen_volledig <- basename(scenario_path)
scenario_naam <- gsub("^DM_Scenario_|^Scenario_|.rds$", "", scen_volledig)

message(paste("Verwerken van soort:", soort, "binnen scenario:", scenario_naam))

df <- read_excel(here::here("data/input/Excel_files/Soorten_bwk_afstanden.xlsx"))
resultaat <- df %>%
  filter(tolower(trimws(Soort)) == soort) %>%
  select(Type, MinOpp_ha, AfstandBiotopen_m, Dispersiecap_m)

oppervlakte_ha <- resultaat$MinOpp_ha[1]
afstand_m      <- resultaat$AfstandBiotopen_m[1]
buffer_m       <- resultaat$Dispersiecap_m[1]

rm(df, resultaat)

area_shape  <- vect(here("data/input/De_Maten.shp"))
master_grid <- rast(here("data/input/Raster_Vlaanderen/Vlaanderen_MasterGrid_10m.tif"))[[1]]

df_namen_sleutel <- read_csv(here("data/input/Excel_files/BWK_Laag_Namen_2025.csv"), show_col_types = FALSE)
gouden_namenlijst <- tolower(trimws(df_namen_sleutel$Laagnaam))

area_shape_proj <- project(area_shape, crs(master_grid))
area_buffer_fix <- buffer(area_shape_proj, width = buffer_m)

message("-> Vertaalraster voor globale/lokale cellen opbouwen via snelle MASK methode...")
id_raster_DM <- crop(master_grid, area_buffer_fix, snap = "near")

globale_id_raster <- master_grid
globale_id_raster <- terra::init(globale_id_raster, fun = "cell")

id_raster_DM_globale_values <- crop(globale_id_raster, area_buffer_fix, snap = "near")
id_raster_DM_masked <- mask(id_raster_DM_globale_values, area_buffer_fix)

message("-> Vertaaltabel bliksemsnel opbouwen via C++ dataframe extractie...")

df_extractie <- as.data.frame(id_raster_DM_masked, cells = TRUE)
vertaal_df <- as.data.table(df_extractie)
setnames(vertaal_df, c(1, 2), c("lokale_id", "globale_id"))

vertaal_df <- vertaal_df[!is.na(globale_id)]
studiegebied_globale_ids <- unique(vertaal_df$globale_id)

values(id_raster_DM) <- NA
template_DM <- terra::rasterize(area_buffer_fix, id_raster_DM, field = 1, background = 0)

grens_web <- sf::st_as_sf(terra::project(area_shape, "EPSG:4326"))

rm(globale_id_raster, id_raster_DM_globale_values, id_raster_DM_masked, df_extractie)
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
  
  tabel_DM <- tabel_gefilterd[cel_id %in% studiegebied_globale_ids]
  tabel_DM_unique <- unique(tabel_DM, by = c("cel_id", "CODE"))
  
  tabel_cel_som <- tabel_DM_unique[, .(Oppervlakte = pmin(sum(BWK_FRAC, na.rm = TRUE), 1.0)), by = .(cel_id)]
  
  r_match_type <- id_raster_DM * NA
  r_opp_type   <- id_raster_DM * NA
  
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

openwater_max <- lijst_matches[["openwater"]]
openwater_opp <- lijst_oppervlaktes[["openwater"]]
rietland_max  <- lijst_matches[["rietland"]]
rietland_opp  <- lijst_oppervlaktes[["rietland"]]
zegge_max     <- lijst_matches[["zegge"]]
zegge_opp     <- lijst_oppervlaktes[["zegge"]]
moeras_max    <- lijst_matches[["moeras"]]
moeras_opp    <- lijst_oppervlaktes[["moeras"]]
grasland_max  <- lijst_matches[["grasland"]]
grasland_opp  <- lijst_oppervlaktes[["grasland"]]
struweel_max  <- lijst_matches[["struweel"]]
struweel_opp  <- lijst_oppervlaktes[["struweel"]]

message("=== STAP 1: Basishabitat bundelen en storing uitsluiten (PARALLEL) ===")

max_lagen <- list(openwater_max, rietland_max, zegge_max, moeras_max, grasland_max, struweel_max)
opp_lagen <- list(openwater_opp, rietland_opp, zegge_opp, moeras_opp, grasland_opp, struweel_opp)

max_lagen <- max_lagen[!sapply(max_lagen, is.null)]
opp_lagen <- opp_lagen[!sapply(opp_lagen, is.null)]

# Spoor A: MAX
biotoop_max <- max_lagen[[1]]
if(length(max_lagen) > 1) {
  for(i in 2:length(max_lagen)) biotoop_max <- terra::cover(biotoop_max, max_lagen[[i]])
}

# Spoor B: OPP
biotoop_opp_som <- template_DM * 0
for(i in 1:length(opp_lagen)) {
  clean_opp <- terra::deepcopy(opp_lagen[[i]])
  clean_opp[is.na(clean_opp)] <- 0
  biotoop_opp_som <- biotoop_opp_som + clean_opp
}
biotoop_opp_raw <- terra::clamp(biotoop_opp_som, upper = 1.0)
biotoop_opp     <- terra::ifel(!is.na(biotoop_max) & biotoop_opp_raw > 0, biotoop_opp_raw, NA)

groenkaart <- terra::rast(here("data/input/ASCI Files/Groenkaart_2021.tif"))
groenkaart_DM <- terra::crop(groenkaart, template_DM, snap = "near")
hooggroen_aligned <- terra::resample(groenkaart_DM, template_DM, method = "near")

r_bos_bin <- terra::ifel(hooggroen_aligned == 1, 1, NA)

tabel_estuarium <- tabel_vlaanderen[CODE == "1130"]
tabel_estuarium_DM <- tabel_estuarium[cel_id %in% studiegebied_globale_ids]
r_estuarium <- id_raster_DM * NA

if(nrow(tabel_estuarium_DM) > 0) {
  tabel_est_mapping <- copy(tabel_estuarium_DM)
  setnames(tabel_est_mapping, "cel_id", "globale_id")
  tabel_finaal_est <- merge(tabel_est_mapping, vertaal_df, by = "globale_id", all.x = TRUE)
  tabel_finaal_est <- tabel_finaal_est[!is.na(lokale_id)]
  if(nrow(tabel_finaal_est) > 0) r_estuarium[tabel_finaal_est$lokale_id] <- 1
}

stoorlaag <- terra::cover(r_bos_bin, r_estuarium)

# Filteren van de hoofdlagen (Parallel)
roerdomp_leefgebied1_max <- terra::mask(biotoop_max, stoorlaag, inverse = TRUE)
roerdomp_leefgebied1_opp <- terra::mask(biotoop_opp, stoorlaag, inverse = TRUE)

rietland_clean_max  <- terra::mask(rietland_max, stoorlaag, inverse = TRUE)
openwater_clean_max <- terra::mask(openwater_max, stoorlaag, inverse = TRUE)

rietland_clean_opp  <- terra::mask(rietland_opp, stoorlaag, inverse = TRUE)
openwater_clean_opp <- terra::mask(openwater_opp, stoorlaag, inverse = TRUE)

rm(biotoop_max, biotoop_opp_som, biotoop_opp_raw, biotoop_opp, groenkaart, groenkaart_DM, hooggroen_aligned, 
   r_bos_bin, r_estuarium, stoorlaag, tabel_estuarium, tabel_estuarium_DM, clean_opp, vertaal_df, lijst_matches, lijst_oppervlaktes, tabel_vlaanderen, resultaten_gegroepeerd, df_nieuw)
gc()

message("=== STAP 2: Open water partijen filteren op minimale grootte (10 ha - PARALLEL) ===")

drempel_10ha_m2 <- 10 * 10000

# Spoor A: MAX
water_clusters_max <- cluster_filter_compleet(
  masker     = openwater_clean_max,
  opp_laag   = openwater_clean_opp,
  drempel_m2 = drempel_10ha_m2,
  dist_m     = 50,
  werkelijk  = FALSE
)
r_grote_plassen_max <- water_clusters_max$raster

if (terra::global(is.na(r_grote_plassen_max), "sum")[[1]] < terra::ncell(r_grote_plassen_max)) {
  dist_tot_water_max <- terra::distance(r_grote_plassen_max)
  mal_nabij_water_max <- terra::ifel(dist_tot_water_max <= 100, 1, NA)
  
  roerdomp_leefgebied2_max <- terra::mask(roerdomp_leefgebied1_max, mal_nabij_water_max)
} else {
  roerdomp_leefgebied2_max <- template_DM * NA
}

# Spoor B: OPP (Parallel & Autonoom)
r_binair_openw_opp <- terra::ifel(!is.na(openwater_clean_opp) & openwater_clean_opp > 0, 1, NA)
water_clusters_opp <- cluster_filter_compleet(
  masker     = r_binair_openw_opp,
  opp_laag   = openwater_clean_opp,
  drempel_m2 = drempel_10ha_m2,
  dist_m     = 50,
  werkelijk  = TRUE
)
r_grote_plassen_opp <- water_clusters_opp$raster

if (terra::global(is.na(r_grote_plassen_opp), "sum")[[1]] < terra::ncell(r_grote_plassen_opp)) {
  dist_tot_water_opp <- terra::distance(!is.na(r_grote_plassen_opp) & r_grote_plassen_opp > 0)
  mal_nabij_water_opp <- terra::ifel(dist_tot_water_opp <= 100, 1, NA)
  
  roerdomp_leefgebied2_opp <- terra::mask(roerdomp_leefgebied1_opp, mal_nabij_water_opp)
} else {
  roerdomp_leefgebied2_opp <- template_DM * NA
}

rm(water_clusters_max, water_clusters_opp, r_grote_plassen_max, r_grote_plassen_opp, 
   r_binair_openw_opp, roerdomp_leefgebied1_max, roerdomp_leefgebied1_opp)
gc()

message("=== STAP 3: Rietkragen / Ecotonen analyseren (Minimaal 500m waterkant - PARALLEL) ===")

drempel_waterkant_m2 <- 5000

# Spoor A: MAX
riet_buffer_max <- terra::buffer(rietland_clean_max, width = 10)
r_waterkant_max <- terra::mask(openwater_clean_max, riet_buffer_max)

waterkant_gefilterd_max <- cluster_filter_compleet(
  masker     = r_waterkant_max,
  opp_laag   = NULL,
  drempel_m2 = drempel_waterkant_m2,
  dist_m     = 50,
  werkelijk  = FALSE
)
r_robuuste_waterkant_max <- waterkant_gefilterd_max$raster

if (terra::global(is.na(r_robuuste_waterkant_max), "sum")[[1]] < terra::ncell(r_robuuste_waterkant_max)) {
  dist_tot_waterkant_max <- terra::distance(r_robuuste_waterkant_max)
  mal_nabij_waterkant_max <- terra::ifel(dist_tot_waterkant_max <= 100, 1, NA)
  
  roerdomp_leefgebied3_max <- terra::mask(roerdomp_leefgebied2_max, mal_nabij_waterkant_max)
} else {
  roerdomp_leefgebied3_max <- template_DM * NA
}

# Spoor B: OPP (Parallel & Autonoom)
r_binair_riet_opp <- terra::ifel(!is.na(rietland_clean_opp) & rietland_clean_opp > 0, 1, NA)
r_binair_open_opp <- terra::ifel(!is.na(openwater_clean_opp) & openwater_clean_opp > 0, 1, NA)

riet_buffer_opp <- terra::buffer(r_binair_riet_opp, width = 10)
r_waterkant_opp <- terra::mask(openwater_clean_opp, riet_buffer_opp)

waterkant_gefilterd_opp <- cluster_filter_compleet(
  masker     = r_binair_open_opp & !is.na(r_waterkant_opp),
  opp_laag   = r_waterkant_opp,
  drempel_m2 = drempel_waterkant_m2,
  dist_m     = 50,
  werkelijk  = TRUE
)
r_robuuste_waterkant_opp <- waterkant_gefilterd_opp$raster

if (terra::global(is.na(r_robuuste_waterkant_opp), "sum")[[1]] < terra::ncell(r_robuuste_waterkant_opp)) {
  dist_tot_waterkant_opp <- terra::distance(!is.na(r_robuuste_waterkant_opp) & r_robuuste_waterkant_opp > 0)
  mal_nabij_waterkant_opp <- terra::ifel(dist_tot_waterkant_opp <= 100, 1, NA)
  
  roerdomp_leefgebied3_opp <- terra::mask(roerdomp_leefgebied2_opp, mal_nabij_waterkant_opp)
} else {
  roerdomp_leefgebied3_opp <- template_DM * NA
}

rm(riet_buffer_max, riet_buffer_opp, r_waterkant_max, r_waterkant_opp, waterkant_gefilterd_max, waterkant_gefilterd_opp,
   r_robuuste_waterkant_max, r_robuuste_waterkant_opp, r_binair_riet_opp, r_binair_open_opp, roerdomp_leefgebied2_max, roerdomp_leefgebied2_opp)
gc()

message("=== STAP 4: Nestkernen filteren (Minimaal 0.5 ha zuiver rietland - PARALLEL) ===")

drempel_05ha_m2 <- 0.5 * 10000

# Spoor A: MAX
riet_kernen_max <- cluster_filter_compleet(
  masker     = rietland_clean_max,
  opp_laag   = rietland_clean_opp,
  drempel_m2 = drempel_05ha_m2,
  dist_m     = 50,
  werkelijk  = FALSE
)
r_riet_kernen_max <- riet_kernen_max$raster

if (terra::global(is.na(r_riet_kernen_max), "sum")[[1]] < terra::ncell(r_riet_kernen_max)) {
  dist_tot_rietkern_max <- terra::distance(r_riet_kernen_max)
  mal_nabij_rietkern_max <- terra::ifel(dist_tot_rietkern_max <= 100, 1, NA)
  
  roerdomp_leefgebied4_max <- terra::mask(roerdomp_leefgebied3_max, mal_nabij_rietkern_max)
} else {
  roerdomp_leefgebied4_max <- template_DM * NA
}

# Spoor B: OPP (Parallel & Autonoom)
r_binair_riet2_opp <- terra::ifel(!is.na(rietland_clean_opp) & rietland_clean_opp > 0, 1, NA)
riet_kernen_opp <- cluster_filter_compleet(
  masker     = r_binair_riet2_opp,
  opp_laag   = rietland_clean_opp,
  drempel_m2 = drempel_05ha_m2,
  dist_m     = 50,
  werkelijk  = TRUE
)
r_riet_kernen_opp <- riet_kernen_opp$raster

if (terra::global(is.na(r_riet_kernen_opp), "sum")[[1]] < terra::ncell(r_riet_kernen_opp)) {
  dist_tot_rietkern_opp <- terra::distance(!is.na(r_riet_kernen_opp) & r_riet_kernen_opp > 0)
  mal_nabij_rietkern_opp <- terra::ifel(dist_tot_rietkern_opp <= 100, 1, NA)
  
  roerdomp_leefgebied4_opp <- terra::mask(roerdomp_leefgebied3_opp, mal_nabij_rietkern_opp)
} else {
  roerdomp_leefgebied4_opp <- template_DM * NA
}

rm(riet_kernen_max, riet_kernen_opp, r_riet_kernen_max, r_riet_kernen_opp, r_binair_riet2_opp, 
   roerdomp_leefgebied3_max, roerdomp_leefgebied3_opp, rietland_clean_max, rietland_clean_opp, openwater_clean_max, openwater_clean_opp)
gc()

message("=== STAP 5: Grote metapopulatie-clustering op 25 hectare (PARALLEL) ===")

drempel_25ha_m2 <- 25 * 10000

# Finaal Spoor 1: Maximale Potentie
finaal_roerdomp_max <- cluster_filter_compleet(
  masker     = roerdomp_leefgebied4_max,
  opp_laag   = roerdomp_leefgebied4_opp,
  drempel_m2 = drempel_25ha_m2,
  dist_m     = 500, 
  werkelijk  = FALSE
)
finaal_leefgebied_max <- finaal_roerdomp_max$raster

# Finaal Spoor 2: Werkelijke Oppervlakte
r_binair_leef4_opp <- terra::ifel(!is.na(roerdomp_leefgebied4_opp) & roerdomp_leefgebied4_opp > 0, 1, NA)
finaal_roerdomp_opp <- cluster_filter_compleet(
  masker     = r_binair_leef4_opp,
  opp_laag   = roerdomp_leefgebied4_opp,
  drempel_m2 = drempel_25ha_m2,
  dist_m     = 500,
  werkelijk  = TRUE
)
finaal_leefgebied_opp <- finaal_roerdomp_opp$raster

# DEFINITIEVE MODELUITGANGEN
final_max <- finaal_leefgebied_max
final_opp <- finaal_leefgebied_opp

if (!all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  cl_opp <- finaal_roerdomp_opp$clusters
} else {
  cl_opp <- template_DM * NA
}

if (!all(is.na(suppressWarnings(terra::minmax(final_max))))) {
  cl_max <- finaal_roerdomp_max$clusters
} else {
  cl_max <- template_DM * NA
}

rm(finaal_roerdomp_max, finaal_roerdomp_opp, roerdomp_leefgebied4_max, roerdomp_leefgebied4_opp, r_binair_leef4_opp)
gc()

# ==============================================================================
# SCHONE EXPORT BIOTOOP EN ANALYTISCH ID-RASTER (VOOR SCRIPT 2 / ARPL)
# ==============================================================================
base_dir <- here::here("data/output/De_Maten/Rasters_Soorten", scenario_naam)

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
  terra::rast(template_DM, vals = NA)
}

# 2. Bepaal Werkelijke Oppervlakte Raster
werkelijk_export_rast <- if (exists("final_opp") && !is.null(final_opp) && !all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  terra::ifel(!is.na(final_opp) & final_opp > 0, 1, NA)
} else {
  terra::rast(template_DM, vals = NA)
}

# 3. Bepaal Analytisch Metacluster ID-raster (EXCLUSIEF OP BASIS VAN WERKELIJKE OPPERVLAKTES)
if (exists("cl_opp") && !is.null(cl_opp) && !all(is.na(suppressWarnings(terra::minmax(cl_opp))))) {
  id_export_rast <- cl_opp
} else if (exists("final_opp") && !is.null(final_opp) && !all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  id_export_rast <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE)
} else {
  id_export_rast <- template_DM * NA
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
  message(paste("    [OK] Geëxporteerd naar scenariomap:", basename(file_path)))
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
