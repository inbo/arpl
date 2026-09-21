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

soort <- "noordsewitsnuitlibel"

# --- DYNAMISCHE SCENARIO PARAMETER CHECK ---
if (exists("SCENARIO_RDS_PAD") && !is.null(SCENARIO_RDS_PAD)) {
  scenario_rds_path <- SCENARIO_RDS_PAD
} else if (exists("params") && !is.null(params$scenario_rds_path)) {
  scenario_rds_path <- params$scenario_rds_path
} else {
  scenario_rds_path <- "data/input/Scenario_rds/MH_Scenario_BWK_2025.rds"
}

p_raw <- gsub("^([.][.]/)+", "", scenario_rds_path)
scenario_path <- here::here(p_raw)

if (!file.exists(scenario_path)) {
  stop(paste("❌ FOUT: Scenario RDS bestand NIET gevonden op:", scenario_path))
}

scen_volledig <- basename(scenario_path)
scenario_naam <- gsub("^MH_Scenario_|^Scenario_|.rds$", "", scen_volledig)

message(paste("Verwerken van soort:", soort, "binnen scenario:", scenario_naam))

df <- read_excel(here::here("data/input/Excel_files/Soorten_bwk_afstanden.xlsx"))
resultaat <- df %>%
  filter(tolower(trimws(Soort)) == soort) %>%
  select(Type, MinOpp_ha, AfstandBiotopen_m, Dispersiecap_m)

buffer_m <- resultaat$Dispersiecap_m[1]
straal_water_m <- 500  # Harde afsnij-afstand rond het water

rm(df, resultaat)

area_shape  <- vect(here("data/input/Mechelse_Heide.shp"))
master_grid <- rast(here("data/input/Raster_Vlaanderen/Vlaanderen_MasterGrid_10m.tif"))[[1]]

df_namen_sleutel <- read_csv(here("data/input/Excel_files/BWK_Laag_Namen_2025.csv"), show_col_types = FALSE)
gouden_namenlijst <- tolower(trimws(df_namen_sleutel$Laagnaam))

area_shape_proj <- project(area_shape, crs(master_grid))
area_buffer_fix <- buffer(area_shape_proj, width = buffer_m)

message("-> Vertaalraster voor globale/lokale cellen opbouwen via snelle MASK methode...")
id_raster_MH <- crop(master_grid, area_buffer_fix, snap = "near")

globale_id_raster <- master_grid
globale_id_raster <- terra::init(globale_id_raster, fun = "cell")

id_raster_MH_globale_values <- crop(globale_id_raster, area_buffer_fix, snap = "near")
id_raster_MH_masked <- mask(id_raster_MH_globale_values, area_buffer_fix)

message("-> Vertaaltabel bliksemsnel opbouwen via C++ dataframe extractie...")

df_extractie <- as.data.frame(id_raster_MH_masked, cells = TRUE)
vertaal_df <- as.data.table(df_extractie)
setnames(vertaal_df, c(1, 2), c("lokale_id", "globale_id"))

vertaal_df <- vertaal_df[!is.na(globale_id)]
studiegebied_globale_ids <- unique(vertaal_df$globale_id)

values(id_raster_MH) <- NA
template_MH <- terra::rasterize(area_buffer_fix, id_raster_MH, field = 1, background = 0)

grens_web <- sf::st_as_sf(terra::project(area_shape, "EPSG:4326"))

rm(globale_id_raster, id_raster_MH_globale_values, id_raster_MH_masked, df_extractie)
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

waterbiotoop_bwk1_max        <- lijst_matches[["waterbiotoop_bwk1"]]
waterbiotoop_bwk1_opp        <- lijst_oppervlaktes[["waterbiotoop_bwk1"]]
waterbiotoop_notbwk_max      <- lijst_matches[["waterbiotoop_notbwk"]]
waterbiotoop_notbwk_opp      <- lijst_oppervlaktes[["waterbiotoop_notbwk"]]
landbiotoop_bwk1_max         <- lijst_matches[["landbiotoop_bwk"]]
landbiotoop_bwk1_opp         <- lijst_oppervlaktes[["landbiotoop_bwk"]]
landbiotoop_bwk_ruim_max     <- lijst_matches[["landbiotoop_bwk_ruim"]]
landbiotoop_bwk_ruim_opp     <- lijst_oppervlaktes[["landbiotoop_bwk_ruim"]]
landbiotoop_bwk_ruim_bos_max <- lijst_matches[["landbiotoop_bwk_ruim_bos"]]
landbiotoop_bwk_ruim_bos_opp <- lijst_oppervlaktes[["landbiotoop_bwk_ruim_bos"]]

rm(vertaal_df, lijst_matches, lijst_oppervlaktes, tabel_vlaanderen, resultaten_gegroepeerd, df_nieuw)
gc()

message("-> Waterbiotoop berekenen (Spoor Max & Opp)...")

waterbiotoop_bwk_basis_max <- terra::ifel(!is.na(waterbiotoop_bwk1_max) & is.na(waterbiotoop_notbwk_max), 1, NA)
waterbiotoop_bwk_basis_opp <- terra::ifel(!is.na(waterbiotoop_bwk1_opp) & waterbiotoop_bwk1_opp > 0 & (is.na(waterbiotoop_notbwk_opp) | waterbiotoop_notbwk_opp == 0), 1, NA)

# --- SPOOR 1: MAXIMALE POTENTIE ---
waterbiotoop1_max_list <- cluster_filter_compleet(
  masker     = waterbiotoop_bwk_basis_max, 
  opp_laag   = waterbiotoop_bwk1_max, 
  drempel_m2 = 200,                  # 0.02 ha = 200 m²
  dist_m     = 50, 
  werkelijk  = FALSE
)
waterbiotoop1_max <- waterbiotoop1_max_list$raster

# --- SPOOR 2: WERKELIJKE OPPERVLAKTE ---
r_binair_water_opp <- terra::ifel(!is.na(waterbiotoop_bwk_basis_opp) & waterbiotoop_bwk_basis_opp > 0, 1, NA)
waterbiotoop1_opp_list <- cluster_filter_compleet(
  masker     = r_binair_water_opp, 
  opp_laag   = waterbiotoop_bwk1_opp, 
  drempel_m2 = 200, 
  dist_m     = 50, 
  werkelijk  = TRUE
)
waterbiotoop1_opp <- waterbiotoop1_opp_list$raster

rm(r_binair_water_opp, waterbiotoop_bwk_basis_max, waterbiotoop_bwk_basis_opp, waterbiotoop1_max_list, waterbiotoop1_opp_list)
gc()

message("-> Landbiotoop nabije omgeving berekenen (Spoor Max & Opp)...")

landbiotoop_bwk_basis_max <- terra::ifel(!is.na(landbiotoop_bwk1_max) & is.na(waterbiotoop_notbwk_max), 1, NA)
landbiotoop_bwk_basis_opp <- terra::ifel(!is.na(landbiotoop_bwk1_opp) & landbiotoop_bwk1_opp > 0 & (is.na(waterbiotoop_notbwk_opp) | waterbiotoop_notbwk_opp == 0), 1, NA)

# --- SPOOR 1: MAXIMALE POTENTIE ---
landbiotoop1_max_list <- cluster_filter_compleet(
  masker     = landbiotoop_bwk_basis_max, 
  opp_laag   = landbiotoop_bwk1_max, 
  drempel_m2 = 50000, 
  dist_m     = 20, 
  werkelijk  = FALSE
)
landbiotoop_nabij_max <- landbiotoop1_max_list$raster

# --- SPOOR 2: WERKELIJKE OPPERVLAKTE ---
r_binair_land_opp <- terra::ifel(!is.na(landbiotoop_bwk_basis_opp) & landbiotoop_bwk_basis_opp > 0, 1, NA)
landbiotoop1_opp_list <- cluster_filter_compleet(
  masker     = r_binair_land_opp, 
  opp_laag   = landbiotoop_bwk1_opp, 
  drempel_m2 = 50000, 
  dist_m     = 20, 
  werkelijk  = TRUE
)
landbiotoop_nabij_opp <- landbiotoop1_opp_list$raster

rm(landbiotoop_bwk_basis_max, landbiotoop_bwk_basis_opp, r_binair_land_opp, landbiotoop1_max_list, landbiotoop1_opp_list)
gc()

message("-> Landbiotoop ruimere omgeving & bos-samenstellingscheck berekenen (Minimaal 40% bos - PARALLEL)...")

# --- SPOOR 1: MAXIMALE POTENTIE ---
landbiotoop_ruim_max_list <- cluster_filter_compleet(
  masker     = landbiotoop_bwk_ruim_max, 
  opp_laag   = landbiotoop_bwk_ruim_max, 
  drempel_m2 = 400000, 
  dist_m     = 50, 
  werkelijk  = FALSE
)
r_ruim_max_src <- landbiotoop_ruim_max_list$raster
r_ruim_max_src[is.na(r_ruim_max_src) | r_ruim_max_src == 0] <- NA

r_bos_max_src <- landbiotoop_bwk_ruim_bos_max
r_bos_max_src[is.na(r_bos_max_src) | r_bos_max_src == 0] <- NA

cl_id_ruim_max <- terra::patches(r_ruim_max_src, directions = 8, zeroAsNA = TRUE)

stats_ruim_totaal_max <- terra::zonal(r_ruim_max_src * 0.01, cl_id_ruim_max, fun = "sum", na.rm = TRUE)
colnames(stats_ruim_totaal_max) <- c("Cluster_ID", "Totaal_Ruim_ha")

bos_in_ruim_max <- terra::mask(r_bos_max_src, cl_id_ruim_max)
stats_bos_in_ruim_max <- terra::zonal(bos_in_ruim_max * 0.01, cl_id_ruim_max, fun = "sum", na.rm = TRUE)
colnames(stats_bos_in_ruim_max) <- c("Cluster_ID", "Bos_ha")

df_samenstelling_max <- merge(stats_ruim_totaal_max, stats_bos_in_ruim_max, by = "Cluster_ID", all.x = TRUE)
df_samenstelling_max$Bos_ha[is.na(df_samenstelling_max$Bos_ha)] <- 0
df_samenstelling_max$Percentage_Bos <- (df_samenstelling_max$Bos_ha / df_samenstelling_max$Totaal_Ruim_ha) * 100

valide_ruim_ids_max <- df_samenstelling_max$Cluster_ID[df_samenstelling_max$Percentage_Bos >= 40]

if(length(valide_ruim_ids_max) > 0) {
  masker_ruim_max <- cl_id_ruim_max %in% valide_ruim_ids_max
  landbiotoop_ruim_gekoppeld_max <- terra::mask(r_ruim_max_src, terra::ifel(masker_ruim_max, 1, NA))
} else {
  landbiotoop_ruim_gekoppeld_max <- template_MH * NA
}

# --- SPOOR 2: WERKELIJKE OPPERVLAKTE ---
r_binair_ruim_opp_init <- terra::ifel(!is.na(landbiotoop_bwk_ruim_opp) & landbiotoop_bwk_ruim_opp > 0, 1, NA)

landbiotoop_ruim_opp_list <- cluster_filter_compleet(
  masker     = r_binair_ruim_opp_init, 
  opp_laag   = landbiotoop_bwk_ruim_opp, 
  drempel_m2 = 400000, 
  dist_m     = 50, 
  werkelijk  = TRUE
)
r_ruim_opp_src <- landbiotoop_ruim_opp_list$raster
r_ruim_opp_bin <- terra::ifel(!is.na(r_ruim_opp_src) & r_ruim_opp_src > 0, 1, NA)

if(!all(is.na(suppressWarnings(terra::minmax(r_ruim_opp_bin))))) {
  cl_id_ruim_opp <- terra::patches(r_ruim_opp_bin, directions = 8, zeroAsNA = TRUE)
  
  stats_ruim_totaal_opp <- terra::zonal(r_ruim_opp_src, cl_id_ruim_opp, fun = "sum", na.rm = TRUE)
  colnames(stats_ruim_totaal_opp) <- c("Cluster_ID", "Totaal_Ruim_ha")
  stats_ruim_totaal_opp$Totaal_Ruim_ha <- stats_ruim_totaal_opp$Totaal_Ruim_ha * 0.01
  
  bos_in_ruim_opp <- terra::mask(landbiotoop_bwk_ruim_bos_opp, cl_id_ruim_opp)
  stats_bos_in_ruim_opp <- terra::zonal(bos_in_ruim_opp, cl_id_ruim_opp, fun = "sum", na.rm = TRUE)
  colnames(stats_bos_in_ruim_opp) <- c("Cluster_ID", "Bos_ha")
  stats_bos_in_ruim_opp$Bos_ha <- stats_bos_in_ruim_opp$Bos_ha * 0.01
  
  df_samenstelling_opp <- merge(stats_ruim_totaal_opp, stats_bos_in_ruim_opp, by = "Cluster_ID", all.x = TRUE)
  df_samenstelling_opp$Bos_ha[is.na(df_samenstelling_opp$Bos_ha)] <- 0
  df_samenstelling_opp$Percentage_Bos <- (df_samenstelling_opp$Bos_ha / df_samenstelling_opp$Totaal_Ruim_ha) * 100
  
  valide_ruim_ids_opp <- df_samenstelling_opp$Cluster_ID[df_samenstelling_opp$Percentage_Bos >= 40]
  
  if(length(valide_ruim_ids_opp) > 0) {
    masker_ruim_opp <- cl_id_ruim_opp %in% valide_ruim_ids_opp
    landbiotoop_ruim_gekoppeld_opp <- terra::mask(r_ruim_opp_src, terra::ifel(masker_ruim_opp, 1, NA))
  } else {
    landbiotoop_ruim_gekoppeld_opp <- template_MH * NA
  }
} else {
  landbiotoop_ruim_gekoppeld_opp <- template_MH * NA
}

rm(landbiotoop_ruim_max_list, r_ruim_max_src, r_bos_max_src, cl_id_ruim_max,
   stats_ruim_totaal_max, bos_in_ruim_max, stats_bos_in_ruim_max, df_samenstelling_max,
   r_binair_ruim_opp_init, landbiotoop_ruim_opp_list, r_ruim_opp_src, r_ruim_opp_bin, cl_id_ruim_opp,
   stats_ruim_totaal_opp, bos_in_ruim_opp, stats_bos_in_ruim_opp, df_samenstelling_opp)
gc()

message("-> Originele ruimtelijke interacties uitvoeren + Water-First 500m afsnijding...")

straal_water_m <- 500

# ==============================================================================
# SPOOR 1: MAXIMALE POTENTIE
# ==============================================================================
water_buffer10_max <- terra::buffer(waterbiotoop1_max, width = 10)
land_buffer10_max  <- terra::buffer(landbiotoop_nabij_max, width = 10)

waterbiotoop2_max <- terra::mask(waterbiotoop1_max, land_buffer10_max)
landbiotoop2_max  <- terra::mask(landbiotoop_nabij_max, water_buffer10_max)

ruim_buffer500_max <- terra::buffer(landbiotoop_ruim_gekoppeld_max, width = 500)

waterbiotoop_finaal_max_ruw <- terra::mask(waterbiotoop2_max, ruim_buffer500_max)
landbiotoop_finaal_max_ruw  <- terra::mask(landbiotoop2_max, ruim_buffer500_max)

if (!all(is.na(suppressWarnings(terra::minmax(waterbiotoop_finaal_max_ruw))))) {
  poly_water_max  <- terra::as.polygons(waterbiotoop_finaal_max_ruw, aggregate = TRUE)
  poly_buffer_max <- terra::buffer(poly_water_max, width = straal_water_m)
  landbiotoop_finaal_max <- terra::mask(landbiotoop_finaal_max_ruw, poly_buffer_max)
} else {
  landbiotoop_finaal_max <- landbiotoop_finaal_max_ruw
}
waterbiotoop_finaal_max <- waterbiotoop_finaal_max_ruw

leefgebied_max <- terra::cover(waterbiotoop_finaal_max, landbiotoop_finaal_max)

# ==============================================================================
# SPOOR 2: WERKELIJKE OPPERVLAKTE (PARALLEL)
# ==============================================================================
water_buffer10_opp <- terra::buffer(!is.na(waterbiotoop1_opp) & waterbiotoop1_opp > 0, width = 10)
land_buffer10_opp  <- terra::buffer(!is.na(landbiotoop_nabij_opp) & landbiotoop_nabij_opp > 0, width = 10)

waterbiotoop2_opp <- terra::mask(waterbiotoop1_opp, land_buffer10_opp)
landbiotoop2_opp  <- terra::mask(landbiotoop_nabij_opp, water_buffer10_opp)

ruim_buffer500_opp <- terra::buffer(!is.na(landbiotoop_ruim_gekoppeld_opp) & landbiotoop_ruim_gekoppeld_opp > 0, width = 500)

waterbiotoop_finaal_opp_ruw <- terra::mask(waterbiotoop2_opp, ruim_buffer500_opp)
landbiotoop_finaal_opp_ruw  <- terra::mask(landbiotoop2_opp, ruim_buffer500_opp)

if (!all(is.na(suppressWarnings(terra::minmax(waterbiotoop_finaal_opp_ruw))))) {
  r_bin_water_opp <- terra::ifel(!is.na(waterbiotoop_finaal_opp_ruw) & waterbiotoop_finaal_opp_ruw > 0, 1, NA)
  poly_water_opp  <- terra::as.polygons(r_bin_water_opp, aggregate = TRUE)
  poly_buffer_opp <- terra::buffer(poly_water_opp, width = straal_water_m)
  landbiotoop_finaal_opp <- terra::mask(landbiotoop_finaal_opp_ruw, poly_buffer_opp)
} else {
  landbiotoop_finaal_opp <- landbiotoop_finaal_opp_ruw
}
waterbiotoop_finaal_opp <- waterbiotoop_finaal_opp_ruw

leefgebied_opp <- terra::cover(waterbiotoop_finaal_opp, landbiotoop_finaal_opp)

# DEFINITIEVE MODELUITGANGEN
final_max <- leefgebied_max
final_opp <- leefgebied_opp

if (!all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  cl_opp <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE)
} else {
  cl_opp <- template_MH * NA
}

if (!all(is.na(suppressWarnings(terra::minmax(final_max))))) {
  cl_max <- terra::patches(final_max, directions = 8, zeroAsNA = TRUE)
} else {
  cl_max <- template_MH * NA
}

rm(water_buffer10_max, land_buffer10_max, waterbiotoop2_max, landbiotoop2_max, ruim_buffer500_max,
   water_buffer10_opp, land_buffer10_opp, waterbiotoop2_opp, landbiotoop2_opp, ruim_buffer500_opp)
gc()

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
potentie_export_rast <- if (exists("final_max") && !is.null(final_max) && !all(is.na(suppressWarnings(terra::minmax(final_max))))) {
  terra::ifel(!is.na(final_max) & final_max > 0, 1, NA)
} else {
  terra::rast(template_MH, vals = NA)
}

# 2. Bepaal Werkelijke Oppervlakte Raster
werkelijk_export_rast <- if (exists("final_opp") && !is.null(final_opp) && !all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  terra::ifel(!is.na(final_opp) & final_opp > 0, 1, NA)
} else {
  terra::rast(template_MH, vals = NA)
}

# 3. Bepaal Analytisch Metacluster ID-raster (EXCLUSIEF OP BASIS VAN WERKELIJKE OPPERVLAKTES)
if (exists("cl_opp") && !is.null(cl_opp) && !all(is.na(suppressWarnings(terra::minmax(cl_opp))))) {
  id_export_rast <- cl_opp
} else if (exists("final_opp") && !is.null(final_opp) && !all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  id_export_rast <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE)
} else {
  id_export_rast <- template_MH * NA
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
  message(paste("    [OK] Geëxporteerd naat scenariomap:", basename(file_path)))
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
