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

cluster_filter_compleet <- function(masker, opp_laag = NULL, drempel_m2, dist_m, werkelijk = FALSE) {
  if (is.null(masker) || terra::global(is.na(masker), "sum")[[1]] == terra::ncell(masker)) {
    return(list(raster = masker * NA, clusters = masker * NA))
  }
  
  r_binair <- terra::ifel(!is.na(masker) & masker > 0, 1, NA)
  if (terra::global(!is.na(r_binair), "sum")[[1]] == 0) {
    return(list(raster = masker * NA, clusters = masker * NA))
  }
  
  p_ruw <- terra::as.polygons(r_binair, aggregate = TRUE)
  if (is.null(p_ruw) || nrow(p_ruw) == 0) {
    return(list(raster = masker * NA, clusters = masker * NA))
  }
  
  p_ruw <- terra::disagg(p_ruw)
  
  if (dist_m > 0) {
    p_buf <- terra::buffer(p_ruw, width = dist_m / 2)
    p_netwerk <- terra::aggregate(p_buf) %>% terra::disagg()
    p_netwerk$Netwerk_ID <- 1:nrow(p_netwerk)
    
    r_net_zones <- terra::rasterize(p_netwerk, r_binair, field = "Netwerk_ID")
    cl_biotoop_only <- terra::mask(r_net_zones, r_binair)
    
    rm(p_buf, p_netwerk, r_net_zones)
  } else {
    cl_biotoop_only <- terra::patches(r_binair, directions = 4, zeroAsNA = TRUE)
  }
  
  if (werkelijk && !is.null(opp_laag)) {
    stats_df <- terra::zonal(opp_laag, cl_biotoop_only, fun = "sum", na.rm = TRUE, as.raster = FALSE)
    colnames(stats_df) <- c("ID", "Waarde")
    stats_df$Area_m2 <- stats_df$Waarde * 100 
  } else {
    f <- terra::freq(cl_biotoop_only)
    stats_df <- data.frame(ID = f$value, Waarde = f$count)
    stats_df$Area_m2 <- stats_df$Waarde * 100 
  }
  
  stats_df <- stats_df[!is.na(stats_df$ID), ]
  if (nrow(stats_df) == 0) return(list(raster = masker * NA, clusters = masker * NA))
  
  voldoet_ids <- stats_df$ID[stats_df$Area_m2 >= drempel_m2]
  if (length(voldoet_ids) == 0) return(list(raster = masker * NA, clusters = masker * NA))
  
  masker_binair <- cl_biotoop_only %in% voldoet_ids
  final_network_mask <- terra::ifel(masker_binair == 1, 1, NA)
  
  r_finaal  <- terra::mask(masker, final_network_mask)
  cl_finaal <- terra::mask(cl_biotoop_only, r_finaal) 
  
  rm(p_ruw, masker_binair, final_network_mask, r_binair)
  gc()
  
  return(list(raster = r_finaal, clusters = cl_finaal))
}

terraOptions(
  memfrac = 0.8,
  tempdir = tempdir(),
  verbose = FALSE
)

soort <- "kwartelkoning"

# --- DYNAMISCHE SCENARIO PARAMETER CHECK ---
if (exists("SCENARIO_RDS_PAD") && !is.null(SCENARIO_RDS_PAD)) {
  scenario_rds_path <- SCENARIO_RDS_PAD
} else if (exists("params") && !is.null(params$scenario_rds_path)) {
  scenario_rds_path <- params$scenario_rds_path
} else {
  scenario_rds_path <- "data/input/Scenario_rds/VS_Scenario_BWK_2025.rds"
}

p_raw <- gsub("^([.][.]/)+", "", scenario_rds_path)
scenario_path <- here::here(p_raw)

if (!file.exists(scenario_path)) {
  stop(paste("❌ FOUT: Scenario RDS bestand NIET gevonden op:", scenario_path))
}

scen_volledig <- basename(scenario_path)
scenario_naam <- gsub("^VS_Scenario_|^Scenario_|.rds$", "", scen_volledig)

message(paste("Verwerken van soort:", soort, "binnen scenario:", scenario_naam))

df <- read_excel(here::here("data/input/Excel_files/Soorten_bwk_afstanden.xlsx"))
resultaat <- df %>%
  filter(tolower(trimws(Soort)) == soort) %>%
  select(Type, MinOpp_ha, AfstandBiotopen_m, Dispersiecap_m)

oppervlakte_ha <- resultaat$MinOpp_ha[1]
afstand_m      <- resultaat$AfstandBiotopen_m[1]
buffer_m       <- resultaat$Dispersiecap_m[1]

rm(df, resultaat)

area_shape  <- vect(here("data/input/Voerstreek.shp"))
master_grid <- rast(here("data/input/Raster_Vlaanderen/Vlaanderen_MasterGrid_10m.tif"))[[1]]

df_namen_sleutel <- read_csv(here("data/input/Excel_files/BWK_Laag_Namen_2025.csv"), show_col_types = FALSE)
gouden_namenlijst <- tolower(trimws(df_namen_sleutel$Laagnaam))

area_shape_proj <- project(area_shape, crs(master_grid))
area_buffer_fix <- buffer(area_shape_proj, width = buffer_m)

message("-> Vertaalraster voor globale/lokale cellen opbouwen via snelle MASK methode...")
id_raster_VS <- crop(master_grid, area_buffer_fix, snap = "near")

globale_id_raster <- master_grid
globale_id_raster <- terra::init(globale_id_raster, fun = "cell")

id_raster_VS_globale_values <- crop(globale_id_raster, area_buffer_fix, snap = "near")
id_raster_VS_masked <- mask(id_raster_VS_globale_values, area_buffer_fix)

message("-> Vertaaltabel bliksemsnel opbouwen via C++ dataframe extractie...")

df_extractie <- as.data.frame(id_raster_VS_masked, cells = TRUE)
vertaal_df <- as.data.table(df_extractie)
setnames(vertaal_df, c(1, 2), c("lokale_id", "globale_id"))

vertaal_df <- vertaal_df[!is.na(globale_id)]
studiegebied_globale_ids <- unique(vertaal_df$globale_id)

values(id_raster_VS) <- NA
template_VS <- terra::rasterize(area_buffer_fix, id_raster_VS, field = 1, background = 0)

rm(globale_id_raster, id_raster_VS_globale_values, id_raster_VS_masked, df_extractie)
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
  
  tabel_VS <- tabel_gefilterd[cel_id %in% studiegebied_globale_ids]
  tabel_VS_unique <- unique(tabel_VS, by = c("cel_id", "CODE"))
  
  tabel_cel_som <- tabel_VS_unique[, .(Oppervlakte = pmin(sum(BWK_FRAC, na.rm = TRUE), 1.0)), by = .(cel_id)]
  
  r_match_type <- id_raster_VS * NA
  r_opp_type   <- id_raster_VS * NA
  
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

grasland_max     <- lijst_matches[["grasland"]]
grasland_opp     <- lijst_oppervlaktes[["grasland"]]
moeras_max       <- lijst_matches[["moeras"]]
moeras_opp       <- lijst_oppervlaktes[["moeras"]]
zeergeschikt_max <- lijst_matches[["zeergeschikt1"]]
zeergeschikt_opp <- lijst_oppervlaktes[["zeergeschikt1"]]

# 1. HOOGGROEN SYNCHRONISEREN AND HARD WEGSCHRIJVEN NAAR SCHIJF
hooggroen_raw         <- rast(here("data/input/ASCI Files/Groenkaart_2021.tif"))
area_buffer_groen_crs <- project(area_buffer_fix, crs(hooggroen_raw))
hooggroen_local_raw   <- crop(hooggroen_raw, area_buffer_groen_crs)

hooggroen_binair  <- terra::classify(hooggroen_local_raw, matrix(c(0.5, 1.5, 1), ncol = 3, byrow = TRUE), others = 0)
hooggroen_10m_raw <- terra::aggregate(hooggroen_binair, fact = 10, fun = "max")

tf_groen_sync <- tempfile(fileext = ".tif")
hooggroen_sync <- project(hooggroen_10m_raw, crs(template_VS)) %>% 
  terra::resample(template_VS, method = "near")
hooggroen_sync <- terra::writeRaster(hooggroen_sync, filename = tf_groen_sync, overwrite = TRUE, datatype = "INT1U", NAflag = 255)

rm(hooggroen_raw, hooggroen_local_raw, hooggroen_binair, hooggroen_10m_raw)
gc()

# 2. TERRILS ISOLEREN (kg%)
tabel_terril <- tabel_vlaanderen[CODE %like% "^kg"]
tabel_terril_VS <- tabel_terril[cel_id %in% studiegebied_globale_ids]

tabel_terril_som <- tabel_terril_VS[, .(Oppervlakte = sum(BWK_FRAC, na.rm = TRUE)), by = .(cel_id)]

r_terril <- template_VS * NA

if(nrow(tabel_terril_som) > 0) {
  setnames(tabel_terril_som, "cel_id", "globale_id")
  tabel_terril_mapping <- merge(tabel_terril_som, vertaal_df, by = "globale_id", all.x = TRUE)
  tabel_terril_mapping <- tabel_terril_mapping[!is.na(lokale_id) & Oppervlakte >= 0.01]
  
  if(nrow(tabel_terril_mapping) > 0) {
    r_terril[tabel_terril_mapping$lokale_id] <- 1
  }
}

# ==============================================================================
# STAP 1: PARALLELLE START
# ==============================================================================
biotoop_initieel_bin <- (grasland_max == 1) | (moeras_max == 1)
biotoop_initieel_bin <- terra::ifel(biotoop_initieel_bin == 1, 1, NA)

biotoop_initieel_frac <- terra::cover(grasland_opp, moeras_opp)

ongeschikt_masker <- (hooggroen_sync == 1) | (!is.na(r_terril) & r_terril == 1)
ongeschikt_masker_bin <- terra::ifel(ongeschikt_masker, 1, NA)

biotoop_fase1_bin  <- terra::mask(biotoop_initieel_bin, ongeschikt_masker_bin, maskvalues = 1)
biotoop_fase1_frac <- terra::mask(biotoop_initieel_frac, ongeschikt_masker_bin, maskvalues = 1)

rm(vertaal_df, tabel_terril, tabel_terril_VS, tabel_terril_som)
gc()

# --- STAP 2: DYNAMISCHE VALLEIGRONDEN & HYDROGRAFIE ---
raw_profiel      <- rast(here("data/input/Raster_Vlaanderen/Vlaanderen_profiel_10m.tif"))
raw_drainage     <- rast(here("data/input/Raster_Vlaanderen/vlaanderen_drainage_n2khab_10m.tif"))
raw_overstroming <- rast(here("data/input/Raster_Vlaanderen/vlaanderen_ovstrg_10m.tif")) 

profiel_sync      <- terra::crop(raw_profiel, template_VS) %>% terra::resample(template_VS, method = "near")
drainage_sync     <- terra::crop(raw_drainage, template_VS) %>% terra::resample(template_VS, method = "near")
overstroming_sync <- terra::crop(raw_overstroming, template_VS) %>% terra::resample(template_VS, method = "near")

gewenste_profielen <- c("p", "p+x")
legenda_profiel    <- terra::cats(profiel_sync)[[1]]
profiel_nummers    <- legenda_profiel$value[tolower(trimws(legenda_profiel$Label)) %in% tolower(gewenste_profielen)]

gewenste_drainages <- c("d", "e", "f", "g", "h", "i", "e-f", "h-i", "e-i") 
legenda_drainage    <- terra::cats(drainage_sync)[[1]]
drainage_nummers    <- legenda_drainage$value[tolower(trimws(legenda_drainage$Label)) %in% tolower(gewenste_drainages)]

overstroming_nummers <- c(1, 2)

mask_profiel      <- profiel_sync %in% profiel_nummers
mask_drainage     <- drainage_sync %in% drainage_nummers
mask_overstroming <- overstroming_sync %in% overstroming_nummers

vallei_geschikt        <- (mask_profiel == 1) | (mask_drainage == 1) | (mask_overstroming == 1)
vallei_geschikt_binair <- terra::ifel(vallei_geschikt == 1, 1, NA)

biotoop_fase2_bin  <- terra::mask(biotoop_fase1_bin, vallei_geschikt_binair)
biotoop_fase2_frac <- terra::mask(biotoop_fase1_frac, vallei_geschikt_binair)

rm(raw_profiel, raw_drainage, raw_overstroming, profiel_sync, drainage_sync, overstroming_sync)

# --- STAP 3: MAJORITY FILTER ---
f_matrix <- terra::focalMat(template_VS, d = 20, type = "circle")
f_matrix[f_matrix > 0] <- 1
drempel <- ceiling(sum(f_matrix) / 2)

binaire_aanwezigheid <- terra::ifel(!is.na(biotoop_fase2_bin) & biotoop_fase2_bin == 1, 1, 0)

temp_focal_file <- tempfile(fileext = ".tif")
buurt_som <- terra::focal(
  binaire_aanwezigheid, 
  w = f_matrix, 
  fun = "sum", 
  na.rm = TRUE,
  filename = temp_focal_file,
  overwrite = TRUE
)

vorm_fase3 <- terra::ifel(buurt_som >= drempel, 1, NA)

biotoop_fase3_bin  <- terra::mask(biotoop_fase2_bin, vorm_fase3)
biotoop_fase3_frac <- terra::mask(biotoop_fase2_frac, vorm_fase3)

rm(binaire_aanwezigheid, buurt_som, vorm_fase3)
if(file.exists(temp_focal_file)) unlink(temp_focal_file)
gc()

# --- STAP 4: MINIMAAL 3 HA FILTER ---
fase4_lijst_bin <- cluster_filter_compleet(
  masker     = biotoop_fase3_bin, 
  opp_laag   = NULL, 
  drempel_m2 = 30000, 
  dist_m     = 50, 
  werkelijk  = FALSE
)
biotoop_fase4_bin <- fase4_lijst_bin$raster

r_binair_fase3_frac <- terra::ifel(!is.na(biotoop_fase3_frac) & biotoop_fase3_frac > 0, 1, NA)

fase4_lijst_frac <- cluster_filter_compleet(
  masker     = r_binair_fase3_frac, 
  opp_laag   = biotoop_fase3_frac, 
  drempel_m2 = 30000, 
  dist_m     = 50, 
  werkelijk  = TRUE
)

masker_fase4_valid <- terra::ifel(!is.na(fase4_lijst_frac$raster), 1, NA)
biotoop_fase4_frac <- terra::mask(biotoop_fase3_frac, masker_fase4_valid)

rm(r_binair_fase3_frac, masker_fase4_valid, fase4_lijst_bin, fase4_lijst_frac)
gc()

# --- STAP 5: HOOGGROEN FILTER ---
r_groen_bin <- terra::ifel(!is.na(hooggroen_sync) & hooggroen_sync > 0, 1, NA)

w_matrix <- matrix(1, nrow = 5, ncol = 5)
temp_groen_file <- tempfile(fileext = ".tif")

groen_som <- terra::focal(
  r_groen_bin, 
  w = w_matrix, 
  fun = "sum", 
  na.rm = TRUE,
  filename = temp_groen_file,
  overwrite = TRUE
)

r_groen_20are <- terra::ifel(groen_som >= 12 & r_groen_bin == 1, 1, NA)

buf_groen <- terra::buffer(r_groen_20are, width = 100)
masker_groen_100m_cluster <- terra::ifel(!is.na(buf_groen) & buf_groen > 0, 1, NA)

biotoop_fase5_bin  <- terra::mask(biotoop_fase4_bin, masker_groen_100m_cluster, inverse = TRUE)
biotoop_fase5_frac <- terra::mask(biotoop_fase4_frac, masker_groen_100m_cluster, inverse = TRUE)

rm(r_groen_bin, groen_som, r_groen_20are, buf_groen, masker_groen_100m_cluster)
if(file.exists(temp_groen_file)) unlink(temp_groen_file)
if(file.exists(tf_groen_sync)) unlink(tf_groen_sync)
gc()

# --- STAP 6: FUNCTIONEEL LEEFGEBIED 10 HA ---
fase6_lijst_bin <- cluster_filter_compleet(
  masker     = biotoop_fase5_bin, 
  opp_laag   = NULL, 
  drempel_m2 = 100000, 
  dist_m     = 50, 
  werkelijk  = FALSE
)
biotoop_fase6_bin <- fase6_lijst_bin$raster
cl_id_fase6_bin   <- fase6_lijst_bin$clusters

r_binair_fase5_frac <- terra::ifel(!is.na(biotoop_fase5_frac) & biotoop_fase5_frac > 0, 1, NA)

fase6_lijst_frac <- cluster_filter_compleet(
  masker     = r_binair_fase5_frac, 
  opp_laag   = biotoop_fase5_frac, 
  drempel_m2 = 100000, 
  dist_m     = 50, 
  werkelijk  = TRUE
)
biotoop_fase6_frac <- fase6_lijst_frac$raster
cl_id_fase6_frac   <- fase6_lijst_frac$clusters

rm(r_binair_fase5_frac, fase6_lijst_bin, fase6_lijst_frac)
gc()

# --- STAP 7: GRASLANDFRACTIE-CONTROLE ---
cl_id_fase6_bin_gefilterd <- terra::mask(cl_id_fase6_bin, biotoop_fase6_bin)

tot_biotoop_bin  <- terra::zonal(biotoop_fase5_bin, cl_id_fase6_bin_gefilterd, fun = "sum", na.rm = TRUE)
tot_grasland_bin <- terra::zonal(grasland_max, cl_id_fase6_bin_gefilterd, fun = "sum", na.rm = TRUE)
colnames(tot_biotoop_bin)  <- c("ID", "Tot_Biotoop")
colnames(tot_grasland_bin) <- c("ID", "Tot_Grasland")

df_bin <- merge(tot_biotoop_bin, tot_grasland_bin, by = "ID", all.x = TRUE)
df_bin$Tot_Grasland[is.na(df_bin$Tot_Grasland)] <- 0
df_bin$Grasland_Pct <- (df_bin$Tot_Grasland / df_bin$Tot_Biotoop) * 100

ids_ok_bin <- df_bin$ID[df_bin$Grasland_Pct >= 75]

if(length(ids_ok_bin) > 0) {
  biotoop_fase7_bin <- terra::ifel(cl_id_fase6_bin_gefilterd %in% ids_ok_bin, biotoop_fase6_bin, NA)
  cl_id_fase7_bin   <- terra::ifel(cl_id_fase6_bin_gefilterd %in% ids_ok_bin, cl_id_fase6_bin_gefilterd, NA)
} else {
  biotoop_fase7_bin <- template_VS * NA
  cl_id_fase7_bin   <- template_VS * NA
}

cl_id_fase6_frac_gefilterd <- terra::mask(cl_id_fase6_frac, biotoop_fase6_frac)

tot_biotoop_frac  <- terra::zonal(biotoop_fase5_frac, cl_id_fase6_frac_gefilterd, fun = "sum", na.rm = TRUE)
tot_grasland_frac <- terra::zonal(grasland_opp, cl_id_fase6_frac_gefilterd, fun = "sum", na.rm = TRUE)
colnames(tot_biotoop_frac)  <- c("ID", "Tot_Biotoop")
colnames(tot_grasland_frac) <- c("ID", "Tot_Grasland")

df_frac <- merge(tot_biotoop_frac, tot_grasland_frac, by = "ID", all.x = TRUE)
df_frac$Tot_Grasland[is.na(df_frac$Tot_Grasland)] <- 0
df_frac$Grasland_Pct <- (df_frac$Tot_Grasland / df_frac$Tot_Biotoop) * 100

ids_ok_frac <- df_frac$ID[df_frac$Grasland_Pct >= 75]

if(length(ids_ok_frac) > 0) {
  biotoop_fase7_frac <- terra::ifel(cl_id_fase6_frac_gefilterd %in% ids_ok_frac, biotoop_fase6_frac, NA)
  cl_id_fase7_frac   <- terra::ifel(cl_id_fase6_frac_gefilterd %in% ids_ok_frac, cl_id_fase6_frac_gefilterd, NA)
} else {
  biotoop_fase7_frac <- template_VS * NA
  cl_id_fase7_frac   <- template_VS * NA
}

# DEFINITIEVE MODELUITGANGEN
final_max <- biotoop_fase7_bin
final_opp <- biotoop_fase7_frac

if (!all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  cl_opp <- cl_id_fase7_frac
} else {
  cl_opp <- template_VS * NA
}

if (!all(is.na(suppressWarnings(terra::minmax(final_max))))) {
  cl_max <- cl_id_fase7_bin
} else {
  cl_max <- template_VS * NA
}

rm(df_bin, df_frac, tot_biotoop_bin, tot_grasland_bin, tot_biotoop_frac, tot_grasland_frac, tabel_vlaanderen, resultaten_gegroepeerd)
gc()

# ==============================================================================
# SCHONE EXPORT BIOTOOP EN ANALYTISCH ID-RASTER (VOOR SCRIPT 2 / ARPL)
# ==============================================================================
base_dir <- here::here("data/output/Voerstreek/Rasters_Soorten", scenario_naam)

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
  terra::rast(template_VS, vals = NA)
}

# 2. Bepaal Werkelijke Oppervlakte Raster
werkelijk_export_rast <- if (exists("final_opp") && !is.null(final_opp) && !all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  terra::ifel(!is.na(final_opp) & final_opp > 0, 1, NA)
} else {
  terra::rast(template_VS, vals = NA)
}

# 3. Bepaal Analytisch Metacluster ID-raster (EXCLUSIEF OP BASIS VAN WERKELIJKE OPPERVLAKTES)
if (exists("cl_opp") && !is.null(cl_opp) && !all(is.na(suppressWarnings(terra::minmax(cl_opp))))) {
  id_export_rast <- cl_opp
} else if (exists("final_opp") && !is.null(final_opp) && !all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  id_export_rast <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE)
} else {
  id_export_rast <- terra::rast(template_VS, vals = NA)
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
