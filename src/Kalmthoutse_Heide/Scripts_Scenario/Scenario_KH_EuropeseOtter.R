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

df <- read_excel(here::here("data/input/Excel_files/Soorten_bwk_afstanden.xlsx"))
soort <- "europeseotter"

# --- DYNAMISCHE SCENARIO PARAMETER CHECK ---
if (exists("SCENARIO_RDS_PAD") && !is.null(SCENARIO_RDS_PAD)) {
  scenario_rds_path <- SCENARIO_RDS_PAD
} else if (exists("params") && !is.null(params$scenario_rds_path)) {
  scenario_rds_path <- params$scenario_rds_path
} else {
  scenario_rds_path <- "data/input/Scenario_rds/KH_Scenario_BWK_2025.rds"
}

p_raw <- gsub("^([.][.]/)+", "", scenario_rds_path)
scenario_path <- here::here(p_raw)

if (!file.exists(scenario_path)) {
  stop(paste("❌ FOUT: Scenario RDS bestand NIET gevonden op:", scenario_path))
}

scen_volledig <- basename(scenario_path)
scenario_naam <- gsub("^KH_Scenario_|^Scenario_|.rds$", "", scen_volledig)

message(paste("Verwerken van soort:", soort, "binnen scenario:", scenario_naam))

resultaat <- df %>%
  filter(tolower(trimws(Soort)) == soort) %>%
  select(Type, MinOpp_ha, AfstandBiotopen_m, Dispersiecap_m)

oppervlakte_ha <- resultaat$MinOpp_ha[1]
afstand_m      <- resultaat$AfstandBiotopen_m[1]
buffer_m       <- resultaat$Dispersiecap_m[1]

rm(df, resultaat)

area_shape  <- vect(here("data/input/Kalmthoutse_Heide.shp"))
master_grid <- rast(here("data/input/Raster_Vlaanderen/Vlaanderen_MasterGrid_10m.tif"))[[1]]

df_namen_sleutel <- read_csv(here("data/input/Excel_files/BWK_Laag_Namen_2025.csv"), show_col_types = FALSE)
gouden_namenlijst <- tolower(trimws(df_namen_sleutel$Laagnaam))

area_shape_proj <- project(area_shape, crs(master_grid))
area_buffer_fix <- buffer(area_shape_proj, width = buffer_m)

message("-> Uitsnede maken en lokale/globale cel-IDs berekenen...")

id_raster_KH <- crop(master_grid, area_buffer_fix, snap = "near")

lokale_coords <- terra::xyFromCell(id_raster_KH, 1:ncell(id_raster_KH))
globale_ids   <- terra::cellFromXY(master_grid, lokale_coords)

vertaal_df <- data.table(
  lokale_id  = 1:ncell(id_raster_KH),
  globale_id = globale_ids
)

studiegebied_globale_ids <- unique(vertaal_df$globale_id[!is.na(vertaal_df$globale_id)])

template_KH <- terra::rasterize(area_buffer_fix, id_raster_KH, field = 1, background = 0)

rm(lokale_coords, globale_ids)
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
  
  tabel_KH  <- tabel_gefilterd[cel_id %in% studiegebied_globale_ids]
  tabel_cel_som <- tabel_KH[, .(Oppervlakte = sum(BWK_FRAC, na.rm = TRUE)), by = .(cel_id)]
  
  if(nrow(tabel_cel_som) > 0) {
    tabel_cel_som[, Oppervlakte := pmin(Oppervlakte, 1)]
    tabel_cel_som[, Match := ifelse(Oppervlakte >= 0.01, 1, 0)]
    
    setnames(tabel_cel_som, "cel_id", "globale_id")
    tabel_finaal_mapping <- merge(tabel_cel_som, vertaal_df, by = "globale_id", all.x = TRUE)
    tabel_finaal_mapping <- tabel_finaal_mapping[!is.na(lokale_id)]
    
    if(nrow(tabel_finaal_mapping) > 0) {
      r_opp_type   <- template_KH * NA
      r_match_type <- template_KH * NA
      
      vec_opp   <- rep(NA_real_, terra::ncell(id_raster_KH))
      vec_match <- rep(NA_real_, terra::ncell(id_raster_KH))
      
      vec_opp[tabel_finaal_mapping$lokale_id] <- tabel_finaal_mapping$Oppervlakte
      
      tabel_matches_clean <- tabel_finaal_mapping[Match == 1]
      if(nrow(tabel_matches_clean) > 0) {
        vec_match[tabel_matches_clean$lokale_id] <- tabel_matches_clean$Match
      }
      
      terra::values(r_opp_type)   <- vec_opp
      terra::values(r_match_type) <- vec_match
      
      lijst_matches[[h_type]]      <- r_match_type
      lijst_oppervlaktes[[h_type]] <- r_opp_type
      
      rm(vec_opp, vec_match)
    }
  }
}

d_max <- lijst_matches[["d"]]; d_opp <- lijst_oppervlaktes[["d"]]
h_max <- lijst_matches[["h"]]; h_opp <- lijst_oppervlaktes[["h"]]
k_max <- lijst_matches[["k"]]; k_opp <- lijst_oppervlaktes[["k"]]
m_max <- lijst_matches[["m"]]; m_opp <- lijst_oppervlaktes[["m"]]
r_max <- lijst_matches[["r"]]; r_opp <- lijst_oppervlaktes[["r"]]
s_max <- lijst_matches[["s"]]; s_opp <- lijst_oppervlaktes[["s"]]
t_max <- lijst_matches[["t"]]; t_opp <- lijst_oppervlaktes[["t"]]
v_max <- lijst_matches[["v"]]; v_opp <- lijst_oppervlaktes[["v"]]

rm(resultaten_gegroepeerd, df_nieuw, lijst_matches, lijst_oppervlaktes, tabel_vlaanderen)
gc()

# ==============================================================================
# OTTER STAP 1: INTEGRAL WATERNETWERK
# ==============================================================================
message("-> Otter Spoor 1: Waterlopen en plassen inlezen en samenvoegen...")

r_huetzon_raw <- terra::rast(here("data/input/Raster_Vlaanderen/vlaanderen_huetzon_10m.tif"))
r_huetzon     <- terra::resample(r_huetzon_raw, id_raster_KH, method = "near") %>% terra::crop(id_raster_KH)

r_waterlopen1 <- r_huetzon == 1

bl_path <- here("data/input/ASCI Files/BlauweLaag_Plassen_06042017_20m_.asc")
if (file.exists(bl_path)) {
  r_bl_raw <- terra::rast(bl_path)
  r_bl     <- terra::resample(r_bl_raw, id_raster_KH, method = "near") %>% terra::crop(id_raster_KH)
  r_plassen <- r_bl == 1
} else {
  r_plassen <- id_raster_KH * NA
}

otter_water_max <- (r_waterlopen1 == 1) | (r_plassen == 1)
otter_water_max <- terra::ifel(otter_water_max == 1, 1, NA)

suppressWarnings(rm(r_huetzon_raw, r_huetzon, r_waterlopen1, r_bl_raw, r_bl, r_plassen))
gc()

# ==============================================================================
# OTTER STAP 2: LANDBIOTOPEN SMELTEN EN ZEVEN OP 1000 HA
# ==============================================================================
message("-> Otter Spoor 2: Alle landbiotopen samenvoegen en opschonen...")

clean_max <- function(r) {
  if (is.null(r)) return(template_KH * NA)
  r[is.na(r) | r == 0] <- NA
  return(r)
}

opp_stack <- c(
  clean_max(d_opp), clean_max(h_opp), clean_max(k_opp), clean_max(m_opp),
  clean_max(r_opp), clean_max(s_opp), clean_max(t_opp), clean_max(v_opp)
)

otter_biotopen_opp <- sum(opp_stack, na.rm = TRUE)
otter_biotopen_opp <- terra::clamp(otter_biotopen_opp, upper = 1.0)
otter_biotopen_max <- terra::ifel(otter_biotopen_opp > 0, 1, NA)

rm(opp_stack)
gc()

# --- SPOOR A: MAXIMALE POTENTIE ---
if (!all(is.na(suppressWarnings(terra::minmax(otter_biotopen_max))))) {
  poly_land_max <- terra::as.polygons(otter_biotopen_max, aggregate = TRUE) %>% terra::disagg()
  poly_buf_max  <- terra::buffer(poly_land_max, width = 50)
  poly_net_max  <- terra::aggregate(poly_buf_max) %>% terra::disagg()
  poly_net_max$Cluster_ID <- 1:nrow(poly_net_max)
  
  poly_land_max$Cluster_ID <- terra::relate(poly_land_max, poly_net_max, "within") %>% 
    apply(1, function(x) { if(any(x)) which(x)[1] else NA })
  
  r_clusters_max <- terra::rasterize(poly_land_max, template_KH, field = "Cluster_ID", touches = FALSE)
  
  freq_max <- terra::freq(r_clusters_max)
  freq_max$Opp_ha <- freq_max$count * 0.01
  
  voldoet_ids_max <- freq_max$value[freq_max$Opp_ha >= 1000]
  
  if (length(voldoet_ids_max) > 0) {
    otter_biotopen1000_max <- terra::mask(otter_biotopen_max, r_clusters_max %in% voldoet_ids_max, maskvalues = FALSE)
  } else {
    otter_biotopen1000_max <- template_KH * NA
  }
  rm(poly_land_max, poly_buf_max, poly_net_max, r_clusters_max, freq_max)
} else {
  otter_biotopen1000_max <- template_KH * NA
}

# --- SPOOR B: WERKELIJKE OPPERVLAKTE ---
r_binair_land_opp <- terra::ifel(!is.na(otter_biotopen_opp) & otter_biotopen_opp > 0, 1, NA)

if (!all(is.na(suppressWarnings(terra::minmax(r_binair_land_opp))))) {
  poly_land_opp <- terra::as.polygons(r_binair_land_opp, aggregate = TRUE) %>% terra::disagg()
  poly_buf_opp  <- terra::buffer(poly_land_opp, width = 50)
  poly_net_opp  <- terra::aggregate(poly_buf_opp) %>% terra::disagg()
  poly_net_opp$Cluster_ID <- 1:nrow(poly_net_opp)
  
  poly_land_opp$Cluster_ID <- terra::relate(poly_land_opp, poly_net_opp, "within") %>% 
    apply(1, function(x) { if(any(x)) which(x)[1] else NA })
  
  r_clusters_opp <- terra::rasterize(poly_land_opp, template_KH, field = "Cluster_ID", touches = FALSE)
  
  zonal_ha_opp <- terra::zonal(otter_biotopen_opp, r_clusters_opp, fun = "sum", na.rm = TRUE)
  colnames(zonal_ha_opp) <- c("Cluster_ID", "Som_Pixels")
  zonal_ha_opp$Opp_ha    <- zonal_ha_opp$Som_Pixels * 0.01
  
  voldoet_ids_opp <- zonal_ha_opp$Cluster_ID[zonal_ha_opp$Opp_ha >= 1000]
  
  if (length(voldoet_ids_opp) > 0) {
    otter_biotopen1000_opp <- terra::mask(otter_biotopen_opp, r_clusters_opp %in% voldoet_ids_opp, maskvalues = FALSE)
  } else {
    otter_biotopen1000_opp <- template_KH * NA
  }
  rm(poly_land_opp, poly_buf_opp, poly_net_opp, r_clusters_opp, zonal_ha_opp)
} else {
  otter_biotopen1000_opp <- template_KH * NA
}

suppressWarnings(rm(otter_biotopen_max, otter_biotopen_opp, r_binair_land_opp))
gc()

# ==============================================================================
# OTTER STAP 3: COMMUTER WATERKOPPELING (< 1000M) EN HABITAT-FINALESERING
# ==============================================================================
message("-> Otter Spoor 3: Wederzijdse afstandsmetingen uitvoeren (max. 1000m)...")

final_max <- template_KH * NA
final_opp <- template_KH * NA

if (exists("otter_water_max") && !all(is.na(suppressWarnings(terra::minmax(otter_water_max))))) {
  
  # --- SPOOR A: MAXIMALE POTENTIE ---
  if (exists("otter_biotopen1000_max") && !all(is.na(suppressWarnings(terra::minmax(otter_biotopen1000_max))))) {
    dist_tot_water_max <- terra::distance(otter_water_max)
    dist_tot_land_max  <- terra::distance(otter_biotopen1000_max)
    
    leefgebied3_max <- otter_biotopen1000_max == 1 & dist_tot_water_max <= 1000
    leefgebied4_max <- otter_water_max == 1 & dist_tot_land_max <= 1000
    
    otter_leefgebied_max <- leefgebied3_max | leefgebied4_max
    final_max <- terra::ifel(otter_leefgebied_max == 1, 1, NA) %>% terra::crop(template_KH)
    
    rm(dist_tot_water_max, dist_tot_land_max, leefgebied3_max, leefgebied4_max, otter_leefgebied_max)
  }
  
  # --- SPOOR B: WERKELIJKE OPPERVLAKTE ---
  r_binair_land1000_opp <- terra::ifel(!is.na(otter_biotopen1000_opp) & otter_biotopen1000_opp > 0, 1, NA)
  
  if (!all(is.na(suppressWarnings(terra::minmax(r_binair_land1000_opp))))) {
    dist_tot_water_opp <- terra::distance(otter_water_max)
    dist_tot_land_opp  <- terra::distance(r_binair_land1000_opp)
    
    leefgebied3_opp_mask <- r_binair_land1000_opp == 1 & dist_tot_water_opp <= 1000
    leefgebied3_opp      <- terra::mask(otter_biotopen1000_opp, leefgebied3_opp_mask)
    
    leefgebied4_opp_mask <- otter_water_max == 1 & dist_tot_land_opp <= 1000
    r_water_opp_full     <- terra::ifel(otter_water_max == 1, 1.0, NA)
    leefgebied4_opp      <- terra::mask(r_water_opp_full, leefgebied4_opp_mask)
    
    l3_clean <- terra::ifel(is.na(leefgebied3_opp), 0, leefgebied3_opp)
    l4_clean <- terra::ifel(is.na(leefgebied4_opp), 0, leefgebied4_opp)
    
    otter_leefgebied_opp_raw <- terra::clamp(l3_clean + l4_clean, upper = 1.0)
    otter_leefgebied_opp     <- terra::ifel(!is.na(otter_leefgebied_opp_raw) & otter_leefgebied_opp_raw > 0, otter_leefgebied_opp_raw, NA)
    
    final_opp <- terra::crop(otter_leefgebied_opp, template_KH)
    
    rm(dist_tot_water_opp, dist_tot_land_opp, leefgebied3_opp_mask, leefgebied3_opp, 
       leefgebied4_opp_mask, r_water_opp_full, l3_clean, l4_clean, otter_leefgebied_opp_raw, otter_leefgebied_opp)
  }
  
  rm(r_binair_land1000_opp)
  gc()
}

if (!all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  cl_opp <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE)
} else {
  cl_opp <- template_KH * NA
}

if (!all(is.na(suppressWarnings(terra::minmax(final_max))))) {
  cl_max <- terra::patches(final_max, directions = 8, zeroAsNA = TRUE)
} else {
  cl_max <- template_KH * NA
}

suppressWarnings(rm(otter_biotopen1000_max, otter_biotopen1000_opp, otter_water_max))
gc()

# ==============================================================================
# SCHONE EXPORT BIOTOOP EN ANALYTISCH ID-RASTER (VOOR SCRIPT 2 / ARPL)
# ==============================================================================
base_dir <- here::here("data/output/Kalmthoutse_Heide/Rasters_Soorten", scenario_naam)

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
  terra::rast(template_KH, vals = NA)
}

# 2. Bepaal Werkelijke Oppervlakte Raster
werkelijk_export_rast <- if (exists("final_opp") && !is.null(final_opp) && !all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  terra::ifel(!is.na(final_opp) & final_opp > 0, 1, NA)
} else {
  terra::rast(template_KH, vals = NA)
}

# 3. Bepaal Analytisch Metacluster ID-raster (EXCLUSIEF OP BASIS VAN WERKELIJKE OPPERVLAKTES)
if (exists("cl_opp") && !is.null(cl_opp) && !all(is.na(suppressWarnings(terra::minmax(cl_opp))))) {
  id_export_rast <- cl_opp
} else if (exists("final_opp") && !is.null(final_opp) && !all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  id_export_rast <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE)
} else {
  id_export_rast <- terra::rast(template_KH, vals = NA)
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
