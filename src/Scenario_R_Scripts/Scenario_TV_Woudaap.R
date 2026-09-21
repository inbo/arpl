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

soort <- "woudaap"

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
buffer_m       <- 50000

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

openwater_max        <- lijst_matches[["openwater"]]
openwater_opp        <- lijst_oppervlaktes[["openwater"]]
voortplanting_bwk_max <- lijst_matches[["voortplanting_bwk"]]
voortplanting_bwk_opp <- lijst_oppervlaktes[["voortplanting_bwk"]]
not_populier_max     <- lijst_matches[["not_populier"]]
not_populier_opp     <- lijst_oppervlaktes[["not_populier"]]

rm(vertaal_df, lijst_matches, lijst_oppervlaktes, tabel_vlaanderen, resultaten_gegroepeerd, df_nieuw)
gc()

# ==============================================================================
# STAP 1: VOORTPLANTINGSHABITAT (>= 10 HA - PARALLEL)
# ==============================================================================
message("-> Spoor 1: Voortplantingsbiotopen zuiveren en clusteren...")

# --- SPOOR A: MAX ---
r_woudaap_voortplanting2_max <- id_raster_TV * NA

if (exists("voortplanting_bwk_max") && !all(is.na(suppressWarnings(terra::minmax(voortplanting_bwk_max))))) {
  r_moeras_schoon_max <- voortplanting_bwk_max
  
  if (exists("not_populier_max") && !all(is.na(suppressWarnings(terra::minmax(not_populier_max))))) {
    r_moeras_schoon_max <- terra::mask(voortplanting_bwk_max, not_populier_max, inverse = TRUE)
  }
  
  if (!all(is.na(suppressWarnings(terra::minmax(r_moeras_schoon_max))))) {
    cl_voortp_max <- cluster_filter_compleet(
      masker     = r_moeras_schoon_max,
      opp_laag   = voortplanting_bwk_opp,
      drempel_m2 = 10 * 10000, 
      dist_m     = 50,          
      werkelijk  = FALSE
    )
    r_woudaap_voortplanting2_max <- cl_voortp_max$raster
  }
}

# --- SPOOR B: OPP ---
r_woudaap_voortplanting2_opp <- id_raster_TV * NA

if (exists("voortplanting_bwk_opp") && !all(is.na(suppressWarnings(terra::minmax(voortplanting_bwk_opp))))) {
  r_moeras_schoon_opp <- voortplanting_bwk_opp
  
  if (exists("not_populier_opp") && !all(is.na(suppressWarnings(terra::minmax(not_populier_opp))))) {
    r_moeras_schoon_opp <- terra::mask(voortplanting_bwk_opp, not_populier_opp, inverse = TRUE)
  }
  
  r_binair_moeras_opp <- terra::ifel(!is.na(r_moeras_schoon_opp) & r_moeras_schoon_opp > 0, 1, NA)
  
  if (!all(is.na(suppressWarnings(terra::minmax(r_binair_moeras_opp))))) {
    cl_voortp_opp <- cluster_filter_compleet(
      masker     = r_binair_moeras_opp,
      opp_laag   = r_moeras_schoon_opp,
      drempel_m2 = 10 * 10000, 
      dist_m     = 50,          
      werkelijk  = TRUE
    )
    r_woudaap_voortplanting2_opp <- cl_voortp_opp$raster
  }
}

suppressWarnings(rm(cl_voortp_max, cl_voortp_opp, r_moeras_schoon_max, r_moeras_schoon_opp, r_binair_moeras_opp))
gc()

# ==============================================================================
# STAP 2: OPEN WATER TOETS (>= 1 HA) EN KOPPELING (<100M - PARALLEL)
# ==============================================================================
message("-> Spoor 2: Open waterpartijen filteren op min. 1 ha...")

# --- SPOOR A: MAX ---
r_woudaap_leefgebied1_max <- id_raster_TV * NA

if (exists("openwater_max") && !all(is.na(suppressWarnings(terra::minmax(openwater_max)))) && 
    !all(is.na(suppressWarnings(terra::minmax(r_woudaap_voortplanting2_max))))) {
  
  cl_water_max <- cluster_filter_compleet(
    masker     = openwater_max,
    opp_laag   = openwater_opp,
    drempel_m2 = 1 * 10000, 
    dist_m     = 10,        
    werkelijk  = FALSE
  )
  r_water_1ha_max <- cl_water_max$raster
  
  if (!all(is.na(suppressWarnings(terra::minmax(r_water_1ha_max))))) {
    dist_tot_water_max <- terra::distance(r_water_1ha_max)
    r_woudaap_leefgebied1_max <- terra::ifel(dist_tot_water_max <= 100 & r_woudaap_voortplanting2_max == 1, 1, NA)
  }
}

# --- SPOOR B: OPP ---
r_woudaap_leefgebied1_opp <- id_raster_TV * NA

if (exists("openwater_opp") && !all(is.na(suppressWarnings(terra::minmax(openwater_opp)))) && 
    !all(is.na(suppressWarnings(terra::minmax(r_woudaap_voortplanting2_opp))))) {
  
  r_binair_water_opp <- terra::ifel(!is.na(openwater_opp) & openwater_opp > 0, 1, NA)
  cl_water_opp <- cluster_filter_compleet(
    masker     = r_binair_water_opp,
    opp_laag   = openwater_opp,
    drempel_m2 = 1 * 10000, 
    dist_m     = 10,        
    werkelijk  = TRUE
  )
  r_water_1ha_opp <- cl_water_opp$raster
  
  if (!all(is.na(suppressWarnings(terra::minmax(r_water_1ha_opp))))) {
    dist_tot_water_opp <- terra::distance(!is.na(r_water_1ha_opp) & r_water_1ha_opp > 0)
    r_woudaap_leefgebied1_opp <- terra::mask(r_woudaap_voortplanting2_opp, terra::ifel(dist_tot_water_opp <= 100, 1, NA))
  }
}

suppressWarnings(rm(cl_water_max, cl_water_opp, r_water_1ha_max, r_water_1ha_opp, r_binair_water_opp, dist_tot_water_max, dist_tot_water_opp))
gc()

# ==============================================================================
# STAP 3: FINALE NETWERKCLUSTERING (FUZZY 100M, MIN 10 HA - PARALLEL)
# ==============================================================================
message("-> 3. Finale netwerkclustering uitvoeren (fuzzy 100m, min. 10 ha)...")

woudaap_finaal_max <- template_TV * NA
woudaap_finaal_opp <- template_TV * NA

r_tot_opp <- terra::cover(voortplanting_bwk_opp, openwater_opp)

# --- SPOOR A: MAX ---
if (!all(is.na(suppressWarnings(terra::minmax(r_woudaap_leefgebied1_max))))) {
  cl_woudaap_finaal_max <- cluster_filter_compleet(
    masker     = r_woudaap_leefgebied1_max,
    opp_laag   = r_tot_opp,
    drempel_m2 = 10 * 10000, 
    dist_m     = 100,        
    werkelijk  = FALSE
  )
  woudaap_finaal_max <- cl_woudaap_finaal_max$raster %>% terra::crop(template_TV)
}

# --- SPOOR B: OPP ---
if (!all(is.na(suppressWarnings(terra::minmax(r_woudaap_leefgebied1_opp))))) {
  r_binair_leef1_opp <- terra::ifel(!is.na(r_woudaap_leefgebied1_opp) & r_woudaap_leefgebied1_opp > 0, 1, NA)
  cl_woudaap_finaal_opp <- cluster_filter_compleet(
    masker     = r_binair_leef1_opp,
    opp_laag   = r_woudaap_leefgebied1_opp,
    drempel_m2 = 10 * 10000, 
    dist_m     = 100,        
    werkelijk  = TRUE
  )
  woudaap_finaal_opp <- cl_woudaap_finaal_opp$raster %>% terra::crop(template_TV)
  rm(r_binair_leef1_opp)
}

# DEFINITIEVE MODELUITGANGEN
final_max <- woudaap_finaal_max
final_opp <- woudaap_finaal_opp

if (!all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  cl_opp <- cl_woudaap_finaal_opp$clusters %>% terra::crop(template_TV)
} else {
  cl_opp <- template_TV * NA
}

if (!all(is.na(suppressWarnings(terra::minmax(final_max))))) {
  cl_max <- cl_woudaap_finaal_max$clusters %>% terra::crop(template_TV)
} else {
  cl_max <- template_TV * NA
}

suppressWarnings(rm(r_tot_opp, cl_woudaap_finaal_max, cl_woudaap_finaal_opp, r_woudaap_leefgebied1_max, r_woudaap_leefgebied1_opp, r_woudaap_voortplanting2_max, r_woudaap_voortplanting2_opp, openwater_max, openwater_opp, voortplanting_bwk_max, voortplanting_bwk_opp, not_populier_max, not_populier_opp))
gc()

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
  id_export_rast <- template_TV * NA
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
