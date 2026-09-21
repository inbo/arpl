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

soort <- "vliegendhert"

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
buffer_m       <- resultaat$Dispersiecap_m[1]

rm(df, resultaat)

area_shape  <- vect(here("data/input/Turnhouts_Vennegebied.shp"))
master_grid <- rast(here("data/input/Raster_Vlaanderen/Vlaanderen_MasterGrid_10m.tif"))[[1]]

df_namen_sleutel <- read_csv(here("data/input/Excel_files/BWK_Laag_Namen_2025.csv"), show_col_types = FALSE)
gouden_namenlijst <- tolower(trimws(df_namen_sleutel$Laagnaam))

area_shape_proj <- project(area_shape, crs(master_grid))
area_buffer_fix <- buffer(area_shape_proj, width = buffer_m)

message("-> Vertaalraster voor globale/lokale cellen opbouwen via mask methode...")
id_raster_VS <- crop(master_grid, area_buffer_fix, snap = "near")

globale_id_raster <- master_grid
globale_id_raster <- terra::init(globale_id_raster, fun = "cell")

id_raster_VS_globale_values <- crop(globale_id_raster, area_buffer_fix, snap = "near")
id_raster_VS_masked <- mask(id_raster_VS_globale_values, area_buffer_fix)

message("-> Vertaaltabel opbouwen via C++ dataframe extractie...")
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
  
  tabel_VS      <- tabel_gefilterd[cel_id %in% studiegebied_globale_ids]
  tabel_cel_som <- tabel_VS[, .(Oppervlakte = sum(BWK_FRAC, na.rm = TRUE)), by = .(cel_id)]
  
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

naaldbos_max <- lijst_matches[["naaldbos"]]
naaldbos_opp <- lijst_oppervlaktes[["naaldbos"]]

rm(vertaal_df, lijst_matches, lijst_oppervlaktes, tabel_vlaanderen, resultaten_gegroepeerd, df_nieuw)
gc()

message("-> Externe GIS-lagen inlezen voor het vliegend hert...")

load_and_sync <- function(pad) {
  if (file.exists(pad)) {
    r_raw <- terra::rast(pad)
    r_sync <- terra::resample(r_raw, id_raster_VS, method = "near") %>% terra::crop(id_raster_VS)
    return(r_sync)
  } else {
    message("⚠️ Waarschuwing: Bestand niet gevonden op pad: ", pad)
    return(id_raster_VS * NA)
  }
}

r_groenkaart <- load_and_sync(here("data/input/ASCI Files/Groenkaart_2021.tif"))

r_groen_hoog_max <- r_groenkaart == 1
r_groen_open_max <- r_groenkaart %in% c(2, 3, 4)

r_bosleeftijd_raw <- load_and_sync(here("data/input/Raster_Vlaanderen/vlaanderen_bosleeftijd_10m.tif"))
r_oud_bos_max      <- r_bosleeftijd_raw %in% c(1, 2)

r_helling_max <- load_and_sync(here("data/input/Raster_Vlaanderen/vliegend_hert_4gr_DTMVLII_10m_comp.tif")) == 1

suppressWarnings(rm(r_groenkaart, r_bosleeftijd_raw))
gc()

message("-> Geschikt loofbos afbakenen met 40% naaldhout-tolerantie...")

r_naaldbos_opp_clean <- terra::ifel(is.na(naaldbos_opp), 0, naaldbos_opp)
r_naaldbos_dominant  <- r_naaldbos_opp_clean > 0.40

r_bossen1_max        <- r_groen_hoog_max == 1 & !r_naaldbos_dominant
r_bossen1_max        <- terra::ifel(r_bossen1_max == 1, 1, NA)

if (!all(is.na(suppressWarnings(terra::minmax(r_oud_bos_max))))) {
  dist_oud_bos_max  <- terra::distance(r_oud_bos_max)
  r_bossen2_max     <- r_bossen1_max == 1 & dist_oud_bos_max <= 1000
} else {
  r_bossen2_max     <- r_bossen1_max
}
r_bossen2_max <- terra::ifel(r_bossen2_max == 1, 1, NA)

if (!all(is.na(suppressWarnings(terra::minmax(r_helling_max))))) {
  dist_helling_max       <- terra::distance(r_helling_max)
  r_bossen_gefilterd_max <- r_bossen2_max == 1 & dist_helling_max <= 1000
} else {
  r_bossen_gefilterd_max <- r_bossen2_max
}
r_bossen_gefilterd_max <- terra::ifel(r_bossen_gefilterd_max == 1, 1, NA)

cl_bos_10ha_max <- cluster_filter_compleet(
  masker     = r_bossen_gefilterd_max,
  opp_laag   = r_bossen_gefilterd_max,
  drempel_m2 = 10 * 10000, 
  dist_m     = 50,
  werkelijk  = FALSE
)
vh_bossen3_max <- cl_bos_10ha_max$raster

r_loof_frac    <- terra::clamp(1.0 - r_naaldbos_opp_clean, lower = 0.0, upper = 1.0)
vh_bossen3_opp <- terra::mask(r_loof_frac, vh_bossen3_max)

suppressWarnings(rm(r_naaldbos_opp_clean, r_naaldbos_dominant, r_bossen1_max, 
                    r_bossen2_max, r_bossen_gefilterd_max, cl_bos_10ha_max))
gc()

message("-> Open groengebieden clusteren (50m fuzzy) en filteren op min. 5 ha...")

cl_open_5ha_max <- cluster_filter_compleet(
  masker     = r_groen_open_max,
  opp_laag   = r_groen_open_max,
  drempel_m2 = 5 * 10000, 
  dist_m     = 50,
  werkelijk  = FALSE
)
vh_open1_max <- cl_open_5ha_max$raster

vh_bosrand2_max <- id_raster_VS * NA

if (!all(is.na(suppressWarnings(terra::minmax(vh_bossen3_max)))) && !all(is.na(suppressWarnings(terra::minmax(vh_open1_max))))) {
  dist_tot_open_max  <- terra::distance(vh_open1_max)
  r_bosrand_rauw_max <- vh_bossen3_max == 1 & dist_tot_open_max <= 15
  r_bosrand_rauw_max <- terra::ifel(r_bosrand_rauw_max == 1, 1, NA)
  
  cl_bosrand_max <- cluster_filter_compleet(
    masker     = r_bosrand_rauw_max,
    opp_laag   = r_bosrand_rauw_max,
    drempel_m2 = 0.25 * 10000, 
    dist_m     = 10,           
    werkelijk  = FALSE
  )
  vh_bosrand2_max <- cl_bosrand_max$raster
}
vh_leefgebied2_max <- vh_bosrand2_max

vh_bosrand2_opp <- terra::mask(vh_bossen3_opp, vh_bosrand2_max)
vh_leefgebied2_opp <- vh_bosrand2_opp

suppressWarnings(rm(cl_open_5ha_max, r_bosrand_rauw_max, cl_bosrand_max))
gc()

message("-> Bosdiepte-buffer van 40m berekenen vanaf de kwalitatieve bosrand...")

if (!all(is.na(suppressWarnings(terra::minmax(vh_leefgebied2_max))))) {
  dist_tot_rand_max <- terra::distance(vh_leefgebied2_max)
  vh_bosrand3_max   <- vh_bossen3_max == 1 & dist_tot_rand_max <= 40
  vh_bosrand3_max   <- terra::ifel(vh_bosrand3_max == 1, 1, NA)
  
  vh_leefgebied_bruto_max <- terra::cover(vh_leefgebied2_max, vh_bosrand3_max)
} else {
  vh_leefgebied_bruto_max <- id_raster_VS * NA
}

cl_finaal_max <- cluster_filter_compleet(
  masker     = vh_leefgebied_bruto_max,
  opp_laag   = vh_leefgebied_bruto_max,
  drempel_m2 = 5 * 10000, 
  dist_m     = 50,
  werkelijk  = FALSE
)
vliegendhert_leefgebied_max <- cl_finaal_max$raster %>% terra::crop(template_VS)
vliegendhert_leefgebied_opp <- terra::mask(vh_bossen3_opp, vliegendhert_leefgebied_max) %>% terra::crop(template_VS)

# DEFINITIEVE MODELUITGANGEN
final_max <- vliegendhert_leefgebied_max
final_opp <- vliegendhert_leefgebied_opp

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

suppressWarnings(rm(vh_bossen3_max, vh_bossen3_opp, vh_open1_max, vh_leefgebied2_max, 
                    vh_leefgebied2_opp, vh_bosrand3_max, vh_leefgebied_bruto_max, cl_finaal_max))
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
