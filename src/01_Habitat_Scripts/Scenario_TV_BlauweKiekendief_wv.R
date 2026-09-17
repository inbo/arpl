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

cluster_filter_compleet <- function(masker, opp_laag, drempel_m2, dist_m, werkelijk = FALSE, straal_m = 0) {
  if (terra::global(is.na(masker), "sum")[[1]] == terra::ncell(masker)) {
    return(list(raster = masker * NA, clusters = masker * NA))
  }
  
  if (exists("id_raster_TV")) {
    masker   <- terra::crop(masker, id_raster_TV, snap = "out")
    opp_laag <- terra::crop(opp_laag, id_raster_TV, snap = "out")
  }
  
  if (straal_m > 0) {
    r_binair <- terra::ifel(!is.na(masker) & masker > 0, 1, NA)
    if (terra::global(is.na(r_binair), "sum")[[1]] == terra::ncell(r_binair)) {
      return(list(raster = masker * NA, clusters = masker * NA))
    }
    
    v_biotoop <- terra::as.polygons(r_binair, aggregate = TRUE)
    v_eroded  <- terra::buffer(v_biotoop, width = -straal_m)
    
    if (length(v_eroded) == 0 || terra::geomtype(v_eroded) == "none") {
      return(list(raster = masker * NA, clusters = masker * NA))
    }
    
    v_dilated <- terra::buffer(v_eroded, width = straal_m)
    r_zuiver  <- terra::rasterize(v_dilated, r_binair, field = 1, background = NA)
    r_masker_werkelijk <- terra::mask(r_zuiver, masker)
    rm(v_biotoop, v_eroded, v_dilated, r_zuiver)
  } else {
    r_masker_werkelijk <- masker
  }
  
  if (terra::global(is.na(r_masker_werkelijk), "sum")[[1]] == terra::ncell(r_masker_werkelijk)) {
    return(list(raster = masker * NA, clusters = masker * NA))
  }
  
  if (dist_m > 0) {
    r_binair_z <- terra::ifel(!is.na(r_masker_werkelijk) & r_masker_werkelijk > 0, 1, NA)
    r_buffered <- terra::buffer(r_binair_z, width = dist_m / 2)
    cl_network <- terra::patches(r_buffered, directions = 4, zeroAsNA = TRUE)
    cl_biotoop_only <- terra::mask(cl_network, r_masker_werkelijk)
  } else {
    cl_network <- terra::patches(r_masker_werkelijk, directions = 8, zeroAsNA = TRUE)
    cl_biotoop_only <- cl_network
  }
  
  if(werkelijk) {
    opp_laag_sub <- terra::crop(opp_laag, cl_biotoop_only)
    stats_df     <- terra::zonal(opp_laag_sub, cl_biotoop_only, fun = "sum", na.rm = TRUE)
    colnames(stats_df) <- c("ID", "Waarde")
    stats_df$Area_m2  <- stats_df$Waarde * 100 
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
  
  r_finaal  <- terra::mask(r_masker_werkelijk, final_network_mask)
  cl_finaal <- terra::mask(cl_biotoop_only, r_finaal) 
  
  r_finaal  <- terra::deepcopy(r_finaal)
  cl_finaal <- terra::deepcopy(cl_finaal)
  
  return(list(raster = r_finaal, clusters = cl_finaal))
}

terraOptions(
  memfrac = 0.8,
  tempdir = tempdir(),
  verbose = FALSE
)

df <- read_excel(here::here("data/input/Excel_files/Soorten_bwk_afstanden.xlsx"))
soort <- "blauwekiekendief"

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

resultaat <- df %>%
  filter(tolower(trimws(Soort)) == soort) %>%
  select(Type, MinOpp_ha, AfstandBiotopen_m, Dispersiecap_m)

oppervlakte_ha <- 50
buffer_m       <- 50000
actieradius_m  <- 50000 

rm(df, resultaat)

area_shape  <- vect(here("data/input/Turnhouts_Vennegebied.shp"))
master_grid <- rast(here("data/input/Raster_Vlaanderen/Vlaanderen_MasterGrid_10m.tif"))[[1]]

df_namen_sleutel  <- read_csv(here("data/input/Excel_files/BWK_Laag_Namen_2025.csv"), show_col_types = FALSE)
gouden_namenlijst <- tolower(trimws(df_namen_sleutel$Laagnaam))

area_shape_proj <- project(area_shape, crs(master_grid))
area_buffer_fix <- buffer(area_shape_proj, width = buffer_m)

message("-> Vertaalraster voor globale/lokale cellen opbouwen...")
id_raster_TV <- crop(master_grid, area_buffer_fix, snap = "out")
values(id_raster_TV) <- 1:terra::ncell(id_raster_TV)

globale_ids_vector <- terra::cells(master_grid, terra::ext(id_raster_TV))

vertaal_df <- data.table(
  lokale_id  = 1:terra::ncell(id_raster_TV),
  globale_id = globale_ids_vector
)

vertaal_df <- vertaal_df[!is.na(globale_id)]
studiegebied_globale_ids <- unique(vertaal_df$globale_id)

template_TV <- terra::rasterize(area_buffer_fix, id_raster_TV, field = 1, background = NA)

rm(globale_ids_vector)
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
  tabel_TV_unique <- unique(tabel_TV, by = c("cel_id", "CODE"))

  tabel_cel_som <- tabel_TV_unique[, .(Oppervlakte = pmin(sum(BWK_FRAC, na.rm = TRUE), 1.0)), by = .(cel_id)]
  
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

foer_raw_max <- lijst_matches[["foerageer"]]
foer_raw_opp <- lijst_oppervlaktes[["foerageer"]]

rm(resultaten_gegroepeerd, df_nieuw, lijst_matches, lijst_oppervlaktes)
gc()

# Landbouwgebruikspercelen (LGP) toevoegen
r_doel_na <- id_raster_TV * NA
r_lgp_vlaanderen <- rast(here("data/input/Raster_Vlaanderen/vlaanderen_lgp_gwscod_2025_10m.tif"))
r_lgp_local <- r_lgp_vlaanderen %>% 
  terra::crop(area_buffer_fix) %>% 
  terra::resample(r_doel_na, method = "near")

kiekendief_winter_gewassen <- c(
  100, 101, 102, 
  201, 202,      
  311, 321,      
  601, 602, 603, 
  800, 801, 810, 
  900, 910, 931, 951 
) 

r_lgp_bin <- terra::ifel(!is.na(r_lgp_local) & r_lgp_local %in% kiekendief_winter_gewassen, 1, NA)

foer_raw_max <- terra::cover(foer_raw_max, r_lgp_bin)
foer_raw_opp <- terra::cover(foer_raw_opp, r_lgp_bin)

# --- BEBOUWINGSFILTER ---
if (!exists("tabel_vlaanderen")) {
  tabel_vlaanderen <- readRDS(scenario_path)
  setDT(tabel_vlaanderen)
  tabel_vlaanderen[, CODE := tolower(trimws(CODE))]
}

bebouwings_lagen <- gouden_namenlijst[grepl("^u[a-z0-9]", gouden_namenlijst) | gouden_namenlijst == "u" | gouden_namenlijst == "gh"]

tabel_bebouw_TV  <- tabel_vlaanderen[CODE %in% bebouwings_lagen & cel_id %in% studiegebied_globale_ids]
tabel_bebouw_som <- tabel_bebouw_TV[, .(Present = ifelse(sum(BWK_FRAC, na.rm = TRUE) >= 0.3, 1, 0)), by = .(cel_id)]

r_bebouwing_raw <- id_raster_TV * NA

if(nrow(tabel_bebouw_som) > 0) {
  setnames(tabel_bebouw_som, "cel_id", "globale_id")
  tabel_bebouw_mapping <- merge(tabel_bebouw_som, vertaal_df, by = "globale_id", all.x = TRUE)
  tabel_bebouw_mapping <- tabel_bebouw_mapping[!is.na(lokale_id) & Present == 1]
  
  if(nrow(tabel_bebouw_mapping) > 0) {
    r_bebouwing_raw[tabel_bebouw_mapping$lokale_id] <- 1
  }
}

r_bebouwing_schoon <- terra::ifel(!is.na(r_bebouwing_raw) & r_bebouwing_raw == 1, 1, NA)

rm(tabel_bebouw_TV, tabel_bebouw_som, tabel_bebouw_mapping, r_bebouwing_raw)
gc()

# ==============================================================================
# STAP 1: VERBOSSINGSFILTER
# ==============================================================================
message("-> Externe groenkaart (fracties) inlezen en afsnijden op studiegebied...")

r_hooggroen_frac <- rast(here("data/input/Raster_Vlaanderen/Groenkaart_2021_Fracties_5Bands_10m.tif"))[[1]] %>% 
  terra::crop(id_raster_TV, snap = "out")

message("-> Verbossingsfilter berekenen: hooggroencomplexen groter dan 20 are uitsluiten...")

drempel_fractie <- 0.3 
r_hooggroen_bin_10m <- terra::ifel(!is.na(r_hooggroen_frac) & r_hooggroen_frac >= drempel_fractie, 1, NA)

w_3x3 <- matrix(1, nrow = 3, ncol = 3)
r_som_3x3 <- terra::focal(r_hooggroen_bin_10m, w = w_3x3, fun = "sum", na.rm = TRUE)
r_hooggroen_schoon_10m <- terra::ifel(r_som_3x3 >= 2 & r_hooggroen_bin_10m == 1, 1, NA)

r_hooggroen_bin_40m <- terra::aggregate(r_hooggroen_schoon_10m, fact = 4, fun = "max", na.rm = TRUE)

verbossing_clusters_40m <- cluster_filter_compleet(
  masker     = r_hooggroen_bin_40m,
  opp_laag   = r_hooggroen_bin_40m,
  drempel_m2 = 2000,                
  dist_m     = 0,                    
  werkelijk  = FALSE      
)

r_bos_10m_back <- terra::disagg(verbossing_clusters_40m$raster, fact = 4) %>% 
  terra::crop(id_raster_TV, snap = "near")

r_bos_20are <- terra::ifel(!is.na(r_bos_10m_back) & r_bos_10m_back > 0, 1, NA)

gc()

# ==============================================================================
# MODEL OP 50 HA OPEN JAAGGEBIED
# ==============================================================================
message("-> Filteren van open jaaggebied (>= 50 ha) zonder dicht bos...")

masker_bebouw_75m <- terra::buffer(r_bebouwing_schoon, width = 75)

r_bos_masker <- terra::ifel(!is.na(r_bos_20are) & r_bos_20are > 0, 1, NA)
r_bebouw_masker <- terra::ifel(!is.na(masker_bebouw_75m) & masker_bebouw_75m > 0, 1, NA)

foer_schoon_max <- terra::ifel(!is.na(r_bos_masker) | !is.na(r_bebouw_masker), NA, foer_raw_max)
foer_schoon_opp <- terra::ifel(!is.na(r_bos_masker) | !is.na(r_bebouw_masker), NA, foer_raw_opp)

res_foer_max <- cluster_filter_compleet(
  masker     = foer_schoon_max,
  opp_laag   = foer_schoon_max,
  drempel_m2 = 500000, 
  dist_m     = 400,    
  werkelijk  = FALSE
)

res_foer_opp <- cluster_filter_compleet(
  masker     = foer_schoon_opp,
  opp_laag   = foer_schoon_opp,
  drempel_m2 = 500000, 
  dist_m     = 400,    
  werkelijk  = TRUE
)

kiekendief_leefgebied_max <- res_foer_max$raster
kiekendief_leefgebied_opp <- res_foer_opp$raster

# ==============================================================================
# DEFINITIEVE TOEWIJSING EXPORT VARIABELEN (CRUCIALE FIX)
# ==============================================================================
final_max <- kiekendief_leefgebied_max
final_opp <- kiekendief_leefgebied_opp

if (!all(is.na(terra::values(final_max, mat=FALSE)))) {
  cl_max <- terra::patches(final_max, directions = 8, zeroAsNA = TRUE)
} else {
  cl_max <- template_TV * NA
}

if (!all(is.na(terra::values(final_opp, mat=FALSE)))) {
  cl_opp <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE)
} else {
  cl_opp <- template_TV * NA
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
potentie_export_rast <- if (exists("final_max") && !is.null(final_max) && !all(is.na(terra::values(final_max, mat=FALSE)))) {
  terra::ifel(!is.na(final_max) & final_max > 0, 1, NA)
} else {
  terra::rast(template_TV, vals = NA)
}

# 2. Bepaal Werkelijke Oppervlakte Raster
werkelijk_export_rast <- if (exists("final_opp") && !is.null(final_opp) && !all(is.na(terra::values(final_opp, mat=FALSE)))) {
  terra::ifel(!is.na(final_opp) & final_opp > 0, 1, NA)
} else {
  terra::rast(template_TV, vals = NA)
}

# 3. Bepaal Analytisch Metacluster ID-raster
if (exists("cl_opp") && !is.null(cl_opp) && !all(is.na(terra::values(cl_opp, mat=FALSE)))) {
  id_export_rast <- cl_opp
} else if (exists("final_opp") && !is.null(final_opp) && !all(is.na(terra::values(final_opp, mat=FALSE)))) {
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
  suffix <- if ((exists("zoek_sleutel") && grepl("_wv$", zoek_sleutel)) || (exists("soort") && grepl("blauwekiekendief", soort, ignore.case = TRUE))) "_wv.tif" else ".tif"
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
suffix <- if ((exists("zoek_sleutel") && grepl("_wv$", zoek_sleutel)) || (exists("soort") && grepl("blauwekiekendief", soort, ignore.case = TRUE))) "_wv.tif" else ".tif"
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
