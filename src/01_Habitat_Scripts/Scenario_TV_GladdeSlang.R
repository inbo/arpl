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

cluster_filter_twee_drempels <- function(masker, opp_laag, min_opp_ha, min_totale_opp_ha, dist_m, werkelijk = FALSE) {
  if (terra::global(is.na(masker), "sum")[[1]] == terra::ncell(masker)) {
    return(list(raster = masker * NA, clusters = masker * NA))
  }
  
  min_opp_m2 <- min_opp_ha * 10000
  min_totale_m2 <- min_totale_opp_ha * 10000
  
  # 1. Eerst individuele clusters identificeren en filteren op MinOpp_ha
  cl_individueel <- terra::patches(masker, directions = 8, zeroAsNA = TRUE)
  
  if(werkelijk) {
    stats_indiv <- terra::zonal(opp_laag, cl_individueel, fun = "sum", na.rm = TRUE)
    colnames(stats_indiv) <- c("ID", "Waarde")
    stats_indiv$Area_m2 <- stats_indiv$Waarde * 100
  } else {
    f_indiv <- terra::freq(cl_individueel)
    stats_indiv <- data.frame(ID = f_indiv$value, Waarde = f_indiv$count)
    stats_indiv$Area_m2 <- stats_indiv$Waarde * 100
  }
  
  voldoet_indiv_ids <- stats_indiv$ID[!is.na(stats_indiv$ID) & stats_indiv$Area_m2 >= min_opp_m2]
  if(length(voldoet_indiv_ids) == 0) return(list(raster = masker * NA, clusters = masker * NA))
  
  masker_gefilterd <- cl_individueel %in% voldoet_indiv_ids
  masker_gefilterd <- terra::ifel(masker_gefilterd == 1, masker, NA)
  
  # 2. Netwerkvorming op basis van Dispersiecap (dist_m)
  if (dist_m > 0) {
    r_binair <- terra::ifel(!is.na(masker_gefilterd) & masker_gefilterd > 0, 1, NA)
    r_buffered <- terra::buffer(r_binair, width = dist_m / 2)
    cl_network <- terra::patches(r_buffered, directions = 4, zeroAsNA = TRUE)
    cl_biotoop_only <- terra::mask(cl_network, masker_gefilterd)
  } else {
    cl_network <- terra::patches(masker_gefilterd, directions = 8, zeroAsNA = TRUE)
    cl_biotoop_only <- cl_network
  }
  
  # 3. Oppervlakte-optelsom per netwerk en filteren op Min_totale_opp
  if(werkelijk) {
    stats_net <- terra::zonal(opp_laag, cl_biotoop_only, fun = "sum", na.rm = TRUE)
    colnames(stats_net) <- c("ID", "Waarde")
    stats_net$Area_m2 <- stats_net$Waarde * 100
  } else {
    f_net <- terra::freq(cl_biotoop_only)
    stats_net <- data.frame(ID = f_net$value, Waarde = f_net$count)
    stats_net$Area_m2 <- stats_net$Waarde * 100
  }
  
  stats_net <- stats_net[!is.na(stats_net$ID), ]
  if(nrow(stats_net) == 0) return(list(raster = masker * NA, clusters = masker * NA))
  
  voldoet_net_ids <- stats_net$ID[stats_net$Area_m2 >= min_totale_m2]
  if(length(voldoet_net_ids) == 0) return(list(raster = masker * NA, clusters = masker * NA))
  
  final_network_mask <- cl_biotoop_only %in% voldoet_net_ids
  final_network_mask <- terra::ifel(final_network_mask == 1, 1, NA)
  
  r_finaal  <- terra::mask(masker_gefilterd, final_network_mask)
  cl_finaal <- terra::mask(cl_biotoop_only, r_finaal)
  
  return(list(raster = r_finaal, clusters = cl_finaal))
}

terraOptions(
  memfrac = 0.8,
  tempdir = tempdir(),
  verbose = FALSE
)

soort <- "gladdeslang"

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

min_oppervlakte_ha <- 20
oppervlakte_ha     <- min_oppervlakte_ha
min_totale_opp_ha  <- 760
afstand_m          <- 0   
buffer_m           <- 500 

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
values(globale_id_raster) <- 1:ncell(globale_id_raster)

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

rm(globale_id_raster, id_raster_TV_globale_values, id_raster_TV_masked, df_extractie)
gc()

# 1. LAAD DE BRON-CROSSWALK/DICTIONARY IN
df_nieuw <- read_csv(here("data/input/Excel_files/Resultaten_Totaal_Samengevoegd.csv"), show_col_types = FALSE)

resultaten_gegroepeerd <- df_nieuw %>%
  mutate(Soort_clean = tolower(trimws(Soort))) %>%
  filter(Soort_clean == soort) %>%
  group_by(Type) %>%
  nest(Data = c(Code, Match))

# 2. LAAD DE SCENARIO TABEL IN
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

heide_max <- lijst_matches[["biotoop_bwk"]]
heide_opp <- lijst_oppervlaktes[["biotoop_bwk"]]
lijn_max  <- lijst_matches[["lijn_bwk"]]
lijn_opp  <- lijst_oppervlaktes[["lijn_bwk"]]

rm(resultaten_gegroepeerd, df_nieuw, lijst_matches, lijst_oppervlaktes, tabel_vlaanderen, vertaal_df)
gc()

gladdeslang_biotoop1_max <- cover(heide_max, lijn_max)
gladdeslang_biotoop1_opp <- cover(heide_opp, lijn_opp)

hooggroen_raw          <- rast(here("data/input/ASCI Files/Groenkaart_2021.tif"))
area_buffer_groen_crs <- project(area_buffer_fix, crs(hooggroen_raw))
hooggroen_local_raw   <- crop(hooggroen_raw, area_buffer_groen_crs)

hooggroen_binair  <- terra::classify(hooggroen_local_raw, matrix(c(0.5, 1.5, 1), ncol = 3, byrow = TRUE), others = 0)
hooggroen_10m_raw <- terra::aggregate(hooggroen_binair, fact = 10, fun = "max")
hooggroen_sync    <- project(hooggroen_10m_raw, crs(id_raster_TV)) %>% terra::resample(id_raster_TV, method = "near")

rand_breedte_m <- 50

# ==============================================================================
# STAP 1 & 2: HABITATOPBOUW (Open habitat + 50m bosrand)
# ==============================================================================
gladdeslang_open_max <- terra::mask(gladdeslang_biotoop1_max, hooggroen_sync, maskvalues = 1)
gladdeslang_open_opp <- terra::mask(gladdeslang_biotoop1_opp, hooggroen_sync, maskvalues = 1)

zoekzone_rand_max <- terra::buffer(gladdeslang_open_max, width = rand_breedte_m)
r_binair_open_opp <- terra::ifel(!is.na(gladdeslang_open_opp) & gladdeslang_open_opp > 0, 1, NA)
zoekzone_rand_opp <- terra::buffer(r_binair_open_opp, width = rand_breedte_m)

gladdeslang_bufferbos_max <- terra::ifel(hooggroen_sync == 1 & zoekzone_rand_max == 1, 1, NA)
gladdeslang_bufferbos_opp <- terra::ifel(hooggroen_sync == 1 & zoekzone_rand_opp == 1, 1, NA)

bruto_biotoop3_max <- terra::cover(gladdeslang_open_max, gladdeslang_bufferbos_max)
bruto_biotoop3_opp <- terra::cover(gladdeslang_open_opp, gladdeslang_bufferbos_opp)

# ==============================================================================
# STAP 3: DUBBELE DREMPELFILTERING (Individueel >= 20 ha EN Netwerk >= 760 ha)
# ==============================================================================
res_max <- cluster_filter_twee_drempels(
  masker            = bruto_biotoop3_max,
  opp_laag          = bruto_biotoop3_opp,
  min_opp_ha        = min_oppervlakte_ha,
  min_totale_opp_ha = min_totale_opp_ha,
  dist_m            = buffer_m, 
  werkelijk         = FALSE
)

r_binair_bruto_opp <- terra::ifel(!is.na(bruto_biotoop3_opp) & bruto_biotoop3_opp > 0, 1, NA)
res_opp <- cluster_filter_twee_drempels(
  masker            = r_binair_bruto_opp,
  opp_laag          = bruto_biotoop3_opp,
  min_opp_ha        = min_oppervlakte_ha,
  min_totale_opp_ha = min_totale_opp_ha,
  dist_m            = buffer_m,
  werkelijk         = TRUE
)

final_max <- res_max$raster
final_opp <- res_opp$raster
cl_max    <- res_max$clusters
cl_opp    <- res_opp$clusters

if (is.null(cl_opp) || all(is.na(terra::values(cl_opp, mat=FALSE)))) {
  if (!is.null(final_opp) && !all(is.na(terra::values(final_opp, mat=FALSE)))) {
    cl_opp <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE)
  } else {
    cl_opp <- template_TV * NA
  }
}

suppressWarnings(
  rm(gladdeslang_open_max, gladdeslang_open_opp, zoekzone_rand_max, zoekzone_rand_opp, 
     gladdeslang_bufferbos_max, gladdeslang_bufferbos_opp, r_binair_open_opp, r_binair_bruto_opp,
     bruto_biotoop3_max, bruto_biotoop3_opp, hooggroen_binair, hooggroen_10m_raw, hooggroen_raw, hooggroen_local_raw)
)
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

# 3. Bepaal Analytisch Metacluster ID-raster (EXCLUSIEF OP BASIS VAN WERKELIJKE OPPERVLAKTES)
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
