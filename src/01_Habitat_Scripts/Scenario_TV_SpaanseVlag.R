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

soort <- "spaansevlag"

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

voortplanting_bwk_max <- lijst_matches[["voortplanting_bwk"]]
voortplanting_bwk_opp <- lijst_oppervlaktes[["voortplanting_bwk"]]
bos_bwk_max          <- lijst_matches[["bos_bwk"]]
bos_bwk_opp          <- lijst_oppervlaktes[["bos_bwk"]]

rm(vertaal_df, lijst_matches, lijst_oppervlaktes, tabel_vlaanderen, resultaten_gegroepeerd, df_nieuw)
gc()

# ==============================================================================
# STAP 1: Voortplantingsgebieden afbakenen (Rupsen - PARALLEL)
# ==============================================================================
message("-> Voortplantingsgebieden rupsen verwerken...")

# --- SPOOR A: MAX ---
voortplanting_geclusterd_max <- cluster_filter_compleet(
  masker     = voortplanting_bwk_max,
  opp_laag   = voortplanting_bwk_opp,
  drempel_m2 = 10000, 
  dist_m     = 50,
  werkelijk  = FALSE
)
spaansevlag_voortplanting2_max <- voortplanting_geclusterd_max$raster

# --- SPOOR B: OPP (Parallel & Autonoom) ---
r_binair_vp_opp <- terra::ifel(!is.na(voortplanting_bwk_opp) & voortplanting_bwk_opp > 0, 1, NA)
voortplanting_geclusterd_opp <- cluster_filter_compleet(
  masker     = r_binair_vp_opp,
  opp_laag   = voortplanting_bwk_opp,
  drempel_m2 = 10000, 
  dist_m     = 50,
  werkelijk  = TRUE
)
spaansevlag_voortplanting2_opp <- voortplanting_geclusterd_opp$raster

rm(voortplanting_geclusterd_max, voortplanting_geclusterd_opp, r_binair_vp_opp)
gc()

# ==============================================================================
# STAP 2: Bossen met bosranden afbakenen (Adulte vlinders - PARALLEL)
# ==============================================================================
message("-> Boshabitats verwerken...")

# --- SPOOR A: MAX ---
bos_geclusterd_max <- cluster_filter_compleet(
  masker     = bos_bwk_max,
  opp_laag   = bos_bwk_opp,
  drempel_m2 = 10000, 
  dist_m     = 50,
  werkelijk  = FALSE
)
spaansevlag_bos_max <- bos_geclusterd_max$raster

# --- SPOOR B: OPP (Parallel & Autonoom) ---
r_binair_bos_opp <- terra::ifel(!is.na(bos_bwk_opp) & bos_bwk_opp > 0, 1, NA)
bos_geclusterd_opp <- cluster_filter_compleet(
  masker     = r_binair_bos_opp,
  opp_laag   = bos_bwk_opp,
  drempel_m2 = 10000, 
  dist_m     = 50,
  werkelijk  = TRUE
)
spaansevlag_bos_opp <- bos_geclusterd_opp$raster

rm(bos_geclusterd_max, bos_geclusterd_opp, r_binair_bos_opp)
gc()

# ==============================================================================
# STAP 3: Scherpe Bosranden (Contactzones) Bepalen (PARALLEL)
# ==============================================================================
message("-> Contactzones scherp identificeren...")

bepaal_bosrand <- function(r_voortplanting, r_bos_cluster) {
  v_bin <- terra::ifel(!is.na(r_voortplanting) & r_voortplanting > 0, 1, 0)
  b_bin <- terra::ifel(!is.na(r_bos_cluster) & r_bos_cluster > 0, 1, 0)
  
  bos_buffer <- terra::focal(b_bin, w = 3, fun = "max", na.rm = TRUE)
  bosrand <- terra::ifel(v_bin == 1 & bos_buffer == 1, 1, NA)
  return(bosrand)
}

# --- SPOOR A: MAX ---
dist_tot_bos_max <- terra::distance(spaansevlag_bos_max)
spaansevlag_bosrand_vpt_max <- terra::ifel(!is.na(spaansevlag_voortplanting2_max) & dist_tot_bos_max <= 10, 1, NA)

dist_tot_vpt_max <- terra::distance(spaansevlag_voortplanting2_max)
spaansevlag_bosrand_bos_max <- terra::ifel(!is.na(spaansevlag_bos_max) & dist_tot_vpt_max <= 10, 1, NA)

spaansevlag_bosrand1_max <- terra::cover(spaansevlag_bosrand_vpt_max, spaansevlag_bosrand_bos_max)

# --- SPOOR B: OPP (Parallel & Autonoom) ---
dist_tot_bos_opp <- terra::distance(!is.na(spaansevlag_bos_opp) & spaansevlag_bos_opp > 0)
spaansevlag_bosrand_vpt_opp <- terra::mask(spaansevlag_voortplanting2_opp, terra::ifel(dist_tot_bos_opp <= 10, 1, NA))

dist_tot_vpt_opp <- terra::distance(!is.na(spaansevlag_voortplanting2_opp) & spaansevlag_voortplanting2_opp > 0)
spaansevlag_bosrand_bos_opp <- terra::mask(spaansevlag_bos_opp, terra::ifel(dist_tot_vpt_opp <= 10, 1, NA))

spaansevlag_bosrand1_opp <- terra::cover(spaansevlag_bosrand_vpt_opp, spaansevlag_bosrand_bos_opp)

suppressWarnings(rm(dist_tot_bos_max, dist_tot_vpt_max, spaansevlag_bosrand_vpt_max, spaansevlag_bosrand_bos_max,
                    dist_tot_bos_opp, dist_tot_vpt_opp, spaansevlag_bosrand_vpt_opp, spaansevlag_bosrand_bos_opp))
gc()

# ==============================================================================
# STAP 4: Echte aaneengesloten lengte filteren (>= 100m - PARALLEL)
# ==============================================================================
message("-> Bosranden filteren op ECHTE aaneengesloten lengte (>= 100m)...")

# --- SPOOR A: MAX ---
r_scherpe_rand_clusters_max <- terra::patches(spaansevlag_bosrand1_max, directions = 8, zeroAsNA = TRUE)
freq_rand_max <- terra::freq(r_scherpe_rand_clusters_max)
freq_rand_max$area_ha <- freq_rand_max$count * 0.01
goedgekeurde_rand_ids_max <- freq_rand_max$value[freq_rand_max$area_ha >= 0.25]

if(length(goedgekeurde_rand_ids_max) > 0) {
  spaansevlag_bosrand2_max <- r_scherpe_rand_clusters_max %in% goedgekeurde_rand_ids_max
  spaansevlag_bosrand2_max <- terra::ifel(spaansevlag_bosrand2_max == 1, 1, NA)
} else {
  spaansevlag_bosrand2_max <- template_TV * NA
}

# --- SPOOR B: OPP (Parallel & Autonoom) ---
r_binair_rand1_opp <- terra::ifel(!is.na(spaansevlag_bosrand1_opp) & spaansevlag_bosrand1_opp > 0, 1, NA)
r_scherpe_rand_clusters_opp <- terra::patches(r_binair_rand1_opp, directions = 8, zeroAsNA = TRUE)

freq_rand_opp <- terra::zonal(spaansevlag_bosrand1_opp, r_scherpe_rand_clusters_opp, fun = "sum", na.rm = TRUE)
colnames(freq_rand_opp) <- c("value", "ha_exact")
freq_rand_opp$ha_exact <- freq_rand_opp$ha_exact * 0.01

goedgekeurde_rand_ids_opp <- freq_rand_opp$value[freq_rand_opp$ha_exact >= 0.25]

if(length(goedgekeurde_rand_ids_opp) > 0) {
  masker_rand_opp <- r_scherpe_rand_clusters_opp %in% goedgekeurde_rand_ids_opp
  spaansevlag_bosrand2_opp <- terra::mask(spaansevlag_bosrand1_opp, terra::ifel(masker_rand_opp, 1, NA))
} else {
  spaansevlag_bosrand2_opp <- template_TV * NA
}

suppressWarnings(rm(r_scherpe_rand_clusters_max, freq_rand_max, goedgekeurde_rand_ids_max,
                    r_binair_rand1_opp, r_scherpe_rand_clusters_opp, freq_rand_opp, goedgekeurde_rand_ids_opp,
                    spaansevlag_bosrand1_max, spaansevlag_bosrand1_opp))
gc()

# ==============================================================================
# STAP 5: Finaal leefgebied via snelle buffer-intersectie (< 100m - PARALLEL)
# ==============================================================================
message("-> Finaal leefgebied berekenen via snelle binaire buffers...")

# --- SPOOR A: MAX ---
bosrand_buffer_100m_max <- terra::buffer(spaansevlag_bosrand2_max, width = 100)

spaansevlag_leefgebied1_max <- terra::mask(spaansevlag_voortplanting2_max, bosrand_buffer_100m_max)
spaansevlag_leefgebied2_max <- terra::mask(spaansevlag_bosrand2_max, bosrand_buffer_100m_max)

spaansevlag_leefgebied_finaal_max <- terra::cover(spaansevlag_leefgebied1_max, spaansevlag_leefgebied2_max)
spaansevlag_leefgebied_finaal_max <- terra::ifel(spaansevlag_leefgebied_finaal_max == 1, 1, NA)

# --- SPOOR B: OPP (Parallel & Autonoom) ---
r_binair_rand2_opp <- terra::ifel(!is.na(spaansevlag_bosrand2_opp) & spaansevlag_bosrand2_opp > 0, 1, NA)
bosrand_buffer_100m_opp <- terra::buffer(r_binair_rand2_opp, width = 100)

spaansevlag_leefgebied1_opp <- terra::mask(spaansevlag_voortplanting2_opp, bosrand_buffer_100m_opp)
spaansevlag_leefgebied2_opp <- terra::mask(spaansevlag_bosrand2_opp, bosrand_buffer_100m_opp)

spaansevlag_leefgebied_finaal_opp <- terra::cover(spaansevlag_leefgebied1_opp, spaansevlag_leefgebied2_opp)

# KOPPELING NAAR HET MODEL
final_max <- spaansevlag_leefgebied_finaal_max %>% terra::crop(template_TV)
final_opp <- spaansevlag_leefgebied_finaal_opp %>% terra::crop(template_TV)

if (!all(is.na(terra::values(final_opp, mat=FALSE)))) {
  cl_opp <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE)
} else {
  cl_opp <- template_TV * NA
}

if (!all(is.na(terra::values(final_max, mat=FALSE)))) {
  cl_max <- terra::patches(final_max, directions = 8, zeroAsNA = TRUE)
} else {
  cl_max <- template_TV * NA
}

suppressWarnings(rm(bosrand_buffer_100m_max, bosrand_buffer_100m_opp, r_binair_rand2_opp,
                    spaansevlag_leefgebied1_max, spaansevlag_leefgebied2_max, spaansevlag_leefgebied1_opp, spaansevlag_leefgebied2_opp,
                    spaansevlag_voortplanting2_max, spaansevlag_voortplanting2_opp, spaansevlag_bos_max, spaansevlag_bos_opp,
                    spaansevlag_bosrand2_max, spaansevlag_bosrand2_opp, spaansevlag_leefgebied_finaal_max, spaansevlag_leefgebied_finaal_opp))
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
