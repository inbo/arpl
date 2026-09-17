library(here)
# CRUCIALE FIX: Dwing R Markdown om te werken vanaf de hoofdmap (arpl/)
# Hierdoor werken relatieve paden met hier() overal hetzelfde.

library(knitr)
library(tidyverse)
library(sf)
library(terra)
library(readxl)
library(tidyterra)
library(leaflet)
library(kableExtra)
library(data.table)

conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::select)
conflicted::conflicts_prefer(terra::intersect)
conflicted::conflicts_prefer(terra::any)

# Exacte hectare-berekening behouden (cellSize voor 100% precisie)
calc_ha_exact <- function(r) {
  if(is.null(r)) return(0)
  if(all(is.na(terra::values(r, mat=FALSE)))) return(0)
  area_raster <- r * terra::cellSize(r, unit = "ha")
  val <- terra::global(area_raster, "sum", na.rm = TRUE)[[1]]
  return(as.numeric(val))
}

# Geoptimaliseerde get_stats helper
get_stats <- function(cat_id, label) {
  cid <- cat_id
  target_mask <- r_status == cid
  f_pix <- terra::freq(target_mask)
  n_pix <- if(nrow(f_pix) > 0) sum(f_pix$count[f_pix$value == 1], na.rm=TRUE) else 0
  area_ha <- (n_pix * 100) / 10000
  
  cl_src <- if(cid %in% c(1, 2)) cl_id_max else cl_id_opp
  cl_zone <- terra::mask(cl_src, target_mask)
  f_cl <- terra::freq(cl_zone)
  n_cl <- if(!is.null(f_cl) && nrow(f_cl) > 0) nrow(f_cl) else 0
  
  return(data.frame(Type = label, Clusters = n_cl, Oppervlakte_ha = round(area_ha, 2)))
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
  # 1. Controleer of het raster leeg is
  if (terra::global(is.na(masker), "sum")[[1]] == terra::ncell(masker)) {
    return(list(raster = masker * NA, clusters = masker * NA))
  }
  
  # ----------------------------------------------------------------------------
  # STAP 1: VEILIGE & SNELLER NETWERKVORMING ZONDER RAM-CRASH
  # ----------------------------------------------------------------------------
  if (dist_m > 0) {
    # Maak binaire kaart (1 = biotoop, NA = rest)
    r_binair <- terra::ifel(!is.na(masker) & masker > 0, 1, NA)
    
    # Buffer de BINAIR kaart (dit kost vrijwel geen geheugen!)
    # Volle afstand dist_m zorgt dat plukjes binnen dist_m gegarandeerd samensmelten
    r_buffered <- terra::buffer(r_binair, width = dist_m / 2)
    
    # Maak unieke Netwerk-ID's op de gebufferde zones
    cl_network <- terra::patches(r_buffered, directions = 4, zeroAsNA = TRUE)
    
    # Snijd de Netwerk-ID's direct terug naar waar de originele biotooppixels liggen
    cl_biotoop_only <- terra::mask(cl_network, masker)
  } else {
    cl_network <- terra::patches(masker, directions = 8, zeroAsNA = TRUE)
    cl_biotoop_only <- cl_network
  }
  
  # ----------------------------------------------------------------------------
  # STAP 2: OPPERVLAKTE-OPTELSOM PER GEKOPPELD NETWERK
  # ----------------------------------------------------------------------------
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
  
  # ----------------------------------------------------------------------------
  # STAP 3: FILTEREN OP TOTALE NETWERK-OPPERVLAKTE >= DREMPEL
  # ----------------------------------------------------------------------------
  voldoet_ids <- stats_df$ID[stats_df$Area_m2 >= drempel_m2]
  if(length(voldoet_ids) == 0) return(list(raster = masker * NA, clusters = masker * NA))
  
  # Behaal alleen de winnende netwerk-ID's
  masker_binair <- cl_biotoop_only %in% voldoet_ids
  final_network_mask <- terra::ifel(masker_binair == 1, 1, NA)
  
  r_finaal  <- terra::mask(masker, final_network_mask)
  cl_finaal <- terra::mask(cl_biotoop_only, r_finaal) 
  
  return(list(raster = r_finaal, clusters = cl_finaal))
}

terraOptions(
  memfrac = 0.8,        # Dwing terra om tot max. 80% van het RAM-geheugen te gebruiken
  tempdir = tempdir(),  # Geef toestemming voor automatische disk-swapping bij zware rasters
  verbose = FALSE
)

# Expert-parameters Hoogveenglanslibel volgens de experttabel van Greet
soort <- "hoogveenglanslibel"

# Scenario pad en naam bepalen

# --- DYNAMISCHE SCENARIO PARAMETER CHECK ---
if (!exists("params") || is.null(params$scenario_rds_path)) {
  scenario_rds_path <- "data/input/Scenario_rds/MH_Scenario_BWK_2025.rds"
} else {
  scenario_rds_path <- params$scenario_rds_path
}

p_raw <- gsub("^([.][.]/)+", "", scenario_rds_path)
scenario_path <- here::here(p_raw)


if (!file.exists(scenario_path)) {
  stop(paste("❌ FOUT: Scenario RDS bestand NIET gevonden op:", scenario_path))
}

scen_volledig <- basename(scenario_path)
scenario_naam <- gsub("^MH_Scenario_|^Scenario_|.rds$", "", scen_volledig)

message(paste("Verwerken van soort:", soort, "binnen scenario:", scenario_naam))

# Parameters uit de tabel van Greet:
buffer_m             <- 5000   # Dispersiecapaciteit voor smear/buffer (5000m)
straal_water_m       <- 500    # Harde afsnij-afstand rond het water (500m)

drempel_water_m2     <- 100    # Min. 0.01 ha water (100 m²)
dist_water_m         <- 50     # 50m overbrugging voor water

drempel_land_m2      <- 50000  # Min. 5 ha land
dist_land_m          <- 10     # 10m overbrugging voor aanpalend land

drempel_ruim_m2      <- 500000 # Min. 50 ha ruim land
dist_ruim_m          <- 50     # 50m overbrugging voor ruim land

area_shape  <- vect(here("data/input/Mechelse_Heide.shp"))
master_grid <- rast(here("data/input/Raster_Vlaanderen/Vlaanderen_MasterGrid_10m.tif"))[[1]]

df_namen_sleutel <- read_csv(here("data/input/Excel_files/BWK_Laag_Namen_2025.csv"), show_col_types = FALSE)
gouden_namenlijst <- tolower(trimws(df_namen_sleutel$Laagnaam))

area_shape_proj <- project(area_shape, crs(master_grid))
area_buffer_fix <- buffer(area_shape_proj, width = buffer_m)

message("-> Vertaalraster voor globale/lokale cellen opbouwen via snelle MASK methode...")
id_raster_MH <- crop(master_grid, area_buffer_fix, snap = "near")

globale_id_raster <- master_grid
values(globale_id_raster) <- 1:ncell(globale_id_raster)

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
    bevat_codes_escaped <- gsub("([\\.\\^\\$\\*\\+\\?\\(\\)\\[\\{\\\\\\|])", "\\\\\\1", bevat_codes)
    bevat_codes_anchored <- paste0("^", bevat_codes_escaped)
    regex_term <- paste0(bevat_codes_anchored, collapse = "|")
    
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

water1_max    <- lijst_matches[["waterbiotoop_bwk1"]]
water1_opp    <- lijst_oppervlaktes[["waterbiotoop_bwk1"]]
not_water_max <- lijst_matches[["waterbiotoop_notbwk"]]
not_water_opp <- lijst_oppervlaktes[["waterbiotoop_notbwk"]]
land1_max     <- lijst_matches[["landbiotoop_bwk1"]]
land1_opp     <- lijst_oppervlaktes[["landbiotoop_bwk1"]]
not_land_max  <- lijst_matches[["landbiotoop_notbwk"]]
not_land_opp  <- lijst_oppervlaktes[["landbiotoop_notbwk"]]
bos_max       <- lijst_matches[["landbiotoop_bos"]]
bos_opp       <- lijst_oppervlaktes[["landbiotoop_bos"]]
land_ruim_max <- lijst_matches[["landbiotoop_bwk_ruim"]]
land_ruim_opp <- lijst_oppervlaktes[["landbiotoop_bwk_ruim"]]

rm(tabel_vlaanderen, vertaal_df, lijst_matches, lijst_oppervlaktes)
gc()

waterbiotoop_max <- (water1_max == 1) & (not_water_max == 0 | is.na(not_water_max))
waterbiotoop_max <- terra::ifel(waterbiotoop_max == 1, 1, NA)
waterbiotoop_opp <- terra::mask(water1_opp, waterbiotoop_max)

# Spoor A: MAX
waterbiotoop_clusters_max <- cluster_filter_compleet(
  masker     = waterbiotoop_max,
  opp_laag   = waterbiotoop_opp,
  drempel_m2 = 100,    # 0.01 ha
  dist_m     = 50,
  werkelijk  = FALSE
)

# Spoor B: OPP (Parallel & Autonoom)
r_binair_water_opp <- terra::ifel(!is.na(waterbiotoop_opp) & waterbiotoop_opp > 0, 1, NA)
waterbiotoop_clusters_opp <- cluster_filter_compleet(
  masker     = r_binair_water_opp,
  opp_laag   = waterbiotoop_opp,
  drempel_m2 = 100, 
  dist_m     = 50,
  werkelijk  = TRUE
)

rm(r_binair_water_opp)
gc()

landbiotoop_max <- (land1_max == 1) & (not_land_max == 0 | is.na(not_land_max))
landbiotoop_max <- terra::ifel(landbiotoop_max == 1, 1, NA)
landbiotoop_opp <- terra::mask(land1_opp, landbiotoop_max)

# Spoor A: MAX
landbiotoop_clusters_max <- cluster_filter_compleet(
  masker     = landbiotoop_max,
  opp_laag   = landbiotoop_opp,
  drempel_m2 = 100000,   # 5 ha
  dist_m     = 10,
  werkelijk  = FALSE
)

# Spoor B: OPP (Parallel & Autonoom)
r_binair_land_opp <- terra::ifel(!is.na(landbiotoop_opp) & landbiotoop_opp > 0, 1, NA)
landbiotoop_clusters_opp <- cluster_filter_compleet(
  masker     = r_binair_land_opp,
  opp_laag   = landbiotoop_opp,
  drempel_m2 = 100000, 
  dist_m     = 10,
  werkelijk  = TRUE
)

rm(r_binair_land_opp)
gc()

# Spoor A: MAX
ruim_land_clusters_max <- cluster_filter_compleet(
  masker     = land_ruim_max,
  opp_laag   = land_ruim_opp,
  drempel_m2 = 1000000,  # 50 ha
  dist_m     = 50,
  werkelijk  = FALSE
)

# Spoor B: OPP (Parallel & Autonoom)
r_binair_ruim_opp <- terra::ifel(!is.na(land_ruim_opp) & land_ruim_opp > 0, 1, NA)
ruim_land_clusters_opp <- cluster_filter_compleet(
  masker     = r_binair_ruim_opp,
  opp_laag   = land_ruim_opp,
  drempel_m2 = 1000000, 
  dist_m     = 50,
  werkelijk  = TRUE
)

rm(r_binair_ruim_opp)
gc()

# ==============================================================================
# GECORRIGEERD: BOSPATCHES EN BOSPERCENTAGES PER CLUSTER (MIN 60% BOS)
# ==============================================================================
message("-> Start samenstellingsanalyse van landbiotoopclusters (min. 60% bos)...")

# --- SPOOR A: MAX ---
r_land_zuiver_max <- landbiotoop_clusters_max$raster
r_land_zuiver_max[is.na(r_land_zuiver_max) | r_land_zuiver_max == 0] <- NA

# Unieke ID's maken per fysieke aaneengesloten landcluster
cl_id_land_max <- terra::patches(r_land_zuiver_max, directions = 8, zeroAsNA = TRUE)

if (!all(is.na(terra::values(cl_id_land_max, mat=FALSE)))) {
  r_binair_patch_max <- terra::ifel(!is.na(cl_id_land_max), 1, NA)
  totale_opp_per_patch_max <- terra::zonal(r_binair_patch_max, cl_id_land_max, fun = "sum", na.rm = TRUE)
  colnames(totale_opp_per_patch_max) <- c("Cluster_ID", "Totaal_Land_ha")
  totale_opp_per_patch_max$Totaal_Land_ha <- totale_opp_per_patch_max$Totaal_Land_ha * 0.01

  bos_binair_max <- terra::ifel(!is.na(bos_max) & bos_max > 0, 1, NA)
  bos_in_land_max <- terra::mask(bos_binair_max, cl_id_land_max)
  stats_bos_in_land_max <- terra::zonal(bos_in_land_max, cl_id_land_max, fun = "sum", na.rm = TRUE)
  colnames(stats_bos_in_land_max) <- c("Cluster_ID", "Bos_ha")
  stats_bos_in_land_max$Bos_ha <- stats_bos_in_land_max$Bos_ha * 0.01

  df_samenstelling_land_max <- merge(totale_opp_per_patch_max, stats_bos_in_land_max, by = "Cluster_ID", all.x = TRUE)
  df_samenstelling_land_max$Bos_ha[is.na(df_samenstelling_land_max$Bos_ha)] <- 0
  df_samenstelling_land_max$Percentage_Bos <- round((df_samenstelling_land_max$Bos_ha / df_samenstelling_land_max$Totaal_Land_ha) * 100, 2)
  df_samenstelling_land_max$Totaal_Land_ha <- round(df_samenstelling_land_max$Totaal_Land_ha, 2)
  df_samenstelling_land_max$Bos_ha <- round(df_samenstelling_land_max$Bos_ha, 2)

  cat("\n=== [SPOOR MAX - HOOGVEENGLANSLIBEL] BOSPATCHES PER CLUSTER ===\n")
  print(head(df_samenstelling_land_max, 20))
  cat("=============================================================\n\n")

  valide_land_ids_max <- df_samenstelling_land_max$Cluster_ID[df_samenstelling_land_max$Percentage_Bos >= 60]

  if(length(valide_land_ids_max) > 0) {
    m_match_max <- terra::match(cl_id_land_max, valide_land_ids_max)
    land_nabij_max <- terra::ifel(!is.na(m_match_max), 1, NA)
  } else {
    land_nabij_max <- template_MH * NA
  }
} else {
  land_nabij_max <- template_MH * NA
  valide_land_ids_max <- c()
}

# --- SPOOR B: OPP ---
r_land_zuiver_opp <- landbiotoop_clusters_opp$raster
r_binair_land_opp <- terra::ifel(!is.na(r_land_zuiver_opp) & r_land_zuiver_opp > 0, 1, NA)

cl_id_land_opp <- terra::patches(r_binair_land_opp, directions = 8, zeroAsNA = TRUE)

if (!all(is.na(terra::values(cl_id_land_opp, mat=FALSE)))) {
  totale_opp_per_patch_opp <- terra::zonal(landbiotoop_opp, cl_id_land_opp, fun = "sum", na.rm = TRUE)
  colnames(totale_opp_per_patch_opp) <- c("Cluster_ID", "Totaal_Land_ha")
  totale_opp_per_patch_opp$Totaal_Land_ha <- totale_opp_per_patch_opp$Totaal_Land_ha * 0.01

  bos_in_land_opp <- terra::mask(bos_opp, cl_id_land_opp)
  stats_bos_in_land_opp <- terra::zonal(bos_in_land_opp, cl_id_land_opp, fun = "sum", na.rm = TRUE)
  colnames(stats_bos_in_land_opp) <- c("Cluster_ID", "Bos_ha")
  stats_bos_in_land_opp$Bos_ha <- stats_bos_in_land_opp$Bos_ha * 0.01

  df_samenstelling_land_opp <- merge(totale_opp_per_patch_opp, stats_bos_in_land_opp, by = "Cluster_ID", all.x = TRUE)
  df_samenstelling_land_opp$Bos_ha[is.na(df_samenstelling_land_opp$Bos_ha)] <- 0
  df_samenstelling_land_opp$Percentage_Bos <- round((df_samenstelling_land_opp$Bos_ha / df_samenstelling_land_opp$Totaal_Land_ha) * 100, 2)
  df_samenstelling_land_opp$Totaal_Land_ha <- round(df_samenstelling_land_opp$Totaal_Land_ha, 2)
  df_samenstelling_land_opp$Bos_ha <- round(df_samenstelling_land_opp$Bos_ha, 2)

  cat("\n=== [SPOOR OPP - HOOGVEENGLANSLIBEL] BOSPATCHES PER CLUSTER ===\n")
  print(head(df_samenstelling_land_opp, 20))
  cat("=============================================================\n\n")

  valide_land_ids_opp <- df_samenstelling_land_opp$Cluster_ID[df_samenstelling_land_opp$Percentage_Bos >= 60]

  if(length(valide_land_ids_opp) > 0) {
    m_match_opp <- terra::match(cl_id_land_opp, valide_land_ids_opp)
    winnende_masker_land_opp <- terra::ifel(!is.na(m_match_opp), 1, NA)
    land_nabij_opp <- terra::mask(landbiotoop_opp, winnende_masker_land_opp)
  } else {
    land_nabij_opp <- template_MH * NA
  }
} else {
  land_nabij_opp <- template_MH * NA
  valide_land_ids_opp <- c()
}

# BEHOUD cl_id_land_max EN cl_id_land_opp VOOR VOLGENDE CHUNK!
rm(r_land_zuiver_max, bos_binair_max, bos_in_land_max, stats_bos_in_land_max,
   r_land_zuiver_opp, r_binair_land_opp, bos_in_land_opp, stats_bos_in_land_opp)
gc()

# ==============================================================================
# WATER-FIRST LANDBIOTOOP FILTERING (HOOGVEENGLANSLIBEL - 500M AFSNIJDING)
# ==============================================================================
message("-> Landbiotoop (min. 60% bos) hard afsnijden op ", straal_water_m, "m van geschikt water...")

# --- SPOOR A: MAXIMALE POTENTIE ---
if (exists("waterbiotoop_clusters_max") && !is.null(waterbiotoop_clusters_max$raster) &&
    !all(is.na(terra::values(waterbiotoop_clusters_max$raster, mat = FALSE)))) {
  
  r_water_src_max <- waterbiotoop_clusters_max$raster
  
  # 1. Maak een vector-buffer van exact 500m rond het geschikte water
  poly_water_max  <- terra::as.polygons(r_water_src_max, aggregate = TRUE)
  poly_buffer_max <- terra::buffer(poly_water_max, width = straal_water_m) # 500m
  
  # 2. Pak het landbiotoop dat voldoet aan de 60% bos-eis (uit de vorige chunk `land_nabij_max`)
  if (exists("land_nabij_max") && !all(is.na(terra::values(land_nabij_max, mat = FALSE)))) {
    
    # 3. SNIJD DE FYSICHE PIXELS HARDAF OP EXACT 500M VAN HET VEN
    landbiotoop2_max <- terra::mask(land_nabij_max, poly_buffer_max)
    waterbiotoop2_max <- r_water_src_max
    
  } else {
    landbiotoop2_max  <- template_MH * NA
    waterbiotoop2_max <- template_MH * NA
  }
} else {
  landbiotoop2_max  <- template_MH * NA
  waterbiotoop2_max <- template_MH * NA
}


# --- SPOOR B: WERKELIJKE OPPERVLAKTE ---
if (exists("waterbiotoop_clusters_opp") && !is.null(waterbiotoop_clusters_opp$raster) &&
    !all(is.na(terra::values(waterbiotoop_clusters_opp$raster, mat = FALSE)))) {
  
  r_water_src_opp <- waterbiotoop_clusters_opp$raster
  r_bin_water_opp <- terra::ifel(!is.na(r_water_src_opp) & r_water_src_opp > 0, 1, NA)
  
  poly_water_opp  <- terra::as.polygons(r_bin_water_opp, aggregate = TRUE)
  poly_buffer_opp <- terra::buffer(poly_water_opp, width = straal_water_m) # 500m
  
  if (exists("land_nabij_opp") && !all(is.na(terra::values(land_nabij_opp, mat = FALSE)))) {
    
    # SNIJD DE WERKELIJKE FRACTIES HARDAF OP EXACT 500M VAN HET VEN
    landbiotoop2_opp <- terra::mask(land_nabij_opp, poly_buffer_opp)
    waterbiotoop2_opp <- r_water_src_opp
    
  } else {
    landbiotoop2_opp  <- template_MH * NA
    waterbiotoop2_opp <- template_MH * NA
  }
} else {
  landbiotoop2_opp  <- template_MH * NA
  waterbiotoop2_opp <- template_MH * NA
}

gc()
message("-> Water-First afsnijding op 500m voltooid voor Hoogveenglanslibel!")

# ==============================================================================
# COMBINATIE 3: RUIME OMGEVING EN FINALE SAMENSMELTING (PARALLEL)
# ==============================================================================
message("-> Ruime omgevingsfactor (50 ha complex binnen 500m) toepassen op water...")

# --- SPOOR A: MAX ---
r_ruim_land_src_max <- ruim_land_clusters_max$raster
r_ruim_land_src_max[is.na(r_ruim_land_src_max) | r_ruim_land_src_max == 0] <- NA
dist_to_ruim_land_max <- terra::distance(r_ruim_land_src_max)

hoogveenglanslibel_waterbiotoop_max <- waterbiotoop2_max == 1 & dist_to_ruim_land_max <= 500
hoogveenglanslibel_waterbiotoop_max <- terra::ifel(hoogveenglanslibel_waterbiotoop_max, 1, NA)
hoogveenglanslibel_landbiotoop_max  <- landbiotoop2_max

# Samenvoegen van gekoppeld water en GEHELE landclusters
hoogveenglanslibel_leefgebied_max <- terra::cover(hoogveenglanslibel_waterbiotoop_max, hoogveenglanslibel_landbiotoop_max)
hoogveenglanslibel_leefgebied_max <- terra::ifel(!is.na(hoogveenglanslibel_leefgebied_max), 1, NA)
hoogveenglanslibel_leefgebied_max <- terra::crop(hoogveenglanslibel_leefgebied_max, area_shape_proj, mask = TRUE)

rm(r_ruim_land_src_max, dist_to_ruim_land_max)

# --- SPOOR B: OPP ---
r_ruim_land_src_opp <- ruim_land_clusters_opp$raster
r_binair_ruim_opp   <- terra::ifel(!is.na(r_ruim_land_src_opp) & r_ruim_land_src_opp > 0, 1, NA)
dist_to_ruim_land_opp <- terra::distance(r_binair_ruim_opp)

m_water_ruim_opp                    <- dist_to_ruim_land_opp <= 500
hoogveenglanslibel_waterbiotoop_opp <- terra::mask(waterbiotoop2_opp, terra::ifel(m_water_ruim_opp, 1, NA))
hoogveenglanslibel_landbiotoop_opp  <- landbiotoop2_opp

w_clean <- terra::ifel(is.na(hoogveenglanslibel_waterbiotoop_opp), 0, hoogveenglanslibel_waterbiotoop_opp)
l_clean <- terra::ifel(is.na(hoogveenglanslibel_landbiotoop_opp), 0, hoogveenglanslibel_landbiotoop_opp)

som_fracties <- w_clean + l_clean
hoogveenglanslibel_leefgebied_opp_raw <- terra::clamp(som_fracties, upper = 1.0)
hoogveenglanslibel_leefgebied_opp     <- terra::ifel(hoogveenglanslibel_leefgebied_opp_raw > 0, hoogveenglanslibel_leefgebied_opp_raw, NA)
hoogveenglanslibel_leefgebied_opp     <- terra::crop(hoogveenglanslibel_leefgebied_opp, area_shape_proj, mask = TRUE)

final_max <- hoogveenglanslibel_leefgebied_max
final_opp <- hoogveenglanslibel_leefgebied_opp

ha_leef_pot  <- calc_ha_exact(final_max)
ha_leef_werk <- calc_ha_exact(final_opp)

cat("\n=========================================================\n")
cat("          EINDPOTENTIE LEEFGEBIED HOOGVEENGLANSLIBEL      \n")
cat("=========================================================\n")
cat(paste("Finaal Leefgebied (Maximale Potentie):   ", round(ha_leef_pot, 2), "ha\n"))
cat(paste("Finaal Leefgebied (Werkelijke Fractie): ", round(ha_leef_werk, 2), "ha\n"))
cat("=========================================================\n\n")

rm(r_ruim_land_src_opp, r_binair_ruim_opp, dist_to_ruim_land_opp, w_clean, l_clean, som_fracties, hoogveenglanslibel_leefgebied_opp_raw)
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
potentie_export_rast <- if (exists("grutto_kaart_A") && !is.null(grutto_kaart_A) && !all(is.na(terra::values(grutto_kaart_A, mat=FALSE)))) {
  terra::ifel(grutto_kaart_A > 0, 1, NA)
} else if (exists("final_max") && !all(is.na(terra::values(final_max, mat=FALSE)))) {
  terra::ifel(!is.na(final_max) & final_max > 0, 1, NA)
} else {
  terra::rast(template_MH, vals = NA)
}

# 2. Bepaal Werkelijke Oppervlakte Raster
werkelijk_export_rast <- if (exists("resB_strikt") && !is.null(resB_strikt) && (!all(is.na(terra::values(resB_strikt$kern, mat=FALSE))) || !all(is.na(terra::values(resB_strikt$bouw, mat=FALSE))))) {
  r_net_totaal_opp <- terra::cover(resB_strikt$kern, resB_strikt$bouw)
  terra::ifel(!is.na(r_net_totaal_opp) & r_net_totaal_opp > 0, 1, NA)
} else if (exists("final_opp") && !all(is.na(terra::values(final_opp, mat=FALSE)))) {
  terra::ifel(!is.na(final_opp) & final_opp > 0, 1, NA)
} else {
  terra::rast(template_MH, vals = NA)
}

# 3. Bepaal Analytisch Metacluster ID-raster
if (exists("resB_strikt") && !is.null(resB_strikt) && (!all(is.na(terra::values(resB_strikt$kern, mat=FALSE))) || !all(is.na(terra::values(resB_strikt$bouw, mat=FALSE))))) {
  r_net_totaal <- terra::cover(resB_strikt$kern, resB_strikt$bouw)
  r_buf        <- terra::buffer(r_net_totaal, width = 100)
  id_export_rast <- terra::mask(terra::patches(r_buf, directions = 8, zeroAsNA = TRUE), r_net_totaal)
} else if (exists("cl_max") && !all(is.na(terra::values(cl_max, mat=FALSE)))) {
  id_export_rast <- cl_max
} else if (exists("cl_opp") && !all(is.na(terra::values(cl_opp, mat=FALSE)))) {
  id_export_rast <- cl_opp
} else {
  id_export_rast <- terra::rast(template_MH, vals = NA)
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

