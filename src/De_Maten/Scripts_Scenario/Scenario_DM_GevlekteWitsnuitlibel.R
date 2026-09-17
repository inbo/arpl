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

# Hardgecodeerde parameters voor Gevlekte witsnuitlibel volgens de experttabel van Greet
soort                <- "gevlektewitsnuitlibel"

# Scenario pad en naam bepalen

# --- DYNAMISCHE SCENARIO PARAMETER CHECK ---
if (!exists("params") || is.null(params$scenario_rds_path)) {
  scenario_rds_path <- "data/input/Scenario_rds/DM_Scenario_BWK_2025.rds"
} else {
  scenario_rds_path <- params$scenario_rds_path
}

p_raw <- gsub("^([.][.]/)+", "", scenario_rds_path)
scenario_path <- here::here(p_raw)


if (!file.exists(scenario_path)) {
  stop(paste("❌ FOUT: Scenario RDS bestand NIET gevonden op:", scenario_path))
}

scen_volledig <- basename(scenario_path)
scenario_naam <- gsub("^DM_Scenario_|^Scenario_|.rds$", "", scen_volledig)

message(paste("Verwerken van soort:", soort, "binnen scenario:", scenario_naam))

# Parameters uit de tabel:
buffer_m             <- 20000  # Dispersiecapaciteit voor smear/buffer
straal_water_m       <- 200    # Afsnij-afstand van land rond het water (200m)

# Waterdrempels (Min. 0.02 ha EN Max. 1 ha)
drempel_water_min_m2 <- 200    # Min. 200 m² (0.02 ha)
drempel_water_max_m2 <- 10000  # Max. 10.000 m² (1.0 ha)
dist_water_m          <- 10     # 10m overbrugging voor water

# Landbiotoop drempels
drempel_land_m2      <- 50000   # Min. 5 ha land
dist_land_m          <- 10     # 10m overbrugging voor aanpalend land

drempel_ruim_m2      <- 300000  # Min. 30 ha ruim land
dist_ruim_m          <- 50     # 50m overbrugging voor ruim land

area_shape  <- vect(here("data/input/De_Maten.shp"))
master_grid <- rast(here("data/input/Raster_Vlaanderen/Vlaanderen_MasterGrid_10m.tif"))[[1]]

df_namen_sleutel <- read_csv(here("data/input/Excel_files/BWK_Laag_Namen_2025.csv"), show_col_types = FALSE)
gouden_namenlijst <- tolower(trimws(df_namen_sleutel$Laagnaam))

area_shape_proj <- project(area_shape, crs(master_grid))
area_buffer_fix <- buffer(area_shape_proj, width = buffer_m)

message("-> Uitsnede maken en lokale/globale cel-IDs berekenen...")

id_raster_DM <- crop(master_grid, area_buffer_fix, snap = "near")
lokale_coords <- terra::xyFromCell(id_raster_DM, 1:ncell(id_raster_DM))
globale_ids   <- terra::cellFromXY(master_grid, lokale_coords)

vertaal_df <- data.table(
  lokale_id  = 1:ncell(id_raster_DM),
  globale_id = globale_ids
)

studiegebied_globale_ids <- unique(vertaal_df$globale_id[!is.na(vertaal_df$globale_id)])

template_DM <- terra::rasterize(area_buffer_fix, id_raster_DM, field = 1, background = 0)
grens_web   <- sf::st_as_sf(terra::project(area_shape, "EPSG:4326"))

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
    bevat_codes_escaped <- gsub("([\\.\\^\\$\\*\\+\\?\\(\\)\\[\\{\\\\\\|])", "\\\\\\1", bevat_codes)
    bevat_codes_anchored <- paste0("^", bevat_codes_escaped)
    regex_term <- paste0(bevat_codes_anchored, collapse = "|")
    
    tabel_gefilterd <- tabel_vlaanderen[CODE %in% exact_codes | grepl(regex_term, CODE)]
  } else {
    tabel_gefilterd <- tabel_vlaanderen[CODE %in% exact_codes]
  }
  
  tabel_DM <- tabel_gefilterd[cel_id %in% studiegebied_globale_ids]
  tabel_DM_unique <- unique(tabel_DM, by = c("cel_id", "CODE"))

  tabel_cel_som <- tabel_DM_unique[, .(
    Oppervlakte = pmin(sum(BWK_FRAC, na.rm = TRUE), 1.0)
  ), by = .(cel_id)]
  
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

voortplanting_bwk_max <- lijst_matches[["voortplanting_bwk"]]
voortplanting_bwk_opp <- lijst_oppervlaktes[["voortplanting_bwk"]]
land_nabij_max        <- lijst_matches[["land_nabij"]]
land_nabij_opp        <- lijst_oppervlaktes[["land_nabij"]]
land_ruim_max         <- lijst_matches[["land_ruim"]]
land_ruim_opp         <- lijst_oppervlaktes[["land_ruim"]]
struweel1_max         <- lijst_matches[["struweel1"]]
struweel1_opp         <- lijst_oppervlaktes[["struweel1"]]

raster_simpel_final   <- id_raster_DM

rm(resultaten_gegroepeerd, df_nieuw, lijst_matches, lijst_oppervlaktes, tabel_vlaanderen)
gc()

message("-> Watervlakken shapefile lokaal inlezen & bundelen met BWK-water...")

watervlakken_v <- vect(here("data/input/ASCI Files/watervlakken2024.shp"))
watervlakken_proj <- project(watervlakken_v, crs(template_DM))

area_buffer_clean  <- makeValid(area_buffer_fix) %>% aggregate()
watervlakken_clean <- makeValid(watervlakken_proj)
watervlakken_DM    <- crop(watervlakken_clean, area_buffer_clean)

# Fuzzy clustering op kleine wateren (<= 200m²)
watervlakken_buffered <- buffer(watervlakken_DM, width = 5)
watervlakken_fuzzy_clusters <- aggregate(watervlakken_buffered, by = NULL)
watervlakken_fuzzy_clusters$fuzzy_area_m2 <- terra::expanse(watervlakken_fuzzy_clusters)

kleine_clusters_v <- watervlakken_fuzzy_clusters[watervlakken_fuzzy_clusters$fuzzy_area_m2 <= 200, ]

if (nrow(kleine_clusters_v) > 0) {
  kleine_watervlakken_finaal <- crop(watervlakken_DM, kleine_clusters_v)
  if (nrow(kleine_watervlakken_finaal) > 0) {
    r_witsnuit_kleinwater        <- terra::rasterize(kleine_watervlakken_finaal, template_DM, field = 1, background = NA)
    r_witsnuit_waterbiotoop1     <- terra::cover(voortplanting_bwk_max, r_witsnuit_kleinwater)
    r_witsnuit_waterbiotoop1_opp <- terra::cover(voortplanting_bwk_opp, r_witsnuit_kleinwater)
  } else {
    r_witsnuit_waterbiotoop1     <- voortplanting_bwk_max
    r_witsnuit_waterbiotoop1_opp <- voortplanting_bwk_opp
  }
} else {
  r_witsnuit_waterbiotoop1     <- voortplanting_bwk_max
  r_witsnuit_waterbiotoop1_opp <- voortplanting_bwk_opp
}

suppressWarnings(
  rm(watervlakken_v, watervlakken_proj, area_buffer_clean, watervlakken_clean, watervlakken_DM, 
     watervlakken_buffered, watervlakken_fuzzy_clusters, kleine_clusters_v, kleine_watervlakken_finaal, r_witsnuit_kleinwater)
)
gc()

message("-> Water filteren: 0.02 ha <= Opp <= 1 ha (10m netwerk) + 3-Buren dichtheidseis...")

# ==============================================================================
# SPOOR A: MAXIMALE POTENTIE WATER
# ==============================================================================
if (!is.null(r_witsnuit_waterbiotoop1) && !all(is.na(terra::values(r_witsnuit_waterbiotoop1, mat = FALSE)))) {
  
  # A. Vorm waterclusters via 10m overbrugging
  res_water_max <- cluster_filter_compleet(
    masker     = r_witsnuit_waterbiotoop1, 
    opp_laag   = r_witsnuit_waterbiotoop1, 
    drempel_m2 = drempel_water_min_m2, # Min 200 m² (0.02 ha)
    dist_m     = dist_water_m,          # 10m overbrugging
    werkelijk  = FALSE
  )
  r_water_cand_max  <- res_water_max$raster
  cl_water_cand_max <- res_water_max$clusters
  
  # B. MAXIMALE BOVENGRENS VAN 1 HA (10.000 m²)
  if (!all(is.na(terra::values(cl_water_cand_max, mat = FALSE)))) {
    f_water_max <- terra::freq(cl_water_cand_max)
    f_water_max$Area_m2 <- f_water_max$count * 100
    
    valide_opp_ids_max <- f_water_max$value[f_water_max$Area_m2 <= drempel_water_max_m2] # <= 1 ha
    
    if (length(valide_opp_ids_max) > 0) {
      r_water_cand_max <- terra::mask(r_water_cand_max, terra::ifel(cl_water_cand_max %in% valide_opp_ids_max, 1, NA))
    } else { r_water_cand_max <- template_DM * NA }
  } else { r_water_cand_max <- template_DM * NA }
  
  # C. ORIGINELE 3-BUREN CHECK (Minstens 3 naburige vennen binnen 200m)
  if (!all(is.na(terra::values(r_water_cand_max, mat = FALSE)))) {
    r_water_cl_max <- terra::patches(r_water_cand_max, directions = 8, zeroAsNA = TRUE)
    
    if (!all(is.na(terra::values(r_water_cl_max, mat = FALSE)))) {
      watervlakken_poly <- terra::as.polygons(r_water_cl_max, dissolve = TRUE)
      
      if (!is.null(watervlakken_poly) && nrow(watervlakken_poly) > 0) {
        centroids <- terra::centroids(watervlakken_poly)
        centroids_buffer <- terra::buffer(centroids, width = 200)
        intersect_matrix <- terra::relate(centroids_buffer, centroids, "intersects")
        cluster_counts   <- rowSums(intersect_matrix)
        
        watervlakken_poly$unieke_buren_count <- cluster_counts
        goedgekeurde_clusters_poly <- watervlakken_poly[watervlakken_poly$unieke_buren_count >= 4, ] # Self + 3 buren
        
        if (nrow(goedgekeurde_clusters_poly) > 0) {
          gevlektewitsnuitlibel_voortplanting_max <- terra::rasterize(goedgekeurde_clusters_poly, template_DM, field = 1, background = NA)
          gevlektewitsnuitlibel_water_finaal_max  <- gevlektewitsnuitlibel_voortplanting_max
        } else {
          gevlektewitsnuitlibel_voortplanting_max <- template_DM * NA
          gevlektewitsnuitlibel_water_finaal_max  <- template_DM * NA
        }
      } else {
        gevlektewitsnuitlibel_voortplanting_max <- template_DM * NA
        gevlektewitsnuitlibel_water_finaal_max  <- template_DM * NA
      }
    } else {
      gevlektewitsnuitlibel_voortplanting_max <- template_DM * NA
      gevlektewitsnuitlibel_water_finaal_max  <- template_DM * NA
    }
  } else {
    gevlektewitsnuitlibel_voortplanting_max <- template_DM * NA
    gevlektewitsnuitlibel_water_finaal_max  <- template_DM * NA
  }

} else {
  gevlektewitsnuitlibel_voortplanting_max <- template_DM * NA
  gevlektewitsnuitlibel_water_finaal_max  <- template_DM * NA
}


# ==============================================================================
# SPOOR B: WERKELIJKE OPPERVLAKTE WATER
# ==============================================================================
if (!is.null(r_witsnuit_waterbiotoop1_opp) && !all(is.na(terra::values(r_witsnuit_waterbiotoop1_opp, mat = FALSE)))) {
  r_binair_water_opp <- terra::ifel(!is.na(r_witsnuit_waterbiotoop1_opp) & r_witsnuit_waterbiotoop1_opp > 0, 1, NA)
  
  res_water_opp <- cluster_filter_compleet(
    masker     = r_binair_water_opp, 
    opp_laag   = r_witsnuit_waterbiotoop1_opp, 
    drempel_m2 = drempel_water_min_m2, 
    dist_m     = dist_water_m, # 10m overbrugging
    werkelijk  = TRUE
  )
  r_water_cand_opp  <- res_water_opp$raster
  cl_water_cand_opp <- res_water_opp$clusters
  
  # BOVENGRENS 1 HA TOEPASSEN OP WERKELIJKE OPPERVLAKTE
  if (!all(is.na(terra::values(cl_water_cand_opp, mat = FALSE)))) {
    stats_water_opp <- terra::zonal(r_water_cand_opp, cl_water_cand_opp, fun = "sum", na.rm = TRUE)
    colnames(stats_water_opp) <- c("ID", "Waarde")
    stats_water_opp$Area_m2 <- stats_water_opp$Waarde * 100
    
    valide_opp_ids_opp <- stats_water_opp$ID[stats_water_opp$Area_m2 <= drempel_water_max_m2]
    
    if (length(valide_opp_ids_opp) > 0) {
      r_water_cand_opp <- terra::mask(r_water_cand_opp, terra::ifel(cl_water_cand_opp %in% valide_opp_ids_opp, 1, NA))
    } else { r_water_cand_opp <- template_DM * NA }
  } else { r_water_cand_opp <- template_DM * NA }
  
  if (!all(is.na(terra::values(gevlektewitsnuitlibel_voortplanting_max, mat = FALSE)))) {
    gevlektewitsnuitlibel_voortplanting_opp <- terra::mask(r_water_cand_opp, gevlektewitsnuitlibel_voortplanting_max)
    gevlektewitsnuitlibel_water_finaal_opp <- gevlektewitsnuitlibel_voortplanting_opp
  } else {
    gevlektewitsnuitlibel_voortplanting_opp <- template_DM * NA
    gevlektewitsnuitlibel_water_finaal_opp <- template_DM * NA
  }
} else {
  gevlektewitsnuitlibel_voortplanting_opp <- template_DM * NA
  gevlektewitsnuitlibel_water_finaal_opp <- template_DM * NA
}

rm(r_water_clusters, r_witsnuit_waterbiotoop1, r_witsnuit_waterbiotoop1_opp)
gc()

# ==============================================================================
# WATER-FIRST LANDBIOTOOP FILTERING (GEVLEKTE WITSNUITLIBEL)
# ==============================================================================
message("-> Landbiotoop filteren (10m & 50m netwerken, 30% bos) & hard afsnijden op ", straal_water_m, "m van water...")

fix_bos_raster <- function(r_bos, template) {
  if (is.null(r_bos) || !inherits(r_bos, "SpatRaster")) {
    return(terra::rast(template, vals = NA))
  }
  return(terra::ifel(!is.na(r_bos) & r_bos > 0, 1, NA))
}

# --- STAP 1: RUIM LANDBIOTOOP FILTEREN (MIN. 30 HA MET 50M OVERBRUGGING) ---
if (exists("land_ruim_max") && !all(is.na(terra::values(land_ruim_max, mat = FALSE)))) {
  res_ruim_max <- cluster_filter_compleet(
    masker     = land_ruim_max,
    opp_laag   = land_ruim_max,
    drempel_m2 = drempel_ruim_m2, # 30 ha
    dist_m     = dist_ruim_m,     # 50m overbrugging
    werkelijk  = FALSE
  )
  land_ruim_ok_max <- res_ruim_max$raster
} else { land_ruim_ok_max <- template_DM * NA }

if (exists("land_ruim_opp") && !all(is.na(terra::values(land_ruim_opp, mat = FALSE)))) {
  r_bin_ruim_opp <- terra::ifel(!is.na(land_ruim_opp) & land_ruim_opp > 0, 1, NA)
  res_ruim_opp <- cluster_filter_compleet(
    masker     = r_bin_ruim_opp,
    opp_laag   = land_ruim_opp,
    drempel_m2 = drempel_ruim_m2, # 30 ha
    dist_m     = dist_ruim_m,     # 50m overbrugging
    werkelijk  = TRUE
  )
  land_ruim_ok_opp <- res_ruim_opp$raster
} else { land_ruim_ok_opp <- template_DM * NA }


# --- STAP 2: SPOOR A - MAXIMALE POTENTIE (200M WATER AFSNIJDING) ---
if (exists("gevlektewitsnuitlibel_voortplanting_max") && !all(is.na(terra::values(gevlektewitsnuitlibel_voortplanting_max, mat = FALSE)))) {
  
  poly_water_max  <- terra::as.polygons(gevlektewitsnuitlibel_voortplanting_max, aggregate = TRUE)
  poly_buffer_max <- terra::buffer(poly_water_max, width = straal_water_m) # 200m!
  
  if (exists("land_nabij_max") && !all(is.na(terra::values(land_nabij_max, mat = FALSE)))) {
    
    # 1. HARDE WATER-FIRST AFSNIJDING OP 200M
    land_afgesneden_max <- terra::mask(land_nabij_max, poly_buffer_max)
    
    # 2. Check 5 ha op het AFGESNEDEN land (met 10m overbrugging)
    res_land_5ha_max <- cluster_filter_compleet(
      masker     = land_afgesneden_max, 
      opp_laag   = land_afgesneden_max, 
      drempel_m2 = drempel_land_m2, # 5 ha
      dist_m     = dist_land_m,     # 10m overbrugging
      werkelijk  = FALSE
    )
    land_5ha_max  <- res_land_5ha_max$raster
    cl_id_5ha_max <- res_land_5ha_max$clusters
    
    # 3. ORIGINELE BOS/STRUWEEL CHECK (Minstens 30% bos/struweel)
    if (!all(is.na(terra::values(land_5ha_max, mat = FALSE)))) {
      r_struweel_max_clean <- fix_bos_raster(struweel1_max, template_DM)
      stats_totaal_max <- terra::zonal(land_5ha_max * 0.01, cl_id_5ha_max, fun = "sum", na.rm = TRUE)
      colnames(stats_totaal_max) <- c("Cluster_ID", "Totaal_ha")
      struweel_in_cl_max <- terra::mask(r_struweel_max_clean * 0.01, cl_id_5ha_max)
      stats_struweel_max <- terra::zonal(struweel_in_cl_max, cl_id_5ha_max, fun = "sum", na.rm = TRUE)
      colnames(stats_struweel_max) <- c("Cluster_ID", "Struweel_ha")
      
      df_check_max <- merge(stats_totaal_max, stats_struweel_max, by = "Cluster_ID", all.x = TRUE)
      df_check_max$Struweel_ha[is.na(df_check_max$Struweel_ha)] <- 0
      df_check_max$Pct_Bos <- (df_check_max$Struweel_ha / df_check_max$Totaal_ha) * 100
      valide_ids_max <- df_check_max$Cluster_ID[df_check_max$Pct_Bos >= 30] # 30% bos/struweel
      
      land_ok_max <- if(length(valide_ids_max) > 0) terra::mask(land_5ha_max, terra::ifel(cl_id_5ha_max %in% valide_ids_max, 1, NA)) else template_DM * NA
    } else { land_ok_max <- template_DM * NA }
    
    # 4. Koppel aan 30ha ruim landbiotoop
    if (!all(is.na(terra::values(land_ok_max, mat = FALSE))) && !all(is.na(terra::values(land_ruim_ok_max, mat = FALSE)))) {
      gevlektewitsnuitlibel_land_finaal_max <- terra::mask(land_ok_max, terra::ifel(!is.na(land_ruim_ok_max) & land_ruim_ok_max > 0, 1, NA))
    } else { gevlektewitsnuitlibel_land_finaal_max <- land_ok_max }
    
  } else { gevlektewitsnuitlibel_land_finaal_max <- template_DM * NA }
} else { gevlektewitsnuitlibel_land_finaal_max <- template_DM * NA }

water_finaal_max  <- if(exists("gevlektewitsnuitlibel_water_finaal_max")) gevlektewitsnuitlibel_water_finaal_max else template_DM * NA
leefgebied_max    <- terra::cover(water_finaal_max, gevlektewitsnuitlibel_land_finaal_max)
finaal_max_binair <- terra::ifel(!is.na(leefgebied_max) & leefgebied_max > 0, 1, NA)
cl_id_max         <- terra::patches(finaal_max_binair, directions = 8, zeroAsNA = TRUE)


# --- STAP 3: SPOOR B - WERKELIJKE OPPERVLAKTE (200M WATER AFSNIJDING) ---
if (exists("gevlektewitsnuitlibel_voortplanting_opp") && !all(is.na(terra::values(gevlektewitsnuitlibel_voortplanting_opp, mat = FALSE)))) {
  r_water_bin_opp <- terra::ifel(!is.na(gevlektewitsnuitlibel_voortplanting_opp) & gevlektewitsnuitlibel_voortplanting_opp > 0, 1, NA)
  poly_water_opp  <- terra::as.polygons(r_water_bin_opp, aggregate = TRUE)
  poly_buffer_opp <- terra::buffer(poly_water_opp, width = straal_water_m)
  
  if (exists("land_nabij_opp") && !all(is.na(terra::values(land_nabij_opp, mat = FALSE)))) {
    
    # 1. HARDE WATER-FIRST AFSNIJDING OP 200M
    land_afgesneden_opp <- terra::mask(land_nabij_opp, poly_buffer_opp)
    r_bin_land_afgesneden <- terra::ifel(!is.na(land_afgesneden_opp) & land_afgesneden_opp > 0, 1, NA)
    
    # 2. Check 5 ha op het AFGESNEDEN land (met 10m overbrugging)
    res_land_5ha_opp <- cluster_filter_compleet(
      masker     = r_bin_land_afgesneden, 
      opp_laag   = land_afgesneden_opp, 
      drempel_m2 = drempel_land_m2, # 5 ha
      dist_m     = dist_land_m,     # 10m overbrugging
      werkelijk  = TRUE
    )
    land_5ha_opp  <- res_land_5ha_opp$raster
    cl_id_5ha_opp <- res_land_5ha_opp$clusters
    
    # 3. ORIGINELE BOS/STRUWEEL CHECK (Minstens 30% bos/struweel)
    if (!all(is.na(terra::values(land_5ha_opp, mat = FALSE)))) {
      r_struweel_opp_clean <- fix_bos_raster(struweel1_opp, template_DM)
      stats_totaal_opp <- terra::zonal(land_5ha_opp, cl_id_5ha_opp, fun = "sum", na.rm = TRUE)
      colnames(stats_totaal_opp) <- c("Cluster_ID", "Totaal_m2")
      struweel_in_cl_opp <- terra::mask(r_struweel_opp_clean, cl_id_5ha_opp)
      stats_struweel_opp <- terra::zonal(struweel_in_cl_opp, cl_id_5ha_opp, fun = "sum", na.rm = TRUE)
      colnames(stats_struweel_opp) <- c("Cluster_ID", "Struweel_m2")
      
      df_check_opp <- merge(stats_totaal_opp, stats_struweel_opp, by = "Cluster_ID", all.x = TRUE)
      df_check_opp$Struweel_m2[is.na(df_check_opp$Struweel_m2)] <- 0
      df_check_opp$Pct_Bos <- (df_check_opp$Struweel_m2 / df_check_opp$Totaal_m2) * 100
      valide_ids_opp <- df_check_opp$Cluster_ID[df_check_opp$Pct_Bos >= 30] # 30% bos/struweel
      
      land_ok_opp <- if(length(valide_ids_opp) > 0) terra::mask(land_5ha_opp, terra::ifel(cl_id_5ha_opp %in% valide_ids_opp, 1, NA)) else template_DM * NA
    } else { land_ok_opp <- template_DM * NA }
    
    # 4. Koppel aan 30ha ruim landbiotoop
    if (!all(is.na(terra::values(land_ok_opp, mat = FALSE))) && !all(is.na(terra::values(land_ruim_ok_opp, mat = FALSE)))) {
      gevlektewitsnuitlibel_land_finaal_opp <- terra::mask(land_ok_opp, terra::ifel(!is.na(land_ruim_ok_opp) & land_ruim_ok_opp > 0, 1, NA))
    } else { gevlektewitsnuitlibel_land_finaal_opp <- land_ok_opp }
    
  } else { gevlektewitsnuitlibel_land_finaal_opp <- template_DM * NA }
} else { gevlektewitsnuitlibel_land_finaal_opp <- template_DM * NA }

water_finaal_opp  <- if(exists("gevlektewitsnuitlibel_voortplanting_opp")) gevlektewitsnuitlibel_voortplanting_opp else template_DM * NA
leefgebied_opp    <- terra::cover(water_finaal_opp, gevlektewitsnuitlibel_land_finaal_opp)
finaal_opp_binair <- terra::ifel(!is.na(leefgebied_opp) & leefgebied_opp > 0, 1, NA)
cl_id_opp         <- terra::patches(finaal_opp_binair, directions = 8, zeroAsNA = TRUE)

message("-> Water-First filtering met alle originele bronnen/checks succesvol doorlopen!")


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
potentie_export_rast <- if (exists("grutto_kaart_A") && !is.null(grutto_kaart_A) && !all(is.na(terra::values(grutto_kaart_A, mat=FALSE)))) {
  terra::ifel(grutto_kaart_A > 0, 1, NA)
} else if (exists("final_max") && !all(is.na(terra::values(final_max, mat=FALSE)))) {
  terra::ifel(!is.na(final_max) & final_max > 0, 1, NA)
} else {
  terra::rast(template_DM, vals = NA)
}

# 2. Bepaal Werkelijke Oppervlakte Raster
werkelijk_export_rast <- if (exists("resB_strikt") && !is.null(resB_strikt) && (!all(is.na(terra::values(resB_strikt$kern, mat=FALSE))) || !all(is.na(terra::values(resB_strikt$bouw, mat=FALSE))))) {
  r_net_totaal_opp <- terra::cover(resB_strikt$kern, resB_strikt$bouw)
  terra::ifel(!is.na(r_net_totaal_opp) & r_net_totaal_opp > 0, 1, NA)
} else if (exists("final_opp") && !all(is.na(terra::values(final_opp, mat=FALSE)))) {
  terra::ifel(!is.na(final_opp) & final_opp > 0, 1, NA)
} else {
  terra::rast(template_DM, vals = NA)
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
  id_export_rast <- terra::rast(template_DM, vals = NA)
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

