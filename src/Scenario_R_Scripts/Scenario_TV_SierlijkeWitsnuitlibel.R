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

soort <- "sierlijkewitsnuitlibel"

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

buffer_m       <- resultaat$Dispersiecap_m[1]
straal_water_m <- 800  # Harde afsnij-afstand van land rond het water (800m)

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

voortplanting_bwk_max        <- lijst_matches[["voortplanting_bwk"]]
voortplanting_bwk_opp        <- lijst_oppervlaktes[["voortplanting_bwk"]]
foerageer0_max                 <- lijst_matches[["foerageer0"]]
foerageer0_notbwk_opp         <- lijst_oppervlaktes[["foerageer0"]]
struweel1_max                  <- lijst_matches[["struweel1"]]
struweel1_opp                  <- lijst_oppervlaktes[["struweel1"]]

rm(tabel_vlaanderen, vertaal_df, lijst_matches, lijst_oppervlaktes, df_nieuw, resultaten_gegroepeerd)
gc()

message("-> Watervlakken shapefile lokaal inlezen en opschonen...")

watervlakken_path <- here("data/input/ASCI Files/watervlakken2024.shp")

if (file.exists(watervlakken_path)) {
  watervlakken_v    <- vect(watervlakken_path)
  watervlakken_proj <- project(watervlakken_v, crs(template_TV))
  
  area_buffer_clean  <- aggregate(makeValid(area_buffer_fix))
  watervlakken_clean <- makeValid(watervlakken_proj)
  watervlakken_TV    <- crop(watervlakken_clean, area_buffer_clean)
  
  watervlakken_buffered       <- buffer(watervlakken_TV, width = 5)
  watervlakken_fuzzy_clusters <- aggregate(watervlakken_buffered, by = NULL)
  watervlakken_fuzzy_clusters$fuzzy_area_m2 <- terra::expanse(watervlakken_fuzzy_clusters)
  
  kleine_clusters_v <- watervlakken_fuzzy_clusters[watervlakken_fuzzy_clusters$fuzzy_area_m2 <= 4000, ]
  
  if (!is.null(kleine_clusters_v) && nrow(kleine_clusters_v) > 0) {
    kleine_watervlakken_finaal <- crop(watervlakken_TV, kleine_clusters_v)
    
    if (!is.null(kleine_watervlakken_finaal) && nrow(kleine_watervlakken_finaal) > 0) {
      r_kleinwater <- terra::rasterize(kleine_watervlakken_finaal, template_TV, field = 1, background = NA)
      
      r_witsnuit_waterbiotoop1     <- terra::cover(voortplanting_bwk_max, r_kleinwater)
      r_witsnuit_waterbiotoop1_opp <- terra::cover(voortplanting_bwk_opp, r_kleinwater)
    } else {
      r_witsnuit_waterbiotoop1     <- voortplanting_bwk_max
      r_witsnuit_waterbiotoop1_opp <- voortplanting_bwk_opp
    }
  } else {
    r_witsnuit_waterbiotoop1     <- voortplanting_bwk_max
    r_witsnuit_waterbiotoop1_opp <- voortplanting_bwk_opp
  }
} else {
  r_witsnuit_waterbiotoop1     <- voortplanting_bwk_max
  r_witsnuit_waterbiotoop1_opp <- voortplanting_bwk_opp
}

suppressWarnings(
  rm(watervlakken_v, watervlakken_proj, area_buffer_clean, watervlakken_clean, watervlakken_TV, 
     watervlakken_buffered, watervlakken_fuzzy_clusters, kleine_clusters_v, kleine_watervlakken_finaal, r_kleinwater)
)
gc()

message("=== STAP 2: Voortplantingsbiotoop filteren op netwerkconnectiviteit (3 buren binnen 500m - PARALLEL) ===")

# --- SPOOR A: MAX ---
r_water_clusters_max <- terra::patches(r_witsnuit_waterbiotoop1, directions = 8, zeroAsNA = TRUE)
witsnuit_water_finaal_max <- template_TV * NA

if (!all(is.na(suppressWarnings(terra::minmax(r_water_clusters_max))))) {
  r_poel_kernen_max <- r_water_clusters_max * NA
  df_cells_max <- terra::as.data.frame(r_water_clusters_max, cells = TRUE)
  colnames(df_cells_max) <- c("cell", "ID")
  
  df_kernen_max <- df_cells_max %>% dplyr::group_by(ID) %>% dplyr::summarise(cell = dplyr::first(cell))
  r_poel_kernen_max[df_kernen_max$cell] <- df_kernen_max$ID
  
  xy_kernen_max <- terra::xyFromCell(r_poel_kernen_max, df_kernen_max$cell)
  n_plassen_max <- nrow(xy_kernen_max)
  
  if (n_plassen_max >= 4) {
    dist_matrix_max <- stats::dist(xy_kernen_max, method = "euclidean") %>% as.matrix()
    cluster_counts_max <- rowSums(dist_matrix_max <= 500)
    
    df_scores_max <- data.table(ID = df_kernen_max$ID, Buren_Count = cluster_counts_max)
    goedgekeurde_ids_max <- df_scores_max[Buren_Count >= 4, ID]
    
    if (length(goedgekeurde_ids_max) > 0) {
      masker_water_bin_max <- r_water_clusters_max %in% goedgekeurde_ids_max
      witsnuit_water_finaal_max <- terra::ifel(masker_water_bin_max == 1, 1, NA)
    }
    rm(dist_matrix_max, df_scores_max, goedgekeurde_ids_max, cluster_counts_max)
  }
  rm(r_poel_kernen_max, df_cells_max, df_kernen_max, xy_kernen_max)
}

cl_water_ids_gefilterd_max <- terra::mask(r_water_clusters_max, witsnuit_water_finaal_max)

# --- SPOOR B: OPP ---
r_binair_water_opp <- terra::ifel(!is.na(r_witsnuit_waterbiotoop1_opp) & r_witsnuit_waterbiotoop1_opp > 0, 1, NA)
r_water_clusters_opp <- terra::patches(r_binair_water_opp, directions = 8, zeroAsNA = TRUE)
witsnuit_water_finaal_opp <- template_TV * NA

if (!all(is.na(suppressWarnings(terra::minmax(r_water_clusters_opp))))) {
  r_poel_kernen_opp <- r_water_clusters_opp * NA
  df_cells_opp <- terra::as.data.frame(r_water_clusters_opp, cells = TRUE)
  colnames(df_cells_opp) <- c("cell", "ID")
  
  df_kernen_opp <- df_cells_opp %>% dplyr::group_by(ID) %>% dplyr::summarise(cell = dplyr::first(cell))
  r_poel_kernen_opp[df_kernen_opp$cell] <- df_kernen_opp$ID
  
  xy_kernen_opp <- terra::xyFromCell(r_poel_kernen_opp, df_kernen_opp$cell)
  n_plassen_opp <- nrow(xy_kernen_opp)
  
  if (n_plassen_opp >= 4) {
    dist_matrix_opp <- stats::dist(xy_kernen_opp, method = "euclidean") %>% as.matrix()
    cluster_counts_opp <- rowSums(dist_matrix_opp <= 500)
    
    df_scores_opp <- data.table(ID = df_kernen_opp$ID, Buren_Count = cluster_counts_opp)
    goedgekeurde_ids_opp <- df_scores_opp[Buren_Count >= 4, ID]
    
    if (length(goedgekeurde_ids_opp) > 0) {
      masker_water_bin_opp <- r_water_clusters_opp %in% goedgekeurde_ids_opp
      witsnuit_water_finaal_opp <- terra::mask(r_witsnuit_waterbiotoop1_opp, terra::ifel(masker_water_bin_opp, 1, NA))
    }
    rm(dist_matrix_opp, df_scores_opp, goedgekeurde_ids_opp, cluster_counts_opp)
  }
  rm(r_poel_kernen_opp, df_cells_opp, df_kernen_opp, xy_kernen_opp)
}

cl_water_ids_gefilterd_opp <- terra::mask(r_water_clusters_opp, terra::ifel(!is.na(witsnuit_water_finaal_opp) & witsnuit_water_finaal_opp > 0, 1, NA))

rm(r_water_clusters_max, r_water_clusters_opp, r_binair_water_opp)
gc()

message("=== STAP 3: Foerageergebieden filteren en struweel-aandeel wegen (Netto 70% Check - PARALLEL) ===")

if (!is.null(foerageer0_max) && !all(is.na(suppressWarnings(terra::minmax(foerageer0_max))))) {
  foerageer_buffered <- terra::buffer(foerageer0_max, width = 25)
  cl_foerageer <- terra::patches(foerageer_buffered, directions = 8, zeroAsNA = TRUE) %>% terra::mask(foerageer0_max)
} else {
  cl_foerageer <- template_TV * NA
}

witsnuit_foerageer_finaal_max <- template_TV * NA
witsnuit_foerageer_finaal_opp <- template_TV * NA

if (!all(is.na(suppressWarnings(terra::minmax(cl_foerageer))))) {
  f_stats <- terra::freq(cl_foerageer)
  df_stats <- data.frame(ID = f_stats$value, Foerageer_ha = f_stats$count * 0.01)
  
  w_50m <- terra::focalMat(cl_foerageer, 50, type = 'circle')
  w_50m[w_50m > 0] <- 1
  cl_foerageer_buf <- terra::focal(cl_foerageer, w = w_50m, fun = "max", na.rm = TRUE)
  
  struweel_in_zone_opp <- terra::mask(struweel1_opp, cl_foerageer_buf)
  zonal_struweel_opp  <- terra::zonal(struweel_in_zone_opp, cl_foerageer_buf, fun = "sum", na.rm = TRUE)
  colnames(zonal_struweel_opp) <- c("ID", "Struweel_Fractie_Sum")
  zonal_struweel_opp$Struweel_ha <- zonal_struweel_opp$Struweel_Fractie_Sum * 0.01
  
  df_criteria <- merge(df_stats, zonal_struweel_opp, by = "ID", all.x = TRUE)
  df_criteria[is.na(df_criteria)] <- 0
  df_criteria$Percentage_Struweel <- (df_criteria$Struweel_ha / df_criteria$Foerageer_ha) * 100
  
  goedgekeurde_foer_ids <- df_criteria$ID[df_criteria$Foerageer_ha >= 50 & df_criteria$Percentage_Struweel >= 70]
  
  if (length(goedgekeurde_foer_ids) > 0) {
    masker_foer_bin <- cl_foerageer %in% goedgekeurde_foer_ids
    witsnuit_foerageer_finaal_max <- terra::mask(struweel1_max, terra::ifel(masker_foer_bin == 1, 1, NA))
    witsnuit_foerageer_finaal_opp <- terra::mask(struweel1_opp, witsnuit_foerageer_finaal_max)
  }
  
  rm(f_stats, df_stats, cl_foerageer_buf, struweel_in_zone_opp, zonal_struweel_opp, df_criteria, goedgekeurde_foer_ids, w_50m)
}

rm(cl_foerageer, foerageer0_max, foerageer0_notbwk_opp, struweel1_max, struweel1_opp)
gc()

message("=== STAP 4: Ruimtelijke relatie tussen water en foerageergebied controleren (800m - PARALLEL) ===")

witsnuit_water_relatie_max <- template_TV * NA
witsnuit_foer_relatie_max  <- template_TV * NA

# --- SPOOR A: MAX ---
if (!all(is.na(suppressWarnings(terra::minmax(witsnuit_water_finaal_max)))) && 
    !all(is.na(suppressWarnings(terra::minmax(witsnuit_foerageer_finaal_max))))) {
  
  water_buffer_800m_max <- terra::distance(witsnuit_water_finaal_max) <= 800
  foer_buffer_800m_max  <- terra::distance(witsnuit_foerageer_finaal_max) <= 800
  
  witsnuit_water_relatie_max <- terra::mask(witsnuit_water_finaal_max, foer_buffer_800m_max, maskvalues = FALSE)
  witsnuit_foer_relatie_max  <- terra::mask(witsnuit_foerageer_finaal_max, water_buffer_800m_max, maskvalues = FALSE)
  
  rm(water_buffer_800m_max, foer_buffer_800m_max)
}

# --- SPOOR B: OPP ---
witsnuit_water_relatie_opp <- template_TV * NA
witsnuit_foer_relatie_opp  <- template_TV * NA

if (!all(is.na(suppressWarnings(terra::minmax(witsnuit_water_finaal_opp)))) && 
    !all(is.na(suppressWarnings(terra::minmax(witsnuit_foerageer_finaal_opp))))) {
  
  r_binair_water_finaal_opp <- terra::ifel(!is.na(witsnuit_water_finaal_opp) & witsnuit_water_finaal_opp > 0, 1, NA)
  r_binair_foer_finaal_opp  <- terra::ifel(!is.na(witsnuit_foerageer_finaal_opp) & witsnuit_foerageer_finaal_opp > 0, 1, NA)
  
  water_buffer_800m_opp <- terra::distance(r_binair_water_finaal_opp) <= 800
  foer_buffer_800m_opp  <- terra::distance(r_binair_foer_finaal_opp) <= 800
  
  witsnuit_water_relatie_opp <- terra::mask(witsnuit_water_finaal_opp, foer_buffer_800m_opp, maskvalues = FALSE)
  witsnuit_foer_relatie_opp  <- terra::mask(witsnuit_foerageer_finaal_opp, water_buffer_800m_opp, maskvalues = FALSE)
  
  rm(r_binair_water_finaal_opp, r_binair_foer_finaal_opp, water_buffer_800m_opp, foer_buffer_800m_opp)
}

rm(witsnuit_water_finaal_max, witsnuit_water_finaal_opp, witsnuit_foerageer_finaal_max, witsnuit_foerageer_finaal_opp)
gc()

message("=== STAP 5: Finaal leefgebied complex samenstellen (35 ha) + Water-First 800m afsnijding ===")

straal_water_m <- 800

# ==============================================================================
# SPOOR A: MAXIMALE POTENTIE
# ==============================================================================
witsnuit_leefgebied1_max <- witsnuit_water_relatie_max | witsnuit_foer_relatie_max
witsnuit_leefgebied1_max <- terra::ifel(witsnuit_leefgebied1_max == 1, 1, NA)

if (!all(is.na(suppressWarnings(terra::minmax(witsnuit_leefgebied1_max))))) {
  leefgebied_buffered_max   <- terra::buffer(witsnuit_leefgebied1_max, width = 25)
  r_leefgebied_clusters_max <- terra::patches(leefgebied_buffered_max, directions = 8)
  r_scherpe_clusters_max    <- terra::mask(r_leefgebied_clusters_max, witsnuit_leefgebied1_max)
  
  cluster_freq_max    <- terra::freq(r_scherpe_clusters_max)
  cluster_freq_max$area_ha <- cluster_freq_max$count * 0.01
  goedgekeurde_ids_max <- cluster_freq_max$value[cluster_freq_max$area_ha >= 35]
  
  if (length(goedgekeurde_ids_max) > 0) {
    leefgebied_35ha_max <- r_scherpe_clusters_max %in% goedgekeurde_ids_max
    leefgebied_35ha_max <- terra::ifel(leefgebied_35ha_max == 1, 1, NA)
    
    water_finaal_max <- terra::mask(witsnuit_water_relatie_max, leefgebied_35ha_max)
    land_finaal_max_ruw <- terra::mask(witsnuit_foer_relatie_max, leefgebied_35ha_max)
  } else {
    water_finaal_max    <- template_TV * NA
    land_finaal_max_ruw <- template_TV * NA
  }
} else {
  water_finaal_max    <- template_TV * NA
  land_finaal_max_ruw <- template_TV * NA
}

if (!all(is.na(suppressWarnings(terra::minmax(water_finaal_max))))) {
  poly_water_max  <- terra::as.polygons(water_finaal_max, aggregate = TRUE)
  poly_buffer_max <- terra::buffer(poly_water_max, width = straal_water_m)
  land_finaal_max <- terra::mask(land_finaal_max_ruw, poly_buffer_max)
} else {
  land_finaal_max <- land_finaal_max_ruw
}

waterbiotoop_finaal_max <- water_finaal_max
landbiotoop_finaal_max  <- land_finaal_max

sierlijkewitsnuitlibel_finaal_max <- terra::cover(waterbiotoop_finaal_max, landbiotoop_finaal_max)

# ==============================================================================
# SPOOR B: WERKELIJKE OPPERVLAKTE (PARALLEL)
# ==============================================================================
w_clean <- terra::ifel(is.na(witsnuit_water_relatie_opp), 0, witsnuit_water_relatie_opp)
f_clean <- terra::ifel(is.na(witsnuit_foer_relatie_opp), 0, witsnuit_foer_relatie_opp)
som_fracties   <- w_clean + f_clean
som_gemedieerd <- terra::clamp(som_fracties, upper = 1.0)

witsnuit_leefgebied1_opp_raw <- terra::ifel(som_gemedieerd > 0, som_gemedieerd, NA)
witsnuit_leefgebied1_opp     <- terra::crop(witsnuit_leefgebied1_opp_raw, template_TV)

r_binair_leef1_opp <- terra::ifel(!is.na(witsnuit_leefgebied1_opp) & witsnuit_leefgebied1_opp > 0, 1, NA)

if (!all(is.na(suppressWarnings(terra::minmax(r_binair_leef1_opp))))) {
  leefgebied_buffered_opp   <- terra::buffer(r_binair_leef1_opp, width = 25)
  r_leefgebied_clusters_opp <- terra::patches(leefgebied_buffered_opp, directions = 8)
  r_scherpe_clusters_opp    <- terra::mask(r_leefgebied_clusters_opp, r_binair_leef1_opp)
  
  cluster_stats_opp <- terra::zonal(witsnuit_leefgebied1_opp, r_scherpe_clusters_opp, fun = "sum", na.rm = TRUE)
  colnames(cluster_stats_opp) <- c("value", "ha_exact")
  goedgekeurde_ids_opp <- cluster_stats_opp$value[cluster_stats_opp$ha_exact >= 35]
  
  if (length(goedgekeurde_ids_opp) > 0) {
    masker_opp_leef <- r_scherpe_clusters_opp %in% goedgekeurde_ids_opp
    leef_35ha_mask <- terra::ifel(masker_opp_leef, 1, NA)
    
    water_finaal_opp    <- terra::mask(witsnuit_water_relatie_opp, leef_35ha_mask)
    land_finaal_opp_ruw <- terra::mask(witsnuit_foer_relatie_opp, leef_35ha_mask)
  } else {
    water_finaal_opp    <- template_TV * NA
    land_finaal_opp_ruw <- template_TV * NA
  }
} else {
  water_finaal_opp    <- template_TV * NA
  land_finaal_opp_ruw <- template_TV * NA
}

if (!all(is.na(suppressWarnings(terra::minmax(water_finaal_opp))))) {
  r_bin_water_opp <- terra::ifel(!is.na(water_finaal_opp) & water_finaal_opp > 0, 1, NA)
  poly_water_opp  <- terra::as.polygons(r_bin_water_opp, aggregate = TRUE)
  poly_buffer_opp <- terra::buffer(poly_water_opp, width = straal_water_m)
  land_finaal_opp <- terra::mask(land_finaal_opp_ruw, poly_buffer_opp)
} else {
  land_finaal_opp <- land_finaal_opp_ruw
}

waterbiotoop_finaal_opp <- water_finaal_opp
landbiotoop_finaal_opp  <- land_finaal_opp

sierlijkewitsnuitlibel_finaal_opp <- terra::cover(waterbiotoop_finaal_opp, landbiotoop_finaal_opp)

# DEFINITIEVE MODELUITGANGEN
final_max <- sierlijkewitsnuitlibel_finaal_max %>% terra::crop(template_TV)
final_opp <- sierlijkewitsnuitlibel_finaal_opp %>% terra::crop(template_TV)

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

rm(witsnuit_water_relatie_max, witsnuit_foer_relatie_max, witsnuit_water_relatie_opp, witsnuit_foer_relatie_opp,
   witsnuit_leefgebied1_max, witsnuit_leefgebied1_opp, w_clean, f_clean, som_fracties, som_gemedieerd, witsnuit_leefgebied1_opp_raw,
   r_binair_leef1_opp)
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
