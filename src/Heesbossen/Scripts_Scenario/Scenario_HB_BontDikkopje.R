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
soort <- "bontdikkopje"

# --- DYNAMISCHE SCENARIO PARAMETER CHECK ---
if (exists("SCENARIO_RDS_PAD") && !is.null(SCENARIO_RDS_PAD)) {
  scenario_rds_path <- SCENARIO_RDS_PAD
} else if (exists("params") && !is.null(params$scenario_rds_path)) {
  scenario_rds_path <- params$scenario_rds_path
} else {
  scenario_rds_path <- "data/input/Scenario_rds/HB_Scenario_BWK_2025.rds"
}

p_raw <- gsub("^([.][.]/)+", "", scenario_rds_path)
scenario_path <- here::here(p_raw)

if (!file.exists(scenario_path)) {
  stop(paste("❌ FOUT: Scenario RDS bestand NIET gevonden op:", scenario_path))
}

scen_volledig <- basename(scenario_path)
scenario_naam <- gsub("^HB_Scenario_|^Scenario_|.rds$", "", scen_volledig)

message(paste("Verwerken van soort:", soort, "binnen scenario:", scenario_naam))

resultaat <- df %>% filter(tolower(trimws(Soort)) == soort) %>% select(Type, MinOpp_ha, AfstandBiotopen_m, Dispersiecap_m)

oppervlakte_ha <- resultaat$MinOpp_ha[1]
afstand_m      <- resultaat$AfstandBiotopen_m[1]
buffer_m       <- resultaat$Dispersiecap_m[1]

rm(df, resultaat)

area_shape  <- vect(here("data/input/Heesbossen.shp"))
master_grid <- rast(here("data/input/Raster_Vlaanderen/Vlaanderen_MasterGrid_10m.tif"))[[1]]

df_namen_sleutel <- read_csv(here("data/input/Excel_files/BWK_Laag_Namen_2025.csv"), show_col_types = FALSE)
gouden_namenlijst <- tolower(trimws(df_namen_sleutel$Laagnaam))

area_shape_proj <- project(area_shape, crs(master_grid))
area_buffer_fix <- buffer(area_shape_proj, width = buffer_m)

message("-> Vertaalraster voor globale/lokale cellen opbouwen...")
id_raster_HB <- crop(master_grid, area_buffer_fix, snap = "near")

globale_id_raster <- master_grid
globale_id_raster <- terra::init(globale_id_raster, fun = "cell")

id_raster_HB_globale_values <- crop(globale_id_raster, area_buffer_fix, snap = "near")
id_raster_HB_masked <- mask(id_raster_HB_globale_values, area_buffer_fix)

message("-> Vertaaltabel bliksemsnel opbouwen via C++ dataframe extractie...")
df_extractie <- as.data.frame(id_raster_HB_masked, cells = TRUE)
vertaal_df <- as.data.table(df_extractie)
setnames(vertaal_df, c(1, 2), c("lokale_id", "globale_id"))

vertaal_df <- vertaal_df[!is.na(globale_id)]
studiegebied_globale_ids <- unique(vertaal_df$globale_id)

values(id_raster_HB) <- NA
template_HB <- terra::rasterize(area_buffer_fix, id_raster_HB, field = 1, background = 0)

rm(globale_id_raster, id_raster_HB_globale_values, id_raster_HB_masked, df_extractie)
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
  
  tabel_HB <- tabel_gefilterd[cel_id %in% studiegebied_globale_ids]
  tabel_HB_unique <- unique(tabel_HB, by = c("cel_id", "CODE"))
  tabel_cel_som <- tabel_HB_unique[, .(Oppervlakte = pmin(sum(BWK_FRAC, na.rm = TRUE), 1.0)), by = .(cel_id)]
  
  r_match_type <- id_raster_HB * NA
  r_opp_type   <- id_raster_HB * NA
  
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
bosranden_max         <- lijst_matches[["bosranden"]]
bosranden_opp         <- lijst_oppervlaktes[["bosranden"]]
water_bwk_max         <- lijst_matches[["water_bwk"]]
water_bwk_opp         <- lijst_oppervlaktes[["water_bwk"]]
foerageer_bwk_max     <- lijst_matches[["foerageer_bwk"]]
foerageer_bwk_opp     <- lijst_oppervlaktes[["foerageer_bwk"]]

rm(tabel_vlaanderen, vertaal_df, lijst_matches, lijst_oppervlaktes)
gc()

# ==============================================================================
# DRAINAGEFILTER
# ==============================================================================
r_drain_raw   <- rast(here("data/input/Raster_Vlaanderen/vlaanderen_drainage_n2khab_10m.tif"))
r_drain_local <- r_drain_raw %>% terra::crop(area_buffer_fix) %>% terra::resample(template_HB, method = "near")

drain_cats            <- terra::cats(r_drain_local)[[1]]
geselecteerde_letters <- c("d", "e", "f", "h", "i", "g", "e-f", "h-i", "e-i")

# Strikt veilige conversie naar character om factor-join crashes te vermijden
drain_labels_char     <- as.character(drain_cats$Label)
bontdikkopje_drain_ids <- drain_cats$value[drain_labels_char %in% geselecteerde_letters]

masker_drainage         <- r_drain_local %in% bontdikkopje_drain_ids
voortplanting_vocht_max <- terra::ifel(masker_drainage == 1, voortplanting_bwk_max, NA)
voortplanting_vocht_opp <- terra::mask(voortplanting_bwk_opp, voortplanting_vocht_max)

rm(r_drain_raw, r_drain_local, masker_drainage)
gc()

# ==============================================================================
# AFSTANDSFILTER TOT OPEN WATER
# ==============================================================================
r_grb_water_raw   <- rast(here("data/input/Raster_Vlaanderen/vlaanderen_grb_water_10m.tif"))
r_grb_water_local <- r_grb_water_raw %>% 
  terra::crop(area_buffer_fix, snap = "near") %>% 
  terra::resample(template_HB, method = "near")

crs(r_grb_water_local) <- crs(template_HB)

masker_grb_water <- r_grb_water_local == 1
masker_grb_water[masker_grb_water == 0] <- NA

if(!is.null(water_bwk_max)) {
  bontdikkopjeWater <- terra::cover(masker_grb_water, water_bwk_max)
} else {
  bontdikkopjeWater <- masker_grb_water
}

message("-> Afstandsmatrix tot water berekenen (max 200m)...")
afstand_tot_water <- terra::distance(bontdikkopjeWater)
zone_binnen_200m  <- terra::ifel(afstand_tot_water <= 200, 1, NA)

if(!is.null(voortplanting_vocht_max)) {
  voortplanting_vocht_max <- terra::mask(voortplanting_vocht_max, zone_binnen_200m)
  voortplanting_vocht_opp <- terra::mask(voortplanting_vocht_opp, zone_binnen_200m)
}

rm(r_grb_water_raw, r_grb_water_local, masker_grb_water, bontdikkopjeWater, afstand_tot_water, zone_binnen_200m)
gc()

# ==============================================================================
# CLUSTER-ANALYSE (VOORTPLANTING & BOSRANDEN)
# ==============================================================================
message("-> VOORTPLANTING: Bosrand-voorwaarde & clusteranalyse opbouwen...")

verwerk_habitat_met_bosrand <- function(r_habitat_max, r_habitat_opp, r_bos_max, r_bos_opp, dist_m = 20, drempel_m2, afstand_m = 0) {
  
  if (!terra::hasValues(r_habitat_max) || terra::global(r_habitat_max, "notNA")$notNA == 0) {
    return(list(max = template_HB * NA, opp = template_HB * NA, clusters_max = template_HB * NA, clusters_opp = template_HB * NA))
  }

  # --- SPOOR A: MAXIMALE POTENTIE ---
  temp_mask_max <- r_habitat_max
  if (afstand_m > 0) temp_mask_max <- terra::buffer(r_habitat_max, width = afstand_m / 2)
  
  cl_habitat_max      <- terra::patches(temp_mask_max, directions = 8, zeroAsNA = TRUE)
  cl_habitat_only_max <- terra::mask(cl_habitat_max, r_habitat_max)
  
  f_hab_max     <- terra::freq(cl_habitat_only_max)
  stats_hab_max <- data.frame(ID = f_hab_max$value, Area_m2 = f_hab_max$count * 100)
  
  voldoet_opp_ids_max <- stats_hab_max$ID[!is.na(stats_hab_max$ID) & stats_hab_max$Area_m2 >= drempel_m2]
  
  voldoet_bos_ids_max <- c()
  if (length(voldoet_opp_ids_max) > 0) {
    cl_habitat_groot_max    <- cl_habitat_only_max %in% voldoet_opp_ids_max
    cl_habitat_groot_id_max <- terra::ifel(cl_habitat_groot_max == 1, cl_habitat_only_max, NA)
    
    if (terra::hasValues(r_bos_max) && terra::global(r_bos_max, "notNA")$notNA > 0) {
      dist_naar_hab_max       <- terra::distance(cl_habitat_groot_id_max)
      bos_contact_pixels_max  <- terra::ifel(dist_naar_hab_max <= 10 & !is.na(r_bos_max), 1, NA)
      
      if (terra::hasValues(bos_contact_pixels_max) && terra::global(bos_contact_pixels_max, "notNA")$notNA > 0) {
        cl_contact_lijnen_max <- terra::patches(bos_contact_pixels_max, directions = 8, zeroAsNA = TRUE)
        f_contact_max         <- terra::freq(cl_contact_lijnen_max)
        stats_contact_max     <- data.frame(ID = f_contact_max$value, Lengte_m = f_contact_max$count * 10)
        raaklijn_50m_ids_max  <- stats_contact_max$ID[!is.na(stats_contact_max$ID) & stats_contact_max$Lengte_m >= 50]
        
        if (length(raaklijn_50m_ids_max) > 0) {
          r_raaklijn_50m_max     <- cl_contact_lijnen_max %in% raaklijn_50m_ids_max
          r_raaklijn_50m_max     <- terra::ifel(r_raaklijn_50m_max == 1, 1, NA)
          hab_raakt_50m_lijn_max <- terra::mask(cl_habitat_groot_id_max, terra::buffer(r_raaklijn_50m_max, width = 5))
          voldoet_bos_ids_max    <- unique(na.omit(terra::values(hab_raakt_50m_lijn_max, mat = FALSE)))
        }
      }
    }
  }
  
  if (length(voldoet_bos_ids_max) > 0) {
    hab_goedgekeurd_mask_max  <- cl_habitat_groot_id_max %in% voldoet_bos_ids_max
    r_hab_goedgekeurd_max     <- terra::ifel(hab_goedgekeurd_mask_max == 1, r_habitat_max, NA)
    
    dist_naar_hab_finaal_max  <- terra::distance(r_hab_goedgekeurd_max)
    bos_binnen_20m_max        <- terra::ifel(dist_naar_hab_finaal_max <= dist_m & !is.na(r_bos_max), 1, NA)
    
    combo_potentieel_max      <- terra::cover(r_hab_goedgekeurd_max, bos_binnen_20m_max)
    cl_combo_max              <- terra::patches(combo_potentieel_max, directions = 8, zeroAsNA = TRUE)
    cl_die_habitat_raken_max  <- terra::mask(cl_combo_max, r_hab_goedgekeurd_max)
    aangesloten_ids_max       <- unique(na.omit(terra::values(cl_die_habitat_raken_max, mat = FALSE)))
    
    combo_finaal_mask_max     <- cl_combo_max %in% aangesloten_ids_max
    complex_raw_max           <- terra::ifel(combo_finaal_mask_max == 1, 1, NA)
    cl_finaal_max             <- terra::mask(cl_combo_max, complex_raw_max)
  } else {
    complex_raw_max <- template_HB * NA
    cl_finaal_max   <- template_HB * NA
  }

  # --- SPOOR B: WERKELIJKE OPPERVLAKTE ---
  r_habitat_bin_opp <- terra::ifel(!is.na(r_habitat_opp) & r_habitat_opp > 0, 1, NA)
  
  if (terra::hasValues(r_habitat_bin_opp) && terra::global(r_habitat_bin_opp, "notNA")$notNA > 0) {
    temp_mask_opp <- r_habitat_bin_opp
    if (afstand_m > 0) temp_mask_opp <- terra::buffer(r_habitat_bin_opp, width = afstand_m / 2)
    
    cl_habitat_opp      <- terra::patches(temp_mask_opp, directions = 8, zeroAsNA = TRUE)
    cl_habitat_only_opp <- terra::mask(cl_habitat_opp, r_habitat_bin_opp)
    
    stats_hab_opp       <- terra::zonal(r_habitat_opp, cl_habitat_only_opp, fun = "sum", na.rm = TRUE)
    colnames(stats_hab_opp) <- c("ID", "Waarde")
    stats_hab_opp$Area_m2   <- stats_hab_opp$Waarde * 100
    
    voldoet_opp_ids_opp <- stats_hab_opp$ID[!is.na(stats_hab_opp$ID) & stats_hab_opp$Area_m2 >= drempel_m2]
    
    voldoet_bos_ids_opp <- c()
    if (length(voldoet_opp_ids_opp) > 0) {
      cl_habitat_groot_opp    <- cl_habitat_only_opp %in% voldoet_opp_ids_opp
      cl_habitat_groot_id_opp <- terra::ifel(cl_habitat_groot_opp == 1, cl_habitat_only_opp, NA)
      
      r_bos_bin_opp <- terra::ifel(!is.na(r_bos_opp) & r_bos_opp > 0, 1, NA)
      if (terra::hasValues(r_bos_bin_opp) && terra::global(r_bos_bin_opp, "notNA")$notNA > 0) {
        dist_naar_hab_opp       <- terra::distance(cl_habitat_groot_id_opp)
        bos_contact_pixels_opp  <- terra::ifel(dist_naar_hab_opp <= 10 & !is.na(r_bos_bin_opp), 1, NA)
        
        if (terra::hasValues(bos_contact_pixels_opp) && terra::global(bos_contact_pixels_opp, "notNA")$notNA > 0) {
          cl_contact_lijnen_opp <- terra::patches(bos_contact_pixels_opp, directions = 8, zeroAsNA = TRUE)
          f_contact_opp         <- terra::freq(cl_contact_lijnen_opp)
          stats_contact_opp     <- data.frame(ID = f_contact_opp$value, Lengte_m = f_contact_opp$count * 10)   
          raaklijn_50m_ids_opp  <- stats_contact_opp$ID[!is.na(stats_contact_opp$ID) & stats_contact_opp$Lengte_m >= 50]
          
          if (length(raaklijn_50m_ids_opp) > 0) {
            r_raaklijn_50m_opp     <- cl_contact_lijnen_opp %in% raaklijn_50m_ids_opp
            r_raaklijn_50m_opp     <- terra::ifel(r_raaklijn_50m_opp == 1, 1, NA)
            hab_raakt_50m_lijn_opp <- terra::mask(cl_habitat_groot_id_opp, terra::buffer(r_raaklijn_50m_opp, width = 5))
            voldoet_bos_ids_opp    <- unique(na.omit(terra::values(hab_raakt_50m_lijn_opp, mat = FALSE)))
          }
        }
      }
    }
    
    if (length(voldoet_bos_ids_opp) > 0) {
      hab_goedgekeurd_mask_opp  <- cl_habitat_groot_id_opp %in% voldoet_bos_ids_opp
      r_hab_goedgekeurd_opp     <- terra::mask(r_habitat_opp, hab_goedgekeurd_mask_opp)
      
      dist_naar_hab_finaal_opp  <- terra::distance(r_hab_goedgekeurd_opp)
      r_bos_bin_opp             <- terra::ifel(!is.na(r_bos_opp) & r_bos_opp > 0, 1, NA)
      bos_binnen_20m_opp        <- terra::ifel(dist_naar_hab_finaal_opp <= dist_m & !is.na(r_bos_bin_opp), 1, NA)
      
      combo_potentieel_opp      <- terra::cover(r_hab_goedgekeurd_opp, terra::mask(r_bos_opp, bos_binnen_20m_opp))
      cl_combo_opp              <- terra::patches(terra::ifel(!is.na(combo_potentieel_opp) & combo_potentieel_opp > 0, 1, NA), directions = 8, zeroAsNA = TRUE)
      cl_die_habitat_raken_opp  <- terra::mask(cl_combo_opp, r_hab_goedgekeurd_opp)
      aangesloten_ids_opp       <- unique(na.omit(terra::values(cl_die_habitat_raken_opp, mat = FALSE)))
      
      combo_finaal_mask_opp     <- cl_combo_opp %in% aangesloten_ids_opp
      complex_raw_opp           <- terra::mask(combo_potentieel_opp, combo_finaal_mask_opp)
      cl_finaal_opp             <- terra::mask(cl_combo_opp, complex_raw_opp)
    } else {
      complex_raw_opp <- template_HB * NA
      cl_finaal_opp   <- template_HB * NA
    }
  } else {
    complex_raw_opp <- template_HB * NA
    cl_finaal_opp   <- template_HB * NA
  }

  return(list(max = complex_raw_max, opp = complex_raw_opp, clusters_max = cl_finaal_max, clusters_opp = cl_finaal_opp))
}

drempel_m2 <- oppervlakte_ha * 10000

vpt_complex <- verwerk_habitat_met_bosrand(
  r_habitat_max = voortplanting_vocht_max,
  r_habitat_opp = voortplanting_vocht_opp,
  r_bos_max     = bosranden_max,
  r_bos_opp     = bosranden_opp,
  dist_m        = 20,
  drempel_m2    = drempel_m2,
  afstand_m     = afstand_m
)

voortplanting_finaal_max <- vpt_complex$max
voortplanting_finaal_opp <- vpt_complex$opp
cl_voortplanting_max     <- vpt_complex$clusters_max
cl_voortplanting_opp     <- vpt_complex$clusters_opp

rm(vpt_complex)
gc()

# ==============================================================================
# CLUSTER-ANALYSE (FOERAGEREN)
# ==============================================================================
message("-> FOERAGEREN: Bosrand-voorwaarde & clusteranalyse opbouwen...")

foe_complex <- verwerk_habitat_met_bosrand(
  r_habitat_max = foerageer_bwk_max,
  r_habitat_opp = foerageer_bwk_opp,
  r_bos_max     = bosranden_max,
  r_bos_opp     = bosranden_opp,
  dist_m        = 20,
  drempel_m2    = drempel_m2,
  afstand_m     = afstand_m
)

res_max_foerageer <- cluster_filter_compleet(foe_complex$max, foe_complex$opp, 
                                              drempel_m2, afstand_m, werkelijk = FALSE)

r_binair_foe_opp  <- terra::ifel(!is.na(foe_complex$opp) & foe_complex$opp > 0, 1, NA)
res_opp_foerageer <- cluster_filter_compleet(r_binair_foe_opp, foe_complex$opp, 
                                              drempel_m2, afstand_m, werkelijk = TRUE)

foerageer_finaal_max <- res_max_foerageer$raster
foerageer_finaal_opp <- res_opp_foerageer$raster
cl_foerageer_max     <- res_max_foerageer$clusters
cl_foerageer_opp     <- res_opp_foerageer$clusters

rm(res_max_foerageer, res_opp_foerageer, foe_complex, r_binair_foe_opp)
gc()

# ==============================================================================
# TOTAAL LEEFGEBIED SAMENSTELLEN
# ==============================================================================
message("-> Totale leefgebiedslagen opbouwen via terra::cover...")

final_max <- terra::cover(voortplanting_finaal_max, foerageer_finaal_max)
final_opp <- terra::cover(voortplanting_finaal_opp, foerageer_finaal_opp)

# Veilige verwerking cl_max
if (!is.null(cl_foerageer_max) && !all(is.na(suppressWarnings(terra::minmax(cl_foerageer_max))))) {
  max_id_voortplanting <- 0
  if (!is.null(cl_voortplanting_max) && !all(is.na(suppressWarnings(terra::minmax(cl_voortplanting_max))))) {
    val_max <- terra::global(cl_voortplanting_max, "max", na.rm = TRUE)[[1]]
    if (!is.na(val_max)) max_id_voortplanting <- val_max
  }
  
  cl_foerageer_uniek_max <- cl_foerageer_max + (max_id_voortplanting + 1000)
  cl_max <- if (!is.null(cl_voortplanting_max)) terra::cover(cl_voortplanting_max, cl_foerageer_uniek_max) else cl_foerageer_uniek_max
} else {
  cl_max <- cl_voortplanting_max
}

# Veilige verwerking cl_opp
if (!is.null(cl_foerageer_opp) && !all(is.na(suppressWarnings(terra::minmax(cl_foerageer_opp))))) {
  max_id_voortplanting_opp <- 0
  if (!is.null(cl_voortplanting_opp) && !all(is.na(suppressWarnings(terra::minmax(cl_voortplanting_opp))))) {
    val_opp <- terra::global(cl_voortplanting_opp, "max", na.rm = TRUE)[[1]]
    if (!is.na(val_opp)) max_id_voortplanting_opp <- val_opp
  }
  
  cl_foerageer_uniek_opp <- cl_foerageer_opp + (max_id_voortplanting_opp + 1000)
  cl_opp <- if (!is.null(cl_voortplanting_opp)) terra::cover(cl_voortplanting_opp, cl_foerageer_uniek_opp) else cl_foerageer_uniek_opp
} else {
  cl_opp <- cl_voortplanting_opp
}

# Failsafe vulling voor ID rasters
if (is.null(cl_max) || all(is.na(terra::values(cl_max, mat=FALSE)))) {
  if (!is.null(final_max) && !all(is.na(suppressWarnings(terra::minmax(final_max))))) {
    cl_max <- terra::patches(final_max, directions = 8, zeroAsNA = TRUE)
  } else {
    cl_max <- template_HB * NA
  }
}

if (is.null(cl_opp) || all(is.na(terra::values(cl_opp, mat=FALSE)))) {
  if (!is.null(final_opp) && !all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
    cl_opp <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE)
  } else {
    cl_opp <- template_HB * NA
  }
}

# ==============================================================================
# SCHONE EXPORT BIOTOOP EN ANALYTISCH ID-RASTER (VOOR SCRIPT 2 / ARPL)
# ==============================================================================
base_dir <- here::here("data/output/Heesbossen/Rasters_Soorten", scenario_naam)

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
  terra::rast(template_HB, vals = NA)
}

# 2. Bepaal Werkelijke Oppervlakte Raster
werkelijk_export_rast <- if (exists("final_opp") && !is.null(final_opp) && !all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  terra::ifel(!is.na(final_opp) & final_opp > 0, 1, NA)
} else {
  terra::rast(template_HB, vals = NA)
}

# 3. Bepaal Analytisch Metacluster ID-raster
if (exists("cl_opp") && !is.null(cl_opp) && !all(is.na(suppressWarnings(terra::minmax(cl_opp))))) {
  id_export_rast <- cl_opp
} else if (exists("final_opp") && !is.null(final_opp) && !all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  id_export_rast <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE)
} else {
  id_export_rast <- terra::rast(template_HB, vals = NA)
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
