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
soort <- "grutto"

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

oppervlakte_ha <- resultaat$MinOpp_ha
afstand_m      <- resultaat$AfstandBiotopen_m[1]
buffer_m       <- resultaat$Dispersiecap_m[1]

rm(df, resultaat)

area_shape  <- vect(here("data/input/Turnhouts_Vennegebied.shp"))
master_grid <- rast(here("data/input/Raster_Vlaanderen/Vlaanderen_MasterGrid_10m.tif"))[[1]]

df_namen_sleutel <- read_csv(here("data/input/Excel_files/BWK_Laag_Namen_2025.csv"), show_col_types = FALSE)
gouden_namenlijst <- tolower(trimws(df_namen_sleutel$Laagnaam))

area_shape_proj <- project(area_shape, crs(master_grid))
area_buffer_fix <- buffer(area_shape_proj, width = buffer_m)

message("-> Vertaalraster voor globale/lokale cellen bliksemsnel opbouwen...")

id_raster_TV <- crop(master_grid, area_buffer_fix, snap = "near")
template_TV_mask <- terra::rasterize(area_buffer_fix, id_raster_TV, field = 1)

lokale_ids  <- terra::cells(template_TV_mask)
coords      <- terra::xyFromCell(template_TV_mask, lokale_ids)
globale_ids <- terra::cellFromXY(master_grid, coords)

vertaal_df <- data.table(
  lokale_id  = lokale_ids,
  globale_id = globale_ids
)[!is.na(globale_id)]

studiegebied_globale_ids <- vertaal_df$globale_id

template_TV <- terra::classify(template_TV_mask, cbind(NA, 0))
values(id_raster_TV) <- NA

rm(template_TV_mask, lokale_ids, coords, globale_ids)
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

cultuur_0a_match <- lijst_matches[["cultuur0a"]]
cultuur_0a_opp   <- lijst_oppervlaktes[["cultuur0a"]]
akker_match      <- lijst_matches[["akker"]]
hooiland_match   <- lijst_matches[["hooiland"]]
hooiland_opp     <- lijst_oppervlaktes[["hooiland"]]
water_bwk_match  <- lijst_matches[["kleinwater_bwk"]]
water_bwk_opp    <- lijst_oppervlaktes[["kleinwater_bwk"]]
kern1_match      <- lijst_matches[["kern1"]]
kern2_match      <- lijst_matches[["kern2"]]
kern1_opp        <- lijst_oppervlaktes[["kern1"]]
kern2_opp        <- lijst_oppervlaktes[["kern2"]]
c1b_match        <- lijst_matches[["cultuurhooi1b"]]
c1b_opp          <- lijst_oppervlaktes[["cultuurhooi1b"]]

raster_simpel_final <- id_raster_TV

rm(resultaten_gegroepeerd, df_nieuw, lijst_matches, lijst_oppervlaktes)
gc()

# --- CULTUURGRASLAND ---
grutto_cultuur_max <- (cultuur_0a_match == 1) & (akker_match == 0 | is.na(akker_match))
grutto_cultuur_max <- terra::ifel(grutto_cultuur_max == 1, 1, NA)
grutto_cultuur_opp <- terra::mask(cultuur_0a_opp, grutto_cultuur_max)

names(grutto_cultuur_max) <- "cultuur_max"
names(grutto_cultuur_opp) <- "cultuur_opp"

# --- HOOILAND ---
grutto_hooiland_max <- hooiland_match
grutto_hooiland_max[grutto_hooiland_max == 0] <- NA
grutto_hooiland_opp <- terra::mask(hooiland_opp, grutto_hooiland_max)

names(grutto_hooiland_max) <- "hooiland_max"
names(grutto_hooiland_opp) <- "hooiland_opp"

# --- WATER (GRB + BWK HYBRID MET RAM-STROOMBEVEILIGING) ---
r_grb_raw <- rast(here("data/input/Raster_Vlaanderen/vlaanderen_grb_water_10m.tif"))

r_grb_crop <- terra::crop(r_grb_raw, area_buffer_fix)
r_grb_sync <- terra::resample(r_grb_crop, raster_simpel_final, method = "near")
rm(r_grb_raw, r_grb_crop)
gc()

grb_water_bin <- terra::ifel(r_grb_sync > 0, 1, NA)
rm(r_grb_sync)
gc()

temp_patchfile <- tempfile(fileext = ".tif")

water_patches <- terra::patches(
  grb_water_bin, 
  directions = 8, 
  zeroAsNA = TRUE, 
  filename = temp_patchfile, 
  overwrite = TRUE
)

water_freq <- terra::freq(water_patches)
klein_water_ids <- water_freq$value[water_freq$count <= 50]

grutto_water_grb_max <- terra::ifel(water_patches %in% klein_water_ids, 1, NA)

terra::crs(grutto_water_grb_max) <- terra::crs(raster_simpel_final)
terra::crs(water_bwk_match)      <- terra::crs(raster_simpel_final)
terra::crs(water_bwk_opp)        <- terra::crs(raster_simpel_final)

grutto_water_max <- terra::cover(grutto_water_grb_max, water_bwk_match)
grutto_water_opp <- terra::ifel(!is.na(grutto_water_grb_max), 0.01, water_bwk_opp)
grutto_water_opp <- terra::mask(grutto_water_opp, grutto_water_max)

names(grutto_water_max) <- "water_max"
names(grutto_water_opp) <- "water_opp"

rm(water_patches, grutto_water_grb_max, grb_water_bin, water_freq)
if(file.exists(temp_patchfile)) unlink(temp_patchfile)
gc()

# --- GROENKAART STAP ---
hooggroen_raw         <- rast(here("data/input/ASCI Files/Groenkaart_2021.tif"))
area_buffer_groen_crs <- project(area_buffer_fix, crs(hooggroen_raw))
hooggroen_local_raw   <- crop(hooggroen_raw, area_buffer_groen_crs)

hooggroen_binair  <- terra::classify(hooggroen_local_raw, matrix(c(0.5, 1.5, 1), ncol = 3, byrow = TRUE), others = 0)
hooggroen_10m_raw <- terra::aggregate(hooggroen_binair, fact = 10, fun = "max")
hooggroen_sync    <- project(hooggroen_10m_raw, crs(raster_simpel_final)) %>% terra::resample(raster_simpel_final, method = "near")

grutto_cultuur_max_f  <- terra::mask(grutto_cultuur_max, hooggroen_sync, maskvalues = 1)
grutto_hooiland_max_f <- terra::mask(grutto_hooiland_max, hooggroen_sync, maskvalues = 1)
grutto_water_max_f    <- terra::mask(grutto_water_max, hooggroen_sync, maskvalues = 1)

grutto_cultuur_opp_f  <- terra::mask(grutto_cultuur_opp, hooggroen_sync, maskvalues = 1)
grutto_hooiland_opp_f <- terra::mask(grutto_hooiland_opp, hooggroen_sync, maskvalues = 1)
grutto_water_opp_f    <- terra::mask(grutto_water_opp, hooggroen_sync, maskvalues = 1)

rm(hooggroen_raw, hooggroen_local_raw, hooggroen_binair, hooggroen_10m_raw)
gc()

# --- KERNEN INRICHTEN ---
grutto_basis_kern_max <- max(kern1_match, kern2_match, na.rm = TRUE)
grutto_basis_kern_opp <- max(kern1_opp, kern2_opp, na.rm = TRUE)

r_drain_raw   <- rast(here("data/input/Raster_Vlaanderen/vlaanderen_drainage_n2khab_10m.tif"))
r_drain_local <- r_drain_raw %>% terra::crop(area_buffer_fix) %>% terra::resample(raster_simpel_final, method = "near")

drain_cats            <- terra::cats(r_drain_local)[[1]]
geselecteerde_letters <- c("d", "e", "f", "h", "i", "g", "e-f", "h-i", "e-i")
grutto_drain_ids      <- drain_cats$value[drain_cats$Label %in% geselecteerde_letters]

masker_drainage  <- r_drain_local %in% grutto_drain_ids
grutto_kern3_max <- terra::ifel(masker_drainage == 1, grutto_basis_kern_max, NA)
grutto_kern3_opp <- terra::mask(grutto_basis_kern_opp, grutto_kern3_max)

rm(r_drain_raw, r_drain_local, masker_drainage)

# --- BEBOUWINGSFILTER ---
if (!exists("tabel_vlaanderen")) {
  tabel_vlaanderen <- readRDS(scenario_path)
  setDT(tabel_vlaanderen)
  tabel_vlaanderen[, CODE := tolower(trimws(CODE))]
}

bebouwings_lagen <- gouden_namenlijst[grepl("^u[a-z0-9]", gouden_namenlijst) | gouden_namenlijst == "u"]
bebouwings_lagen <- bebouwings_lagen[bebouwings_lagen != "ulm"]

tabel_vlaanderen_bebouw <- tabel_vlaanderen
tabel_bebouw_TV <- tabel_vlaanderen_bebouw[CODE %in% bebouwings_lagen & cel_id %in% studiegebied_globale_ids]
tabel_bebouw_som <- tabel_bebouw_TV[, .(Present = ifelse(sum(BWK_FRAC, na.rm = TRUE) > 0, 1, 0)), by = .(cel_id)]

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

res_bebouw_groot <- cluster_filter_compleet(
  masker     = r_bebouwing_schoon,
  opp_laag   = r_bebouwing_schoon,
  drempel_m2 = 100000, 
  dist_m     = 100,    
  werkelijk  = FALSE
)

r_bebouw_groot <- res_bebouw_groot$raster

if (!is.null(r_bebouw_groot) && !terra::global(is.na(r_bebouw_groot), "sum")[[1]] == terra::ncell(r_bebouw_groot)) {
  buf_bebouw <- terra::buffer(r_bebouw_groot, width = 100)
  masker_bebouw_100m_groot <- terra::ifel(!is.na(buf_bebouw) & buf_bebouw > 0, 1, NA)
  
  grutto_kern_bebouw_max <- terra::mask(grutto_kern3_max, masker_bebouw_100m_groot, inverse = TRUE)
  grutto_kern_bebouw_opp <- terra::mask(grutto_kern3_opp, masker_bebouw_100m_groot, inverse = TRUE)
} else {
  grutto_kern_bebouw_max <- grutto_kern3_max
  grutto_kern_bebouw_opp <- grutto_kern3_opp
  masker_bebouw_100m_groot <- id_raster_TV * NA
}

rm(tabel_vlaanderen_bebouw, tabel_bebouw_TV, tabel_bebouw_som, tabel_bebouw_mapping, 
   r_bebouwing_raw, r_bebouwing_schoon, res_bebouw_groot, r_bebouw_groot)
gc()

# --- HOOGGROEN FILTER ---
r_groen_bin <- terra::ifel(!is.na(hooggroen_sync) & hooggroen_sync > 0, 1, NA)

w_matrix <- matrix(1, nrow = 5, ncol = 5)
temp_groen_file <- tempfile(fileext = ".tif")

groen_som <- terra::focal(
  r_groen_bin, 
  w = w_matrix, 
  fun = "sum", 
  na.rm = TRUE,
  filename = temp_groen_file,
  overwrite = TRUE
)

r_groen_20are <- terra::ifel(groen_som >= 12 & r_groen_bin == 1, 1, NA)
buf_groen <- terra::buffer(r_groen_20are, width = 50)
masker_groen_50m_cluster <- terra::ifel(!is.na(buf_groen) & buf_groen > 0, 1, NA)

grutto_kern_finaal_max <- terra::mask(grutto_kern_bebouw_max, masker_groen_50m_cluster, inverse = TRUE)
grutto_kern_finaal_opp <- terra::mask(grutto_kern_bebouw_opp, masker_groen_50m_cluster, inverse = TRUE)

rm(r_groen_bin, groen_som, r_groen_20are, hooggroen_sync, buf_groen)
gc()

# --- CLUSTER-TABEL MET MAXIMALE POTENTIE ---
drempel_kern_m2 <- 10 * 10000 

res_kern_max <- cluster_filter_compleet(
  masker     = grutto_kern_finaal_max,
  opp_laag   = grutto_kern_finaal_max,
  drempel_m2 = drempel_kern_m2,
  dist_m     = 50,
  werkelijk  = FALSE
)
grutto_kern_10ha_max <- res_kern_max$raster

r_binair_finaal_opp <- terra::ifel(!is.na(grutto_kern_finaal_opp) & grutto_kern_finaal_opp > 0, 1, NA)
res_kern_opp <- cluster_filter_compleet(
  masker     = r_binair_finaal_opp,
  opp_laag   = grutto_kern_finaal_opp,
  drempel_m2 = drempel_kern_m2,
  dist_m     = 50,
  werkelijk  = TRUE
)
grutto_kern_10ha_opp <- res_kern_opp$raster

gc()

# --- SAMENVOEGEN BOUWSTENEN ---
grutto_bouwsteen_totaal_max <- terra::cover(grutto_cultuur_max, grutto_hooiland_max)
grutto_bouwsteen_totaal_max <- terra::cover(grutto_bouwsteen_totaal_max, grutto_water_max)
names(grutto_bouwsteen_totaal_max) <- "bouwsteen_totaal_max"

grutto_bouwsteen_totaal_opp <- terra::cover(grutto_cultuur_opp, grutto_hooiland_opp)
grutto_bouwsteen_totaal_opp <- terra::cover(grutto_bouwsteen_totaal_opp, grutto_water_opp)
names(grutto_bouwsteen_totaal_opp) <- "bouwsteen_totaal_opp"

grutto_bouwsteen_totaal_max <- terra::ifel(grutto_bouwsteen_totaal_max == 0, NA, grutto_bouwsteen_totaal_max)
grutto_bouwsteen_totaal_opp <- terra::ifel(grutto_bouwsteen_totaal_opp == 0, NA, grutto_bouwsteen_totaal_opp)

# --- SAMENVOEGEN BIJKOMEND HABITAT ---
grutto_1b_match_sync <- c1b_match
grutto_1b_opp_sync   <- c1b_opp

grutto_1b_max <- grutto_1b_match_sync
grutto_1b_max[grutto_1b_max == 0] <- NA
grutto_1b_opp <- terra::mask(grutto_1b_opp_sync, grutto_1b_max)

names(grutto_1b_max) <- "1b_max"
names(grutto_1b_opp) <- "1b_opp"

grutto_basis_bouw_1b_max <- terra::cover(grutto_bouwsteen_totaal_max, grutto_1b_max)
grutto_cultuurhooi1_opp  <- terra::cover(grutto_bouwsteen_totaal_opp, grutto_1b_opp)

filter_akker <- terra::ifel(!is.na(akker_match) & akker_match == 1, 1, NA)

if(exists("masker_bebouw_100m_groot") && !all(is.na(suppressWarnings(terra::minmax(masker_bebouw_100m_groot))))) {
  filter_bebouw <- terra::ifel(!is.na(masker_bebouw_100m_groot) & masker_bebouw_100m_groot > 0, 1, NA)
} else {
  filter_bebouw <- id_raster_TV * NA
}

if(exists("masker_groen_50m_cluster") && !all(is.na(suppressWarnings(terra::minmax(masker_groen_50m_cluster))))) {
  filter_groen <- terra::ifel(!is.na(masker_groen_50m_cluster) & masker_groen_50m_cluster > 0, 1, NA)
} else {
  filter_groen <- id_raster_TV * NA
}

grutto_c1_m1_max     <- terra::mask(grutto_basis_bouw_1b_max, filter_akker, inverse = TRUE)
grutto_c1_m2_max     <- terra::mask(grutto_c1_m1_max, filter_bebouw, inverse = TRUE)
grutto_c1_finaal_max <- terra::mask(grutto_c1_m2_max, filter_groen, inverse = TRUE)

grutto_c1_m1_opp     <- terra::mask(grutto_cultuurhooi1_opp, filter_akker, inverse = TRUE)
grutto_c1_m2_opp     <- terra::mask(grutto_c1_m1_opp, filter_bebouw, inverse = TRUE)
grutto_c1_finaal_opp <- terra::mask(grutto_c1_m2_opp, filter_groen, inverse = TRUE)

grutto_scenA_finaal_max <- grutto_c1_finaal_max
grutto_scenB_finaal_opp <- grutto_c1_finaal_opp

rm(grutto_basis_bouw_1b_max, grutto_c1_m1_max, grutto_c1_m2_max, grutto_c1_m1_opp, grutto_c1_m2_opp, 
   filter_akker, filter_bebouw, filter_groen, grutto_cultuurhooi1_opp, grutto_bouwsteen_totaal_max, grutto_bouwsteen_totaal_opp)
gc()

# --- WEDERZIJDSE KOPPELING (KERN & BOUWSTEEN) ---
voer_koppeling_en_drempel_uit_strikt <- function(kern_rast, bouw_rast, afstand = 200, min_kern_ha = 10, min_totaal_ha = 25) {
  if (is.null(kern_rast) || terra::global(is.na(kern_rast), "sum")[[1]] == terra::ncell(kern_rast)) return(NULL)
  
  r_totaal_bin <- terra::cover(
    terra::ifel(!is.na(kern_rast) & kern_rast > 0, 1, NA),
    terra::ifel(!is.na(bouw_rast) & bouw_rast > 0, 1, NA)
  )
  
  if (terra::global(is.na(r_totaal_bin), "sum")[[1]] == terra::ncell(r_totaal_bin)) return(NULL)
  
  r_buf <- terra::buffer(r_totaal_bin, width = afstand / 2)
  cl_netwerk <- terra::patches(r_buf, directions = 8, zeroAsNA = TRUE)
  cl_netwerk_biotoop <- terra::mask(cl_netwerk, r_totaal_bin)
  
  r_comb_val <- terra::cover(kern_rast, bouw_rast)
  stats_netwerk <- terra::zonal(r_comb_val * 0.01, cl_netwerk_biotoop, fun = "sum", na.rm = TRUE)
  colnames(stats_netwerk) <- c("netwerk_id", "totaal_ha")
  
  cl_kernen_only <- terra::mask(cl_netwerk_biotoop, kern_rast)
  stats_kernen <- terra::zonal(kern_rast * 0.01, cl_kernen_only, fun = "sum", na.rm = TRUE)
  colnames(stats_kernen) <- c("netwerk_id", "kern_ha")
  
  valide_netwerken <- stats_netwerk %>%
    inner_join(stats_kernen, by = "netwerk_id") %>%
    filter(totaal_ha >= min_totaal_ha & kern_ha >= min_kern_ha) %>%
    pull(netwerk_id)
  
  if (length(valide_netwerken) == 0) return(NULL)
  
  m_valide <- cl_netwerk_biotoop %in% valide_netwerken
  m_valide <- terra::ifel(m_valide == 1, 1, NA)
  
  return(list(
    kern = terra::mask(kern_rast, m_valide),
    bouw = terra::mask(bouw_rast, m_valide)
  ))
}

resA_strikt <- voer_koppeling_en_drempel_uit_strikt(grutto_kern_10ha_max, grutto_scenA_finaal_max, 200, 10, 25)
resB_strikt <- voer_koppeling_en_drempel_uit_strikt(grutto_kern_10ha_opp, grutto_scenB_finaal_opp, 200, 10, 25)

grutto_kaart_A <- maak_kaart_laag(resA_strikt)

# ==============================================================================
# DEFINITIEVE TOEWIJSING AAN EXPORT VARIABELEN (CRUCIALE FIX)
# ==============================================================================
if (!is.null(resA_strikt)) {
  final_max <- terra::cover(resA_strikt$kern, resA_strikt$bouw)
} else {
  final_max <- template_TV * NA
}

if (!is.null(resB_strikt)) {
  final_opp <- terra::cover(resB_strikt$kern, resB_strikt$bouw)
} else {
  final_opp <- template_TV * NA
}

if (!all(is.na(suppressWarnings(terra::minmax(final_max))))) {
  cl_max <- terra::patches(final_max, directions = 8, zeroAsNA = TRUE)
} else {
  cl_max <- template_TV * NA
}

if (!all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
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

# 3. Bepaal Analytisch Metacluster ID-raster
if (exists("cl_opp") && !is.null(cl_opp) && !all(is.na(suppressWarnings(terra::minmax(cl_opp))))) {
  id_export_rast <- cl_opp
} else if (exists("final_opp") && !is.null(final_opp) && !all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
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
