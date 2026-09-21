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

soort <- "grotemodderkruiper"

# --- DYNAMISCHE SCENARIO PARAMETER CHECK ---
if (exists("SCENARIO_RDS_PAD") && !is.null(SCENARIO_RDS_PAD)) {
  scenario_rds_path <- SCENARIO_RDS_PAD
} else if (exists("params") && !is.null(params$scenario_rds_path)) {
  scenario_rds_path <- params$scenario_rds_path
} else {
  scenario_rds_path <- "data/input/Scenario_rds/DM_Scenario_BWK_2025.rds"
}

p_raw <- gsub("^([.][.]/)+", "", scenario_rds_path)
scenario_path <- here::here(p_raw)

if (!file.exists(scenario_path)) {
  stop(paste("❌ FOUT: Scenario RDS bestand NIET gevonden op:", scenario_path))
}

scen_volledig <- basename(scenario_path)
scenario_naam <- gsub("^DM_Scenario_|^Scenario_|.rds$", "", scen_volledig)

message(paste("Verwerken van soort:", soort, "binnen scenario:", scenario_naam))

df <- read_excel(here::here("data/input/Excel_files/Soorten_bwk_afstanden.xlsx"))
resultaat <- df %>%
  filter(tolower(trimws(Soort)) == soort) %>%
  select(Type, MinOpp_ha, AfstandBiotopen_m, Dispersiecap_m)

oppervlakte_ha <- resultaat$MinOpp_ha[1]
afstand_m      <- resultaat$AfstandBiotopen_m[1]
buffer_m       <- resultaat$Dispersiecap_m[1]

rm(df, resultaat)

area_shape  <- vect(here("data/input/De_Maten.shp"))
master_grid <- rast(here("data/input/Raster_Vlaanderen/Vlaanderen_MasterGrid_10m.tif"))[[1]]

df_namen_sleutel <- read_csv(here("data/input/Excel_files/BWK_Laag_Namen_2025.csv"), show_col_types = FALSE)
gouden_namenlijst <- tolower(trimws(df_namen_sleutel$Laagnaam))

area_shape_proj <- project(area_shape, crs(master_grid))
area_buffer_fix <- buffer(area_shape_proj, width = buffer_m)

message("-> Vertaalraster voor globale/lokale cellen opbouwen via snelle MASK methode...")
id_raster_DM <- crop(master_grid, area_buffer_fix, snap = "near")

globale_id_raster <- master_grid
globale_id_raster <- terra::init(globale_id_raster, fun = "cell")

id_raster_DM_globale_values <- crop(globale_id_raster, area_buffer_fix, snap = "near")
id_raster_DM_masked <- mask(id_raster_DM_globale_values, area_buffer_fix)

message("-> Vertaaltabel bliksemsnel opbouwen via C++ dataframe extractie...")

df_extractie <- as.data.frame(id_raster_DM_masked, cells = TRUE)
vertaal_df <- as.data.table(df_extractie)
setnames(vertaal_df, c(1, 2), c("lokale_id", "globale_id"))

vertaal_df <- vertaal_df[!is.na(globale_id)]
studiegebied_globale_ids <- unique(vertaal_df$globale_id)

values(id_raster_DM) <- NA
template_DM <- terra::rasterize(area_buffer_fix, id_raster_DM, field = 1, background = 0)

rm(globale_id_raster, id_raster_DM_globale_values, id_raster_DM_masked, df_extractie)
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
  tabel_cel_som <- tabel_DM_unique[, .(Oppervlakte = pmin(sum(BWK_FRAC, na.rm = TRUE), 1.0)), by = .(cel_id)]
  
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

bwk_max     <- lijst_matches[["bwk"]]
bwk_opp     <- lijst_oppervlaktes[["bwk"]]
not_bwk_max <- lijst_matches[["not_bwk"]]
not_bwk_opp <- lijst_oppervlaktes[["not_bwk"]]

rm(tabel_vlaanderen, vertaal_df, lijst_matches, lijst_oppervlaktes)
gc()

# ==============================================================================
# 7. HYDROGRAFISCHE DATA (HUET-ZONES) INLADEN EN SYNCHRONISEREN
# ==============================================================================
message("-> Huet-zone waterlopen inlezen...")
r_huetzon_vlaanderen <- rast(here("data/input/Raster_Vlaanderen/vlaanderen_huetzon_10m.tif"))

r_huetzon_local <- r_huetzon_vlaanderen %>% 
  terra::crop(area_buffer_fix) %>% 
  terra::resample(template_DM, method = "near")

legende_huet  <- terra::cats(r_huetzon_local)[[1]]
target_labels <- c("P1", "P2")

if(!is.null(legende_huet) && is.data.frame(legende_huet) && "Label" %in% colnames(legende_huet)) {
  valide_ids <- legende_huet$value[legende_huet$Label %in% target_labels]
  message(paste("-> Gekoppelde Huet-zone ID's op basis van labels:", paste(valide_ids, collapse = ", ")))
} else {
  valide_ids <- 5:10 
  message("⚠️ Waarschuwing: Rasterlegende-structuur wijkt af, terugvallen op harde ID's 5:10.")
}

r_waterlopen_bin <- r_huetzon_local %in% valide_ids
r_waterlopen_bin <- terra::ifel(r_waterlopen_bin == 1, 1, NA)

# ==============================================================================
# 7b. PLASSEN INLADEN EN NORMEREN (EXCL. NOT_BWK / AP)
# ==============================================================================
message("-> Plassen raster inlezen en filteren op not_bwk...")
r_plassen_vlaanderen <- rast(here("data/input/Raster_Vlaanderen/vlaanderen_watervlakken_2024_10m.tif"))

r_plassen_local <- r_plassen_vlaanderen %>% 
  terra::crop(area_buffer_fix) %>% 
  terra::resample(template_DM, method = "near")

r_plassen_bin <- terra::ifel(r_plassen_local == 1, 1, NA)

if (exists("not_bwk_max") && !is.null(not_bwk_max)) {
  r_plassen_geschikt_max <- terra::ifel(!is.na(not_bwk_max) & not_bwk_max > 0, NA, r_plassen_bin)
  r_not_bwk_zero <- terra::ifel(is.na(not_bwk_opp), 0, not_bwk_opp)
  r_rest_opp <- terra::clamp(1 - r_not_bwk_zero, lower = 0, upper = 1)
  r_plassen_geschikt_opp <- terra::ifel(!is.na(r_plassen_bin) & r_rest_opp > 0, r_rest_opp, NA)
} else {
  r_plassen_geschikt_max <- r_plassen_bin
  r_plassen_geschikt_opp <- r_plassen_bin
}

# ==============================================================================
# 7c. DRAINAGEKAART INLADEN EN FILTEREN OP NATTE GRONDEN
# ==============================================================================
message("-> Drainagekaart inlezen en natte bodems filteren op categorieën...")

r_drain_raw   <- rast(here("data/input/Raster_Vlaanderen/vlaanderen_ovstrg_10m.tif"))
r_drain_local <- r_drain_raw %>% 
  terra::crop(area_buffer_fix) %>% 
  terra::resample(template_DM, method = "near")

valide_drain_ids <- c(1,2)

r_drain_bin <- r_drain_local %in% valide_drain_ids
r_natte_gronden_bin <- terra::ifel(r_drain_bin == 1, 1, NA)

# ==============================================================================
# 8. PARALLELLE FUSIE: SPOOR A (MAX) EN SPOOR B (OPP) MET DRAINAGEKAART
# ==============================================================================
message("-> Parallelle fusie uitvoeren (BWK + Plassen + Waterlopen x Natte Drainagegronden)...")

leefgebied1_max <- terra::cover(bwk_max, r_plassen_geschikt_max)
leefgebied1_max <- terra::cover(leefgebied1_max, r_waterlopen_bin)

r_waterlopen_opp <- terra::ifel(r_waterlopen_bin == 1, 1, NA)
leefgebied1_opp  <- terra::cover(bwk_opp, r_plassen_geschikt_opp)
leefgebied1_opp  <- terra::cover(leefgebied1_opp, r_waterlopen_opp)

grotemodderkruiper_voortplanting_max <- terra::mask(leefgebied1_max, r_natte_gronden_bin)
grotemodderkruiper_voortplanting_opp <- terra::mask(leefgebied1_opp, r_natte_gronden_bin)

names(grotemodderkruiper_voortplanting_max) <- "Match_Max"
names(grotemodderkruiper_voortplanting_opp) <- "Oppervlakte_Real"

alle_matches      <- grotemodderkruiper_voortplanting_max
alle_oppervlaktes <- grotemodderkruiper_voortplanting_opp

rm(r_plassen_vlaanderen, r_plassen_local, r_plassen_bin, r_plassen_geschikt_max, r_plassen_geschikt_opp,
   r_drain_raw, r_drain_local, r_drain_bin, r_natte_gronden_bin, leefgebied1_max, leefgebied1_opp)
gc()

# ==============================================================================
# STAP 11: RUIMTELIJKE CLUSTERING & OPPERVLAKTEFILTER
# ==============================================================================
message("-> Start ruimtelijke clustering (10m drempel)...")

w_matrix <- matrix(1, nrow = 3, ncol = 3) 

message("-> Verwerken Spoor A (Theorie / Maximale Potentie)...")
r_bridge_max <- terra::focal(alle_matches, w = w_matrix, fun = "max", na.rm = TRUE)
cl_id_max    <- terra::patches(r_bridge_max, directions = 8, zeroAsNA = TRUE)

habitat_pixels_max <- alle_matches * 0.01
stats_max <- terra::zonal(habitat_pixels_max, cl_id_max, fun = "sum", na.rm = TRUE)
colnames(stats_max) <- c("Cluster_ID", "Potentie_ha")

valide_cluster_ids_max <- stats_max$Cluster_ID[stats_max$Potentie_ha >= 0.04]

if(length(valide_cluster_ids_max) > 0) {
  masker_max <- terra::ifel(cl_id_max %in% valide_cluster_ids_max, 1, NA)
  grotemodderkruiper_leefgebied_max <- terra::mask(alle_matches, masker_max)
} else {
  grotemodderkruiper_leefgebied_max <- terra::rast(template_DM, vals = NA)
}
names(grotemodderkruiper_leefgebied_max) <- "Leefgebied_Max"

message("-> Verwerken Spoor B (Werkelijkheid / Reële Fracties)...")
r_binair_opp <- terra::ifel(!is.na(alle_oppervlaktes) & alle_oppervlaktes > 0, 1, NA)

r_bridge_opp <- terra::focal(r_binair_opp, w = w_matrix, fun = "max", na.rm = TRUE)
cl_id_opp    <- terra::patches(r_bridge_opp, directions = 8, zeroAsNA = TRUE)

habitat_pixels_opp <- alle_oppervlaktes * 0.01
stats_opp <- terra::zonal(habitat_pixels_opp, cl_id_opp, fun = "sum", na.rm = TRUE)
colnames(stats_opp) <- c("Cluster_ID", "Werkelijke_ha")

valide_cluster_ids_opp <- stats_opp$Cluster_ID[stats_opp$Werkelijke_ha >= 0.04]

if(length(valide_cluster_ids_opp) > 0) {
  masker_opp <- cl_id_opp %in% valide_cluster_ids_opp
  grotemodderkruiper_leefgebied_opp <- terra::mask(alle_oppervlaktes, terra::ifel(masker_opp, 1, NA))
} else {
  grotemodderkruiper_leefgebied_opp <- terra::rast(template_DM, vals = NA)
}
names(grotemodderkruiper_leefgebied_opp) <- "Leefgebied_Opp"

rm(r_bridge_max, r_bridge_opp, cl_id_max, cl_id_opp, masker_max, masker_opp, stats_max, stats_opp)
gc()

# ==============================================================================
# DEFINITIEVE TOEWIJSING EXPORT VARIABELEN
# ==============================================================================
final_max <- grotemodderkruiper_leefgebied_max
final_opp <- grotemodderkruiper_leefgebied_opp

if (!all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  cl_opp <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE)
} else {
  cl_opp <- template_DM * NA
}

if (!all(is.na(suppressWarnings(terra::minmax(final_max))))) {
  cl_max <- terra::patches(final_max, directions = 8, zeroAsNA = TRUE)
} else {
  cl_max <- template_DM * NA
}

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
potentie_export_rast <- if (exists("final_max") && !is.null(final_max) && !all(is.na(suppressWarnings(terra::minmax(final_max))))) {
  terra::ifel(!is.na(final_max) & final_max > 0, 1, NA)
} else {
  terra::rast(template_DM, vals = NA)
}

# 2. Bepaal Werkelijke Oppervlakte Raster
werkelijk_export_rast <- if (exists("final_opp") && !is.null(final_opp) && !all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  terra::ifel(!is.na(final_opp) & final_opp > 0, 1, NA)
} else {
  terra::rast(template_DM, vals = NA)
}

# 3. Bepaal Analytisch Metacluster ID-raster (EXCLUSIEF OP BASIS VAN WERKELIJKE OPPERVLAKTES)
if (exists("cl_opp") && !is.null(cl_opp) && !all(is.na(suppressWarnings(terra::minmax(cl_opp))))) {
  id_export_rast <- cl_opp
} else if (exists("final_opp") && !is.null(final_opp) && !all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  id_export_rast <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE)
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
