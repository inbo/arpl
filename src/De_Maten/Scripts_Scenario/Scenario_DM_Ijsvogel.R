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

soort <- "ijsvogel"

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

grens_web <- sf::st_as_sf(terra::project(area_shape, "EPSG:4326"))

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

bwk_max <- lijst_matches[["bwk"]]
bwk_opp <- lijst_oppervlaktes[["bwk"]]

rm(resultaten_gegroepeerd, df_nieuw, lijst_matches, lijst_oppervlaktes)
gc()

# Spoor A: MAX
bwk_clusters_max <- cluster_filter_compleet(
  masker     = bwk_max,
  opp_laag   = bwk_opp,
  drempel_m2 = 50000,   # 5 ha
  dist_m     = 50,
  werkelijk  = FALSE
)

# Spoor B: OPP
r_binair_bwk_opp <- terra::ifel(!is.na(bwk_opp) & bwk_opp > 0, 1, NA)
bwk_clusters_opp <- cluster_filter_compleet(
  masker     = r_binair_bwk_opp,
  opp_laag   = bwk_opp,
  drempel_m2 = 50000,   # 5 ha
  dist_m     = 50,
  werkelijk  = TRUE
)

rm(r_binair_bwk_opp)
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
target_labels <- c("Bg", "BgK", "Bk1", "Bk2", "BkK1", "BkK2", "KabPKunst", "Kan", "KanPKunst", "Kunst", "Mlz", "O1brak", "O1o", "O2zout", "P1", "P2", "Rg", "Rk", "Rzg")

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
# 8. PARALLELLE FUSIE
# ==============================================================================
message("-> Strikt gescheiden parallelle fusie uitvoeren...")

ijsvogel_voortplanting_max <- terra::cover(bwk_clusters_max$raster, r_waterlopen_bin)
names(ijsvogel_voortplanting_max) <- "Match_Max"

r_waterlopen_opp <- terra::ifel(r_waterlopen_bin == 1, 0.01, NA)
ijsvogel_voortplanting_opp <- terra::cover(bwk_clusters_opp$raster, r_waterlopen_opp)
names(ijsvogel_voortplanting_opp) <- "Oppervlakte_Real"

writeRaster(ijsvogel_voortplanting_max, "temp_ijsvogel_fusie_max.tif", overwrite = TRUE)
writeRaster(ijsvogel_voortplanting_opp, "temp_ijsvogel_fusie_opp.tif", overwrite = TRUE)

rm(ijsvogel_voortplanting_max, ijsvogel_voortplanting_opp, bwk_clusters_max, bwk_clusters_opp, 
   tabel_vlaanderen, tabel_gefilterd, tabel_DM, tabel_cel_som, vertaal_df,
   r_huetzon_vlaanderen, r_huetzon_local, r_waterlopen_bin, r_waterlopen_opp)
gc()

alle_matches      <- rast("temp_ijsvogel_fusie_max.tif")
alle_oppervlaktes <- rast("temp_ijsvogel_fusie_opp.tif")
names(alle_matches)      <- "Match_Max"
names(alle_oppervlaktes) <- "Oppervlakte_Real"

# ==============================================================================
# 9. GEOPTIMALISEERDE OEVERANALYSE
# ==============================================================================
message("-> Start geoptimaliseerde oeveranalyse (exacte lengteberekening)...")

blauwe_laag_raw <- rast(here("data/input/Raster_Vlaanderen/vlaanderen_watervlakken_2024_10m.tif"))
r_blauwe_laag   <- terra::resample(blauwe_laag_raw, template_DM, method = "near")

r_water <- r_blauwe_laag == 1
r_land  <- r_blauwe_laag == 0

ijsvogel_waterkant <- terra::boundaries(r_water, classes = FALSE, inner = TRUE)
ijsvogel_landkant  <- terra::boundaries(r_land, classes = FALSE, inner = TRUE)

ijsvogel_oever <- (ijsvogel_waterkant == 1) | (ijsvogel_landkant == 1)
ijsvogel_oever <- terra::ifel(ijsvogel_oever == 1, 1, NA)

rm(r_water, r_land, ijsvogel_waterkant, ijsvogel_landkant)

temp_oever_mask        <- terra::buffer(ijsvogel_oever, width = 25)
ijsvogel_oever_cluster <- terra::patches(temp_oever_mask, directions = 8, zeroAsNA = TRUE) %>% 
  terra::mask(ijsvogel_oever)

oever_polys <- terra::as.polygons(ijsvogel_oever_cluster, dissolve = TRUE)
oever_polys$lengte_m <- terra::perim(oever_polys) / 2

valide_oever_clusters <- oever_polys$patches[oever_polys$lengte_m >= 1000]

if(length(valide_oever_clusters) > 0) {
  ijsvogel_oever1 <- ijsvogel_oever_cluster %in% valide_oever_clusters
  ijsvogel_oever1 <- terra::ifel(ijsvogel_oever1 == 1, 1, NA)
} else {
  ijsvogel_oever1 <- template_DM * NA
  warning("⚠️ Geen enkele oever voldoet aan de minimale lengte van 1 km!")
}

dist_to_oever  <- terra::distance(ijsvogel_oever1)
oever_zone_20m <- terra::ifel(dist_to_oever <= 20, 1, NA)

ijsvogel_leefgebied_max_raw <- terra::mask(alle_matches, oever_zone_20m)
ijsvogel_leefgebied_opp_raw <- terra::mask(alle_oppervlaktes, oever_zone_20m)

writeRaster(ijsvogel_leefgebied_max_raw, "temp_ijsvogel_leefgebied_max.tif", overwrite = TRUE)
writeRaster(ijsvogel_leefgebied_opp_raw, "temp_ijsvogel_leefgebied_opp.tif", overwrite = TRUE)

rm(alle_matches, alle_oppervlaktes, blauwe_laag_raw, r_blauwe_laag, temp_oever_mask, 
   ijsvogel_oever_cluster, oever_polys, valide_oever_clusters, dist_to_oever, oever_zone_20m,
   ijsvogel_leefgebied_max_raw, ijsvogel_leefgebied_opp_raw)
gc()

ijsvogel_leefgebied_max <- rast("temp_ijsvogel_leefgebied_max.tif")
ijsvogel_leefgebied_opp <- rast("temp_ijsvogel_leefgebied_opp.tif")

names(ijsvogel_leefgebied_max) <- "Leefgebied_Max"
names(ijsvogel_leefgebied_opp) <- "Leefgebied_Real"

# DEFINITIEVE MODELUITGANGEN
final_max <- ijsvogel_leefgebied_max
final_opp <- ijsvogel_leefgebied_opp

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

# --- EINDOPRUIMING TIJDELIJKE CACHE BESTANDEN ---
temp_files <- c(
  "temp_ijsvogel_fusie_max.tif", 
  "temp_ijsvogel_fusie_opp.tif",
  "temp_ijsvogel_leefgebied_max.tif", 
  "temp_ijsvogel_leefgebied_opp.tif"
)
file.remove(temp_files[file.exists(temp_files)])

suppressWarnings(
  rm(potentie_export_rast, werkelijk_export_rast, id_export_rast, export_config, export_rast, export_id_rast)
)
gc()

message(paste("🏁 SCENARIO EXPORT VOLLEDIG AFGEROND VOOR:", toupper(soort)))
