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

soort <- "zwartespecht"

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

bwk_max <- lijst_matches[["bwk"]]
bwk_opp <- lijst_oppervlaktes[["bwk"]]

rm(vertaal_df, lijst_matches, lijst_oppervlaktes, tabel_vlaanderen, resultaten_gegroepeerd, df_nieuw)
gc()

# ==============================================================================
# STAP 1: BASISBIOTOOP EN HOOGGROEN COUPLING (100M FUZZY - PARALLEL)
# ==============================================================================
message("-> Groenkaart inlezen en synchroniseren voor de bosstructuur...")

r_groen_raw <- rast(here("data/input/ASCI Files/Groenkaart_2021.tif"))
r_groen     <- terra::resample(r_groen_raw, id_raster_TV, method = "near") %>% terra::crop(id_raster_TV)

zwartespecht_groen_ok <- r_groen == 1

# --- SPOOR A: MAX ---
zwartespecht_basis_max <- (bwk_max == 1) | zwartespecht_groen_ok
zwartespecht_basis_max <- terra::ifel(zwartespecht_basis_max == 1, 1, NA)

zwartespecht_clusters_max <- cluster_filter_compleet(
  masker     = zwartespecht_basis_max,
  opp_laag   = zwartespecht_basis_max,
  drempel_m2 = 0,   
  dist_m     = 100,        
  werkelijk  = FALSE      
)
r_bos_patches_max <- zwartespecht_clusters_max$clusters

# --- SPOOR B: OPP (Parallel & Autonoom) ---
r_groen_opp_clean <- terra::ifel(!is.na(zwartespecht_groen_ok) & zwartespecht_groen_ok == 1, 1.0, 0)
b_clean           <- terra::ifel(is.na(bwk_opp), 0, bwk_opp)
bos_som_raw       <- b_clean + r_groen_opp_clean
zwartespecht_basis_opp <- terra::clamp(bos_som_raw, upper = 1.0)
zwartespecht_basis_opp <- terra::ifel(zwartespecht_basis_opp > 0, zwartespecht_basis_opp, NA)

r_binair_basis_opp <- terra::ifel(!is.na(zwartespecht_basis_opp) & zwartespecht_basis_opp > 0, 1, NA)

zwartespecht_clusters_opp <- cluster_filter_compleet(
  masker     = r_binair_basis_opp,
  opp_laag   = zwartespecht_basis_opp,
  drempel_m2 = 0,   
  dist_m     = 100,        
  werkelijk  = TRUE       
)
r_bos_patches_opp <- zwartespecht_clusters_opp$clusters

rm(r_groen_raw, r_groen, zwartespecht_groen_ok, r_groen_opp_clean, b_clean, bos_som_raw, 
   r_binair_basis_opp, zwartespecht_clusters_max, zwartespecht_clusters_opp)
gc()

# ==============================================================================
# STAP 2: CLASSIFICATIE BOSGROOTTE (>300 HA VS. 50-300 HA - PARALLEL)
# ==============================================================================
r_bos_groot_300ha_max <- template_TV * NA
r_bos_klein_50ha_max  <- template_TV * NA

if (!all(is.na(suppressWarnings(terra::minmax(r_bos_patches_max))))) {
  stats_ha_max <- terra::freq(r_bos_patches_max)
  stats_ha_max$Grootte_ha <- stats_ha_max$count * 0.01
  
  ids_300ha_max    <- stats_ha_max$value[stats_ha_max$Grootte_ha >= 300]
  ids_50_300ha_max <- stats_ha_max$value[stats_ha_max$Grootte_ha >= 50 & stats_ha_max$Grootte_ha < 300]
  
  if (length(ids_300ha_max) > 0)    r_bos_groot_300ha_max <- terra::ifel(r_bos_patches_max %in% ids_300ha_max, 1, NA)
  if (length(ids_50_300ha_max) > 0) r_bos_klein_50ha_max  <- terra::ifel(r_bos_patches_max %in% ids_50_300ha_max, 1, NA)
  
  rm(stats_ha_max)
}

r_bos_groot_300ha_opp <- template_TV * NA
r_bos_klein_50ha_opp  <- template_TV * NA

if (!all(is.na(suppressWarnings(terra::minmax(r_bos_patches_opp))))) {
  stats_ha_opp <- terra::zonal(zwartespecht_basis_opp, r_bos_patches_opp, fun = "sum", na.rm = TRUE)
  colnames(stats_ha_opp) <- c("ID", "Grootte_ha")
  stats_ha_opp$Grootte_ha <- stats_ha_opp$Grootte_ha * 0.01
  
  ids_300ha_opp    <- stats_ha_opp$ID[stats_ha_opp$Grootte_ha >= 300]
  ids_50_300ha_opp <- stats_ha_opp$ID[stats_ha_opp$Grootte_ha >= 50 & stats_ha_opp$Grootte_ha < 300]
  
  if (length(ids_300ha_opp) > 0)    r_bos_groot_300ha_opp <- terra::mask(zwartespecht_basis_opp, r_bos_patches_opp %in% ids_300ha_opp)
  if (length(ids_50_300ha_opp) > 0) r_bos_klein_50ha_opp  <- terra::mask(zwartespecht_basis_opp, r_bos_patches_opp %in% ids_50_300ha_opp)
  
  rm(stats_ha_opp)
}
gc()

# ==============================================================================
# STAP 3: RECONSTRUCTIE NETWERK (4 KM BUFFER INTERACTIE - PARALLEL)
# ==============================================================================
message("-> Netwerkrelaties controleren binnen een actieradius van 4 km...")

# --- SPOOR A: MAX ---
r_bos_klein_goedgekeurd_max <- template_TV * NA

if (!all(is.na(suppressWarnings(terra::minmax(r_bos_klein_50ha_max))))) {
  buffer_4km_klein_max <- terra::buffer(r_bos_klein_50ha_max, width = 4000)
  
  if (!all(is.na(suppressWarnings(terra::minmax(r_bos_groot_300ha_max))))) {
    buffer_4km_groot_max <- terra::buffer(r_bos_groot_300ha_max, width = 4000)
    netwerk_masker_max <- !is.na(buffer_4km_groot_max) | !is.na(buffer_4km_klein_max)
    suppressWarnings(rm(buffer_4km_groot_max))
  } else {
    netwerk_masker_max <- !is.na(buffer_4km_klein_max)
  }
  
  r_bos_klein_goedgekeurd_max <- terra::mask(r_bos_klein_50ha_max, netwerk_masker_max)
  suppressWarnings(rm(buffer_4km_klein_max, netwerk_masker_max))
}

leefgebied1_max <- !is.na(r_bos_groot_300ha_max) | !is.na(r_bos_klein_goedgekeurd_max)
leefgebied1_max <- terra::ifel(leefgebied1_max == 1, 1, NA)

# --- SPOOR B: OPP ---
r_bos_klein_goedgekeurd_opp <- template_TV * NA

if (!all(is.na(suppressWarnings(terra::minmax(r_bos_klein_50ha_opp))))) {
  r_binair_klein_opp  <- terra::ifel(!is.na(r_bos_klein_50ha_opp) & r_bos_klein_50ha_opp > 0, 1, NA)
  buffer_4km_klein_opp <- terra::buffer(r_binair_klein_opp, width = 4000)
  
  if (!all(is.na(suppressWarnings(terra::minmax(r_bos_groot_300ha_opp))))) {
    r_binair_groot_opp  <- terra::ifel(!is.na(r_bos_groot_300ha_opp) & r_bos_groot_300ha_opp > 0, 1, NA)
    buffer_4km_groot_opp <- terra::buffer(r_binair_groot_opp, width = 4000)
    netwerk_masker_opp   <- !is.na(buffer_4km_groot_opp) | !is.na(buffer_4km_klein_opp)
    suppressWarnings(rm(r_binair_groot_opp, buffer_4km_groot_opp))
  } else {
    netwerk_masker_opp <- !is.na(buffer_4km_klein_opp)
  }
  
  r_bos_klein_goedgekeurd_opp <- terra::mask(r_bos_klein_50ha_opp, netwerk_masker_opp)
  suppressWarnings(rm(r_binair_klein_opp, buffer_4km_klein_opp, netwerk_masker_opp))
}

g_clean <- terra::ifel(is.na(r_bos_groot_300ha_opp), 0, r_bos_groot_300ha_opp)
k_clean <- terra::ifel(is.na(r_bos_klein_goedgekeurd_opp), 0, r_bos_klein_goedgekeurd_opp)
som_raw <- g_clean + k_clean
som_cl  <- terra::clamp(som_raw, upper = 1.0)
leefgebied1_opp <- terra::ifel(som_cl > 0, som_cl, NA)

rm(r_bos_groot_300ha_max, r_bos_klein_50ha_max, r_bos_klein_goedgekeurd_max,
   r_bos_groot_300ha_opp, r_bos_klein_50ha_opp, r_bos_klein_goedgekeurd_opp, g_clean, k_clean, som_raw, som_cl)
gc()

# ==============================================================================
# STAP 4: FINALE COMPACTHEIDS-CHECK EN SAMENSMELTING (PARALLEL)
# ==============================================================================
message("-> Finale 100m clustering en 100 ha minimumoppervlakte check...")

# --- SPOOR A: MAX ---
finaal_leefgebied_clusters_max <- cluster_filter_compleet(
  masker     = leefgebied1_max,
  opp_laag   = leefgebied1_max,
  drempel_m2 = 500000, 
  dist_m     = 100,     
  werkelijk  = FALSE    
)
zwartespecht_leefgebied_max <- terra::crop(finaal_leefgebied_clusters_max$raster, template_TV)

# --- SPOOR B: OPP ---
r_binair_leef1_opp <- terra::ifel(!is.na(leefgebied1_opp) & leefgebied1_opp > 0, 1, NA)

finaal_leefgebied_clusters_opp <- cluster_filter_compleet(
  masker     = r_binair_leef1_opp,
  opp_laag   = leefgebied1_opp,
  drempel_m2 = 500000, 
  dist_m     = 100,     
  werkelijk  = TRUE    
)
zwartespecht_leefgebied_opp <- terra::crop(finaal_leefgebied_clusters_opp$raster, template_TV)

# DEFINITIEVE MODELUITGANGEN
final_max <- zwartespecht_leefgebied_max
final_opp <- zwartespecht_leefgebied_opp

if (!all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  cl_opp <- finaal_leefgebied_clusters_opp$clusters
} else {
  cl_opp <- template_TV * NA
}

if (!all(is.na(suppressWarnings(terra::minmax(final_max))))) {
  cl_max <- finaal_leefgebied_clusters_max$clusters
} else {
  cl_max <- template_TV * NA
}

rm(zwartespecht_basis_max, zwartespecht_basis_opp, r_bos_patches_max, r_bos_patches_opp, 
   leefgebied1_max, leefgebied1_opp, r_binair_leef1_opp, finaal_leefgebied_clusters_max, finaal_leefgebied_clusters_opp)
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
