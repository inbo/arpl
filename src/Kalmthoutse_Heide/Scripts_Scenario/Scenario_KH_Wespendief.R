library(here)
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
soort <- "wespendief"

# --- DYNAMISCHE SCENARIO PARAMETER CHECK ---
if (exists("SCENARIO_RDS_PAD") && !is.null(SCENARIO_RDS_PAD)) {
  scenario_rds_path <- SCENARIO_RDS_PAD
} else if (exists("params") && !is.null(params$scenario_rds_path)) {
  scenario_rds_path <- params$scenario_rds_path
} else {
  scenario_rds_path <- "data/input/Scenario_rds/KH_Scenario_BWK_2025.rds"
}

p_raw <- gsub("^([.][.]/)+", "", scenario_rds_path)
scenario_path <- here::here(p_raw)

if (!file.exists(scenario_path)) {
  stop(paste("❌ FOUT: Scenario RDS bestand NIET gevonden op:", scenario_path))
}

scen_volledig <- basename(scenario_path)
scenario_naam <- gsub("^KH_Scenario_|^Scenario_|.rds$", "", scen_volledig)

message(paste("Verwerken van soort:", soort, "binnen scenario:", scenario_naam))

resultaat <- df %>%
  filter(tolower(trimws(Soort)) == soort) %>%
  select(Type, MinOpp_ha, AfstandBiotopen_m, Dispersiecap_m)

oppervlakte_ha <- resultaat$MinOpp_ha[1]
afstand_m      <- resultaat$AfstandBiotopen_m[1]
buffer_m       <- resultaat$Dispersiecap_m[1]

rm(df, resultaat)

area_shape  <- vect(here("data/input/Kalmthoutse_Heide.shp"))
master_grid <- rast(here("data/input/Raster_Vlaanderen/Vlaanderen_MasterGrid_10m.tif"))[[1]]

df_namen_sleutel <- read_csv(here("data/input/Excel_files/BWK_Laag_Namen_2025.csv"), show_col_types = FALSE)
gouden_namenlijst <- tolower(trimws(df_namen_sleutel$Laagnaam))

area_shape_proj <- project(area_shape, crs(master_grid))
area_buffer_fix <- buffer(area_shape_proj, width = buffer_m)

message("-> Vertaalraster voor globale/lokale cellen opbouwen via snelle MASK methode...")
id_raster_KH <- crop(master_grid, area_buffer_fix, snap = "near")

globale_id_raster <- master_grid
globale_id_raster <- terra::init(globale_id_raster, fun = "cell")

id_raster_KH_globale_values <- crop(globale_id_raster, area_buffer_fix, snap = "near")
id_raster_KH_masked <- mask(id_raster_KH_globale_values, area_buffer_fix)

message("-> Vertaaltabel bliksemsnel opbouwen via C++ dataframe extractie...")

df_extractie <- as.data.frame(id_raster_KH_masked, cells = TRUE)
vertaal_df <- as.data.table(df_extractie)
setnames(vertaal_df, c(1, 2), c("lokale_id", "globale_id"))

vertaal_df <- vertaal_df[!is.na(globale_id)]
studiegebied_globale_ids <- unique(vertaal_df$globale_id)

values(id_raster_KH) <- NA
template_KH <- terra::rasterize(area_buffer_fix, id_raster_KH, field = 1, background = 0)

grens_web <- sf::st_as_sf(terra::project(area_shape, "EPSG:4326"))

rm(globale_id_raster, id_raster_KH_globale_values, id_raster_KH_masked, df_extractie)
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
  
  tabel_KH <- tabel_gefilterd[cel_id %in% studiegebied_globale_ids]
  tabel_KH_unique <- unique(tabel_KH, by = c("cel_id", "CODE"))

  tabel_cel_som <- tabel_KH_unique[, .(Oppervlakte = pmin(sum(BWK_FRAC, na.rm = TRUE), 1.0)), by = .(cel_id)]
  
  r_match_type <- id_raster_KH * NA
  r_opp_type   <- id_raster_KH * NA
  
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
foerageer_open_max    <- lijst_matches[["foerageer_open"]]
foerageer_open_opp    <- lijst_oppervlaktes[["foerageer_open"]]

rm(vertaal_df, lijst_matches, lijst_oppervlaktes)
gc()

# ==============================================================================
# STAP 1: VOORTPLANTINGSGEBIED (MINIMAAL 50 HA, 100M CLUSTERING)
# ==============================================================================
message("-> Groenkaart inlezen en synchroniseren als fallback voor oude bossen...")

r_groen_raw <- rast(here("data/input/ASCI Files/Groenkaart_2021.tif"))
r_groen     <- terra::resample(r_groen_raw, id_raster_KH, method = "near") %>% terra::crop(id_raster_KH)

wespendief_hooggroen_ok <- r_groen == 1

# --- SPOOR A: MAX ---
wespendief_voort_basis_max <- (voortplanting_bwk_max == 1) & wespendief_hooggroen_ok
wespendief_voort_basis_max <- terra::ifel(wespendief_voort_basis_max == 1, 1, NA)

voort_clusters_max <- cluster_filter_compleet(
  masker     = wespendief_voort_basis_max,
  opp_laag   = voortplanting_bwk_opp,
  drempel_m2 = 500000,    
  dist_m     = 100,       
  werkelijk  = FALSE      
)
r_wespendief_bos_max <- voort_clusters_max$raster

# --- SPOOR B: OPP ---
wespendief_voort_basis_opp_raw <- terra::mask(voortplanting_bwk_opp, wespendief_hooggroen_ok)
wespendief_voort_basis_opp     <- terra::ifel(!is.na(wespendief_voort_basis_opp_raw) & wespendief_voort_basis_opp_raw > 0, wespendief_voort_basis_opp_raw, NA)

r_binair_voort_opp <- terra::ifel(!is.na(wespendief_voort_basis_opp) & wespendief_voort_basis_opp > 0, 1, NA)
voort_clusters_opp <- cluster_filter_compleet(
  masker     = r_binair_voort_opp,
  opp_laag   = wespendief_voort_basis_opp,
  drempel_m2 = 500000,    
  dist_m     = 100,       
  werkelijk  = TRUE      
)
r_wespendief_bos_opp <- voort_clusters_opp$raster

rm(r_groen_raw, wespendief_hooggroen_ok, wespendief_voort_basis_max, wespendief_voort_basis_opp_raw, 
   wespendief_voort_basis_opp, r_binair_voort_opp, voort_clusters_max, voort_clusters_opp)
gc()

# ==============================================================================
# STAP 2: OPEN FOERAGEERGEBIED (MINIMAAL 50 HA, 500M CLUSTERING)
# ==============================================================================
message("-> Open foerageergebieden groeperen (500m fuzzy) en filteren op min. 50 ha...")

open_foer_clusters_max <- cluster_filter_compleet(
  masker     = foerageer_open_max,
  opp_laag   = foerageer_open_opp,
  drempel_m2 = 500000,    
  dist_m     = 500,       
  werkelijk  = FALSE      
)
r_wespendief_open_max <- open_foer_clusters_max$raster

r_binair_open_opp <- terra::ifel(!is.na(foerageer_open_opp) & foerageer_open_opp > 0, 1, NA)
open_foer_clusters_opp <- cluster_filter_compleet(
  masker     = r_binair_open_opp,
  opp_laag   = foerageer_open_opp,
  drempel_m2 = 500000,    
  dist_m     = 500,       
  werkelijk  = TRUE      
)
r_wespendief_open_opp <- open_foer_clusters_opp$raster

rm(r_binair_open_opp, open_foer_clusters_max, open_foer_clusters_opp)
gc()

# ==============================================================================
# STAP 3: OPTIMAAL FOERAGEERCOMPLEX (MINIMAAL 250 HA, 500M)
# ==============================================================================
message("-> Bos en open gebieden samenvoegen tot gecombineerd foerageerlandschap...")

foerageer_basis_max <- !is.na(r_wespendief_bos_max) | !is.na(r_wespendief_open_max)
foerageer_basis_max <- terra::ifel(foerageer_basis_max == 1, 1, NA)

foerageer_complex_clusters_max <- cluster_filter_compleet(
  masker     = foerageer_basis_max,
  opp_laag   = foerageer_basis_max,
  drempel_m2 = 2500000, 
  dist_m     = 500,     
  werkelijk  = FALSE    
)
r_foerageer2_max <- foerageer_complex_clusters_max$raster

b_clean <- terra::ifel(is.na(r_wespendief_bos_opp), 0, r_wespendief_bos_opp)
o_clean <- terra::ifel(is.na(r_wespendief_open_opp), 0, r_wespendief_open_opp)
foer_som_raw <- b_clean + o_clean
foerageer_basis_opp <- terra::clamp(foer_som_raw, upper = 1.0)
foerageer_basis_opp <- terra::ifel(foerageer_basis_opp > 0, foerageer_basis_opp, NA)

r_binair_foer_basis_opp <- terra::ifel(!is.na(foerageer_basis_opp) & foerageer_basis_opp > 0, 1, NA)

foerageer_complex_clusters_opp <- cluster_filter_compleet(
  masker     = r_binair_foer_basis_opp,
  opp_laag   = foerageer_basis_opp,
  drempel_m2 = 2500000, 
  dist_m     = 500,     
  werkelijk  = TRUE    
)
r_foerageer2_opp <- foerageer_complex_clusters_opp$raster

rm(foerageer_basis_max, b_clean, o_clean, foer_som_raw, foerageer_basis_opp, r_binair_foer_basis_opp, 
   foerageer_complex_clusters_max, foerageer_complex_clusters_opp)
gc()

# ==============================================================================
# STAP 4: ACTIERADIUS CHECK (5 KM) EN FINALE BROEDBOS SELECTIE
# ==============================================================================
message("-> Actieradiuskoppeling (5000m) berekenen tussen broedbos en foerageercomplex...")

if (!all(is.na(suppressWarnings(terra::minmax(r_foerageer2_max))))) {
  buffer_5000m_foer_max <- terra::buffer(r_foerageer2_max, width = 5000)
  wespendief_leefgebied_max <- terra::mask(r_wespendief_bos_max, buffer_5000m_foer_max) %>% terra::crop(template_KH)
  suppressWarnings(rm(buffer_5000m_foer_max))
} else {
  wespendief_leefgebied_max <- template_KH * NA
}

if (!all(is.na(suppressWarnings(terra::minmax(r_foerageer2_opp))))) {
  r_binair_foer2_opp <- terra::ifel(!is.na(r_foerageer2_opp) & r_foerageer2_opp > 0, 1, NA)
  buffer_5000m_foer_opp <- terra::buffer(r_foerageer2_opp, width = 5000)
  
  wespendief_leefgebied_opp <- terra::mask(r_wespendief_bos_opp, buffer_5000m_foer_opp) %>% terra::crop(template_KH)
  suppressWarnings(rm(r_binair_foer2_opp, buffer_5000m_foer_opp))
} else {
  wespendief_leefgebied_opp <- template_KH * NA
}

rm(r_wespendief_bos_max, r_wespendief_bos_opp, r_wespendief_open_max, r_wespendief_open_opp, 
   r_foerageer2_max, r_foerageer2_opp)
gc()

# ==============================================================================
# 5. ESSENTIËLE KOPPELING MET SCHONE EXPORT VARIABELEN (FIX VOOR LEEG RASTER)
# ==============================================================================
final_max <- wespendief_leefgebied_max
final_opp <- wespendief_leefgebied_opp

# Bepaal analytische netwerk-ID raster
if (!all(is.na(suppressWarnings(terra::minmax(final_max))))) {
  cl_max <- terra::patches(final_max, directions = 8, zeroAsNA = TRUE)
} else {
  cl_max <- template_KH * NA
}

if (!all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  cl_opp <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE)
} else {
  cl_opp <- template_KH * NA
}

# ==============================================================================
# SCHONE EXPORT BIOTOOP EN ANALYTISCH ID-RASTER (VOOR SCRIPT 2 / ARPL)
# ==============================================================================
base_dir <- here::here("data/output/Kalmthoutse_Heide/Rasters_Soorten", scenario_naam)

folders <- list(
  potentie  = file.path(base_dir, "01_Maximale_Potentie"),
  werkelijk = file.path(base_dir, "02_Werkelijke_Oppervlaktes"),
  id_raster = file.path(base_dir, "00_ID_Rasters")
)
purrr::walk(folders, ~if (!dir.exists(.x)) dir.create(.x, showWarnings = FALSE, recursive = TRUE))

# 1. Bepaal Maximale Potentie Raster
potentie_export_rast <- if (exists("final_max") && !all(is.na(suppressWarnings(terra::minmax(final_max))))) {
  terra::ifel(!is.na(final_max) & final_max > 0, 1, NA)
} else {
  terra::rast(template_KH, vals = NA)
}

# 2. Bepaal Werkelijke Oppervlakte Raster
werkelijk_export_rast <- if (exists("final_opp") && !all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  terra::ifel(!is.na(final_opp) & final_opp > 0, 1, NA)
} else {
  terra::rast(template_KH, vals = NA)
}

# 3. Bepaal Analytisch Metacluster ID-raster
if (exists("cl_opp") && !is.null(cl_opp) && !all(is.na(suppressWarnings(terra::minmax(cl_opp))))) {
  id_export_rast <- cl_opp
} else if (exists("final_opp") && !is.null(final_opp) && !all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  id_export_rast <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE)
} else {
  id_export_rast <- terra::rast(template_KH, vals = NA)
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
