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
conflicted::conflicts_prefer(terra::intersect)
conflicted::conflicts_prefer(terra::any)

calc_ha_exact <- function(r) {
  if(is.null(r)) return(0)
  if(all(is.na(terra::values(r, mat=FALSE)))) return(0)
  area_raster <- r * terra::cellSize(r, unit = "ha")
  val <- terra::global(area_raster, "sum", na.rm = TRUE)[[1]]
  return(as.numeric(val))
}

cluster_filter_compleet <- function(masker, opp_laag, drempel_m2, dist_m, werkelijk = FALSE, straal_m = 0) {
  if (terra::global(is.na(masker), "sum")[[1]] == terra::ncell(masker)) {
    return(list(raster = masker * NA, clusters = masker * NA))
  }
  
  if (exists("id_raster_KH")) {
    masker   <- terra::crop(masker, id_raster_KH, snap = "out")
    opp_laag <- terra::crop(opp_laag, id_raster_KH, snap = "out")
  }
  
  if (straal_m > 0) {
    r_binair <- terra::ifel(!is.na(masker) & masker > 0, 1, NA)
    
    if (terra::global(is.na(r_binair), "sum")[[1]] == terra::ncell(r_binair)) {
      return(list(raster = masker * NA, clusters = masker * NA))
    }
    
    v_biotoop <- terra::as.polygons(r_binair, aggregate = TRUE)
    v_eroded  <- terra::buffer(v_biotoop, width = -straal_m)
    
    if (length(v_eroded) == 0 || terra::geomtype(v_eroded) == "none") {
      return(list(raster = masker * NA, clusters = masker * NA))
    }
    
    v_dilated <- terra::buffer(v_eroded, width = straal_m)
    
    r_zuiver <- terra::rasterize(v_dilated, r_binair, field = 1, background = NA)
    r_masker_werkelijk <- terra::mask(r_zuiver, masker)
    
    rm(v_biotoop, v_eroded, v_dilated, r_zuiver)
  } else {
    r_masker_werkelijk <- masker
  }
  
  if (terra::global(is.na(r_masker_werkelijk), "sum")[[1]] == terra::ncell(r_masker_werkelijk)) {
    return(list(raster = masker * NA, clusters = masker * NA))
  }
  
  if (dist_m > 0) {
    r_binair_z <- terra::ifel(!is.na(r_masker_werkelijk) & r_masker_werkelijk > 0, 1, NA)
    r_buffered <- terra::buffer(r_binair_z, width = dist_m / 2)
    cl_network <- terra::patches(r_buffered, directions = 4, zeroAsNA = TRUE)
    cl_biotoop_only <- terra::mask(cl_network, r_masker_werkelijk)
  } else {
    cl_network <- terra::patches(r_masker_werkelijk, directions = 8, zeroAsNA = TRUE)
    cl_biotoop_only <- cl_network
  }
  
  if(werkelijk) {
    opp_laag_sub <- terra::crop(opp_laag, cl_biotoop_only)
    stats_df     <- terra::zonal(opp_laag_sub, cl_biotoop_only, fun = "sum", na.rm = TRUE)
    colnames(stats_df) <- c("ID", "Waarde")
    stats_df$Area_m2  <- stats_df$Waarde * 100 
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
  
  r_finaal  <- terra::mask(r_masker_werkelijk, final_network_mask)
  cl_finaal <- terra::mask(cl_biotoop_only, r_finaal) 
  
  r_finaal  <- terra::deepcopy(r_finaal)
  cl_finaal <- terra::deepcopy(cl_finaal)
  
  return(list(raster = r_finaal, clusters = cl_finaal))
}

terraOptions(
  memfrac = 0.8,
  tempdir = tempdir(),
  verbose = FALSE
)

df <- read_excel(here::here("data/input/Excel_files/Soorten_bwk_afstanden.xlsx"))
soort <- "bruinekiekendief"

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

message("-> Vertaalraster voor globale/lokale cellen opbouwen...")
id_raster_KH <- crop(master_grid, area_buffer_fix, snap = "near")

globale_id_raster <- master_grid
globale_id_raster <- terra::init(globale_id_raster, fun = "cell")
id_raster_KH_globale_values <- crop(globale_id_raster, area_buffer_fix, snap = "near")

id_raster_KH_masked <- mask(id_raster_KH_globale_values, area_buffer_fix)

lokale_ids  <- cells(id_raster_KH_masked) 
globale_ids <- id_raster_KH_masked[lokale_ids][[1]]

vertaal_df <- data.table(
  lokale_id  = lokale_ids,
  globale_id = globale_ids
)

vertaal_df <- vertaal_df[!is.na(globale_id)]
studiegebied_globale_ids <- unique(vertaal_df$globale_id)

values(id_raster_KH) <- NA
template_KH <- terra::rasterize(area_buffer_fix, id_raster_KH, field = 1, background = 0)

rm(globale_id_raster, id_raster_KH_globale_values, id_raster_KH_masked, lokale_ids, globale_ids)
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
bos_max               <- lijst_matches[["bos"]]
bos_opp               <- lijst_oppervlaktes[["bos"]]
foerageer_max         <- lijst_matches[["foerageer"]]
foerageer_opp         <- lijst_oppervlaktes[["foerageer"]]

rm(tabel_vlaanderen, vertaal_df, lijst_matches, lijst_oppervlaktes)
gc()

# ==============================================================================
# STAP 2: VOORTPLANTINGSGEBIED (BWK + LGP)
# ==============================================================================
r_lgp_KH <- crop(rast(here("data/input/Raster_Vlaanderen/vlaanderen_lgp_gwscod_2025_10m.tif")), id_raster_KH, snap = "near")
kiekendief_codes <- c(311, 321, 331)

bruinekiekendief_lgp <- terra::ifel(r_lgp_KH %in% kiekendief_codes, 1, NA)

bruinekiekendief_voortplanting_max <- terra::cover(bruinekiekendief_lgp, voortplanting_bwk_max)
bruinekiekendief_voortplanting_opp <- terra::cover(bruinekiekendief_lgp, voortplanting_bwk_opp)

rm(r_lgp_KH, bruinekiekendief_lgp)

bruinekiekendief_vp_zonderbos_max <- bruinekiekendief_voortplanting_max
bruinekiekendief_vp_zonderbos_opp <- bruinekiekendief_voortplanting_opp

# ==============================================================================
# STAP 3: VOORTPLANTINGSCLUSTERS (LINIE-ZUIVERING 40M + 100M FUZZY + MIN. 10 HA)
# ==============================================================================
message("-> Voortplantingsgebieden zuiveren (min. 40m breed), clusteren (100m) en filteren op min. 10 ha...")

drempel_vp_m2 <- 10 * 10000

res_vp_max <- cluster_filter_compleet(
  masker     = bruinekiekendief_vp_zonderbos_max,
  opp_laag   = bruinekiekendief_vp_zonderbos_max,
  drempel_m2 = drempel_vp_m2,
  dist_m     = 100,
  werkelijk  = FALSE,
  straal_m   = 20
)
bruinekiekendief_voortplanting_final_max <- res_vp_max$raster

r_binair_vp_opp <- terra::ifel(!is.na(bruinekiekendief_vp_zonderbos_opp) & bruinekiekendief_vp_zonderbos_opp > 0, 1, NA)
res_vp_opp <- cluster_filter_compleet(
  masker     = r_binair_vp_opp,
  opp_laag   = bruinekiekendief_vp_zonderbos_opp,
  drempel_m2 = drempel_vp_m2,
  dist_m     = 100,
  werkelijk  = TRUE,
  straal_m   = 20
)
bruinekiekendief_voortplanting_final_opp <- res_vp_opp$raster

rm(bruinekiekendief_vp_zonderbos_max, bruinekiekendief_vp_zonderbos_opp, r_binair_vp_opp)
gc()

# ==============================================================================
# STAP 4: FOERAGEERGEBIED CLUSTEREN (100M FUZZY, MIN. 100 HA)
# ==============================================================================
message("-> Foerageergebieden clusteren (100m fuzzy) en filteren op min. 100 ha...")

drempel_foer_m2 <- 100 * 10000

res_foer_max <- cluster_filter_compleet(
  masker     = foerageer_max,
  opp_laag   = foerageer_max,
  drempel_m2 = drempel_foer_m2,
  dist_m     = 100,
  werkelijk  = FALSE,
  straal_m   = 0
)
bruinekiekendief_foerageer_final_max <- res_foer_max$raster

r_binair_foer_opp <- terra::ifel(!is.na(foerageer_opp) & foerageer_opp > 0, 1, NA)
res_foer_opp <- cluster_filter_compleet(
  masker     = r_binair_foer_opp,
  opp_laag   = foerageer_opp,
  drempel_m2 = drempel_foer_m2,
  dist_m     = 100,
  werkelijk  = TRUE,
  straal_m   = 0
)
bruinekiekendief_foerageer_final_opp <- res_foer_opp$raster

rm(foerageer_max, foerageer_opp, r_binair_foer_opp)
gc()

# ==============================================================================
# STAP 5: FINALE LEEFGEBIED SYNTHESE (5 KM OMGEVINGSTOETSING)
# ==============================================================================
message("-> Starten finale leefgebied synthese (5 km omgevingstoetsing)...")

# --- 5.1: THEORETISCH MAXIMUM (MAX-SPOOR) ---
r_patches_vp_max   <- terra::patches(bruinekiekendief_voortplanting_final_max, directions = 8, zeroAsNA = TRUE)
r_patches_foer_max <- terra::patches(bruinekiekendief_foerageer_final_max, directions = 8, zeroAsNA = TRUE)

if (!all(is.na(suppressWarnings(terra::minmax(r_patches_vp_max)))) && !all(is.na(suppressWarnings(terra::minmax(r_patches_foer_max))))) {
  dist_foer_max <- terra::distance(r_patches_foer_max)
  vp_dist_max <- terra::mask(dist_foer_max, r_patches_vp_max)
  
  min_dist_df_max <- terra::zonal(vp_dist_max, r_patches_vp_max, fun = "min", na.rm = TRUE)
  colnames(min_dist_df_max) <- c("patch_id", "dist")
  
  vp_ids_valide_max <- min_dist_df_max$patch_id[min_dist_df_max$dist <= 5000]
  bruinekiekendief_leefgebied1_max <- terra::ifel(r_patches_vp_max %in% vp_ids_valide_max, 1, NA)
  
  rm(dist_foer_max, vp_dist_max, min_dist_df_max)
} else {
  bruinekiekendief_leefgebied1_max <- template_KH * NA
}

# --- 5.2: REALISTISCHE OPPERVLAKTE (OPP-SPOOR) ---
r_binair_vp_opp   <- terra::ifel(!is.na(bruinekiekendief_voortplanting_final_opp) & bruinekiekendief_voortplanting_final_opp > 0, 1, NA)
r_binair_foer_opp <- terra::ifel(!is.na(bruinekiekendief_foerageer_final_opp) & bruinekiekendief_foerageer_final_opp > 0, 1, NA)

if (!all(is.na(suppressWarnings(terra::minmax(r_binair_vp_opp)))) && !all(is.na(suppressWarnings(terra::minmax(r_binair_foer_opp))))) {
  r_patches_vp_opp   <- terra::patches(r_binair_vp_opp, directions = 8, zeroAsNA = TRUE)
  r_patches_foer_opp <- terra::patches(r_binair_foer_opp, directions = 8, zeroAsNA = TRUE)
  
  dist_foer_opp <- terra::distance(r_patches_foer_opp)
  vp_dist_opp <- terra::mask(dist_foer_opp, r_patches_vp_opp)
  
  min_dist_df_opp <- terra::zonal(vp_dist_opp, r_patches_vp_opp, fun = "min", na.rm = TRUE)
  colnames(min_dist_df_opp) <- c("patch_id", "dist")
  
  vp_ids_valide_opp <- min_dist_df_opp$patch_id[min_dist_df_opp$dist <= 5000]
  m_vp_opp <- r_patches_vp_opp %in% vp_ids_valide_opp
  
  bruinekiekendief_leefgebied1_opp <- terra::mask(bruinekiekendief_voortplanting_final_opp, terra::ifel(m_vp_opp, 1, NA))
  
  rm(r_patches_vp_opp, r_patches_foer_opp, dist_foer_opp, vp_dist_opp, min_dist_df_opp)
} else {
  bruinekiekendief_leefgebied1_opp <- template_KH * NA
}

rm(r_patches_vp_max, r_patches_foer_max, r_binair_vp_opp, r_binair_foer_opp)
gc()

# ==============================================================================
# DEFINITIEVE TOEWIJSING EXPORT VARIABELEN (CRUCIALE FIX)
# ==============================================================================
final_max <- bruinekiekendief_leefgebied1_max
final_opp <- bruinekiekendief_leefgebied1_opp

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
potentie_export_rast <- if (exists("final_max") && !is.null(final_max) && !all(is.na(suppressWarnings(terra::minmax(final_max))))) {
  terra::ifel(!is.na(final_max) & final_max > 0, 1, NA)
} else {
  terra::rast(template_KH, vals = NA)
}

# 2. Bepaal Werkelijke Oppervlakte Raster
werkelijk_export_rast <- if (exists("final_opp") && !is.null(final_opp) && !all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
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
