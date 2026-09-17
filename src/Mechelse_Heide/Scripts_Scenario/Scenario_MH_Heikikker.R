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

soort          <- "heikikker"

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

# Expert-parameters Heikikker
opp_water_ha   <- 0.025  # MinOpp waterbiotoop (0,025 ha = 250 m²)
opp_land_ha    <- 10.0   # MinOpp landbiotoop per eiland (10 ha)
buffer_m       <- 2000   # Dispersiecapaciteit (2000 m)
min_totale_ha  <- 148    # Minimale totale netwerkoppervlakte (148 ha)
min_plassen    <- 5      # Minimaal aantal geschikte plassen in het netwerk
afstand_m      <- 0      # Afstand land tot water (0 m)

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
    bevat_codes_escaped  <- gsub("([\\.\\^\\$\\*\\+\\?\\(\\)\\[\\{\\\\\\|])", "\\\\\\1", bevat_codes)
    bevat_codes_anchored <- paste0("^", bevat_codes_escaped)
    regex_term           <- paste0(bevat_codes_anchored, collapse = "|")
    
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

waterbiotoop_bwk_max <- lijst_matches[["waterbiotoop0"]]
waterbiotoop_bwk_opp <- lijst_oppervlaktes[["waterbiotoop0"]]
landbiotoop_bwk_max  <- lijst_matches[["landbiotoop1"]]
landbiotoop_bwk_opp  <- lijst_oppervlaktes[["landbiotoop1"]]

rm(tabel_vlaanderen, vertaal_df, lijst_matches, lijst_oppervlaktes)
gc()

message("-> Landbiotoop filteren op minimale eiland-oppervlakte (>= 10 ha) voor MAX en OPP parallel...")

drempel_land_m2 <- opp_land_ha * 10000 # 10 ha = 100.000 m²

# --- SPOOR MAX ---
if (!all(is.na(terra::values(landbiotoop_bwk_max, mat = FALSE)))) {
  land_patches_max <- terra::patches(landbiotoop_bwk_max, directions = 8, zeroAsNA = TRUE)
  land_freq_max    <- as.data.frame(terra::freq(land_patches_max))
  ids_land_ok_max  <- land_freq_max$value[(land_freq_max$count * 100) >= drempel_land_m2]
  
  if (length(ids_land_ok_max) > 0) {
    heikikker_land_gefilterd_max <- terra::mask(landbiotoop_bwk_max, land_patches_max %in% ids_land_ok_max)
  } else {
    heikikker_land_gefilterd_max <- template_MH * NA
  }
} else {
  heikikker_land_gefilterd_max <- template_MH * NA
}

# --- SPOOR OPP (Autonoom op basis van werkelijke bedekking) ---
if (!all(is.na(terra::values(landbiotoop_bwk_opp, mat = FALSE)))) {
  r_binair_land_opp <- terra::ifel(!is.na(landbiotoop_bwk_opp) & landbiotoop_bwk_opp > 0, 1, NA)
  land_patches_opp  <- terra::patches(r_binair_land_opp, directions = 8, zeroAsNA = TRUE)
  
  # Zonal som op reële fracties per eiland
  zonal_land_opp <- terra::zonal(landbiotoop_bwk_opp, land_patches_opp, fun = "sum", na.rm = TRUE)
  colnames(zonal_land_opp) <- c("ID", "Werk_ha")
  ids_land_ok_opp <- zonal_land_opp$ID[(zonal_land_opp$Werk_ha * 100) >= drempel_land_m2]
  
  if (length(ids_land_ok_opp) > 0) {
    heikikker_land_gefilterd_opp <- terra::mask(landbiotoop_bwk_opp, land_patches_opp %in% ids_land_ok_opp)
  } else {
    heikikker_land_gefilterd_opp <- template_MH * NA
  }
} else {
  heikikker_land_gefilterd_opp <- template_MH * NA
}

message("-> Ruimtelijke koppeling land- en waterbiotoop (AfstandBiotopen = 0 m) voor MAX en OPP parallel...")

# --- SPOOR MAX ---
if (!all(is.na(terra::values(heikikker_land_gefilterd_max, mat = FALSE))) && 
    !all(is.na(terra::values(r_heikikker_waterbiotoop1, mat = FALSE)))) {
  
  r_land_mask_max  <- heikikker_land_gefilterd_max
  r_water_mask_max <- r_heikikker_waterbiotoop1
  
  hk_land_gekoppeld_max  <- terra::mask(heikikker_land_gefilterd_max, r_water_mask_max)
  hk_water_gekoppeld_max <- terra::mask(r_heikikker_waterbiotoop1, r_land_mask_max)
  
  heikikker_leefgebied1_max <- hk_land_gekoppeld_max | hk_water_gekoppeld_max
  heikikker_leefgebied1_max <- terra::ifel(heikikker_leefgebied1_max == 1, 1, NA)
} else {
  heikikker_leefgebied1_max <- template_MH * NA
}

# --- SPOOR OPP (Volledig autonoom) ---
if (!all(is.na(terra::values(heikikker_land_gefilterd_opp, mat = FALSE))) && 
    !all(is.na(terra::values(r_heikikker_waterbiotoop1_opp, mat = FALSE)))) {
  
  r_land_mask_opp  <- terra::ifel(!is.na(heikikker_land_gefilterd_opp) & heikikker_land_gefilterd_opp > 0, 1, NA)
  r_water_mask_opp <- terra::ifel(!is.na(r_heikikker_waterbiotoop1_opp) & r_heikikker_waterbiotoop1_opp > 0, 1, NA)
  
  hk_land_gekoppeld_opp  <- terra::mask(heikikker_land_gefilterd_opp, r_water_mask_opp)
  hk_water_gekoppeld_opp <- terra::mask(r_heikikker_waterbiotoop1_opp, r_land_mask_opp)
  
  water_opp_clean <- terra::ifel(is.na(hk_water_gekoppeld_opp), 0, hk_water_gekoppeld_opp)
  land_opp_clean  <- terra::ifel(is.na(hk_land_gekoppeld_opp), 0, hk_land_gekoppeld_opp)

  som_opp <- terra::clamp(water_opp_clean + land_opp_clean, upper = 1.0)
  heikikker_leefgebied1_opp <- terra::ifel(som_opp > 0, som_opp, NA)
} else {
  heikikker_leefgebied1_opp <- template_MH * NA
}

message("-> Netwerkvorming binnen 2000m en controle op >= 148 ha EN >= 5 plassen voor MAX en OPP parallel...")

drempel_netwerk_m2 <- min_totale_ha * 10000 # 148 ha = 1.480.000 m²

# ==============================================================================
# SPOOR MAX (AUTONOOM)
# ==============================================================================
res_netwerk_max <- cluster_filter_compleet(
  masker     = heikikker_leefgebied1_max,
  opp_laag   = heikikker_leefgebied1_max,
  drempel_m2 = drempel_netwerk_m2,
  dist_m     = buffer_m,
  werkelijk  = FALSE
)

cl_netwerk_max <- res_netwerk_max$clusters
r_netwerk_max  <- res_netwerk_max$raster

if (!all(is.na(terra::values(cl_netwerk_max, mat = FALSE)))) {
  netwerk_poly_max <- terra::as.polygons(cl_netwerk_max, dissolve = TRUE) %>% terra::makeValid()
  watervlakken_in_netwerk_max <- terra::crop(watervlakken_clean, netwerk_poly_max)
  
  if (!is.null(watervlakken_in_netwerk_max) && nrow(watervlakken_in_netwerk_max) > 0) {
    rel_max <- matrix(terra::is.related(netwerk_poly_max, watervlakken_in_netwerk_max, "intersects"), 
                      nrow = nrow(netwerk_poly_max), ncol = nrow(watervlakken_in_netwerk_max))
    netwerk_poly_max$Aantal_Plassen <- rowSums(rel_max)
    id_col_max <- names(netwerk_poly_max)[1]
    valid_ids_max <- as.numeric(netwerk_poly_max[[id_col_max]][netwerk_poly_max$Aantal_Plassen >= min_plassen, 1])
  } else {
    valid_ids_max <- c()
  }
  
  if (length(na.omit(valid_ids_max)) > 0) {
    heikikker_leefgebied_max <- terra::mask(r_netwerk_max, cl_netwerk_max %in% valid_ids_max)
    cl_max                   <- terra::mask(cl_netwerk_max, cl_netwerk_max %in% valid_ids_max)
  } else {
    heikikker_leefgebied_max <- template_MH * NA
    cl_max                   <- template_MH * NA
  }
} else {
  heikikker_leefgebied_max <- template_MH * NA
  cl_max                   <- template_MH * NA
}

# ==============================================================================
# SPOOR OPP (AUTONOOM OP BASIS VAN WERKELIJKE HECTARES)
# ==============================================================================
r_binair_leef1_opp <- terra::ifel(!is.na(heikikker_leefgebied1_opp) & heikikker_leefgebied1_opp > 0, 1, NA)

res_netwerk_opp <- cluster_filter_compleet(
  masker     = r_binair_leef1_opp,
  opp_laag   = heikikker_leefgebied1_opp,
  drempel_m2 = drempel_netwerk_m2,
  dist_m     = buffer_m,
  werkelijk  = TRUE
)

cl_netwerk_opp <- res_netwerk_opp$clusters
r_netwerk_opp  <- res_netwerk_opp$raster

if (!all(is.na(terra::values(cl_netwerk_opp, mat = FALSE)))) {
  netwerk_poly_opp <- terra::as.polygons(cl_netwerk_opp, dissolve = TRUE) %>% terra::makeValid()
  watervlakken_in_netwerk_opp <- terra::crop(watervlakken_clean, netwerk_poly_opp)
  
  if (!is.null(watervlakken_in_netwerk_opp) && nrow(watervlakken_in_netwerk_opp) > 0) {
    rel_opp <- matrix(terra::is.related(netwerk_poly_opp, watervlakken_in_netwerk_opp, "intersects"), 
                      nrow = nrow(netwerk_poly_opp), ncol = nrow(watervlakken_in_netwerk_opp))
    netwerk_poly_opp$Aantal_Plassen <- rowSums(rel_opp)
    id_col_opp <- names(netwerk_poly_opp)[1]
    valid_ids_opp <- as.numeric(netwerk_poly_opp[[id_col_opp]][netwerk_poly_opp$Aantal_Plassen >= min_plassen, 1])
  } else {
    valid_ids_opp <- c()
  }
  
  if (length(na.omit(valid_ids_opp)) > 0) {
    heikikker_leefgebied_opp <- terra::mask(heikikker_leefgebied1_opp, cl_netwerk_opp %in% valid_ids_opp)
    cl_opp                   <- terra::mask(cl_netwerk_opp, cl_netwerk_opp %in% valid_ids_opp)
  } else {
    heikikker_leefgebied_opp <- template_MH * NA
    cl_opp                   <- template_MH * NA
  }
} else {
  heikikker_leefgebied_opp <- template_MH * NA
  cl_opp                   <- template_MH * NA
}

final_max <- heikikker_leefgebied_max
final_opp <- heikikker_leefgebied_opp


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

