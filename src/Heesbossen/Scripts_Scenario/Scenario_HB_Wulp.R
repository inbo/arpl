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
conflicted::conflicts_prefer(dplyr::first)
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

df <- read_excel(here::here("data/input/Excel_files/Soorten_bwk_afstanden.xlsx"))
soort <- "wulp"

# Scenario pad en naam bepalen

# --- DYNAMISCHE SCENARIO PARAMETER CHECK ---
if (!exists("params") || is.null(params$scenario_rds_path)) {
  scenario_rds_path <- "data/input/Scenario_rds/HB_Scenario_BWK_2025.rds"
} else {
  scenario_rds_path <- params$scenario_rds_path
}

p_raw <- gsub("^([.][.]/)+", "", scenario_rds_path)
scenario_path <- here::here(p_raw)


if (!file.exists(scenario_path)) {
  stop(paste("❌ FOUT: Scenario RDS bestand NIET gevonden op:", scenario_path))
}

scen_volledig <- basename(scenario_path)
scenario_naam <- gsub("^HB_Scenario_|^Scenario_|.rds$", "", scen_volledig)

message(paste("Verwerken van soort:", soort, "binnen scenario:", scenario_naam))

resultaat <- df %>%
  filter(tolower(trimws(Soort)) == soort) %>%
  select(Type, MinOpp_ha, AfstandBiotopen_m, Dispersiecap_m)

# Variabelen definiëren
oppervlakte_ha <- resultaat$MinOpp_ha[1]
afstand_m      <- resultaat$AfstandBiotopen_m[1]
buffer_m       <- resultaat$Dispersiecap_m[1]

print(resultaat)
rm(df, resultaat)

area_shape  <- vect(here("data/input/Heesbossen.shp"))
master_grid <- rast(here("data/input/Raster_Vlaanderen/Vlaanderen_MasterGrid_10m.tif"))[[1]]

df_namen_sleutel <- read_csv(here("data/input/Excel_files/BWK_Laag_Namen_2025.csv"), show_col_types = FALSE)
gouden_namenlijst <- tolower(trimws(df_namen_sleutel$Laagnaam))

area_shape_proj <- project(area_shape, crs(master_grid))
area_buffer_fix <- buffer(area_shape_proj, width = buffer_m)

message("-> Vertaalraster voor globale/lokale cellen opbouwen via snelle MASK methode...")
id_raster_HB <- crop(master_grid, area_buffer_fix, snap = "near")

globale_id_raster <- master_grid
values(globale_id_raster) <- 1:ncell(globale_id_raster)

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

cat("Gecorrigeerd aantal pixels in template_HB: ", sum(terra::values(template_HB) == 1, na.rm=TRUE), "\n")

grens_web <- sf::st_as_sf(terra::project(area_shape, "EPSG:4326"))

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
    bevat_codes_escaped  <- gsub("([\\.\\^\\$\\*\\+\\?\\(\\)\\[\\{\\\\\\|])", "\\\\\\1", bevat_codes)
    bevat_codes_anchored <- paste0("^", bevat_codes_escaped)
    regex_term           <- paste0(bevat_codes_anchored, collapse = "|")
    
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

biotoop_grasland_nat_max <- lijst_matches[["biotoop_grasland_nat"]]
biotoop_grasland_nat_opp <- lijst_oppervlaktes[["biotoop_grasland_nat"]]
biotoop_overig_max       <- lijst_matches[["biotoop_overig"]]
biotoop_overig_opp       <- lijst_oppervlaktes[["biotoop_overig"]]

tabel_akkers <- tabel_vlaanderen[grepl("^b", CODE) & cel_id %in% studiegebied_globale_ids]
tabel_akkers_som <- tabel_akkers[, .(Oppervlakte = sum(BWK_FRAC, na.rm = TRUE)), by = .(cel_id)]

r_akkers_max <- id_raster_HB * NA
if(nrow(tabel_akkers_som) > 0) {
  setnames(tabel_akkers_som, "cel_id", "globale_id")
  mapping_akkers <- merge(tabel_akkers_som, vertaal_df, by = "globale_id", all.x = TRUE)[!is.na(lokale_id)]
  r_akkers_max[mapping_akkers$lokale_id] <- ifelse(mapping_akkers$Oppervlakte >= 0.01, 1, 0)
}
rm(tabel_akkers, tabel_akkers_som, mapping_akkers, tabel_vlaanderen)

rm(vertaal_df, lijst_matches, lijst_oppervlaktes)
gc()

# ==============================================================================
# STAP 1: OMGEVINGSLAAG SYNCHRONISATIE & VERBOSSINGSFILTER (>20 ARE HOOGGROEN)
# ==============================================================================
message("-> Externe omgevingslagen en groenkaart inlezen...")

r_ovstr_raw      <- rast(here("data/input/Raster_Vlaanderen/vlaanderen_ovstrg_10m.tif"))
r_drainage_raw   <- rast(here("data/input/Raster_Vlaanderen/vlaanderen_drainage_10m.tif"))
r_profiel_raw    <- rast(here("data/input/Raster_Vlaanderen/Vlaanderen_profiel_10m.tif"))
r_groenkaart_raw <- rast(here("data/input/ASCI Files/Groenkaart_2021.tif"))

# Grid-synchronisatie (C++ resample naar id_raster_HB grid)
r_ovstr      <- terra::resample(r_ovstr_raw, id_raster_HB, method = "near") %>% terra::crop(id_raster_HB)
r_drainage   <- terra::resample(r_drainage_raw, id_raster_HB, method = "near") %>% terra::crop(id_raster_HB)
r_profiel    <- terra::resample(r_profiel_raw, id_raster_HB, method = "near") %>% terra::crop(id_raster_HB)
r_groenkaart <- terra::resample(r_groenkaart_raw, id_raster_HB, method = "near") %>% terra::crop(id_raster_HB)

message("-> Verbossingsfilter berekenen: hooggroencomplexen groter dan 20 are uitsluiten...")
# Oude model: Hooggroen (waarde 1) wegknippen indien groter dan 20 are (2000 m²)
# We voeren een interne binaire clustering uit op de bomenlaag met dist_m = 0
r_hooggroen_bin <- r_groenkaart == 1
r_hooggroen_bin <- terra::ifel(r_hooggroen_bin == 1, 1, NA)

verbossing_clusters <- cluster_filter_compleet(
  masker       = r_hooggroen_bin,
  opp_laag     = r_hooggroen_bin, # Dummy opp laag
  drempel_m2   = 2000,   # 20 are = 2000 m²
  dist_m       = 0,       
  werkelijk    = FALSE      
)

# Dit binaire raster markeert de storende bosvlakken
r_verbossing_masker <- !is.na(verbossing_clusters$raster)

# ==============================================================================
# STAP 2: VALLEIFILTER & AKKER-EXCLUSIE (VOLLEDIG PARALLEL)
# ==============================================================================
message("-> Valleicriteria configureren op basis van biologische klassen...")

gewenste_drainage_letters <- c("d", "e", "f", "g", "h", "i", "e-f", "h-i", "e-i")
gewenste_profiel_types    <- c("p", "p+x")
gewenste_overstroming     <- c(1, 2)

drain_cats <- terra::cats(r_drainage)[[1]]
ids_drainage_nat <- drain_cats$value[drain_cats$Label %in% gewenste_drainage_letters]

prof_cats <- terra::cats(r_profiel)[[1]]
ids_profiel_val  <- prof_cats$value[prof_cats$Label %in% gewenste_profiel_types]

drainage_nat_ok   <- r_drainage %in% ids_drainage_nat
profiel_vallei_ok <- r_profiel %in% ids_profiel_val
overstroming_ok   <- r_ovstr %in% gewenste_overstroming

# Binaire vallei-indicator
vallei_indicator <- overstroming_ok | drainage_nat_ok | profiel_vallei_ok
vallei_masker    <- terra::ifel(vallei_indicator, 1, NA)

message("-> Graslanden opsplitsen conform GDX-sets...")
r_akkers_masker <- !is.na(r_akkers_max) & r_akkers_max == 1

# --- SPOOR A: MAX ---
r_grasland_gecorrigeerd_max <- terra::ifel(biotoop_grasland_nat_max == 1 & r_akkers_masker, NA, biotoop_grasland_nat_max)
r_grasland_vallei_max       <- terra::ifel(r_grasland_gecorrigeerd_max == 1 & !is.na(vallei_masker), 1, NA)

r_biotoop_geschikt_max <- !is.na(r_grasland_vallei_max) | (biotoop_overig_max == 1)
r_biotoop_finaal_max   <- terra::ifel(r_biotoop_geschikt_max == 1 & r_verbossing_masker == 0, 1, NA)

# --- SPOOR B: OPP (VEILIGE MASKERING) ---
r_grasland_gecorrigeerd_opp <- terra::ifel(r_akkers_masker, NA, biotoop_grasland_nat_opp)
r_grasland_vallei_opp       <- terra::ifel(!is.na(vallei_masker), r_grasland_gecorrigeerd_opp, NA)

g_clean <- terra::ifel(is.na(r_grasland_vallei_opp), 0, r_grasland_vallei_opp)
o_clean <- terra::ifel(is.na(biotoop_overig_opp), 0, biotoop_overig_opp)

som_raw <- g_clean + o_clean
som_cl  <- terra::clamp(som_raw, upper = 1.0)

r_biotoop_finaal_opp_raw <- terra::ifel(som_cl > 0, som_cl, NA)
r_biotoop_finaal_opp     <- terra::ifel(r_verbossing_masker == 1, NA, r_biotoop_finaal_opp_raw)

gc()

# ==============================================================================
# STAP 3: RUIMTELIJKE CLUSTERING (100M FUZZY, PARALLEL EN AUTONOOM)
# ==============================================================================
message("-> Landschappelijke clustering (100m fuzzy) op beide sporen parallel...")

# --- SPOOR A: MAX (Filtert op binaire footprint >= 100 ha) ---
wulp_clusters_max <- cluster_filter_compleet(
  masker     = r_biotoop_finaal_max,
  opp_laag   = r_biotoop_finaal_opp,
  drempel_m2 = 1000000, 
  dist_m     = 100,     
  werkelijk  = FALSE    
)
r_wulp_cluster_ids_max <- wulp_clusters_max$clusters

# --- SPOOR B: OPP (Filtert op werkelijk gesommeerde fracties >= 100 ha) ---
r_binair_biotoop_opp <- terra::ifel(!is.na(r_biotoop_finaal_opp) & r_biotoop_finaal_opp > 0, 1, NA)

wulp_clusters_opp <- cluster_filter_compleet(
  masker     = r_binair_biotoop_opp,
  opp_laag   = r_biotoop_finaal_opp,
  drempel_m2 = 1000000, 
  dist_m     = 100,     
  werkelijk  = TRUE     
)
r_wulp_cluster_ids_opp <- wulp_clusters_opp$clusters

gc()

# ==============================================================================
# STAP 4: FINALE SAMENSMELTING (PARALLEL)
# ==============================================================================
message("-> Finaal leefgebied vastleggen en synchroniseren met template...")

# Spoor A: Binaire footprint op basis van Spoor A clusters
wulp_leefgebied_max <- terra::ifel(!is.na(r_wulp_cluster_ids_max), 1, NA) %>% terra::crop(template_HB)

# Spoor B: Werkelijke habitat-fracties op basis van Spoor B clusters
wulp_leefgebied_opp <- terra::mask(r_biotoop_finaal_opp, !is.na(r_wulp_cluster_ids_opp)) %>% terra::crop(template_HB)

cat("Definitief leefgebied Wulp MAX (ha):", round(calc_ha_exact(wulp_leefgebied_max), 2), "\n")
cat("Definitief leefgebied Wulp OPP (ha):", round(calc_ha_exact(wulp_leefgebied_opp), 2), "\n")

gc()


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
potentie_export_rast <- if (exists("grutto_kaart_A") && !is.null(grutto_kaart_A) && !all(is.na(terra::values(grutto_kaart_A, mat=FALSE)))) {
  terra::ifel(grutto_kaart_A > 0, 1, NA)
} else if (exists("final_max") && !all(is.na(terra::values(final_max, mat=FALSE)))) {
  terra::ifel(!is.na(final_max) & final_max > 0, 1, NA)
} else {
  terra::rast(template_HB, vals = NA)
}

# 2. Bepaal Werkelijke Oppervlakte Raster
werkelijk_export_rast <- if (exists("resB_strikt") && !is.null(resB_strikt) && (!all(is.na(terra::values(resB_strikt$kern, mat=FALSE))) || !all(is.na(terra::values(resB_strikt$bouw, mat=FALSE))))) {
  r_net_totaal_opp <- terra::cover(resB_strikt$kern, resB_strikt$bouw)
  terra::ifel(!is.na(r_net_totaal_opp) & r_net_totaal_opp > 0, 1, NA)
} else if (exists("final_opp") && !all(is.na(terra::values(final_opp, mat=FALSE)))) {
  terra::ifel(!is.na(final_opp) & final_opp > 0, 1, NA)
} else {
  terra::rast(template_HB, vals = NA)
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

