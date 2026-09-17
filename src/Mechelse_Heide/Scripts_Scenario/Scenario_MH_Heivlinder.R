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

soort <- "heivlinder"

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

df <- read_excel(here::here("data/input/Excel_files/Soorten_bwk_afstanden.xlsx"))
resultaat <- df %>% filter(tolower(trimws(Soort)) == soort) %>% select(Type, MinOpp_ha, AfstandBiotopen_m, Dispersiecap_m)

oppervlakte_ha <- resultaat$MinOpp_ha[1]
afstand_m      <- resultaat$AfstandBiotopen_m[1]
buffer_m       <- resultaat$Dispersiecap_m[1]

print(resultaat)
rm(df, resultaat)

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
    bevat_codes_escaped <- gsub("([\\.\\^\\$\\*\\+\\?\\(\\)\\[\\{\\\\\\|])", "\\\\\\1", bevat_codes)
    bevat_codes_anchored <- paste0("^", bevat_codes_escaped)
    regex_term <- paste0(bevat_codes_anchored, collapse = "|")
    
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

heide_max          <- lijst_matches[["bwk_heide"]]
heide_opp          <- lijst_oppervlaktes[["bwk_heide"]]
duin_max           <- lijst_matches[["bwk_duin"]]
duin_opp           <- lijst_oppervlaktes[["bwk_duin"]]
grasland_max       <- lijst_matches[["bwk_grasland"]]
grasland_opp       <- lijst_oppervlaktes[["bwk_grasland"]]
groeve_terril_max  <- lijst_matches[["bwk_groeve_terril"]]
groeve_terril_opp  <- lijst_oppervlaktes[["bwk_groeve_terril"]]
zand_max           <- lijst_matches[["zand"]]
zand_opp           <- lijst_oppervlaktes[["zand"]]
beschutting_max    <- lijst_matches[["bwk_beschutting"]]
beschutting_opp    <- lijst_oppervlaktes[["bwk_beschutting"]]

rm(tabel_vlaanderen, vertaal_df, lijst_matches, lijst_oppervlaktes)
gc()

heivlinder_biotoop_max <- heide_max | duin_max | grasland_max | groeve_terril_max
heivlinder_biotoop_opp <- terra::app(
  terra::sds(heide_opp, duin_opp, grasland_opp, groeve_terril_opp), 
  fun = function(x) pmin(sum(x, na.rm = TRUE), 1)
)

# Spoor A: MAX
heivlinder_clusters_max <- cluster_filter_compleet(
  masker    = heivlinder_biotoop_max,
  opp_laag  = heivlinder_biotoop_opp,
  drempel_m2 = 100000,   # 10 ha
  dist_m    = 50,
  werkelijk = FALSE
)

# Spoor B: OPP (Parallel & Autonoom)
r_binair_biotoop_opp <- terra::ifel(!is.na(heivlinder_biotoop_opp) & heivlinder_biotoop_opp > 0, 1, NA)
heivlinder_clusters_opp <- cluster_filter_compleet(
  masker    = r_binair_biotoop_opp,
  opp_laag  = heivlinder_biotoop_opp,
  drempel_m2 = 100000,   # 10 ha
  dist_m    = 50,
  werkelijk = TRUE
)

rm(r_binair_biotoop_opp)
gc()

# Spoor A: MAX
zand_clusters_max <- cluster_filter_compleet(
  masker    = zand_max,
  opp_laag  = zand_opp,
  drempel_m2 = 5000,     # 0.5 ha
  dist_m    = 100,
  werkelijk = FALSE
)

# Spoor B: OPP (Parallel & Autonoom)
r_binair_zand_opp <- terra::ifel(!is.na(zand_opp) & zand_opp > 0, 1, NA)
zand_clusters_opp <- cluster_filter_compleet(
  masker    = r_binair_zand_opp,
  opp_laag  = zand_opp,
  drempel_m2 = 5000,     # 0.5 ha
  dist_m    = 100,
  werkelijk = TRUE
)

rm(r_binair_zand_opp)
gc()

# ==============================================================================
# 1. VERWERK GROENKAART (MET BEHOUD VAN DE MAXIMALE FRACTIE NAAR 10M)
# ==============================================================================
hooggroen_raw <- rast(here("data/input/ASCI Files/Groenkaart_2021.tif"))

message("-> Groenkaart synchroniseren met 10m MasterGrid (Max-waarde behouden)...")

# We projecteren en resamplen direct naar het 10m studiegebied (id_raster_MH).
# Door 'method = "max"' te gebruiken, zal terra bij het aggregeren/resamplen 
# van de subpixels altijd de hoogste waarde (bijv. 0.6) selecteren voor de 10m cel.
groenkaart_10m_synchronized <- project(hooggroen_raw, id_raster_MH, method = "max")


# ==============================================================================
# 2. COMBINEER BWK-BESCHUTTING EN GROENKAART (CELL-BY-CELL MAXIMAAL)
# ==============================================================================
message("-> Beschuttingslagen combineren voor de heivlinder (Max van BWK vs Groenkaart)...")

# We zetten eventuele NA's tijdelijk om naar 0 om een zuivere wiskundige vergelijking te kunnen maken
groen_0 <- terra::ifel(is.na(groenkaart_10m_synchronized), 0, groenkaart_10m_synchronized)
bwk_max_0 <- terra::ifel(is.na(beschutting_max), 0, beschutting_max)
bwk_opp_0 <- terra::ifel(is.na(beschutting_opp), 0, beschutting_opp)

# MAX-spoor: Binaire aanwezigheid (als een van beide > 0 is, of specifiek als een van beide 1 is)
# Omdat je Groenkaart nu fracties bevat, is het logischer om te kijken of de waarde groter is dan een drempel (bijv. 0)
heivlinder_beschutting_max <- terra::ifel(bwk_max_0 == 1 | groen_0 > 0, 1, NA)

# OPP-spoor: Neem cel-per-cel de hoogste waarde (als Groenkaart 0.6 is en BWK 0.4, wordt het 0.6)
heivlinder_beschutting_opp <- terra::app(
  terra::sds(bwk_opp_0, groen_0),
  fun = function(x) {
    res <- max(x, na.rm = TRUE)
    return(pmin(res, 1)) # Veiligheidsmarge: nooit boven de 100% (1) gaan
  }
)

# Zet cellen met waarde 0 weer netjes terug naar NA voor de consistentie in je script
heivlinder_beschutting_opp <- terra::ifel(heivlinder_beschutting_opp == 0, NA, heivlinder_beschutting_opp)

# Grote schoonmaak van zware objecten
rm(hooggroen_raw, groenkaart_10m_synchronized, groen_0, bwk_max_0, bwk_opp_0)
gc()

# Spoor A: MAX
beschutting_clusters_max <- cluster_filter_compleet(
  masker    = heivlinder_beschutting_max,
  opp_laag  = heivlinder_beschutting_opp,
  drempel_m2 = 400,      # 0.04 ha
  dist_m    = 100,
  werkelijk = FALSE
)

# Spoor B: OPP (Parallel & Autonoom)
r_binair_beschutting_opp <- terra::ifel(!is.na(heivlinder_beschutting_opp) & heivlinder_beschutting_opp > 0, 1, NA)
beschutting_clusters_opp <- cluster_filter_compleet(
  masker    = r_binair_beschutting_opp,
  opp_laag  = heivlinder_beschutting_opp,
  drempel_m2 = 400,      # 0.04 ha
  dist_m    = 100,
  werkelijk = TRUE
)

rm(r_binair_beschutting_opp)
gc()

# ==============================================================================
# 1. BASISHABITAT OPBOUWEN (HEIDE + CORRESPONDERENDE 50M BOSRAND)
# ==============================================================================
message("-> Stap 1: Basis open biotoop en aangrenzende 50m bosranden combineren...")

# --- SPOOR A: MAX ---
biotoop_max <- heivlinder_clusters_max$raster
beschutting_ruw_max <- terra::ifel(!is.na(beschutting_clusters_max$raster) & beschutting_clusters_max$raster > 0, 1, NA)

dist_tot_biotoop_max <- terra::distance(biotoop_max)
bosrand_masker_max   <- terra::ifel(beschutting_ruw_max == 1 & dist_tot_biotoop_max > 0 & dist_tot_biotoop_max <= 50, 1, NA)
bosrand_direct_max   <- bosrand_masker_max

foerageer_base_max   <- terra::ifel(!is.na(biotoop_max) | !is.na(bosrand_direct_max), 1, NA)
voortplanting_base_max <- zand_clusters_max$raster

# --- SPOOR B: OPP (Parallel & Autonoom) ---
biotoop_opp <- heivlinder_clusters_opp$raster
beschutting_ruw_opp <- beschutting_clusters_opp$raster

r_binair_biotoop_opp <- terra::ifel(!is.na(biotoop_opp) & biotoop_opp > 0, 1, NA)
dist_tot_biotoop_opp <- terra::distance(r_binair_biotoop_opp)

bosrand_masker_opp   <- terra::ifel(!is.na(beschutting_ruw_opp) & beschutting_ruw_opp > 0 & dist_tot_biotoop_opp > 0 & dist_tot_biotoop_opp <= 50, 1, NA)
bosrand_direct_opp   <- terra::mask(beschutting_ruw_opp, bosrand_masker_opp)

foerageer_base_opp <- terra::app(
  terra::sds(
    terra::ifel(is.na(biotoop_opp), 0, biotoop_opp),
    terra::ifel(is.na(bosrand_direct_opp), 0, bosrand_direct_opp)
  ),
  fun = function(x) pmin(sum(x, na.rm = TRUE), 1)
)
foerageer_base_opp <- terra::ifel(foerageer_base_opp == 0, NA, foerageer_base_opp)
voortplanting_base_opp <- zand_clusters_opp$raster

# ==============================================================================
# 2. AFSTANDSFILTER SPOOR A (MAX - 500m)
# ==============================================================================
message("-> Stap 2: Wederzijdse afstandsfilters berekenen voor Maximaal Potentie (max 500m)...")

if (!all(is.na(terra::values(voortplanting_base_max, mat=FALSE))) && 
    !all(is.na(terra::values(foerageer_base_max, mat=FALSE)))) {
  
  dist_naar_voortplanting_max <- terra::distance(voortplanting_base_max)
  zone_voortplanting_500m_max <- terra::ifel(dist_naar_voortplanting_max <= 500, 1, NA)
  
  dist_naar_foerageer_max <- terra::distance(foerageer_base_max)
  zone_foerageer_500m_max <- terra::ifel(dist_naar_foerageer_max <= 500, 1, NA)
  
  foerageer_finaal_max     <- terra::mask(foerageer_base_max, zone_voortplanting_500m_max)
  voortplanting_finaal_max <- terra::mask(voortplanting_base_max, zone_foerageer_500m_max)
  
  rm(dist_naar_voortplanting_max, zone_voortplanting_500m_max, dist_naar_foerageer_max, zone_foerageer_500m_max)
} else {
  foerageer_finaal_max     <- template_MH * NA
  voortplanting_finaal_max <- template_MH * NA
}

# ==============================================================================
# 3. AFSTANDSFILTER SPOOR B (OPP - 500m Parallel)
# ==============================================================================
message("-> Stap 3: Wederzijdse afstandsfilters berekenen voor Werkelijke Oppervlakte (max 500m)...")

if (!all(is.na(terra::values(voortplanting_base_opp, mat=FALSE))) && 
    !all(is.na(terra::values(foerageer_base_opp, mat=FALSE)))) {
  
  r_binair_vpt_base_opp  <- terra::ifel(!is.na(voortplanting_base_opp) & voortplanting_base_opp > 0, 1, NA)
  r_binair_foer_base_opp <- terra::ifel(!is.na(foerageer_base_opp) & foerageer_base_opp > 0, 1, NA)
  
  dist_naar_voortplanting_opp <- terra::distance(r_binair_vpt_base_opp)
  zone_voortplanting_500m_opp <- terra::ifel(dist_naar_voortplanting_opp <= 500, 1, NA)
  
  dist_naar_foerageer_opp <- terra::distance(r_binair_foer_base_opp)
  zone_foerageer_500m_opp <- terra::ifel(dist_naar_foerageer_opp <= 500, 1, NA)
  
  foerageer_finaal_opp     <- terra::mask(foerageer_base_opp, zone_voortplanting_500m_opp)
  voortplanting_finaal_opp <- terra::mask(voortplanting_base_opp, zone_foerageer_500m_opp)
  
  rm(dist_naar_voortplanting_opp, zone_voortplanting_500m_opp, dist_naar_foerageer_opp, zone_foerageer_500m_opp, r_binair_vpt_base_opp, r_binair_foer_base_opp)
} else {
  foerageer_finaal_opp     <- template_MH * NA
  voortplanting_finaal_opp <- template_MH * NA
}

# ==============================================================================
# 4. FINALE COMBINATIE & SYNCHRONISATIE DEELLAGEN
# ==============================================================================
message("-> Stap 4: Eindresultaat samenstellen...")

heivlinder_leefgebied_metzand_max <- terra::ifel(!is.na(foerageer_finaal_max) | !is.na(voortplanting_finaal_max), 1, NA)

heivlinder_leefgebied_metzand_opp <- terra::app(
  terra::sds(
    terra::ifel(is.na(foerageer_finaal_opp), 0, foerageer_finaal_opp),
    terra::ifel(is.na(voortplanting_finaal_opp), 0, voortplanting_finaal_opp)
  ),
  fun = function(x) pmin(sum(x, na.rm = TRUE), 1)
)
heivlinder_leefgebied_metzand_opp <- terra::ifel(heivlinder_leefgebied_metzand_opp == 0, NA, heivlinder_leefgebied_metzand_opp)

deellaag_heide_opp        <- terra::mask(biotoop_opp, foerageer_finaal_opp)
deellaag_bosrand_opp      <- terra::mask(bosrand_direct_opp, foerageer_finaal_opp)
deellaag_zand_opp         <- voortplanting_finaal_opp
deellaag_beschutting_500m <- bosrand_direct_opp

rm(dist_tot_biotoop_max, dist_tot_biotoop_opp, bosrand_masker_max, bosrand_masker_opp, bosrand_direct_max, bosrand_direct_opp, 
   foerageer_base_max, foerageer_base_opp, voortplanting_base_max, voortplanting_base_opp,
   foerageer_finaal_max, voortplanting_finaal_max, foerageer_finaal_opp, voortplanting_finaal_opp)
gc()


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

