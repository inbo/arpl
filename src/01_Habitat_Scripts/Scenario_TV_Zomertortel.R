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
soort <- "zomertortel"

# Scenario pad en naam bepalen

# --- DYNAMISCHE SCENARIO PARAMETER CHECK ---
if (!exists("params") || is.null(params$scenario_rds_path)) {
  scenario_rds_path <- "data/input/Scenario_rds/TV_Scenario_BWK_2025.rds"
} else {
  scenario_rds_path <- params$scenario_rds_path
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

# Variabelen definiëren
oppervlakte_ha <- resultaat$MinOpp_ha[1]
afstand_m      <- resultaat$AfstandBiotopen_m[1]
buffer_m       <- resultaat$Dispersiecap_m[1]

foerageer_afstand_m  <- 300
drinkwater_afstand_m <- 500

print(resultaat)
rm(df, resultaat)

# Gebied voorbereiden (Turnhouts Vennegebied met waterdichte 10km buffer)
area_shape  <- vect(here("data/input/Turnhouts_Vennegebied.shp"))
master_grid <- rast(here("data/input/Raster_Vlaanderen/Vlaanderen_MasterGrid_10m.tif"))[[1]]

df_namen_sleutel <- read_csv(here("data/input/Excel_files/BWK_Laag_Namen_2025.csv"), show_col_types = FALSE)
gouden_namenlijst <- tolower(trimws(df_namen_sleutel$Laagnaam))

area_shape_proj <- project(area_shape, crs(master_grid))

# Dwing een ruime buffer van 10.000m (10 km) af voor de Bounding Box
buffer_raster_m <- max(buffer_m, 10000) 
area_buffer_fix <- buffer(area_shape_proj, width = buffer_raster_m)

message("-> Vertaalraster voor globale/lokale cellen opbouwen INCLUSIEF 10KM BUFFER...")

# FIX 1: Gebruik snap = "out" zodat GEEN ENKELE randpixel buiten de bounding box valt!
id_raster_TV <- crop(master_grid, area_buffer_fix, snap = "out")

# FIX 2: Geef id_raster_TV een unieke opeenvolgende index voor lokale cellen (1..N)
values(id_raster_TV) <- 1:terra::ncell(id_raster_TV)

# FIX 3: Haal exact alle globale cel-IDs op voor de volledige omvang van id_raster_TV
globale_ids_vector <- terra::cells(master_grid, terra::ext(id_raster_TV))

# FIX 4: Bouw de vertaaltabel waterdicht op
vertaal_df <- data.table(
  lokale_id  = 1:terra::ncell(id_raster_TV),
  globale_id = globale_ids_vector
)

vertaal_df <- vertaal_df[!is.na(globale_id)]
studiegebied_globale_ids <- unique(vertaal_df$globale_id)

message("Aantal actieve globale cellen in vertaal_df: ", length(studiegebied_globale_ids))

# Template opbouwen voor geografische exports
template_TV <- terra::rasterize(area_buffer_fix, id_raster_TV, field = 1, background = NA)

grens_web <- sf::st_as_sf(terra::project(area_shape, "EPSG:4326"))

rm(globale_ids_vector)
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

voortplanting_max <- lijst_matches[["voortplanting"]]
voortplanting_opp <- lijst_oppervlaktes[["voortplanting"]]
water_max         <- lijst_matches[["water"]]
water_opp         <- lijst_oppervlaktes[["water"]]
foerageer_max     <- lijst_matches[["foerageer"]]
foerageer_opp     <- lijst_oppervlaktes[["foerageer"]]

raster_simpel_final <- id_raster_TV

rm(lijst_matches, lijst_oppervlaktes)
gc()

# ==============================================================================
# STAP 2: INTEGRATIE EXTERNE KAARTEN
# ==============================================================================
message("-> Aanvullende externe kaarten inlezen...")
r_doel_na <- id_raster_TV * NA

# 1. LANDBOUWPERCELEN (LGP)
r_lgp_vlaanderen <- rast(here("data/input/Raster_Vlaanderen/vlaanderen_lgp_gwscod_2025_10m.tif"))
r_lgp_local <- r_lgp_vlaanderen %>% 
  terra::crop(area_buffer_fix) %>% 
  terra::resample(r_doel_na, method = "near")

zomer_gewascodes <- c(201, 202, 311, 321, 601, 800) 
r_lgp_bin <- terra::ifel(!is.na(r_lgp_local) & r_lgp_local %in% zomer_gewascodes, 1, NA)
r_lgp_opp <- terra::ifel(!is.na(r_lgp_local) & r_lgp_local %in% zomer_gewascodes, 0.01, NA)

foerageer_totaal_max <- terra::cover(foerageer_max, r_lgp_bin)
foerageer_totaal_opp <- terra::cover(foerageer_opp, r_lgp_opp)

# 2. WATERVLAKKEN EN HUET-ZONES
r_watervlakken <- rast(here("data/input/Raster_Vlaanderen/vlaanderen_watervlakken_2024_10m.tif")) %>%
  terra::crop(area_buffer_fix) %>% 
  terra::resample(r_doel_na, method = "near")

r_watervlakken_bin <- terra::ifel(!is.na(r_watervlakken) & r_watervlakken > 0, 1, NA)
r_watervlakken_opp <- terra::ifel(!is.na(r_watervlakken) & r_watervlakken > 0, 0.01, NA)

r_huetzon <- rast(here("data/input/Raster_Vlaanderen/vlaanderen_huetzon_10m.tif")) %>%
  terra::crop(area_buffer_fix) %>% 
  terra::resample(r_doel_na, method = "near")

r_huet_bin <- terra::ifel(!is.na(r_huetzon) & r_huetzon > 0, 1, NA)
r_huet_opp <- terra::ifel(!is.na(r_huetzon) & r_huetzon > 0, 0.01, NA)

water_ext_bin <- terra::cover(r_watervlakken_bin, r_huet_bin)
water_ext_opp <- terra::cover(r_watervlakken_opp, r_huet_opp)

water_totaal_max <- terra::cover(water_max, water_ext_bin)
water_totaal_opp <- terra::cover(water_opp, water_ext_opp)

# ==============================================================================
# OPBOUWEN VAN DE NORMALE BOS- EN BEBOUWINGSFILTER (OPTIES 2 & 4)
# ==============================================================================
message("-> Opbouwen van de Groenkaart (Bosfilter) en Dorpskern-bebouwingsfilter...")

# A. GROENKAART / BOSFILTER (OPTIE 2)
hooggroen_raw         <- rast(here("data/input/ASCI Files/Groenkaart_2021.tif"))
area_buffer_groen_crs <- project(area_buffer_fix, crs(hooggroen_raw))
hooggroen_local_raw   <- crop(hooggroen_raw, area_buffer_groen_crs)

hooggroen_binair  <- terra::classify(hooggroen_local_raw, matrix(c(0.5, 1.5, 1), ncol = 3, byrow = TRUE), others = 0)
hooggroen_10m_raw <- terra::aggregate(hooggroen_binair, fact = 10, fun = "max")
hooggroen_sync    <- project(hooggroen_10m_raw, crs(raster_simpel_final)) %>% terra::resample(raster_simpel_final, method = "near")

r_groen_bin <- terra::ifel(!is.na(hooggroen_sync) & hooggroen_sync > 0, 1, NA)

w_matrix <- matrix(1, nrow = 5, ncol = 5)
temp_groen_file <- tempfile(fileext = ".tif")

groen_som <- terra::focal(r_groen_bin, w = w_matrix, fun = "sum", na.rm = TRUE, filename = temp_groen_file, overwrite = TRUE)

# Dichte boskernen (> 20 are)
r_bos_20are <- terra::ifel(groen_som >= 12 & r_groen_bin == 1, 1, NA)

rm(hooggroen_raw, hooggroen_local_raw, hooggroen_binair, hooggroen_10m_raw, r_groen_bin, groen_som)

# B. GROTE DORPSKERNEN BEBOUWINGSFILTER (OPTIE 4)
bebouwings_lagen <- gouden_namenlijst[grepl("^u[a-z0-9]", gouden_namenlijst) | gouden_namenlijst == "u"]
bebouwings_lagen <- bebouwings_lagen[bebouwings_lagen != "ulm"]

tabel_bebouw_TV  <- tabel_vlaanderen[CODE %in% bebouwings_lagen & cel_id %in% studiegebied_globale_ids]
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

# Dorpskernen >= 10 ha met 100m verstoringszone
res_bebouw_groot <- cluster_filter_compleet(
  masker     = r_bebouwing_schoon,
  opp_laag   = r_bebouwing_schoon,
  drempel_m2 = 100000, 
  dist_m     = 100,
  werkelijk  = FALSE
)

r_bebouw_groot <- res_bebouw_groot$raster

if (!is.null(r_bebouw_groot) && !terra::global(is.na(r_bebouw_groot), "sum")[[1]] == terra::ncell(r_bebouw_groot)) {
  masker_bebouw_100m_groot <- terra::buffer(r_bebouw_groot, width = 100)
  masker_bebouw_100m_groot <- terra::ifel(!is.na(masker_bebouw_100m_groot) & masker_bebouw_100m_groot > 0, 1, NA)
} else {
  masker_bebouw_100m_groot <- id_raster_TV * NA
}

rm(tabel_bebouw_TV, tabel_bebouw_som, tabel_bebouw_mapping, r_bebouwing_raw, res_bebouw_groot, r_bebouw_groot, tabel_vlaanderen)
gc()

# ==============================================================================
# STAP 3: RUIMTELIJKE AFSTANDSFILTERS (PARALLEL SPOOR A EN SPOOR B)
# ==============================================================================
message("-> Toepassen van expertfilters voor Zomertortel op Spoor A en Spoor B...")

# 1. OPTIE 2 & 4: BROEDHABITAT SCHOONMAKEN (Boskernen & Dorpskernen uitsluiten)
voortplanting_schoon_max <- terra::mask(voortplanting_max, r_bos_20are, inverse = TRUE)
voortplanting_schoon_opp <- terra::mask(voortplanting_opp, r_bos_20are, inverse = TRUE)

if (exists("masker_bebouw_100m_groot") && !all(is.na(terra::values(masker_bebouw_100m_groot, mat=FALSE)))) {
  voortplanting_schoon_max <- terra::mask(voortplanting_schoon_max, masker_bebouw_100m_groot, inverse = TRUE)
  voortplanting_schoon_opp <- terra::mask(voortplanting_schoon_opp, masker_bebouw_100m_groot, inverse = TRUE)
}

# 2. OPTIE 1: CLUSTERING VOORTPLANTING (300 m² = meidoornhagen & oeverstruweel)
res_broed_max <- cluster_filter_compleet(
  masker     = voortplanting_schoon_max,
  opp_laag   = voortplanting_schoon_max,
  drempel_m2 = 300,  # 300 m² (0.03 ha) drempel
  dist_m     = 50,   # 50m overbrugging
  werkelijk  = FALSE
)
broed_gefilterd_max <- res_broed_max$raster

r_binair_broed_opp <- terra::ifel(!is.na(voortplanting_schoon_opp) & voortplanting_schoon_opp > 0, 1, NA)
res_broed_opp <- cluster_filter_compleet(
  masker     = r_binair_broed_opp,
  opp_laag   = voortplanting_schoon_opp,
  drempel_m2 = 300,  
  dist_m     = 50,    
  werkelijk  = TRUE
)
broed_gefilterd_opp <- res_broed_opp$raster

# 3. CLUSTERING FOERAGEERGEBIED (Minimaal 1 ha / 10.000 m²)
res_foerageer_max <- cluster_filter_compleet(
  masker     = foerageer_totaal_max,
  opp_laag   = foerageer_totaal_max,
  drempel_m2 = 10000, # 1 ha drempel
  dist_m     = 50,     
  werkelijk  = FALSE
)
foerageer_gefilterd_max <- res_foerageer_max$raster

r_binair_foerageer_opp <- terra::ifel(!is.na(foerageer_totaal_opp) & foerageer_totaal_opp > 0, 1, NA)
res_foerageer_opp <- cluster_filter_compleet(
  masker     = r_binair_foerageer_opp,
  opp_laag   = foerageer_totaal_opp,
  drempel_m2 = 10000, 
  dist_m     = 50,     
  werkelijk  = TRUE
)
foerageer_gefilterd_opp <- res_foerageer_opp$raster

water_gefilterd_max <- water_totaal_max
water_gefilterd_opp <- water_totaal_opp

# ==============================================================================
# STAP 4: OMGEVINGSBUFFERS EN KOPPELINGEN (PARALLEL BEREKEND)
# ==============================================================================
message("-> Omgevingsbuffers berekenen: 300m naar voedsel, 500m naar drinkwater...")

alleen_actieve_pixels <- function(r) {
  if (is.null(r) || all(is.na(terra::values(r, mat=FALSE)))) return(r * NA)
  terra::ifel(!is.na(r) & r > 0, 1, NA)
}

# --- Spoor A: Buffers voor Maximale Potentie ---
foer_max_clean <- alleen_actieve_pixels(foerageer_gefilterd_max)
wat_max_clean  <- alleen_actieve_pixels(water_gefilterd_max)

r_foerageer_buf_max  <- terra::buffer(foer_max_clean, width = foerageer_afstand_m) # 300m
r_foerageer_mask_max <- alleen_actieve_pixels(r_foerageer_buf_max)

r_water_buf_max      <- terra::buffer(wat_max_clean, width = drinkwater_afstand_m)  # 500m
r_water_mask_max     <- alleen_actieve_pixels(r_water_buf_max)

# Omgeving is geschikt als er én voedsel (300m) én drinkwater (500m) aanwezig is
r_omgeving_ok_max    <- terra::mask(r_foerageer_mask_max, r_water_mask_max)
r_omgeving_ok_max    <- alleen_actieve_pixels(r_omgeving_ok_max)


# --- Spoor B: Buffers voor Werkelijke Oppervlakte ---
foer_opp_clean <- alleen_actieve_pixels(foerageer_gefilterd_opp)
wat_opp_clean  <- alleen_actieve_pixels(water_gefilterd_opp)

r_foerageer_buf_opp  <- terra::buffer(foer_opp_clean, width = foerageer_afstand_m)
r_foerageer_mask_opp <- alleen_actieve_pixels(r_foerageer_buf_opp)

r_water_buf_opp      <- terra::buffer(wat_opp_clean, width = drinkwater_afstand_m)
r_water_mask_opp     <- alleen_actieve_pixels(r_water_buf_opp)

r_omgeving_ok_opp    <- terra::mask(r_foerageer_mask_opp, r_water_mask_opp)
r_omgeving_ok_opp    <- alleen_actieve_pixels(r_omgeving_ok_opp)

# ==============================================================================
# STAP 5: FINALE LEEFGEBIEDEN OPBOUWEN (SPOOR A EN B)
# ==============================================================================
message("-> Finaal leefgebied berekenen per spoor...")

# Spoor A: Binaire leefgebiedenkaart
zomertortel_leefgebied_max <- terra::mask(broed_gefilterd_max, r_omgeving_ok_max)
names(zomertortel_leefgebied_max) <- "Leefgebied_Max"

# Spoor B: Echte bedekkingsfracties voor exacte hectare-rapportage
zomertortel_leefgebied_opp <- terra::mask(broed_gefilterd_opp, r_omgeving_ok_opp)
names(zomertortel_leefgebied_opp) <- "Leefgebied_Opp"

alle_matches      <- zomertortel_leefgebied_max
alle_oppervlaktes <- zomertortel_leefgebied_opp


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
potentie_export_rast <- if (exists("grutto_kaart_A") && !is.null(grutto_kaart_A) && !all(is.na(terra::values(grutto_kaart_A, mat=FALSE)))) {
  terra::ifel(grutto_kaart_A > 0, 1, NA)
} else if (exists("final_max") && !all(is.na(terra::values(final_max, mat=FALSE)))) {
  terra::ifel(!is.na(final_max) & final_max > 0, 1, NA)
} else {
  terra::rast(template_TV, vals = NA)
}

# 2. Bepaal Werkelijke Oppervlakte Raster
werkelijk_export_rast <- if (exists("resB_strikt") && !is.null(resB_strikt) && (!all(is.na(terra::values(resB_strikt$kern, mat=FALSE))) || !all(is.na(terra::values(resB_strikt$bouw, mat=FALSE))))) {
  r_net_totaal_opp <- terra::cover(resB_strikt$kern, resB_strikt$bouw)
  terra::ifel(!is.na(r_net_totaal_opp) & r_net_totaal_opp > 0, 1, NA)
} else if (exists("final_opp") && !all(is.na(terra::values(final_opp, mat=FALSE)))) {
  terra::ifel(!is.na(final_opp) & final_opp > 0, 1, NA)
} else {
  terra::rast(template_TV, vals = NA)
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

