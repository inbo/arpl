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

df <- read_excel(here::here("data/input/Excel_files/Soorten_bwk_afstanden.xlsx"))
soort <- "grutto"

# Scenario pad en naam bepalen

# --- DYNAMISCHE SCENARIO PARAMETER CHECK ---
if (!exists("params") || is.null(params$scenario_rds_path)) {
  scenario_rds_path <- "data/input/Scenario_rds/KH_Scenario_BWK_2025.rds"
} else {
  scenario_rds_path <- params$scenario_rds_path
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

# Variabelen definiëren
oppervlakte_ha <- resultaat$MinOpp_ha
afstand_m      <- resultaat$AfstandBiotopen_m[1]
buffer_m       <- resultaat$Dispersiecap_m[1]

print(resultaat)
rm(df, resultaat)

area_shape  <- vect(here("data/input/Kalmthoutse_Heide.shp"))
master_grid <- rast(here("data/input/Raster_Vlaanderen/Vlaanderen_MasterGrid_10m.tif"))[[1]]

df_namen_sleutel <- read_csv(here("data/input/Excel_files/BWK_Laag_Namen_2025.csv"), show_col_types = FALSE)
gouden_namenlijst <- tolower(trimws(df_namen_sleutel$Laagnaam))

area_shape_proj <- project(area_shape, crs(master_grid))
area_buffer_fix <- buffer(area_shape_proj, width = buffer_m)

message("-> Vertaalraster voor globale/lokale cellen bliksemsnel opbouwen...")

# 1. Maak de lokale uitsnede van het raster
id_raster_KH <- crop(master_grid, area_buffer_fix, snap = "near")

# 2. Rasterize de buffer-vector naar de uitsnede (vele malen sneller dan terra::extract!)
template_KH_mask <- terra::rasterize(area_buffer_fix, id_raster_KH, field = 1)

# 3. Wiskundige index-conversie (lokale cel-index -> coordinaten -> globale cel-index)
lokale_ids  <- terra::cells(template_KH_mask)
coords      <- terra::xyFromCell(template_KH_mask, lokale_ids)
globale_ids <- terra::cellFromXY(master_grid, coords)

# 4. Bouw de vertaaltabel (data.table) op
vertaal_df <- data.table(
  lokale_id  = lokale_ids,
  globale_id = globale_ids
)[!is.na(globale_id)]

studiegebied_globale_ids <- vertaal_df$globale_id

# 5. Templates & grenzen instellen voor vervolgstappen
template_KH <- terra::classify(template_KH_mask, cbind(NA, 0))
values(id_raster_KH) <- NA

grens_web <- sf::st_as_sf(terra::project(area_shape, "EPSG:4326"))

# Schoonmaken van tijdelijke variabelen
rm(template_KH_mask, lokale_ids, coords, globale_ids)
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

raster_simpel_final <- id_raster_KH

rm(tabel_vlaanderen, vertaal_df, lijst_matches, lijst_oppervlaktes)
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

# 1. SNEL UITSNIJDEN: Eerst crop/projectie uitvoeren op het bronraster
r_grb_crop <- terra::crop(r_grb_raw, area_buffer_fix)
r_grb_sync <- terra::resample(r_grb_crop, raster_simpel_final, method = "near")
rm(r_grb_raw, r_grb_crop)
gc()

grb_water_bin <- terra::ifel(r_grb_sync > 0, 1, NA)
rm(r_grb_sync)
gc()

# 2. RAM-VEILIGE PATCHES: Gebruik een tijdelijk bestand op de schijf om std::bad_alloc te voorkomen
temp_patchfile <- tempfile(fileext = ".tif")

water_patches <- terra::patches(
  grb_water_bin, 
  directions = 8, 
  zeroAsNA = TRUE, 
  filename = temp_patchfile, # Dwingt terra om op de schijf te verwerken i.p.v. in het RAM
  overwrite = TRUE
)

# 3. Bereken de frequenties
water_freq <- terra::freq(water_patches)

# Klein water is <= 50 pixels (<= 0.5 ha)
klein_water_ids <- water_freq$value[water_freq$count <= 50]

grutto_water_grb_max <- terra::ifel(water_patches %in% klein_water_ids, 1, NA)

# CRS Gelijktrekken
terra::crs(grutto_water_grb_max) <- terra::crs(raster_simpel_final)
terra::crs(water_bwk_match)      <- terra::crs(raster_simpel_final)
terra::crs(water_bwk_opp)        <- terra::crs(raster_simpel_final)

# Combineren met BWK water
grutto_water_max <- terra::cover(grutto_water_grb_max, water_bwk_match)
grutto_water_opp <- terra::ifel(!is.na(grutto_water_grb_max), 0.01, water_bwk_opp)
grutto_water_opp <- terra::mask(grutto_water_opp, grutto_water_max)

names(grutto_water_max) <- "water_max"
names(grutto_water_opp) <- "water_opp"

# Schoonmaken
rm(water_patches, grutto_water_grb_max, grb_water_bin, water_freq)
if(file.exists(temp_patchfile)) unlink(temp_patchfile)
gc()

# --- AUDIT: BASIS BOUWSTENEN ---
bouwstenen_tabel <- data.frame(
  Bouwsteen = c("Cultuurgrasland", "Hooiland", "Klein Water", "TOTAAL"),
  Potentie_ha = c(
    calc_ha_exact(grutto_cultuur_max),
    calc_ha_exact(grutto_hooiland_max),
    calc_ha_exact(grutto_water_max),
    calc_ha_exact(max(grutto_cultuur_max, grutto_hooiland_max, grutto_water_max, na.rm=TRUE))
  ),
  Werkelijk_ha = c(
    calc_ha_exact(grutto_cultuur_opp),
    calc_ha_exact(grutto_hooiland_opp),
    calc_ha_exact(grutto_water_opp),
    calc_ha_exact(max(grutto_cultuur_opp, grutto_hooiland_opp, grutto_water_opp, na.rm=TRUE))
  )
)
print(bouwstenen_tabel)

# --- GEOPTIMALISEERDE GROENKAART STAP ---
hooggroen_raw         <- rast(here("data/input/ASCI Files/Groenkaart_2021.tif"))
area_buffer_groen_crs <- project(area_buffer_fix, crs(hooggroen_raw))
hooggroen_local_raw   <- crop(hooggroen_raw, area_buffer_groen_crs)

hooggroen_binair  <- terra::classify(hooggroen_local_raw, matrix(c(0.5, 1.5, 1), ncol = 3, byrow = TRUE), others = 0)
hooggroen_10m_raw <- terra::aggregate(hooggroen_binair, fact = 10, fun = "max")
hooggroen_sync    <- project(hooggroen_10m_raw, crs(raster_simpel_final)) %>% terra::resample(raster_simpel_final, method = "near")

# Toepassen masker (Uitsluiten waar hooggroen == 1)
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

r_drain_raw   <- rast(here("data/input/Raster_Vlaanderen/vlaanderen_drainage_10m.tif"))
r_drain_local <- r_drain_raw %>% terra::crop(area_buffer_fix) %>% terra::resample(raster_simpel_final, method = "near")

drain_cats            <- terra::cats(r_drain_local)[[1]]
geselecteerde_letters <- c("d", "e", "f", "h", "i", "g", "e-f", "h-i", "e-i")
grutto_drain_ids      <- drain_cats$value[drain_cats$Label %in% geselecteerde_letters]

masker_drainage  <- r_drain_local %in% grutto_drain_ids
grutto_kern3_max <- terra::ifel(masker_drainage == 1, grutto_basis_kern_max, NA)
grutto_kern3_opp <- terra::mask(grutto_basis_kern_opp, grutto_kern3_max)

rm(r_drain_raw, r_drain_local, masker_drainage)

# --- STAP 4: HERSTELDE BEBOUWINGSFILTER (SCHOON & BINAIR) ---

# 1. Bepaal alle relevante bebouwingslagen uit de BWK-crosswalk (codes gestart met 'u', zonder 'ulm')
bebouwings_lagen <- gouden_namenlijst[grepl("^u[a-z0-9]", gouden_namenlijst) | gouden_namenlijst == "u"]
bebouwings_lagen <- bebouwings_lagen[bebouwings_lagen != "ulm"]

# 2. Filter de regionale database op bebouwing binnen het studiegebied
tabel_vlaanderen_bebouw <- tabel_vlaanderen
tabel_bebouw_KH <- tabel_vlaanderen_bebouw[CODE %in% bebouwings_lagen & cel_id %in% studiegebied_globale_ids]
tabel_bebouw_som <- tabel_bebouw_KH[, .(Present = ifelse(sum(BWK_FRAC, na.rm = TRUE) > 0, 1, 0)), by = .(cel_id)]

r_bebouwing_raw <- id_raster_KH * NA

if(nrow(tabel_bebouw_som) > 0) {
  setnames(tabel_bebouw_som, "cel_id", "globale_id")
  tabel_bebouw_mapping <- merge(tabel_bebouw_som, vertaal_df, by = "globale_id", all.x = TRUE)
  tabel_bebouw_mapping <- tabel_bebouw_mapping[!is.na(lokale_id) & Present == 1]
  
  if(nrow(tabel_bebouw_mapping) > 0) {
    r_bebouwing_raw[tabel_bebouw_mapping$lokale_id] <- 1
  }
}

# 3. HARDE BINAIR-SCHOONMAAK (Voorkomt dat NA's of achtergrondwaarden meegebuffered worden)
r_bebouwing_schoon <- terra::ifel(!is.na(r_bebouwing_raw) & r_bebouwing_raw == 1, 1, NA)

# 4. CLUSTEREN EN FILTEREN OP GROTE DORPSKERNEN (>= 10 ha na 100m netwerkvorming)
res_bebouw_groot <- cluster_filter_compleet(
  masker     = r_bebouwing_schoon,
  opp_laag   = r_bebouwing_schoon,
  drempel_m2 = 100000, # 10 ha
  dist_m     = 100,    # 100m koppelafstand
  werkelijk  = FALSE
)

r_bebouw_groot <- res_bebouw_groot$raster

# 5. UITSENIJDEN VAN DE 100M VERSTORINGSZONE RONDOM ENKEL DE GROTE CLUSTERS
if (!is.null(r_bebouw_groot) && !terra::global(is.na(r_bebouw_groot), "sum")[[1]] == terra::ncell(r_bebouw_groot)) {
  
  # Buffer alleen de goedgekeurde dorpskernen met 100m
  buf_bebouw <- terra::buffer(r_bebouw_groot, width = 100)
  
  # Binariseer de verstoringszone (1 = uitsluiten, NA = geschikt leefgebied)
  masker_bebouw_100m_groot <- terra::ifel(!is.na(buf_bebouw) & buf_bebouw > 0, 1, NA)
  
  # Snijd de 100m verstoringszone uit het kernbiotoop (inverse = TRUE verwijdert de 1-zones)
  grutto_kern_bebouw_max <- terra::mask(grutto_kern3_max, masker_bebouw_100m_groot, inverse = TRUE)
  grutto_kern_bebouw_opp <- terra::mask(grutto_kern3_opp, masker_bebouw_100m_groot, inverse = TRUE)
  
} else {
  # Indien er geen clusters >= 10 ha zijn, treedt er geen verlies op
  grutto_kern_bebouw_max <- grutto_kern3_max
  grutto_kern_bebouw_opp <- grutto_kern3_opp
  masker_bebouw_100m_groot <- id_raster_KH * NA
}

# 6. SCHOONMAAK VAN HET GEHEUGEN (gouden_namenlijst blijft behouden!)
rm(tabel_vlaanderen_bebouw, tabel_bebouw_KH, tabel_bebouw_som, tabel_bebouw_mapping, 
   r_bebouwing_raw, r_bebouwing_schoon, res_bebouw_groot, r_bebouw_groot, tabel_vlaanderen)
gc()

# --- STAP 5: HERSTELDE EN SCHONE HOOGGROEN FILTER ---
# 1. Binaire basiskaart van hooggroen (1 = groen, NA = rest)
r_groen_bin <- terra::ifel(!is.na(hooggroen_sync) & hooggroen_sync > 0, 1, NA)

# 2. MORFOLOGISCHE DENSITY CHECK (5x5 venster = 25 cellen = ~25 are)
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

# 3. Behoud boskernen (minstens 12 buren in 5x5 venster = >= ~20 are)
r_groen_20are <- terra::ifel(groen_som >= 12 & r_groen_bin == 1, 1, NA)

# 4. MAAK DE 50M ULEITSLUITINGSZONE RONDOM BOSKERNEN
buf_groen <- terra::buffer(r_groen_20are, width = 50)

# Binariseer het bosmasker strak: 1 = verstoringszone, NA = geschikt leefgebied
masker_groen_50m_cluster <- terra::ifel(!is.na(buf_groen) & buf_groen > 0, 1, NA)

# 5. SNIJD DE 50M BOSZONE UIT HET KERNBIOTOOP (inverse = TRUE verwijdert de 1-zones)
grutto_kern_finaal_max <- terra::mask(grutto_kern_bebouw_max, masker_groen_50m_cluster, inverse = TRUE)
grutto_kern_finaal_opp <- terra::mask(grutto_kern_bebouw_opp, masker_groen_50m_cluster, inverse = TRUE)

# Opruimen van tijdelijke bestanden en geheugen
rm(r_groen_bin, groen_som, r_groen_20are, hooggroen_sync, buf_groen)
gc()

# --- STAP 6.1 (GEFIKST): CLUSTER-TABEL MET MAXIMALE POTENTIE ---

drempel_kern_m2 <- 10 * 10000 # 10 ha

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

# Finale cijfers ophalen
opp_max_finaal <- terra::global(grutto_kern_10ha_max * 0.01, "sum", na.rm = TRUE)[1,1]
opp_wer_finaal <- terra::global(grutto_kern_10ha_opp * 0.01, "sum", na.rm = TRUE)[1,1]

cat("\n--- FINALE ANALYSE AFGEROND (STRIKT GESCHEIDEN) ---\n")
cat("Totaal ha in maximale potentie kaart:", round(opp_max_finaal, 2), "ha\n")
cat("Totaal ha in werkelijke oppervlakte kaart:  ", round(opp_wer_finaal, 2), "ha\n")

# --- STAP 7: SAMENVOEGEN BOUWSTENEN (GECORRIGEERDE RASTER MAX LOGICA) ---

# 1. MAXIMALE POTENTIE (Binair)
grutto_bouwsteen_totaal_max <- terra::cover(grutto_cultuur_max, grutto_hooiland_max)
grutto_bouwsteen_totaal_max <- terra::cover(grutto_bouwsteen_totaal_max, grutto_water_max)
names(grutto_bouwsteen_totaal_max) <- "bouwsteen_totaal_max"

# 2. WERKELIJKE OPPERVLAKTE (Conservatieve samenvoeging van fracties)
grutto_bouwsteen_totaal_opp <- terra::cover(grutto_cultuur_opp, grutto_hooiland_opp)
grutto_bouwsteen_totaal_opp <- terra::cover(grutto_bouwsteen_totaal_opp, grutto_water_opp)
names(grutto_bouwsteen_totaal_opp) <- "bouwsteen_totaal_opp"

# 3. Opschonen: pixels zonder habitat op NA zetten
grutto_bouwsteen_totaal_max <- terra::ifel(grutto_bouwsteen_totaal_max == 0, NA, grutto_bouwsteen_totaal_max)
grutto_bouwsteen_totaal_opp <- terra::ifel(grutto_bouwsteen_totaal_opp == 0, NA, grutto_bouwsteen_totaal_opp)

ha_bouw_start <- as.numeric(terra::global(grutto_bouwsteen_totaal_opp, "sum", na.rm=TRUE)[1,1]) * 0.01
print(paste("Startoppervlakte gecombineerde bouwstenen:", round(ha_bouw_start, 2), "ha"))

# --- STAP 8, 9 & 10: SAMENVOEGEN BIJKOMEND HABITAT EN HIERSTELDE FILTERS ---

# 1. Cultuurhooi1b synchroniseren
grutto_1b_match_sync <- c1b_match
grutto_1b_opp_sync   <- c1b_opp

grutto_1b_max <- grutto_1b_match_sync
grutto_1b_max[grutto_1b_max == 0] <- NA
grutto_1b_opp <- terra::mask(grutto_1b_opp_sync, grutto_1b_max)

names(grutto_1b_max) <- "1b_max"
names(grutto_1b_opp) <- "1b_opp"

# 2. Samenvoegen van alle bijkomende bouwstenen
grutto_basis_bouw_1b_max <- terra::cover(grutto_bouwsteen_totaal_max, grutto_1b_max)
grutto_cultuurhooi1_opp  <- terra::cover(grutto_bouwsteen_totaal_opp, grutto_1b_opp)

# 3. HARDE BINAIR-SCHOONMAAK VAN DE 3 UITSLUITINGSMASKERS (1 = uitsluiten, NA = geschikt)
# A. Akkerfilter
filter_akker <- terra::ifel(!is.na(akker_match) & akker_match == 1, 1, NA)

# B. Bebouwingsfilter (hergebruik het masker uit Stap 4 of maak het strak aan)
if(exists("masker_bebouw_100m_groot") && !all(is.na(terra::values(masker_bebouw_100m_groot, mat=FALSE)))) {
  filter_bebouw <- terra::ifel(!is.na(masker_bebouw_100m_groot) & masker_bebouw_100m_groot > 0, 1, NA)
} else {
  filter_bebouw <- id_raster_KH * NA
}

# C. Hooggroenfilter (hergebruik het masker uit Stap 5 of maak het strak aan)
if(exists("masker_groen_50m_cluster") && !all(is.na(terra::values(masker_groen_50m_cluster, mat=FALSE)))) {
  filter_groen <- terra::ifel(!is.na(masker_groen_50m_cluster) & masker_groen_50m_cluster > 0, 1, NA)
} else {
  filter_groen <- id_raster_KH * NA
}

# 4. STAPSGEWIJS UITSENIJDEN MET INVERSE = TRUE
# Maximaal Potentieel (Spoor A)
grutto_c1_m1_max     <- terra::mask(grutto_basis_bouw_1b_max, filter_akker, inverse = TRUE)
grutto_c1_m2_max     <- terra::mask(grutto_c1_m1_max, filter_bebouw, inverse = TRUE)
grutto_c1_finaal_max <- terra::mask(grutto_c1_m2_max, filter_groen, inverse = TRUE)

# Werkelijke Oppervlakte (Spoor B)
grutto_c1_m1_opp     <- terra::mask(grutto_cultuurhooi1_opp, filter_akker, inverse = TRUE)
grutto_c1_m2_opp     <- terra::mask(grutto_c1_m1_opp, filter_bebouw, inverse = TRUE)
grutto_c1_finaal_opp <- terra::mask(grutto_c1_m2_opp, filter_groen, inverse = TRUE)

grutto_scenA_finaal_max <- grutto_c1_finaal_max
grutto_scenB_finaal_opp <- grutto_c1_finaal_opp

# 5. RESULTATEN BEREKENEN EN CONTROLEREN
ha_finaal_clean_maxpot <- calc_ha_exact(grutto_c1_finaal_max)
ha_finaal_clean_opp    <- calc_ha_exact(grutto_c1_finaal_opp)

print(paste("Oppervlakte grutto_cultuurhooi1 max potentie (Zonder akker/bebouw/hooggroen):", round(ha_finaal_clean_maxpot, 2), "ha"))
print(paste("Oppervlakte grutto_cultuurhooi1 werkelijke oppervlaktes (Zonder akker/bebouw/hooggroen):", round(ha_finaal_clean_opp, 2), "ha"))

# Opruimen
rm(grutto_basis_bouw_1b_max, grutto_c1_m1_max, grutto_c1_m2_max, grutto_c1_m1_opp, grutto_c1_m2_opp, 
   filter_akker, filter_bebouw, filter_groen, grutto_cultuurhooi1_opp, grutto_bouwsteen_totaal_max, grutto_bouwsteen_totaal_opp)
gc()

# --- STAP 10.5: RUIMTELIJKE CONTEXT & CATEGORISATIE (VOORBEREIDING) ---

# De 25 ha filter is hier weggehaald, omdat we deze pas NA de 200m koppeling kunnen berekenen.
grutto_scenA_finaal_max <- grutto_c1_finaal_max
grutto_scenB_finaal_opp <- grutto_c1_finaal_opp

# Randvoorwaarden voor r_status bepaling
grens_merged  <- terra::aggregate(area_shape_proj)
buffer_lijn   <- terra::buffer(grens_merged, width = buffer_m)
binnen_masker <- terra::rasterize(grens_merged, template_KH, field = 1, background = NA)

if (!all(is.na(terra::values(grutto_kern_10ha_max, mat=FALSE)))) {
  is_inside <- !is.na(binnen_masker) 
  is_loss   <- !is.na(grutto_kern_10ha_max) & is.na(grutto_kern_10ha_opp)
  is_kept   <- !is.na(grutto_kern_10ha_opp)

  r_status <- terra::ifel(is_loss & is_inside, 1, NA)       
  r_status <- terra::cover(r_status, terra::ifel(is_loss & !is_inside, 2, NA)) 
  r_status <- terra::cover(r_status, terra::ifel(is_kept & is_inside, 3, NA))  
  r_status <- terra::cover(r_status, terra::ifel(is_kept & !is_inside, 4, NA)) 
} else {
  r_status <- terra::rast(template_KH, vals = NA)
  message("Let op: Geen geschikte clusters gevonden voor deze soort.")
}

# --- STAP 11: WEDERZIJDSE KOPPELING MET INTEGRALE FRACTIE-DREMPELS (STRIKT PER RASTER) ---

voer_koppeling_en_drempel_uit_strikt <- function(kern_rast, bouw_rast, afstand = 200, min_kern_ha = 10, min_totaal_ha = 25) {
  if (is.null(kern_rast) || terra::global(is.na(kern_rast), "sum")[[1]] == terra::ncell(kern_rast)) return(NULL)
  
  # 1. Maak binaire kaarten
  r_totaal_bin <- terra::cover(
    terra::ifel(!is.na(kern_rast) & kern_rast > 0, 1, NA),
    terra::ifel(!is.na(bouw_rast) & bouw_rast > 0, 1, NA)
  )
  
  if (terra::global(is.na(r_totaal_bin), "sum")[[1]] == terra::ncell(r_totaal_bin)) return(NULL)
  
  # 2. Netwerkvorming via 100m buffer (200m totale overbrugging tussen 2 ruimtelijke elementen)
  r_buf <- terra::buffer(r_totaal_bin, width = afstand / 2)
  cl_netwerk <- terra::patches(r_buf, directions = 8, zeroAsNA = TRUE)
  cl_netwerk_biotoop <- terra::mask(cl_netwerk, r_totaal_bin)
  
  # 3. Bereken het TOTALE netto-oppervlak per netwerk (Spoor A = cellen, Spoor B = fracties)
  r_comb_val <- terra::cover(kern_rast, bouw_rast)
  stats_netwerk <- terra::zonal(r_comb_val * 0.01, cl_netwerk_biotoop, fun = "sum", na.rm = TRUE)
  colnames(stats_netwerk) <- c("netwerk_id", "totaal_ha")
  
  # 4. Filter op netto ankerkern (Check of netwerk een kern bevat die >= min_kern_ha is)
  cl_kernen_only <- terra::mask(cl_netwerk_biotoop, kern_rast)
  stats_kernen <- terra::zonal(kern_rast * 0.01, cl_kernen_only, fun = "sum", na.rm = TRUE)
  colnames(stats_kernen) <- c("netwerk_id", "kern_ha")
  
  valide_netwerken <- stats_netwerk %>%
    inner_join(stats_kernen, by = "netwerk_id") %>%
    filter(totaal_ha >= min_totaal_ha & kern_ha >= min_kern_ha) %>%
    pull(netwerk_id)
  
  if (length(valide_netwerken) == 0) return(NULL)
  
  # 5. Maskeer en retourneer
  m_valide <- cl_netwerk_biotoop %in% valide_netwerken
  m_valide <- terra::ifel(m_valide == 1, 1, NA)
  
  return(list(
    kern = terra::mask(kern_rast, m_valide),
    bouw = terra::mask(bouw_rast, m_valide)
  ))
}

# --- UITVOERING: ELK RASTER VOLLEDIG AUTONOOM GEFILTERD ---
# Spoor A (Theorie): Checkt op binaire pixels (min_kern=10ha, min_totaal=25ha)
resA_strikt <- voer_koppeling_en_drempel_uit_strikt(grutto_kern_10ha_max, grutto_scenA_finaal_max, 200, 10, 25)

# Spoor B (Werkelijkheid): Checkt op PURE FRACTIE-WAARDEN (min_kern=10ha, min_totaal=25ha netto!)
resB_strikt <- voer_koppeling_en_drempel_uit_strikt(grutto_kern_10ha_opp, grutto_scenB_finaal_opp, 200, 10, 25)

# Ververs de kaart-raster variabele voor de visualisatie
grutto_kaart_A <- maak_kaart_laag(resA_strikt)

# --- STAP 12: FUNCTIES & POLYGONEN VOOR DE KAART ---

# GEFIXT: De functie herkent nu of het om Spoor A (max) of Spoor B (opp) gaat voor zuivere popups
prepare_cluster_polys_v2 <- function(cl_raster, habitat_raster, is_opp = FALSE) {
  if(is.null(cl_raster) || all(is.na(terra::values(cl_raster, mat=FALSE)))) return(NULL)
  
  polys <- terra::as.polygons(cl_raster, aggregate = TRUE) %>% 
    terra::disagg() %>% 
    sf::st_as_sf() %>% 
    sf::st_transform(4326)
  colnames(polys)[1] <- "ID"
  
  # Bereken de exacte habitat-hectares per cluster via zonal sum
  z <- terra::zonal(habitat_raster * 0.01, cl_raster, fun = "sum", na.rm = TRUE)
  colnames(z) <- c("ID", "Habitat_ha")
  
  info <- as.data.frame(z) %>% mutate(
    Popup = paste0("<strong>Netwerk ID:</strong> ", ID, 
                   "<br><strong>Netto Kernoppervlakte:</strong> ", round(Habitat_ha, 2), " ha")
  )
  return(left_join(polys, info, by = "ID"))
}

# --- UPDATE VOOR KAARTLAGEN (GEBASEERD OP DE INTEGRALE NETWERKEN) ---
if (!is.null(resB_strikt) && !all(is.na(terra::values(resB_strikt$kern, mat=FALSE)))) {
  # We maken on-the-fly de finale, unieke cluster-IDs van Spoor B (de werkelijkheid)
  cl_kernen_finaal <- terra::patches(resB_strikt$kern, directions = 8, zeroAsNA = TRUE)
  cl_bouw_finaal   <- terra::patches(resB_strikt$bouw, directions = 8, zeroAsNA = TRUE)
  
  poly_kernen_all <- prepare_cluster_polys_v2(cl_kernen_finaal, resB_strikt$kern, is_opp = TRUE)
  poly_bouw_all   <- prepare_cluster_polys_v2(cl_bouw_finaal, resB_strikt$bouw, is_opp = TRUE)
} else {
  poly_kernen_all <- NULL
  poly_bouw_all   <- NULL
}

# Finale Gekoppelde Selectie (De felle ingekleurde vlakken op de kaart)
# GEFIXT: We bouwen de kaart op basis van final_opp_export zodat de fracties/netwerken kloppen met de tabel!
if(!is.null(resB_strikt) && !all(is.na(terra::values(resB_strikt$kern, mat=FALSE)))) {
  r_kaart_opp <- terra::cover(
    terra::ifel(!is.na(resB_strikt$kern), 1, NA),
    terra::ifel(!is.na(resB_strikt$bouw), 2, NA)
  )
  
  poly_finaal <- terra::as.polygons(r_kaart_opp, aggregate = TRUE) %>% 
    terra::disagg() %>% 
    sf::st_as_sf() %>% 
    sf::st_transform(4326)
  
  colnames(poly_finaal)[1] <- "Type_ID"
  poly_finaal$Label <- ifelse(poly_finaal$Type_ID == 1, "Geselecteerde Kern", "Geselecteerde Bouwsteen")
} else {
  poly_finaal <- NULL
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
potentie_export_rast <- if (exists("grutto_kaart_A") && !is.null(grutto_kaart_A) && !all(is.na(terra::values(grutto_kaart_A, mat=FALSE)))) {
  terra::ifel(grutto_kaart_A > 0, 1, NA)
} else if (exists("final_max") && !all(is.na(terra::values(final_max, mat=FALSE)))) {
  terra::ifel(!is.na(final_max) & final_max > 0, 1, NA)
} else {
  terra::rast(template_KH, vals = NA)
}

# 2. Bepaal Werkelijke Oppervlakte Raster
werkelijk_export_rast <- if (exists("resB_strikt") && !is.null(resB_strikt) && (!all(is.na(terra::values(resB_strikt$kern, mat=FALSE))) || !all(is.na(terra::values(resB_strikt$bouw, mat=FALSE))))) {
  r_net_totaal_opp <- terra::cover(resB_strikt$kern, resB_strikt$bouw)
  terra::ifel(!is.na(r_net_totaal_opp) & r_net_totaal_opp > 0, 1, NA)
} else if (exists("final_opp") && !all(is.na(terra::values(final_opp, mat=FALSE)))) {
  terra::ifel(!is.na(final_opp) & final_opp > 0, 1, NA)
} else {
  terra::rast(template_KH, vals = NA)
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

