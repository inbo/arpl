library(here)
# CRUCIALE FIX: Dwing R Markdown om te werken vanaf de hoofdmap

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

# ==============================================================================
# HOOFDFUNCTIE: CLUSTEREN + OPTIONELE LINIE-ZUIVERING (MORFOLOGISCHE OPENING)
# ==============================================================================
cluster_filter_compleet <- function(masker, opp_laag, drempel_m2, dist_m, werkelijk = FALSE, straal_m = 0) {
  # 1. Controleer of het raster leeg is
  if (terra::global(is.na(masker), "sum")[[1]] == terra::ncell(masker)) {
    return(list(raster = masker * NA, clusters = masker * NA))
  }
  
  # ----------------------------------------------------------------------------
  # GEHEUGEN-OPTIMALISATIE: Snijd rasters bij tot het actieve studiegebied
  # ----------------------------------------------------------------------------
  if (exists("id_raster_TV")) {
    masker   <- terra::crop(masker, id_raster_TV, snap = "out")
    opp_laag <- terra::crop(opp_laag, id_raster_TV, snap = "out")
  }
  
  # ----------------------------------------------------------------------------
  # STAP 1: LINIE-ZUIVERING VIA VECTORIËLE INWAARTSE BUFFER
  # ----------------------------------------------------------------------------
  if (straal_m > 0) {
    r_binair <- terra::ifel(!is.na(masker) & masker > 0, 1, NA)
    
    if (terra::global(is.na(r_binair), "sum")[[1]] == terra::ncell(r_binair)) {
      return(list(raster = masker * NA, clusters = masker * NA))
    }
    
    # Vectoriële inwaartse buffer
    v_biotoop <- terra::as.polygons(r_binair, aggregate = TRUE)
    v_eroded  <- terra::buffer(v_biotoop, width = -straal_m)
    
    if (length(v_eroded) == 0 || terra::geomtype(v_eroded) == "none") {
      return(list(raster = masker * NA, clusters = masker * NA))
    }
    
    # Uitwaartse buffer om randen te herstellen
    v_dilated <- terra::buffer(v_eroded, width = straal_m)
    
    # Terugzetten naar raster
    r_zuiver <- terra::rasterize(v_dilated, r_binair, field = 1, background = NA)
    r_masker_werkelijk <- terra::mask(r_zuiver, masker)
    
    rm(v_biotoop, v_eroded, v_dilated, r_zuiver)
  } else {
    r_masker_werkelijk <- masker
  }
  
  if (terra::global(is.na(r_masker_werkelijk), "sum")[[1]] == terra::ncell(r_masker_werkelijk)) {
    return(list(raster = masker * NA, clusters = masker * NA))
  }
  
  # ----------------------------------------------------------------------------
  # STAP 2: NETWERKVORMING / CLUSTERING
  # ----------------------------------------------------------------------------
  if (dist_m > 0) {
    r_binair_z <- terra::ifel(!is.na(r_masker_werkelijk) & r_masker_werkelijk > 0, 1, NA)
    r_buffered <- terra::buffer(r_binair_z, width = dist_m / 2)
    cl_network <- terra::patches(r_buffered, directions = 4, zeroAsNA = TRUE)
    cl_biotoop_only <- terra::mask(cl_network, r_masker_werkelijk)
  } else {
    cl_network <- terra::patches(r_masker_werkelijk, directions = 8, zeroAsNA = TRUE)
    cl_biotoop_only <- cl_network
  }
  
  # ----------------------------------------------------------------------------
  # STAP 3: OPPERVLAKTE-OPTELSOM PER GEKOPPELD NETWERK
  # ----------------------------------------------------------------------------
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
  
  # ----------------------------------------------------------------------------
  # STAP 4: FILTEREN OP TOTALE NETWERK-OPPERVLAKTE >= DREMPEL
  # ----------------------------------------------------------------------------
  voldoet_ids <- stats_df$ID[stats_df$Area_m2 >= drempel_m2]
  if(length(voldoet_ids) == 0) return(list(raster = masker * NA, clusters = masker * NA))
  
  # ----------------------------------------------------------------------------
  # STAP 5: FINALE RASTER OPBOUW & IN-MEMORY FIX
  # ----------------------------------------------------------------------------
  masker_binair <- cl_biotoop_only %in% voldoet_ids
  final_network_mask <- terra::ifel(masker_binair == 1, 1, NA)
  
  r_finaal  <- terra::mask(r_masker_werkelijk, final_network_mask)
  cl_finaal <- terra::mask(cl_biotoop_only, r_finaal) 
  
  # HOOFDFIX: Dwing R om de uitkomst fysiek in de RAM te laden
  # Hierdoor is de afhankelijkheid van tijdelijke .tif bestanden op schijf verdwenen
  r_finaal  <- terra::deepcopy(r_finaal)
  cl_finaal <- terra::deepcopy(cl_finaal)
  
  return(list(raster = r_finaal, clusters = cl_finaal))
}

terraOptions(
  memfrac = 0.8,         # Dwing terra om tot max. 80% van het RAM-geheugen te gebruiken
  tempdir = tempdir(),   # Geef toestemming voor automatische disk-swapping bij zware rasters
  verbose = FALSE
)

df <- read_excel(here::here("data/input/Excel_files/Soorten_bwk_afstanden.xlsx"))
soort <- "bruinekiekendief"

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

print(resultaat)
rm(df, resultaat)

area_shape  <- vect(here("data/input/Turnhouts_Vennegebied.shp"))
master_grid <- rast(here("data/input/Raster_Vlaanderen/Vlaanderen_MasterGrid_10m.tif"))[[1]]

df_namen_sleutel <- read_csv(here("data/input/Excel_files/BWK_Laag_Namen_2025.csv"), show_col_types = FALSE)
gouden_namenlijst <- tolower(trimws(df_namen_sleutel$Laagnaam))

area_shape_proj <- project(area_shape, crs(master_grid))
area_buffer_fix <- buffer(area_shape_proj, width = buffer_m)

message("-> Vertaalraster voor globale/lokale cellen opbouwen...")
id_raster_TV <- crop(master_grid, area_buffer_fix, snap = "near")

globale_id_raster <- master_grid
values(globale_id_raster) <- 1:ncell(globale_id_raster)
id_raster_TV_globale_values <- crop(globale_id_raster, area_buffer_fix, snap = "near")

id_raster_TV_masked <- mask(id_raster_TV_globale_values, area_buffer_fix)

lokale_ids  <- cells(id_raster_TV_masked) 
globale_ids <- id_raster_TV_masked[lokale_ids][[1]]

vertaal_df <- data.table(
  lokale_id  = lokale_ids,
  globale_id = globale_ids
)

vertaal_df <- vertaal_df[!is.na(globale_id)]
studiegebied_globale_ids <- unique(vertaal_df$globale_id)

values(id_raster_TV) <- NA
template_TV <- terra::rasterize(area_buffer_fix, id_raster_TV, field = 1, background = 0)
grens_web   <- sf::st_as_sf(terra::project(area_shape, "EPSG:4326"))

rm(globale_id_raster, id_raster_TV_globale_values, id_raster_TV_masked, lokale_ids, globale_ids)
gc()

# Biotoopfiltering (BWK-codes)
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

voortplanting_bwk_max <- lijst_matches[["voortplanting_bwk"]]
voortplanting_bwk_opp <- lijst_oppervlaktes[["voortplanting_bwk"]]
bos_max               <- lijst_matches[["bos"]]
bos_opp               <- lijst_oppervlaktes[["bos"]]
foerageer_max         <- lijst_matches[["foerageer"]]
foerageer_opp         <- lijst_oppervlaktes[["foerageer"]]

rm(tabel_vlaanderen, vertaal_df, lijst_matches, lijst_oppervlaktes)
gc()

# 1. Laad het gecropte coderaster in
r_lgp_TV <- crop(rast(here("data/input/Raster_Vlaanderen/vlaanderen_lgp_gwscod_2025_10m.tif")), id_raster_TV, snap = "near")

# 2. Definieer de EXACTE gewascodes
kiekendief_codes <- c(311, 321, 331)

# 3. Filter het raster
bruinekiekendief_lgp <- terra::ifel(r_lgp_TV %in% kiekendief_codes, 1, NA)

# 4. Combineer met BWK
bruinekiekendief_voortplanting_max <- terra::cover(bruinekiekendief_lgp, voortplanting_bwk_max)
bruinekiekendief_voortplanting_opp <- terra::cover(bruinekiekendief_lgp, voortplanting_bwk_opp)

rm(r_lgp_TV, bruinekiekendief_lgp)

# Bosmaskering overgeslagen conform instellingen
bruinekiekendief_vp_zonderbos_max <- bruinekiekendief_voortplanting_max
bruinekiekendief_vp_zonderbos_opp <- bruinekiekendief_voortplanting_opp

# ==============================================================================
# STAP 3: VOORTPLANTINGSCLUSTERS (LINIE-ZUIVERING 40M + 100M FUZZY + MIN. 10 HA)
# ==============================================================================
message("-> Voortplantingsgebieden zuiveren (min. 40m breed), clusteren (100m) en filteren op min. 10 ha...")

drempel_vp_m2 <- 10 * 10000 # 10 ha

# Spoor A (Max) - Met linie-zuivering (straal_m = 20 -> minstens 40m fysieke breedte)
res_vp_max <- cluster_filter_compleet(
  masker     = bruinekiekendief_vp_zonderbos_max,
  opp_laag   = bruinekiekendief_vp_zonderbos_max,
  drempel_m2 = drempel_vp_m2,
  dist_m     = 100,
  werkelijk  = FALSE,
  straal_m   = 20
)
bruinekiekendief_voortplanting_final_max <- res_vp_max$raster

# Spoor B (Opp) - Met linie-zuivering op de ruimtelijke structuur
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
# STAP 4: FOERAGEERGEBIED CLUSTEREN (100M FUZZY, MIN. 100 HA - ZONDER LINIE-ZUIVERING)
# ==============================================================================
message("-> Foerageergebieden clusteren (100m fuzzy) en filteren op min. 100 ha...")

drempel_foer_m2 <- 100 * 10000 # 100 ha

# Spoor A (Max)
res_foer_max <- cluster_filter_compleet(
  masker     = foerageer_max,
  opp_laag   = foerageer_max,
  drempel_m2 = drempel_foer_m2,
  dist_m     = 100,
  werkelijk  = FALSE,
  straal_m   = 0 # Geen linie-zuivering voor foerageergebied
)
bruinekiekendief_foerageer_final_max <- res_foer_max$raster

# Spoor B (Opp)
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
# STAP 5: FINALE LEEFGEBIED SYNTHESE (Foerageergebied als omgevingsfilter)
# GEHEUGEN-GEOPTIMALISEERD: Geen st_distance matrix meer om std::bad_alloc te voorkomen
# ==============================================================================
message("-> Starten finale leefgebied synthese (5 km omgevingstoetsing)...")

# --- 5.1: THEORETISCH MAXIMUM (MAX-SPOOR) ---
r_patches_vp_max   <- terra::patches(bruinekiekendief_voortplanting_final_max, directions = 8, zeroAsNA = TRUE)
r_patches_foer_max <- terra::patches(bruinekiekendief_foerageer_final_max, directions = 8, zeroAsNA = TRUE)

if (!all(is.na(terra::values(r_patches_vp_max, mat = FALSE))) && !all(is.na(terra::values(r_patches_foer_max, mat = FALSE)))) {
  
  # 1. Bereken een raster met de afstand tot de dichtstbijzijnde foerageerpatch
  dist_foer_max <- terra::distance(r_patches_foer_max)
  
  # 2. Snijd de afstandskaart af op de voortplantingspatches
  vp_dist_max <- terra::mask(dist_foer_max, r_patches_vp_max)
  
  # 3. Bepaal per voortplantingspatch (ID) de minimale afstand tot foerageergebied
  min_dist_df_max <- terra::zonal(vp_dist_max, r_patches_vp_max, fun = "min", na.rm = TRUE)
  colnames(min_dist_df_max) <- c("patch_id", "dist")
  
  # 4. Selecteer de ID's die binnen 5000m liggen
  vp_ids_valide_max <- min_dist_df_max$patch_id[min_dist_df_max$dist <= 5000]
  
  bruinekiekendief_leefgebied1_max <- terra::ifel(r_patches_vp_max %in% vp_ids_valide_max, 1, NA)
  
  rm(dist_foer_max, vp_dist_max, min_dist_df_max)
} else {
  bruinekiekendief_leefgebied1_max <- template_TV * NA
}

# --- 5.2: REALISTISCHE OPPERVLAKTE (OPP-SPOOR) ---
r_binair_vp_opp   <- terra::ifel(!is.na(bruinekiekendief_voortplanting_final_opp) & bruinekiekendief_voortplanting_final_opp > 0, 1, NA)
r_binair_foer_opp <- terra::ifel(!is.na(bruinekiekendief_foerageer_final_opp) & bruinekiekendief_foerageer_final_opp > 0, 1, NA)

if (!all(is.na(terra::values(r_binair_vp_opp, mat = FALSE))) && !all(is.na(terra::values(r_binair_foer_opp, mat = FALSE)))) {
  
  r_patches_vp_opp   <- terra::patches(r_binair_vp_opp, directions = 8, zeroAsNA = TRUE)
  r_patches_foer_opp <- terra::patches(r_binair_foer_opp, directions = 8, zeroAsNA = TRUE)
  
  # 1. Bereken afstandsgrid tot foerageergebied
  dist_foer_opp <- terra::distance(r_patches_foer_opp)
  
  # 2. Maskeer met voortplantingspatches
  vp_dist_opp <- terra::mask(dist_foer_opp, r_patches_vp_opp)
  
  # 3. Vind minimale afstand per patch_id
  min_dist_df_opp <- terra::zonal(vp_dist_opp, r_patches_vp_opp, fun = "min", na.rm = TRUE)
  colnames(min_dist_df_opp) <- c("patch_id", "dist")
  
  # 4. Selecteer geldige ID's
  vp_ids_valide_opp <- min_dist_df_opp$patch_id[min_dist_df_opp$dist <= 5000]
  
  m_vp_opp <- r_patches_vp_opp %in% vp_ids_valide_opp
  bruinekiekendief_leefgebied1_opp <- terra::mask(bruinekiekendief_voortplanting_final_opp, terra::ifel(m_vp_opp, 1, NA))
  
  rm(r_patches_vp_opp, r_patches_foer_opp, dist_foer_opp, vp_dist_opp, min_dist_df_opp)
} else {
  bruinekiekendief_leefgebied1_opp <- template_TV * NA
}

rm(r_patches_vp_max, r_patches_foer_max, r_binair_vp_opp, r_binair_foer_opp)
gc()

totaal_ha_max <- calc_ha_exact(bruinekiekendief_leefgebied1_max)
totaal_ha_opp <- calc_ha_exact(bruinekiekendief_leefgebied1_opp)

message("========================================================================")
message("-> FINALE GEVALIDEERDE BROEDLOCATIES BRUINE KIEKENDIEF VOLTOOID:")
message("-> [Spoor MAX] Totaal goedgekeurd broedgebied (Bruto): ", round(totaal_ha_max, 2), " ha")
message("-> [Spoor OPP] Totaal goedgekeurd broedgebied (Netto): ", round(totaal_ha_opp, 2), " ha")
message("========================================================================")


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

