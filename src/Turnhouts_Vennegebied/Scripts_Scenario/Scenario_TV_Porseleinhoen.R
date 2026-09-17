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

cluster_filter_compleet <- function(masker, opp_laag = NULL, drempel_m2, dist_m, werkelijk = FALSE) {
  if (is.null(masker) || terra::global(is.na(masker), "sum")[[1]] == terra::ncell(masker)) {
    return(list(raster = masker * NA, clusters = masker * NA))
  }
  
  r_binair <- terra::ifel(!is.na(masker) & masker > 0, 1, NA)
  
  if (dist_m > 0) {
    r_buffered <- terra::buffer(r_binair, width = dist_m / 2)
    cl_network <- terra::patches(r_buffered, directions = 4, zeroAsNA = TRUE)
    cl_biotoop_only <- terra::mask(cl_network, masker)
    rm(r_buffered, cl_network)
  } else {
    cl_network <- terra::patches(masker, directions = 8, zeroAsNA = TRUE)
    cl_biotoop_only <- cl_network
  }
  
  if (werkelijk && !is.null(opp_laag)) {
    stats_df <- terra::zonal(opp_laag, cl_biotoop_only, fun = "sum", na.rm = TRUE)
    colnames(stats_df) <- c("ID", "Waarde")
    stats_df$Area_m2 <- stats_df$Waarde * 100 
  } else {
    f <- terra::freq(cl_biotoop_only)
    stats_df <- data.frame(ID = f$value, Waarde = f$count)
    stats_df$Area_m2 <- stats_df$Waarde * 100 
  }
  
  stats_df <- stats_df[!is.na(stats_df$ID), ]
  if (nrow(stats_df) == 0) return(list(raster = masker * NA, clusters = masker * NA))
  
  voldoet_ids <- stats_df$ID[stats_df$Area_m2 >= drempel_m2]
  if (length(voldoet_ids) == 0) return(list(raster = masker * NA, clusters = masker * NA))
  
  masker_binair <- cl_biotoop_only %in% voldoet_ids
  final_network_mask <- terra::ifel(masker_binair == 1, 1, NA)
  
  r_finaal  <- terra::mask(masker, final_network_mask)
  cl_finaal <- terra::mask(cl_biotoop_only, r_finaal) 
  
  rm(masker_binair, final_network_mask, r_binair)
  gc()
  
  return(list(raster = r_finaal, clusters = cl_finaal))
}

terraOptions(
  memfrac = 0.8,        # Dwing terra om tot max. 80% van het RAM-geheugen te gebruiken
  tempdir = tempdir(),  # Geef toestemming voor automatische disk-swapping bij zware rasters
  verbose = FALSE
)

df <- read_excel(here::here("data/input/Excel_files/Soorten_bwk_afstanden.xlsx"))
soort <- "porseleinhoen"

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
buffer_m       <- 50000 #resultaat$Dispersiecap_m[1]

print(resultaat)
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
values(globale_id_raster) <- 1:ncell(globale_id_raster)

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

cat("Gecorrigeerd aantal pixels in template_TV: ", sum(terra::values(template_TV) == 1, na.rm=TRUE), "\n")

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

grasland_max <- lijst_matches[["grasland"]]
grasland_opp <- lijst_oppervlaktes[["grasland"]]
moeras_max   <- lijst_matches[["moeras"]]
moeras_opp   <- lijst_oppervlaktes[["moeras"]]

rm(vertaal_df, lijst_matches, lijst_oppervlaktes)
gc()

message("=== STAP 1: Biotopen combineren en uitsluitingszones toepassen ===")

# --- 1. COMBINEER GRASLAND EN MOERAS ---
biotoop_max <- terra::cover(grasland_max, moeras_max)

biotoop_opp_clean_gras <- terra::ifel(is.na(grasland_opp), 0, grasland_opp)
biotoop_opp_clean_moer <- terra::ifel(is.na(moeras_opp), 0, moeras_opp)
biotoop_opp_som <- biotoop_opp_clean_gras + biotoop_opp_clean_moer
biotoop_opp_raw <- terra::clamp(biotoop_opp_som, upper = 1.0)
biotoop_opp <- terra::ifel(!is.na(biotoop_max) & biotoop_opp_raw > 0, biotoop_opp_raw, NA)

# --- 2. HOOGGROEN EN TERRILS INLEZEN & ALIGNEREN ---
groenkaart <- terra::rast(here("data/input/ASCI Files/Groenkaart_2021.tif"))
groenkaart_TV <- terra::crop(groenkaart, template_TV, snap = "near")
hooggroen_aligned <- terra::resample(groenkaart_TV, template_TV, method = "near")

tabel_terril <- tabel_vlaanderen[CODE %like% "^kg"]
tabel_terril_TV <- tabel_terril[cel_id %in% studiegebied_globale_ids]

r_terril <- id_raster_TV * NA

if(nrow(tabel_terril_TV) > 0) {
  setnames(tabel_terril_TV, "cel_id", "globale_id")
  tabel_terril_mapping <- merge(tabel_terril_TV, vertaal_df, by = "globale_id", all.x = TRUE)
  tabel_terril_mapping <- tabel_terril_mapping[!is.na(lokale_id)]
  
  if(nrow(tabel_terril_mapping) > 0) {
    r_terril[tabel_terril_mapping$lokale_id] <- 1
  }
}

hooggroen_binair <- terra::ifel(hooggroen_aligned == 1, 1, NA)
stoorlaag <- terra::cover(hooggroen_binair, r_terril)

# --- 4. FILTEREN VAN DE BIOTOPEN (PARALLEL) ---
porseleinhoen_geschikt2_max <- terra::mask(biotoop_max, stoorlaag, inverse = TRUE)
porseleinhoen_geschikt2_opp <- terra::mask(biotoop_opp, stoorlaag, inverse = TRUE)

rm(biotoop_max, biotoop_opp_clean_gras, biotoop_opp_clean_moer, biotoop_opp_som, biotoop_opp_raw, biotoop_opp, 
   groenkaart, groenkaart_TV, hooggroen_binair, r_terril, stoorlaag, tabel_terril, tabel_terril_TV, tabel_terril_mapping)
gc()

message("=== STAP 2: Beperken tot natte gronden via Drainagekaart (PARALLEL) ===")

r_drain_raw   <- terra::rast(here("data/input/Raster_Vlaanderen/vlaanderen_drainage_10m.tif"))
r_drain_local <- r_drain_raw %>% 
                   terra::crop(area_buffer_fix, snap = "near") %>% 
                   terra::resample(template_TV, method = "near")

drain_cats            <- terra::cats(r_drain_local)[[1]]
geselecteerde_letters <- c("d", "e", "f", "h", "i", "g", "e-f", "h-i", "e-i")

porseleinhoen_drain_ids <- drain_cats$value[drain_cats$Label %in% geselecteerde_letters]

masker_drainage <- r_drain_local %in% porseleinhoen_drain_ids
mal_drainage    <- terra::ifel(masker_drainage == 1, 1, NA)

porseleinhoen_geschikt3_max <- terra::mask(porseleinhoen_geschikt2_max, mal_drainage)
porseleinhoen_geschikt3_opp <- terra::mask(porseleinhoen_geschikt2_opp, mal_drainage)

cat("   [Drainage Check] Oppervlakte na drainagefilter MAX:", calc_ha_exact(porseleinhoen_geschikt3_max), "ha\n")
cat("   [Drainage Check] Oppervlakte na drainagefilter OPP:", calc_ha_exact(porseleinhoen_geschikt3_opp), "ha\n")

rm(r_drain_raw, r_drain_local, masker_drainage, mal_drainage, drain_cats)
gc()

message("=== STAP 3: Eerste clustering en filtering op 3 hectare (PARALLEL) ===")

drempel_3ha_m2 <- 3 * 10000

# Spoor 1: Maximale Potentie (Binaire check)
cluster_3ha_max <- cluster_filter_compleet(
  masker     = porseleinhoen_geschikt3_max,
  opp_laag   = porseleinhoen_geschikt3_opp,
  drempel_m2 = drempel_3ha_m2,
  dist_m     = 50,
  werkelijk  = FALSE
)
porseleinhoen_geschikt4_max <- cluster_3ha_max$raster

# Spoor 2: Werkelijke Oppervlakte (Fractionele check - Parallel & Autonoom)
r_binair_geschikt3_opp <- terra::ifel(!is.na(porseleinhoen_geschikt3_opp) & porseleinhoen_geschikt3_opp > 0, 1, NA)
cluster_3ha_opp <- cluster_filter_compleet(
  masker     = r_binair_geschikt3_opp,
  opp_laag   = porseleinhoen_geschikt3_opp,
  drempel_m2 = drempel_3ha_m2,
  dist_m     = 50,
  werkelijk  = TRUE
)
porseleinhoen_geschikt4_opp <- cluster_3ha_opp$raster

rm(cluster_3ha_max, cluster_3ha_opp, r_binair_geschikt3_opp)
gc()

message("=== STAP 4: Open zicht filter & re-evaluatie (20 are bosdrempel - VECTOR SPEEDUP) ===")

# 1. Inlezen en croppen
groenkaart <- terra::rast(here("data/input/ASCI Files/Groenkaart_2021.tif"))
groenkaart_TV <- terra::crop(groenkaart, template_TV, snap = "near")
hooggroen_aligned <- terra::resample(groenkaart_TV, template_TV, method = "near")

# 2. Zet direct om naar polygonen (enkel celwaarden == 1)
# Dit vervangt de trage focal() + patches() stap
r_groen_bin <- terra::ifel(hooggroen_aligned == 1, 1, NA)
bos_polys <- terra::as.polygons(r_groen_bin, dissolve = TRUE) %>% 
  sf::st_as_sf()

porseleinhoen_geschikt5_max <- template_TV * NA
porseleinhoen_geschikt5_opp <- template_TV * NA

if (nrow(bos_polys) > 0) {
  # 3. Opsplitsen van multipolygonen naar losse polygonen & filter op >= 2000 m² (20 are)
  bos_polys_single <- sf::st_cast(bos_polys, "POLYGON")
  bos_polys_single$opp_m2 <- as.numeric(sf::st_area(bos_polys_single))
  
  bos_groot_sf <- bos_polys_single %>% 
    dplyr::filter(opp_m2 >= 2000)
  
  if (nrow(bos_groot_sf) > 0) {
    # 4. Maak een buffer van 100 meter rondom de grote boscomplexen
    bos_buffer_100m <- sf::st_buffer(bos_groot_sf, dist = 100) %>% 
      sf::st_union()
    
    # 5. Rasteriseer de buffer direct als maskeerafbeelding (1 binnen buffer, NA buiten)
    bos_buffer_vect <- terra::vect(bos_buffer_100m)
    masker_groen_100m_cluster <- terra::rasterize(bos_buffer_vect, template_TV, field = 1, background = NA)
    
    # Bewaar 'bos_groot' als SpatVector t.b.v. de latere Leaflet inspectiekaart
    bos_groot <- terra::vect(bos_groot_sf)
    
    # 6. Maskeer de geschikte biotopen (inverse = TRUE verwijdert alles binnen 100m van bos)
    biotoop_fase5_kaal_bin  <- terra::mask(porseleinhoen_geschikt4_max, masker_groen_100m_cluster, inverse = TRUE)
    biotoop_fase5_kaal_frac <- terra::mask(porseleinhoen_geschikt4_opp, masker_groen_100m_cluster, inverse = TRUE)
    
    ha_kaal_check_bin  <- calc_ha_exact(biotoop_fase5_kaal_bin)
    ha_kaal_check_frac <- calc_ha_exact(biotoop_fase5_kaal_frac)
    
    message("   [Diagnose] Ruwe open ruimte buiten de bosranden MAX: ", round(ha_kaal_check_bin, 2), " ha")
    message("   [Diagnose] Ruwe open ruimte buiten de bosranden OPP: ", round(ha_kaal_check_frac, 2), " ha")
    
    if (ha_kaal_check_bin > 0) {
      porseleinhoen_geschikt5_max <- cluster_filter_compleet(biotoop_fase5_kaal_bin, NULL, 30000, 50, FALSE)$raster
    }
    
    if (ha_kaal_check_frac > 0) {
      r_binair_kaal_frac <- terra::ifel(!is.na(biotoop_fase5_kaal_frac) & biotoop_fase5_kaal_frac > 0, 1, NA)
      porseleinhoen_geschikt5_opp <- cluster_filter_compleet(r_binair_kaal_frac, opp_laag = biotoop_fase5_kaal_frac, 30000, 50, TRUE)$raster
      rm(r_binair_kaal_frac)
    }
  } else {
    message("Opmerking: Geen storende boscomplexen van >= 20 are gevonden. Volledig open zicht!")
    porseleinhoen_geschikt5_max <- porseleinhoen_geschikt4_max
    porseleinhoen_geschikt5_opp <- porseleinhoen_geschikt4_opp
  }
} else {
  message("Opmerking: Geen bos op de groenkaart gevonden binnen de uitsnede.")
  porseleinhoen_geschikt5_max <- porseleinhoen_geschikt4_max
  porseleinhoen_geschikt5_opp <- porseleinhoen_geschikt4_opp
}

# Opruimen van tijdelijke objecten
rm(r_groen_bin, bos_polys, groenkaart, groenkaart_TV, hooggroen_aligned)
if (exists("bos_polys_single")) rm(bos_polys_single, bos_groot_sf, bos_buffer_100m, bos_buffer_vect, masker_groen_100m_cluster)
if (exists("biotoop_fase5_kaal_bin")) rm(biotoop_fase5_kaal_bin, biotoop_fase5_kaal_frac)
gc()

# ==============================================================================
# FINALE STAP 5: FINALE CLUSTERING & EINDDREMPEL (10 ha - PARALLEL)
# ==============================================================================
message("=== FINALE STAP: Filtering van het leefgebiedcomplex op 10 hectare ===")

drempel_10ha_m2 <- 10 * 10000

# Finaal Spoor 1: Maximale Potentie
finaal_10ha_max <- cluster_filter_compleet(
  masker     = porseleinhoen_geschikt5_max,
  opp_laag   = porseleinhoen_geschikt5_opp,
  drempel_m2 = drempel_10ha_m2,
  dist_m     = 50,
  werkelijk  = FALSE
)
finaal_leefgebied_max <- finaal_10ha_max$raster

# Finaal Spoor 2: Werkelijke Oppervlakte (Parallel & Autonoom)
r_binair_geschikt5_opp <- terra::ifel(!is.na(porseleinhoen_geschikt5_opp) & porseleinhoen_geschikt5_opp > 0, 1, NA)
finaal_10ha_opp <- cluster_filter_compleet(
  masker     = r_binair_geschikt5_opp,
  opp_laag   = porseleinhoen_geschikt5_opp,
  drempel_m2 = drempel_10ha_m2,
  dist_m     = 50,
  werkelijk  = TRUE
)
finaal_leefgebied_opp <- finaal_10ha_opp$raster

paapje_finaal_max <- finaal_leefgebied_max
paapje_finaal_opp <- finaal_leefgebied_opp

cat("\n==================================================\n")
cat("    GEVOELIG GEOCORRIGEERD HABITAT PORSELEINHOEN    \n")
cat("==================================================\n")
cat("Finaal Leefgebied Spoor 1 (Max Potentie):  ", calc_ha_exact(finaal_leefgebied_max), "ha\n")
cat("Finaal Leefgebied Spoor 2 (Werkelijke Opp): ", calc_ha_exact(finaal_leefgebied_opp), "ha\n")
cat("==================================================\n")

rm(finaal_10ha_max, finaal_10ha_opp, r_binair_geschikt5_opp)
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

