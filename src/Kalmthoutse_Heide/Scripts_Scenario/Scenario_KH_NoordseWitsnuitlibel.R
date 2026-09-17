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
soort <- "noordsewitsnuitlibel"

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
buffer_m <- resultaat$Dispersiecap_m[1]
straal_water_m <- 500  # Harde afsnij-afstand rond het water

print(resultaat)
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
values(globale_id_raster) <- 1:ncell(globale_id_raster)

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

cat("Gecorrigeerd aantal pixels in template_KH: ", sum(terra::values(template_KH) == 1, na.rm=TRUE), "\n")

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

waterbiotoop_bwk1_max        <- lijst_matches[["waterbiotoop_bwk1"]]
waterbiotoop_bwk1_opp        <- lijst_oppervlaktes[["waterbiotoop_bwk1"]]
waterbiotoop_notbwk_max      <- lijst_matches[["waterbiotoop_notbwk"]]
waterbiotoop_notbwk_opp      <- lijst_oppervlaktes[["waterbiotoop_notbwk"]]
landbiotoop_bwk1_max         <- lijst_matches[["landbiotoop_bwk"]]
landbiotoop_bwk1_opp         <- lijst_oppervlaktes[["landbiotoop_bwk"]]
landbiotoop_bwk_ruim_max     <- lijst_matches[["landbiotoop_bwk_ruim"]]
landbiotoop_bwk_ruim_opp     <- lijst_oppervlaktes[["landbiotoop_bwk_ruim"]]
landbiotoop_bwk_ruim_bos_max <- lijst_matches[["landbiotoop_bwk_ruim_bos"]]
landbiotoop_bwk_ruim_bos_opp <- lijst_oppervlaktes[["landbiotoop_bwk_ruim_bos"]]

rm(tabel_vlaanderen, vertaal_df, lijst_matches, lijst_oppervlaktes)
gc()

message("-> Waterbiotoop berekenen (Spoor Max & Opp)...")

waterbiotoop_bwk_basis_max <- terra::ifel(!is.na(waterbiotoop_bwk1_max) & is.na(waterbiotoop_notbwk_max), 1, NA)
waterbiotoop_bwk_basis_opp <- terra::ifel(!is.na(waterbiotoop_bwk1_opp) & waterbiotoop_bwk1_opp > 0 & (is.na(waterbiotoop_notbwk_opp) | waterbiotoop_notbwk_opp == 0), 1, NA)

# --- SPOOR 1: MAXIMALE POTENTIE ---
waterbiotoop1_max_list <- cluster_filter_compleet(
  masker     = waterbiotoop_bwk_basis_max, 
  opp_laag   = waterbiotoop_bwk1_max, 
  drempel_m2 = 200,                  # 0.02 ha = 200 m²
  dist_m     = 50, 
  werkelijk  = FALSE
)
waterbiotoop1_max <- waterbiotoop1_max_list$raster

# --- SPOOR 2: WERKELIJKE OPPERVLAKTE (PARALLEL & AUTONOOM) ---
r_binair_water_opp <- terra::ifel(!is.na(waterbiotoop_bwk_basis_opp) & waterbiotoop_bwk_basis_opp > 0, 1, NA)
waterbiotoop1_opp_list <- cluster_filter_compleet(
  masker     = r_binair_water_opp, 
  opp_laag   = waterbiotoop_bwk1_opp, 
  drempel_m2 = 200, 
  dist_m     = 50, 
  werkelijk  = TRUE
)
waterbiotoop1_opp <- waterbiotoop1_opp_list$raster

cl_id_max <- waterbiotoop1_max_list$clusters
cl_id_opp <- waterbiotoop1_opp_list$clusters

rm(r_binair_water_opp, waterbiotoop_bwk_basis_max, waterbiotoop_bwk_basis_opp)
gc()

message("-> Landbiotoop nabije omgeving berekenen (Spoor Max & Opp)...")

landbiotoop_bwk_basis_max <- terra::ifel(!is.na(landbiotoop_bwk1_max) & is.na(waterbiotoop_notbwk_max), 1, NA)
landbiotoop_bwk_basis_opp <- terra::ifel(!is.na(landbiotoop_bwk1_opp) & landbiotoop_bwk1_opp > 0 & (is.na(waterbiotoop_notbwk_opp) | waterbiotoop_notbwk_opp == 0), 1, NA)

# --- SPOOR 1: MAXIMALE POTENTIE ---
landbiotoop1_max_list <- cluster_filter_compleet(
  masker     = landbiotoop_bwk_basis_max, 
  opp_laag   = landbiotoop_bwk1_max, 
  drempel_m2 = 50000, 
  dist_m     = 20, 
  werkelijk  = FALSE
)
landbiotoop_nabij_max <- landbiotoop1_max_list$raster

# --- SPOOR 2: WERKELIJKE OPPERVLAKTE (Parallel & Autonoom) ---
r_binair_land_opp <- terra::ifel(!is.na(landbiotoop_bwk_basis_opp) & landbiotoop_bwk_basis_opp > 0, 1, NA)
landbiotoop1_opp_list <- cluster_filter_compleet(
  masker     = r_binair_land_opp, 
  opp_laag   = landbiotoop_bwk1_opp, 
  drempel_m2 = 50000, 
  dist_m     = 20, 
  werkelijk  = TRUE
)
landbiotoop_nabij_opp <- landbiotoop1_opp_list$raster

rm(landbiotoop_bwk_basis_max, landbiotoop_bwk_basis_opp, r_binair_land_opp, landbiotoop1_max_list, landbiotoop1_opp_list)
gc()

message("-> Landbiotoop ruimere omgeving & bos-samenstellingscheck berekenen (Minimaal 40% bos - PARALLEL)...")

# ==============================================================================
# --- SPOOR 1: MAXIMALE POTENTIE ---
# ==============================================================================
landbiotoop_ruim_max_list <- cluster_filter_compleet(
  masker     = landbiotoop_bwk_ruim_max, 
  opp_laag   = landbiotoop_bwk_ruim_max, 
  drempel_m2 = 400000, 
  dist_m     = 50, 
  werkelijk  = FALSE
)
r_ruim_max_src <- landbiotoop_ruim_max_list$raster
r_ruim_max_src[is.na(r_ruim_max_src) | r_ruim_max_src == 0] <- NA

r_bos_max_src <- landbiotoop_bwk_ruim_bos_max
r_bos_max_src[is.na(r_bos_max_src) | r_bos_max_src == 0] <- NA

cl_id_ruim_max <- terra::patches(r_ruim_max_src, directions = 8, zeroAsNA = TRUE)

stats_ruim_totaal_max <- terra::zonal(r_ruim_max_src * 0.01, cl_id_ruim_max, fun = "sum", na.rm = TRUE)
colnames(stats_ruim_totaal_max) <- c("Cluster_ID", "Totaal_Ruim_ha")

bos_in_ruim_max <- terra::mask(r_bos_max_src, cl_id_ruim_max)
stats_bos_in_ruim_max <- terra::zonal(bos_in_ruim_max * 0.01, cl_id_ruim_max, fun = "sum", na.rm = TRUE)
colnames(stats_bos_in_ruim_max) <- c("Cluster_ID", "Bos_ha")

df_samenstelling_max <- merge(stats_ruim_totaal_max, stats_bos_in_ruim_max, by = "Cluster_ID", all.x = TRUE)
df_samenstelling_max$Bos_ha[is.na(df_samenstelling_max$Bos_ha)] <- 0
df_samenstelling_max$Percentage_Bos <- (df_samenstelling_max$Bos_ha / df_samenstelling_max$Totaal_Ruim_ha) * 100

valide_ruim_ids_max <- df_samenstelling_max$Cluster_ID[df_samenstelling_max$Percentage_Bos >= 40]

if(length(valide_ruim_ids_max) > 0) {
  masker_ruim_max <- cl_id_ruim_max %in% valide_ruim_ids_max
  landbiotoop_ruim_gekoppeld_max <- terra::mask(r_ruim_max_src, terra::ifel(masker_ruim_max, 1, NA))
} else {
  message("⚠️ Waarschuwing: Geen enkel ruim landbiotoopcluster voldoet aan de 40% bos-eis in Spoor A.")
  landbiotoop_ruim_gekoppeld_max <- template_KH * NA
}

# ==============================================================================
# --- SPOOR 2: WERKELIJKE OPPERVLAKTE (AUTONOOM EN PARALLEL) ---
# ==============================================================================
r_binair_ruim_opp_init <- terra::ifel(!is.na(landbiotoop_bwk_ruim_opp) & landbiotoop_bwk_ruim_opp > 0, 1, NA)

landbiotoop_ruim_opp_list <- cluster_filter_compleet(
  masker     = r_binair_ruim_opp_init, 
  opp_laag   = landbiotoop_bwk_ruim_opp, 
  drempel_m2 = 400000, 
  dist_m     = 50, 
  werkelijk  = TRUE
)
r_ruim_opp_src <- landbiotoop_ruim_opp_list$raster
r_ruim_opp_bin <- terra::ifel(!is.na(r_ruim_opp_src) & r_ruim_opp_src > 0, 1, NA)

if(!all(is.na(terra::values(r_ruim_opp_bin, mat=FALSE)))) {
  cl_id_ruim_opp <- terra::patches(r_ruim_opp_bin, directions = 8, zeroAsNA = TRUE)
  
  stats_ruim_totaal_opp <- terra::zonal(r_ruim_opp_src, cl_id_ruim_opp, fun = "sum", na.rm = TRUE)
  colnames(stats_ruim_totaal_opp) <- c("Cluster_ID", "Totaal_Ruim_ha")
  stats_ruim_totaal_opp$Totaal_Ruim_ha <- stats_ruim_totaal_opp$Totaal_Ruim_ha * 0.01
  
  bos_in_ruim_opp <- terra::mask(landbiotoop_bwk_ruim_bos_opp, cl_id_ruim_opp)
  stats_bos_in_ruim_opp <- terra::zonal(bos_in_ruim_opp, cl_id_ruim_opp, fun = "sum", na.rm = TRUE)
  colnames(stats_bos_in_ruim_opp) <- c("Cluster_ID", "Bos_ha")
  stats_bos_in_ruim_opp$Bos_ha <- stats_bos_in_ruim_opp$Bos_ha * 0.01
  
  df_samenstelling_opp <- merge(stats_ruim_totaal_opp, stats_bos_in_ruim_opp, by = "Cluster_ID", all.x = TRUE)
  df_samenstelling_opp$Bos_ha[is.na(df_samenstelling_opp$Bos_ha)] <- 0
  df_samenstelling_opp$Percentage_Bos <- (df_samenstelling_opp$Bos_ha / df_samenstelling_opp$Totaal_Ruim_ha) * 100
  
  valide_ruim_ids_opp <- df_samenstelling_opp$Cluster_ID[df_samenstelling_opp$Percentage_Bos >= 40]
  
  if(length(valide_ruim_ids_opp) > 0) {
    masker_ruim_opp <- cl_id_ruim_opp %in% valide_ruim_ids_opp
    landbiotoop_ruim_gekoppeld_opp <- terra::mask(r_ruim_opp_src, terra::ifel(masker_ruim_opp, 1, NA))
  } else {
    landbiotoop_ruim_gekoppeld_opp <- template_KH * NA
  }
} else {
  landbiotoop_ruim_gekoppeld_opp <- template_KH * NA
}

message("   - Aantal geëvalueerde ruime landclusters (MAX): ", nrow(df_samenstelling_max))
message("   - Aantal clusters dat voldoet (>= 40% bos MAX): ", length(valide_ruim_ids_max))

rm(landbiotoop_ruim_max_list, r_ruim_max_src, r_bos_max_src, cl_id_ruim_max,
   stats_ruim_totaal_max, bos_in_ruim_max, stats_bos_in_ruim_max, df_samenstelling_max,
   r_binair_ruim_opp_init, landbiotoop_ruim_opp_list, r_ruim_opp_src, r_ruim_opp_bin, cl_id_ruim_opp,
   stats_ruim_totaal_opp, bos_in_ruim_opp, stats_bos_in_ruim_opp, df_samenstelling_opp)
gc()

message("-> Originele ruimtelijke interacties uitvoeren + Water-First 500m afsnijding...")

straal_water_m <- 500

# ==============================================================================
# SPOOR 1: MAXIMALE POTENTIE
# ==============================================================================
water_buffer10_max <- terra::buffer(waterbiotoop1_max, width = 10)
land_buffer10_max  <- terra::buffer(landbiotoop_nabij_max, width = 10)

waterbiotoop2_max <- terra::mask(waterbiotoop1_max, land_buffer10_max)
landbiotoop2_max  <- terra::mask(landbiotoop_nabij_max, water_buffer10_max)

ruim_buffer500_max <- terra::buffer(landbiotoop_ruim_gekoppeld_max, width = 500)

waterbiotoop_finaal_max_ruw <- terra::mask(waterbiotoop2_max, ruim_buffer500_max)
landbiotoop_finaal_max_ruw  <- terra::mask(landbiotoop2_max, ruim_buffer500_max)

# HARDE WATER-FIRST AFSNIJDING OP 500M (EXACT ZOALS MAANWATERJUFFER)
if (!all(is.na(terra::values(waterbiotoop_finaal_max_ruw, mat = FALSE)))) {
  poly_water_max  <- terra::as.polygons(waterbiotoop_finaal_max_ruw, aggregate = TRUE)
  poly_buffer_max <- terra::buffer(poly_water_max, width = straal_water_m)
  landbiotoop_finaal_max <- terra::mask(landbiotoop_finaal_max_ruw, poly_buffer_max)
} else {
  landbiotoop_finaal_max <- landbiotoop_finaal_max_ruw
}
waterbiotoop_finaal_max <- waterbiotoop_finaal_max_ruw

leefgebied_max    <- terra::cover(waterbiotoop_finaal_max, landbiotoop_finaal_max)
finaal_max_binair <- terra::ifel(!is.na(leefgebied_max) & leefgebied_max > 0, 1, NA)
cl_id_max         <- terra::patches(finaal_max_binair, directions = 8, zeroAsNA = TRUE)


# ==============================================================================
# SPOOR 2: WERKELIJKE OPPERVLAKTE (PARALLEL)
# ==============================================================================
water_buffer10_opp <- terra::buffer(!is.na(waterbiotoop1_opp) & waterbiotoop1_opp > 0, width = 10)
land_buffer10_opp  <- terra::buffer(!is.na(landbiotoop_nabij_opp) & landbiotoop_nabij_opp > 0, width = 10)

waterbiotoop2_opp <- terra::mask(waterbiotoop1_opp, land_buffer10_opp)
landbiotoop2_opp  <- terra::mask(landbiotoop_nabij_opp, water_buffer10_opp)

ruim_buffer500_opp <- terra::buffer(!is.na(landbiotoop_ruim_gekoppeld_opp) & landbiotoop_ruim_gekoppeld_opp > 0, width = 500)

waterbiotoop_finaal_opp_ruw <- terra::mask(waterbiotoop2_opp, ruim_buffer500_opp)
landbiotoop_finaal_opp_ruw  <- terra::mask(landbiotoop2_opp, ruim_buffer500_opp)

# HARDE WATER-FIRST AFSNIJDING OP 500M (EXACT ZOALS MAANWATERJUFFER)
if (!all(is.na(terra::values(waterbiotoop_finaal_opp_ruw, mat = FALSE)))) {
  r_bin_water_opp <- terra::ifel(!is.na(waterbiotoop_finaal_opp_ruw) & waterbiotoop_finaal_opp_ruw > 0, 1, NA)
  poly_water_opp  <- terra::as.polygons(r_bin_water_opp, aggregate = TRUE)
  poly_buffer_opp <- terra::buffer(poly_water_opp, width = straal_water_m)
  landbiotoop_finaal_opp <- terra::mask(landbiotoop_finaal_opp_ruw, poly_buffer_opp)
} else {
  landbiotoop_finaal_opp <- landbiotoop_finaal_opp_ruw
}
waterbiotoop_finaal_opp <- waterbiotoop_finaal_opp_ruw

leefgebied_opp    <- terra::cover(waterbiotoop_finaal_opp, landbiotoop_finaal_opp)
finaal_opp_binair <- terra::ifel(!is.na(leefgebied_opp) & leefgebied_opp > 0, 1, NA)
cl_id_opp         <- terra::patches(finaal_opp_binair, directions = 8, zeroAsNA = TRUE)

rm(water_buffer10_max, land_buffer10_max, waterbiotoop2_max, landbiotoop2_max, ruim_buffer500_max,
   water_buffer10_opp, land_buffer10_opp, waterbiotoop2_opp, landbiotoop2_opp, ruim_buffer500_opp)
gc()

message("-> Modellering voltooid! Landbiotoop exact afgesneden binnen 500m van goedgekeurd water.")


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

