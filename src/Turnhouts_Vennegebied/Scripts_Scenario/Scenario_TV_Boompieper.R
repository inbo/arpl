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

# Dwing terra om heel behoudend te zijn met RAM (gebruik max 2 GB RAM, rest op schijf)
terra::terraOptions(memfrac = 0.6, memmin = 1)

cluster_filter_compleet <- function(masker, opp_laag, drempel_m2, dist_m, werkelijk = FALSE) {
  # 1. Controleer of het raster leeg is
  if (terra::global(is.na(masker), "sum")[[1]] == terra::ncell(masker)) {
    return(list(raster = masker * NA, clusters = masker * NA))
  }
  
  # Tijdelijk bestand aanmaken voor geheugensafety
  tmp_buf <- tempfile(pattern = "cl_buf_", fileext = ".tif")
  on.exit(unlink(tmp_buf), add = TRUE)
  
  # ----------------------------------------------------------------------------
  # STAP 1: VEILIGE NETWERKVORMING VIA SCHIJF
  # ----------------------------------------------------------------------------
  if (dist_m > 0) {
    r_binair <- terra::ifel(!is.na(masker) & masker > 0, 1, NA)
    
    # Schrijf de buffer naar schijf om RAM-crash te voorkomen
    r_buffered <- terra::buffer(
      r_binair, 
      width = dist_m, 
      filename = tmp_buf, 
      overwrite = TRUE, 
      gdal = c("COMPRESS=LZW")
    )
    
    cl_network <- terra::patches(r_buffered, directions = 4, zeroAsNA = TRUE)
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

soort <- "boompieper"

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

df <- read_excel(here::here("data/input/Excel_files/Soorten_bwk_afstanden.xlsx"))
resultaat <- df %>%
  filter(tolower(trimws(Soort)) == soort) %>%
  select(Type, MinOpp_ha, AfstandBiotopen_m, Dispersiecap_m)

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

message("-> Vertaalraster voor globale/lokale cellen bliksemsnel opbouwen...")

# 1. Maak de lokale uitsnede van het raster
id_raster_TV <- crop(master_grid, area_buffer_fix, snap = "near")

# 2. Rasterize de buffer-vector naar de uitsnede (vele malen sneller dan terra::extract!)
template_TV_mask <- terra::rasterize(area_buffer_fix, id_raster_TV, field = 1)

# 3. Wiskundige index-conversie (lokale cel-index -> coordinaten -> globale cel-index)
lokale_ids  <- terra::cells(template_TV_mask)
coords      <- terra::xyFromCell(template_TV_mask, lokale_ids)
globale_ids <- terra::cellFromXY(master_grid, coords)

# 4. Bouw de vertaaltabel (data.table) op
vertaal_df <- data.table(
  lokale_id  = lokale_ids,
  globale_id = globale_ids
)[!is.na(globale_id)]

studiegebied_globale_ids <- vertaal_df$globale_id

# 5. Templates & grenzen instellen voor vervolgstappen
template_TV <- terra::classify(template_TV_mask, cbind(NA, 0))
values(id_raster_TV) <- NA

grens_web <- sf::st_as_sf(terra::project(area_shape, "EPSG:4326"))

# Schoonmaken van tijdelijke variabelen
rm(template_TV_mask, lokale_ids, coords, globale_ids)
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

bwk_max        <- lijst_matches[["bwk"]]
bwk_opp        <- lijst_oppervlaktes[["bwk"]]
bwk_bossen_max <- lijst_matches[["bwk_bossen"]]
bwk_bossen_opp <- lijst_oppervlaktes[["bwk_bossen"]]
heide_max      <- lijst_matches[["heide"]]
heide_opp      <- lijst_oppervlaktes[["heide"]]
naaldbos_max   <- lijst_matches[["naaldbos"]]
naaldbos_opp   <- lijst_oppervlaktes[["naaldbos"]]

raster_simpel_final <- id_raster_TV

rm(tabel_vlaanderen, vertaal_df, lijst_matches, lijst_oppervlaktes)
gc()

# ==============================================================================
# STAP 1: ONAFHANKELIJK CLUSTEREN OP 5 HA (50M AFSTAND) - REGENWULP LOGICA
# ==============================================================================
message("-> Eerst onafhankelijk clusteren van open biotopen en bos op 5 ha...")

# --- 1. A: Binaire (MAX) lagen clusteren op 5 ha / 50m ---
open_cluster_5ha_max <- cluster_filter_compleet(
  masker     = bwk_max, 
  opp_laag   = bwk_opp,
  drempel_m2 = 5 * 10000,
  dist_m     = 50, 
  werkelijk  = FALSE
)$raster
gc()

bos_cluster_5ha_max <- cluster_filter_compleet(
  masker     = bwk_bossen_max, 
  opp_laag   = bwk_bossen_opp,
  drempel_m2 = 5 * 10000,
  dist_m     = 50, 
  werkelijk  = FALSE
)$raster
gc()

# --- 1. B: Oppervlaktefractie (OPP) lagen (Spoor B) ---
# Net als bij regenwulp: maak eerst een binaire maskerkaart voor de buffer!
r_binair_open_opp <- terra::ifel(!is.na(bwk_opp) & bwk_opp > 0, 1, NA)
open_cluster_5ha_opp <- cluster_filter_compleet(
  masker     = r_binair_open_opp,
  opp_laag   = bwk_opp, 
  drempel_m2 = 5 * 10000,
  dist_m     = 50, 
  werkelijk  = TRUE
)$raster
rm(r_binair_open_opp)
gc()

r_binair_bos_opp <- terra::ifel(!is.na(bwk_bossen_opp) & bwk_bossen_opp > 0, 1, NA)
bos_cluster_5ha_opp <- cluster_filter_compleet(
  masker     = r_binair_bos_opp,
  opp_laag   = bwk_bossen_opp, 
  drempel_m2 = 5 * 10000,
  dist_m     = 50, 
  werkelijk  = TRUE
)$raster
rm(r_binair_bos_opp)
gc()

# Schrijf de zware 5ha basisclusters direct weg naar schijf (RAM leegmaken!)
writeRaster(open_cluster_5ha_max, "temp_open_max.tif", overwrite = TRUE)
writeRaster(bos_cluster_5ha_max,  "temp_bos_max.tif",  overwrite = TRUE)
writeRaster(open_cluster_5ha_opp, "temp_open_opp.tif", overwrite = TRUE)
writeRaster(bos_cluster_5ha_opp, "temp_bos_opp.tif", overwrite = TRUE)

rm(open_cluster_5ha_max, bos_cluster_5ha_max, open_cluster_5ha_opp, bos_cluster_5ha_opp)
gc()

# Herlaad als superlichte bestandspoointers
open_cluster_5ha_max <- rast("temp_open_max.tif")
bos_cluster_5ha_max  <- rast("temp_bos_max.tif")
open_cluster_5ha_opp <- rast("temp_open_opp.tif")
bos_cluster_5ha_opp  <- rast("temp_bos_opp.tif")

# ==============================================================================
# STAP 2: 100M BOSRAND UITSNIJDEN VAN DE 5 HA BOSCLUSTERS (BLIKSEMSNEL EN STABIEL)
# ==============================================================================
message("-> 100m bosranden uitsnijden van de grote bosclusters via afstandskaart...")

# 1. Bepaal de afstand (in meters) van elke bos-cel tot de buitenrand/niet-bos
dist_tot_rand <- terra::distance(!is.na(bos_cluster_5ha_max))

# 2. De 100m bosrand zijn alle bos-pixels die binnen 100 meter van de rand liggen
bos_100m_5ha_max <- terra::mask(bos_cluster_5ha_max, dist_tot_rand <= 100, maskvalue = FALSE)

# 3. Koppel direct aan de OPP-laag
bos_100m_5ha_opp <- terra::mask(bos_cluster_5ha_opp, bos_100m_5ha_max)

rm(dist_tot_rand)
gc()

# ==============================================================================
# STAP 3: CONTACTZONE (20M) & EDGEN TUSSEN OPEN BIOTOOP EN BOSRAND
# ==============================================================================
message("-> Edgen: Alleen bosranden behouden die grenzen aan open biotoopclusters...")

# 20m komt overeen met een 3x3 matrix (1 cel in elke richting rondom)
w_edge_20m <- matrix(1, nrow = 3, ncol = 3)

# GECORRIGEERD: tempfile() in plaats van tmp_file()
f_buf_open_tif <- tempfile(pattern = "buf_open_", fileext = ".tif")
buf_open_max <- terra::focal(
  terra::ifel(!is.na(open_cluster_5ha_max), 1, NA), 
  w = w_edge_20m, fun = "max", na.rm = TRUE,
  filename = f_buf_open_tif, overwrite = TRUE, gdal = c("COMPRESS=LZW")
)

# GECORRIGEERD: tempfile() in plaats van tmp_file()
f_buf_bos_tif <- tempfile(pattern = "buf_bos_", fileext = ".tif")
buf_bos_max <- terra::focal(
  terra::ifel(!is.na(bos_100m_5ha_max), 1, NA), 
  w = w_edge_20m, fun = "max", na.rm = TRUE,
  filename = f_buf_bos_tif, overwrite = TRUE, gdal = c("COMPRESS=LZW")
)

contactzone_groot_max <- terra::ifel(!is.na(buf_open_max) & !is.na(buf_bos_max), 1, NA)

bosrand_gekoppeld_max <- terra::mask(bos_100m_5ha_max, contactzone_groot_max)
bosrand_gekoppeld_opp <- terra::mask(bos_100m_5ha_opp, contactzone_groot_max)

bwk_cluster_max <- terra::cover(open_cluster_5ha_max, bosrand_gekoppeld_max)
bwk_cluster_opp <- terra::cover(open_cluster_5ha_opp, bosrand_gekoppeld_opp)

# ==============================================================================
# STAP 4: RANDZONES BEPALEN VOOR KLEIN LEEFGEBIED (HEIDE + NAALDBOS)
# ==============================================================================
message("-> Randzones (edges) bepalen tussen heide en naaldbos voor klein leefgebied...")

if(!is.null(naaldbos_max) && !is.null(heide_max) && !all(is.na(terra::minmax(naaldbos_max)))) {
  # GECORRIGEERD: tempfile() in plaats van tmp_file()
  f_buf_naald_tif <- tempfile(pattern = "buf_naald_", fileext = ".tif")
  buf_naaldbos <- terra::focal(
    terra::ifel(naaldbos_max > 0, 1, NA), 
    w = w_edge_20m, fun = "max", na.rm = TRUE,
    filename = f_buf_naald_tif, overwrite = TRUE, gdal = c("COMPRESS=LZW")
  )
  
  # GECORRIGEERD: tempfile() in plaats van tmp_file()
  f_buf_heide_tif <- tempfile(pattern = "buf_heide_", fileext = ".tif")
  buf_heide <- terra::focal(
    terra::ifel(heide_max > 0, 1, NA), 
    w = w_edge_20m, fun = "max", na.rm = TRUE,
    filename = f_buf_heide_tif, overwrite = TRUE, gdal = c("COMPRESS=LZW")
  )
  
  contactzone_klein <- terra::ifel(!is.na(buf_naaldbos) & !is.na(buf_heide), 1, NA)
  
  heide_naaldbos_max <- terra::mask(heide_max, contactzone_klein)
  heide_naaldbos_opp <- terra::mask(heide_opp, contactzone_klein)
} else {
  heide_naaldbos_max <- id_raster_TV * NA
  heide_naaldbos_opp <- id_raster_TV * NA
}

gc()

# ==============================================================================
# STAP 5: CLUSTEREN EN FILTEREN OP 2 HECTARE (KLEIN LEEFGEBIED)
# ==============================================================================
message("-> Clusteren en filteren van heide-bosranden op minimale oppervlakte van 2 ha...")

bwk_cluster_klein_max_res <- cluster_filter_compleet(heide_naaldbos_max, heide_naaldbos_opp, dist_m = 50, drempel_m2 = 2 * 10000, werkelijk = FALSE)
gc()

bwk_cluster_klein_opp_res <- cluster_filter_compleet(heide_naaldbos_max, heide_naaldbos_opp, dist_m = 50, drempel_m2 = 2 * 10000, werkelijk = TRUE)
gc()

bwk_cluster_klein_max <- bwk_cluster_klein_max_res$raster
bwk_cluster_klein_opp <- bwk_cluster_klein_opp_res$raster

# Alleen de specifieke focal-bestanden verwijderen, NIET de algemene terra tmpFiles!
tijdelijke_bestanden <- c(
  if(exists("f_buf_open_tif")) f_buf_open_tif,
  if(exists("f_buf_bos_tif")) f_buf_bos_tif,
  if(exists("f_buf_naald_tif")) f_buf_naald_tif,
  if(exists("f_buf_heide_tif")) f_buf_heide_tif
)
unlink(tijdelijke_bestanden)

rm(
  heide_naaldbos_max, heide_naaldbos_opp, open_cluster_5ha_max, bos_cluster_5ha_max, 
  bos_100m_5ha_max, bosrand_gekoppeld_max, buf_open_max, buf_bos_max, tijdelijke_bestanden
)
gc()

# ==============================================================================
# STAP 6: COMBINEREN (GROOT OR KLEIN) - MET TYPE-ID'S VOOR KAART
# ==============================================================================
message("-> Leefgebieden combineren (Groot OF Klein)...")

# --- Combinatie voor MAX (Binair/Cluster IDs) ---
r_groot_binaire_status <- terra::ifel(!is.na(bwk_cluster_max), 1, NA)
r_klein_binaire_status <- terra::ifel(!is.na(bwk_cluster_klein_max), 2, NA)

leefgebied_finaal_max_binair <- terra::cover(r_groot_binaire_status, r_klein_binaire_status)
leefgebied_finaal_max        <- terra::cover(bwk_cluster_max, bwk_cluster_klein_max)

# --- Combinatie voor OPP (Parallel & Autonoom) ---
r_groot_opp_status <- terra::ifel(!is.na(bwk_cluster_opp) & bwk_cluster_opp > 0, bwk_cluster_opp, NA)
r_klein_opp_status <- terra::ifel(!is.na(bwk_cluster_klein_opp) & bwk_cluster_klein_opp > 0, bwk_cluster_klein_opp, NA)

leefgebied_finaal_opp        <- terra::cover(r_groot_opp_status, r_klein_opp_status)
leefgebied_finaal_opp_binair <- terra::ifel(!is.na(leefgebied_finaal_opp) & leefgebied_finaal_opp > 0, 1, NA)

# Pas NU opruimen als de combinatie klaar is!
rm(bwk_cluster_max, bwk_cluster_klein_max, bwk_cluster_opp, bwk_cluster_klein_opp, 
   r_groot_binaire_status, r_klein_binaire_status, r_groot_opp_status, r_klein_opp_status)
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

