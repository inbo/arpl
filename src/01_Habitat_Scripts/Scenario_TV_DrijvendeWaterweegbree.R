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
soort <- "drijvendewaterweegbree"

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
buffer_masker <- terra::rasterize(area_buffer_fix, id_raster_TV, field = 1)

lokale_cellen <- terra::cells(buffer_masker)
coordinaten <- terra::xyFromCell(id_raster_TV, lokale_cellen)
studiegebied_globale_ids <- terra::cellFromXY(master_grid, coordinaten)

vertaal_df <- data.table(
  lokale_id = lokale_cellen,
  globale_id = studiegebied_globale_ids
)

template_TV <- terra::rasterize(area_buffer_fix, id_raster_TV, field = 1, background = 0)
grens_web <- sf::st_as_sf(terra::project(area_shape, "EPSG:4326"))

rm(buffer_masker, lokale_cellen, coordinaten)
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

r_bwk_max <- lijst_matches[["voortplanting_bwk"]]
r_bwk_opp <- lijst_oppervlaktes[["voortplanting_bwk"]]

rm(tabel_vlaanderen, vertaal_df, lijst_matches, lijst_oppervlaktes)
gc()

# ==============================================================================
# 7. HYDROGRAFISCHE DATA (HUET-ZONES) INLADEN EN SYNCHRONISEREN (IF-CHECK FIX)
# ==============================================================================
message("-> Huet-zone waterlopen inlezen...")
r_huetzon_vlaanderen <- rast(here("data/input/Raster_Vlaanderen/vlaanderen_huetzon_10m.tif"))

r_huetzon_local <- r_huetzon_vlaanderen %>% 
  terra::crop(area_buffer_fix) %>% 
  terra::resample(template_TV, method = "near")

# Haal de legende-ID's op (dit print nu netjes jouw tabel van 38 rijen)
legende_huet  <- terra::cats(r_huetzon_local)[[1]]
target_labels <- c("BkK1", "BkK2", "KabPKunst", "Kan", "KanPKunst", "Kunst")

# GEFIXT: We checken of het een data.frame is in plaats van is.na() te gebruiken
if(!is.null(legende_huet) && is.data.frame(legende_huet) && "Label" %in% colnames(legende_huet)) {
  # Dynamische match: welke ID's horen bij jouw 6 doelklassen?
  valide_ids <- legende_huet$value[legende_huet$Label %in% target_labels]
  message(paste("-> Gekoppelde Huet-zone ID's op basis van labels:", paste(valide_ids, collapse = ", ")))
} else {
  # Failsafe terugval op basis van jouw handmatige controle (5:10)
  valide_ids <- 5:10 
  message("⚠️ Waarschuwing: Rasterlegende-structuur wijk af, terugvallen op harde ID's 5:10.")
}

# Genereer het binaire waterlopen-raster (1/NA)
r_waterlopen_bin <- r_huetzon_local %in% valide_ids
r_waterlopen_bin <- terra::ifel(r_waterlopen_bin == 1, 1, NA)

# ==============================================================================
# 8. PARALLELLE FUSIE: SPOOR A (MAX) EN SPOOR B (OPP) STRIKT GESCHEIDEN
# ==============================================================================
message("-> Strikt gescheiden parallelle fusie uitvoeren...")

# --- SPOOR A: MAXIMALE POTENTIE (Binair spoor) ---
# Een pixel is geschikt als hij in de BWK óf in de waterlopen zit
drijvendewaterweegbree_voortplanting_max <- terra::cover(r_bwk_max, r_waterlopen_bin)
names(drijvendewaterweegbree_voortplanting_max) <- "Match_Max"

# --- SPOOR B: WERKELIJKE OPPERVLAKTE (Fractie spoor) ---
# De waterlopen krijgen in dit reële spoor de maximale fractie (0.01 ha per pixel)
r_waterlopen_opp <- terra::ifel(r_waterlopen_bin == 1, 0.01, NA)
drijvendewaterweegbree_voortplanting_opp <- terra::cover(r_bwk_opp, r_waterlopen_opp)
names(drijvendewaterweegbree_voortplanting_opp) <- "Oppervlakte_Real"

# ==============================================================================
# BACKWARDS COMPATIBILITY FOR FOLLOWER CHUNKS
# ==============================================================================
# We overschrijven de centrale namen zodat je kaart/tabel chunks naadloos doorlopen
alle_matches      <- drijvendewaterweegbree_voortplanting_max
alle_oppervlaktes <- drijvendewaterweegbree_voortplanting_opp

# Schoonmaak van alle tijdelijke objecten om RAM vrij te maken
rm(tabel_vlaanderen, tabel_gefilterd, tabel_TV, tabel_cel_som, vertaal_df,
   r_bwk_max, r_bwk_opp, r_huetzon_vlaanderen, r_huetzon_local, r_waterlopen_bin, r_waterlopen_opp)
gc()

# ==============================================================================
# STAP 11: RUIMTELIJKE CLUSTERING & OPPERVLAKTEFILTER VIA UNIVERSELE FUNCTIE
# ==============================================================================
message("-> Start ruimtelijke clustering (20m overbrugging & 0.16 ha drempel)...")

drempel_plant_m2 <- oppervlakte_ha * 10000 # 0.16 ha = 1600 m²

# --- SPOOR A: MAXIMALE POTENTIE (Binair) ---
res_plant_max <- cluster_filter_compleet(
  masker     = alle_matches,
  opp_laag   = alle_matches,
  drempel_m2 = drempel_plant_m2,
  dist_m     = 20, # 20m overbruggingsafstand
  werkelijk  = FALSE
)

drijvendewaterweegbree_leefgebied_max <- res_plant_max$raster
names(drijvendewaterweegbree_leefgebied_max) <- "Leefgebied_Max"


# --- SPOOR B: WERKELIJKE OPPERVLAKTE (Fracties) ---
r_binair_opp <- terra::ifel(!is.na(alle_oppervlaktes) & alle_oppervlaktes > 0, 1, NA)

res_plant_opp <- cluster_filter_compleet(
  masker     = r_binair_opp,
  opp_laag   = alle_oppervlaktes,
  drempel_m2 = drempel_plant_m2,
  dist_m     = 20,
  werkelijk  = TRUE
)

drijvendewaterweegbree_leefgebied_opp <- res_plant_opp$raster
names(drijvendewaterweegbree_leefgebied_opp) <- "Leefgebied_Opp"

# ==============================================================================
# DIAGNOSE EN OUTPUT PRINTS (GESYNCHRONISEERD MET DE FINALE TABEL)
# ==============================================================================
cat("\n--- AUDIT: EFFECT VAN DE MINIMALE CLUSTERDREMPEL (", oppervlakte_ha, "ha) ---\n")

# Spoor A berekent pure cellen (elke cel is 0.01 ha)
ha_max_voor <- calc_ha_exact(alle_matches)
ha_max_na   <- calc_ha_exact(drijvendewaterweegbree_leefgebied_max)
cat("Spoor A (Max Potentie) VOOR filter:", round(ha_max_voor, 2), "ha\n")
cat("Spoor A (Max Potentie) NA filter:  ", round(ha_max_na, 2), "ha\n\n")

# Spoor B berekent de som van de reële fracties uit het cluster-resultaat
ha_opp_voor <- calc_ha_exact(alle_oppervlaktes)
ha_opp_na   <- calc_ha_exact(drijvendewaterweegbree_leefgebied_opp)

cat("Spoor B (Werkelijke Opp) VOOR filter:", round(ha_opp_voor, 2), "ha\n")
cat("Spoor B (Werkelijke Opp) NA filter:  ", round(ha_opp_na, 2), "ha\n")

# Opruimen van tijdelijke objecten om RAM vrij te maken
rm(res_plant_max, res_plant_opp)
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

