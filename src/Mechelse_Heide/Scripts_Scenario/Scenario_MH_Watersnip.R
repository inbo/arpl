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
soort <- "watersnip"

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

resultaat <- df %>%
  filter(tolower(trimws(Soort)) == soort) %>%
  select(Type, MinOpp_ha, AfstandBiotopen_m, Dispersiecap_m)

# Variabelen definiëren
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

cat("Gecorrigeerd aantal pixels in template_MH: ", sum(terra::values(template_MH) == 1, na.rm=TRUE), "\n")

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

biotopen_max <- lijst_matches[["biotopen"]]
biotopen_opp <- lijst_oppervlaktes[["biotopen"]]

rm(vertaal_df, lijst_matches, lijst_oppervlaktes)
gc()

# ==============================================================================
# STAP 1: OMGEVINGSFILTERS & PRIMAIRE CLUSTERING (PARALLEL)
# ==============================================================================
message("-> Drainage en groen/waterkaarten inlezen en synchroniseren...")

r_drainage_raw <- rast(here("data/input/Raster_Vlaanderen/vlaanderen_drainage_10m.tif"))
r_groen_raw    <- rast(here("data/input/ASCI Files/Groenkaart_2021.tif"))
r_grb_raw      <- rast(here("data/input/Raster_Vlaanderen/vlaanderen_grb_water_10m.tif"))

r_drainage <- terra::resample(r_drainage_raw, id_raster_MH, method = "near") %>% terra::crop(id_raster_MH)
r_groen    <- terra::resample(r_groen_raw, id_raster_MH, method = "near") %>% terra::crop(id_raster_MH)
r_grb      <- terra::resample(r_grb_raw, id_raster_MH, method = "near") %>% terra::crop(id_raster_MH)

drain_cats            <- terra::cats(r_drainage)[[1]]
geselecteerde_letters <- c("d", "e", "f", "h", "i", "g", "e-f", "h-i", "e-i")
watersnip_drain_ids   <- drain_cats$value[drain_cats$Label %in% geselecteerde_letters]
watersnip_drainage_ok <- r_drainage %in% watersnip_drain_ids

# Spoor A: MAX
watersnip_basis_max <- (biotopen_max == 1) & watersnip_drainage_ok
watersnip_basis_max <- terra::ifel(watersnip_basis_max == 1, 1, NA)

# Spoor B: OPP
watersnip_basis_opp <- terra::mask(biotopen_opp, watersnip_basis_max)

rm(r_drainage_raw, r_groen_raw, r_grb_raw, watersnip_drainage_ok)
gc()

# ==============================================================================
# STAP 2: PRIMAIRE CLUSTERING (10 HA, 100M BUFFER - PARALLEL)
# ==============================================================================
message("-> Primaire clustering (100m fuzzy koppeling) en filteren op min. 10 ha...")

# Spoor A: MAX
watersnip_clusters_max <- cluster_filter_compleet(
  masker     = watersnip_basis_max,
  opp_laag   = watersnip_basis_opp,
  drempel_m2 = 100000,   
  dist_m     = 100,       
  werkelijk  = FALSE      
)
r_cluster_ids_max <- watersnip_clusters_max$clusters

# Spoor B: OPP (Parallel & Autonoom)
r_binair_basis_opp <- terra::ifel(!is.na(watersnip_basis_opp) & watersnip_basis_opp > 0, 1, NA)
watersnip_clusters_opp <- cluster_filter_compleet(
  masker     = r_binair_basis_opp,
  opp_laag   = watersnip_basis_opp,
  drempel_m2 = 100000,   
  dist_m     = 100,       
  werkelijk  = TRUE      
)
r_cluster_ids_opp <- watersnip_clusters_opp$clusters

rm(r_binair_basis_opp)
gc()

# ==============================================================================
# STAP 3: LANDSCHAPSKWALITEIT CHECKS PER UNIEKE CLUSTER (PARALLEL)
# ==============================================================================
ids_kwaliteit_ok_max <- c()
ids_kwaliteit_ok_opp <- c()

# --- SPOOR A: MAX ---
if (!all(is.na(terra::values(r_cluster_ids_max, mat=FALSE)))) {
  area_ha_raster_max <- r_cluster_ids_max * terra::cellSize(r_cluster_ids_max, unit = "ha")
  totale_opp_max <- terra::zonal(area_ha_raster_max, r_cluster_ids_max, fun = "sum", na.rm = TRUE)
  colnames(totale_opp_max) <- c("Cluster_ID", "Totaal_ha")
  
  bos_masker_max <- terra::mask(r_groen == 1, r_cluster_ids_max)
  bos_ha_raster_max <- bos_masker_max * terra::cellSize(bos_masker_max, unit = "ha")
  bos_max <- terra::zonal(bos_ha_raster_max, r_cluster_ids_max, fun = "sum", na.rm = TRUE)
  colnames(bos_max) <- c("Cluster_ID", "Bos_ha")
  
  water_masker_max <- terra::mask(r_grb == 1, r_cluster_ids_max)
  water_ha_raster_max <- water_masker_max * terra::cellSize(water_masker_max, unit = "ha")
  water_max <- terra::zonal(water_ha_raster_max, r_cluster_ids_max, fun = "sum", na.rm = TRUE)
  colnames(water_max) <- c("Cluster_ID", "Water_ha")
  
  stats_max <- merge(totale_opp_max, bos_max, by = "Cluster_ID", all.x = TRUE) %>%
    merge(water_max, by = "Cluster_ID", all.x = TRUE)
  stats_max[is.na(stats_max)] <- 0
  stats_max$Pct_Bos   <- (stats_max$Bos_ha / stats_max$Totaal_ha) * 100
  stats_max$Pct_Water <- (stats_max$Water_ha / stats_max$Totaal_ha) * 100
  
  ids_kwaliteit_ok_max <- stats_max$Cluster_ID[stats_max$Pct_Bos < 75 & stats_max$Pct_Water < 90]
  rm(area_ha_raster_max, totale_opp_max, bos_masker_max, bos_ha_raster_max, bos_max, water_masker_max, water_ha_raster_max, water_max, stats_max)
}

# --- SPOOR B: OPP (Parallel & Autonoom) ---
if (!all(is.na(terra::values(r_cluster_ids_opp, mat=FALSE)))) {
  totale_opp_opp <- terra::zonal(watersnip_basis_opp, r_cluster_ids_opp, fun = "sum", na.rm = TRUE)
  colnames(totale_opp_opp) <- c("Cluster_ID", "Totaal_ha")
  totale_opp_opp$Totaal_ha <- totale_opp_opp$Totaal_ha * 0.01
  
  bos_in_zone_opp <- terra::mask(terra::ifel(r_groen == 1, 1.0, NA), r_cluster_ids_opp)
  bos_opp <- terra::zonal(bos_in_zone_opp, r_cluster_ids_opp, fun = "sum", na.rm = TRUE)
  colnames(bos_opp) <- c("Cluster_ID", "Bos_ha")
  bos_opp$Bos_ha <- bos_opp$Bos_ha * 0.01
  
  water_in_zone_opp <- terra::mask(terra::ifel(r_grb == 1, 1.0, NA), r_cluster_ids_opp)
  water_opp <- terra::zonal(water_in_zone_opp, r_cluster_ids_opp, fun = "sum", na.rm = TRUE)
  colnames(water_opp) <- c("Cluster_ID", "Water_ha")
  water_opp$Water_ha <- water_opp$Water_ha * 0.01
  
  stats_opp <- merge(totale_opp_opp, bos_opp, by = "Cluster_ID", all.x = TRUE) %>%
    merge(water_opp, by = "Cluster_ID", all.x = TRUE)
  stats_opp[is.na(stats_opp)] <- 0
  stats_opp$Pct_Bos   <- (stats_opp$Bos_ha / stats_opp$Totaal_ha) * 100
  stats_opp$Pct_Water <- (stats_opp$Water_ha / stats_opp$Totaal_ha) * 100
  
  ids_kwaliteit_ok_opp <- stats_opp$Cluster_ID[stats_opp$Pct_Bos < 75 & stats_opp$Pct_Water < 90]
  rm(totale_opp_opp, bos_in_zone_opp, bos_opp, water_in_zone_opp, water_opp, stats_opp)
}
gc()

# ==============================================================================
# STAP 4: TERRESTRISCHE KERN CHECK & FINALE SAMENVOEGING (PARALLEL)
# ==============================================================================
watersnip_leefgebied_max <- template_MH * NA
watersnip_leefgebied_opp <- template_MH * NA

# --- SPOOR A: MAX ---
if (length(ids_kwaliteit_ok_max) > 0) {
  r_terr_basis_max <- watersnip_basis_max & (r_grb == 0 | is.na(r_grb))
  r_terr_basis_max <- terra::ifel(r_terr_basis_max == 1, 1, NA)
  
  terr_patches_max <- cluster_filter_compleet(
    masker     = r_terr_basis_max,
    opp_laag   = watersnip_basis_opp,
    drempel_m2 = 100000,   
    dist_m     = 0,        
    werkelijk  = FALSE      
  )
  r_terr_patch_ids_max <- terr_patches_max$clusters
  
  if (!all(is.na(terra::values(r_terr_patch_ids_max, mat=FALSE)))) {
    df_overlap_max <- as.data.frame(c(r_cluster_ids_max, r_terr_patch_ids_max), na.rm = TRUE)
    if (nrow(df_overlap_max) > 0) {
      colnames(df_overlap_max) <- c("Cluster_ID", "Terr_Patch_ID")
      ids_finaal_max <- unique(df_overlap_max$Cluster_ID[df_overlap_max$Cluster_ID %in% ids_kwaliteit_ok_max])
      
      if (length(ids_finaal_max) > 0) {
        watersnip_leefgebied_max <- terra::ifel(r_cluster_ids_max %in% ids_finaal_max, 1, NA)
      }
    }
    rm(df_overlap_max)
  }
}

# --- SPOOR B: OPP (Parallel & Autonoom) ---
if (length(ids_kwaliteit_ok_opp) > 0) {
  r_terr_basis_opp_raw <- terra::mask(watersnip_basis_opp, (r_grb == 0 | is.na(r_grb)))
  r_terr_basis_opp     <- terra::ifel(!is.na(r_terr_basis_opp_raw) & r_terr_basis_opp_raw > 0, r_terr_basis_opp_raw, NA)
  
  r_binair_terr_opp <- terra::ifel(!is.na(r_terr_basis_opp) & r_terr_basis_opp > 0, 1, NA)
  terr_patches_opp <- cluster_filter_compleet(
    masker     = r_binair_terr_opp,
    opp_laag   = r_terr_basis_opp,
    drempel_m2 = 100000,   
    dist_m     = 0,        
    werkelijk  = TRUE      
  )
  r_terr_patch_ids_opp <- terr_patches_opp$clusters
  
  if (!all(is.na(terra::values(r_terr_patch_ids_opp, mat=FALSE)))) {
    df_overlap_opp <- as.data.frame(c(r_cluster_ids_opp, r_terr_patch_ids_opp), na.rm = TRUE)
    if (nrow(df_overlap_opp) > 0) {
      colnames(df_overlap_opp) <- c("Cluster_ID", "Terr_Patch_ID")
      ids_finaal_opp <- unique(df_overlap_opp$Cluster_ID[df_overlap_opp$Cluster_ID %in% ids_kwaliteit_ok_opp])
      
      if (length(ids_finaal_opp) > 0) {
        watersnip_leefgebied_opp <- terra::mask(watersnip_clusters_opp$raster, r_cluster_ids_opp %in% ids_finaal_opp)
      }
    }
    rm(df_overlap_opp, r_binair_terr_opp)
  }
}

watersnip_leefgebied_max <- terra::crop(watersnip_leefgebied_max, template_MH)
watersnip_leefgebied_opp <- terra::crop(watersnip_leefgebied_opp, template_MH)

cat("Definitief leefgebied Watersnip MAX (ha):", round(calc_ha_exact(watersnip_leefgebied_max), 2), "\n")
cat("Definitief leefgebied Watersnip OPP (ha):", round(calc_ha_exact(watersnip_leefgebied_opp), 2), "\n")

rm(r_drainage, r_ecoregio, r_groen, r_grb, watersnip_basis_max, watersnip_basis_opp, 
   watersnip_clusters_max, watersnip_clusters_opp, r_cluster_ids_max, r_cluster_ids_opp, 
   ids_kwaliteit_ok_max, ids_kwaliteit_ok_opp)
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

