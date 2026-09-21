library(here)
library(knitr)
library(tidyverse)
library(sf)
library(terra)
library(readxl)
library(tidyterra)
library(data.table)

conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::select)
conflicted::conflicts_prefer(terra::intersect)

calc_ha_exact <- function(r) {
  if(is.null(r) || all(is.na(suppressWarnings(terra::minmax(r))))) return(0)
  area_raster <- r * terra::cellSize(r, unit = "ha")
  val <- terra::global(area_raster, "sum", na.rm = TRUE)[[1]]
  return(as.numeric(val))
}

cluster_filter_compleet <- function(masker, opp_laag, drempel_m2, dist_m, werkelijk = FALSE) {
  if (all(is.na(suppressWarnings(terra::minmax(masker))))) {
    return(list(raster = masker * NA, clusters = masker * NA))
  }
  
  if (dist_m > 0) {
    r_binair <- terra::ifel(!is.na(masker) & masker > 0, 1, NA)
    r_buffered <- terra::buffer(r_binair, width = dist_m / 2)
    cl_network <- terra::patches(r_buffered, directions = 4, zeroAsNA = TRUE)
    cl_biotoop_only <- terra::mask(cl_network, masker)
  } else {
    cl_network <- terra::patches(masker, directions = 8, zeroAsNA = TRUE)
    cl_biotoop_only <- cl_network
  }
  
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
  
  voldoet_ids <- stats_df$ID[stats_df$Area_m2 >= drempel_m2]
  if(length(voldoet_ids) == 0) return(list(raster = masker * NA, clusters = masker * NA))
  
  masker_binair <- cl_biotoop_only %in% voldoet_ids
  final_network_mask <- terra::ifel(masker_binair == 1, 1, NA)
  
  r_finaal  <- terra::mask(masker, final_network_mask)
  cl_finaal <- terra::mask(cl_biotoop_only, r_finaal) 
  
  return(list(raster = r_finaal, clusters = cl_finaal))
}

terraOptions(
  memfrac = 0.5,
  tempdir = tempdir(),
  verbose = FALSE
)

# ------------------------------------------------------------------------------
# ROBUUSTE DYNAMISCHE SOORT PARAMETER CHECK (INCLUSIEF ALLE VARIANTEN)
# ------------------------------------------------------------------------------
if (exists("huidige_soort") && !is.null(huidige_soort)) {
  soort <- tolower(trimws(huidige_soort))
} else if (exists("SOORT_INVOER") && !is.null(SOORT_INVOER)) {
  soort <- tolower(trimws(SOORT_INVOER))
} else if (exists("soort_invoer") && !is.null(soort_invoer)) {
  soort <- tolower(trimws(soort_invoer))
} else if (exists("SOORT") && !is.null(SOORT)) {
  soort <- tolower(trimws(SOORT))
} else if (exists("params") && !is.null(params$soort_invoer)) {
  soort <- tolower(trimws(params$soort_invoer))
} else {
  soort <- "kleintasjeskruid"
}

# --- DYNAMISCHE SCENARIO PARAMETER CHECK ---
if (exists("SCENARIO_RDS_PAD") && !is.null(SCENARIO_RDS_PAD)) {
  scenario_rds_path <- SCENARIO_RDS_PAD
} else if (exists("params") && !is.null(params$scenario_rds_path)) {
  scenario_rds_path <- params$scenario_rds_path
} else {
  scenario_rds_path <- "data/input/Scenario_rds/VS_Scenario_BWK_2025.rds"
}

p_raw <- gsub("^([.][.]/)+", "", scenario_rds_path)
scenario_path <- here::here(p_raw)

if (!file.exists(scenario_path)) {
  stop(paste("❌ FOUT: Scenario RDS bestand NIET gevonden op:", scenario_path))
}

scen_volledig <- basename(scenario_path)
scenario_naam <- gsub("^VS_Scenario_|^Scenario_|.rds$", "", scen_volledig)

message(paste("Verwerken van soort:", soort, "binnen scenario:", scenario_naam))

df <- read_excel(here::here("data/input/Excel_files/Soorten_bwk_afstanden.xlsx"))
resultaat <- df %>%
  mutate(Soort_clean = tolower(trimws(Soort))) %>% 
  filter(Soort_clean == soort) %>%
  select(Type, MinOpp_ha, AfstandBiotopen_m, Dispersiecap_m)

if(nrow(resultaat) == 0) {
  warning("⚠️ Soort '", soort, "' niet aangetroffen in Soorten_bwk_afstanden.xlsx. Standaardwaarden worden gebruikt.")
  type           <- "landbiotoop_bwk"
  oppervlakte_ha <- 0.1
  afstand_m      <- 0
  buffer_m       <- 500
} else {
  type           <- resultaat$Type[1]
  oppervlakte_ha <- resultaat$MinOpp_ha[1]
  afstand_m      <- resultaat$AfstandBiotopen_m[1]
  buffer_m       <- resultaat$Dispersiecap_m[1]
}

rm(df, resultaat)

area_shape  <- vect(here("data/input/Voerstreek.shp"))
master_grid <- rast(here("data/input/Raster_Vlaanderen/Vlaanderen_MasterGrid_10m.tif"))[[1]]

area_shape_proj <- project(area_shape, crs(master_grid))
area_buffer_fix <- buffer(area_shape_proj, width = buffer_m)

message("-> Vertaalraster voor globale/lokale cellen opbouwen via snelle MASK methode...")
id_raster_VS <- crop(master_grid, area_buffer_fix, snap = "near")

globale_id_raster <- master_grid
globale_id_raster <- terra::init(globale_id_raster, fun = "cell")

id_raster_VS_globale_values <- crop(globale_id_raster, area_buffer_fix, snap = "near")
id_raster_VS_masked <- mask(id_raster_VS_globale_values, area_buffer_fix)

message("-> Vertaaltabel bliksemsnel opbouwen via C++ dataframe extractie...")

df_extractie <- terra::as.data.frame(id_raster_VS_masked, cells = TRUE, na.rm = TRUE)
vertaal_df <- as.data.table(df_extractie)
setnames(vertaal_df, c(1, 2), c("lokale_id", "globale_id"))

vertaal_df <- vertaal_df[!is.na(globale_id)]
studiegebied_globale_ids <- unique(vertaal_df$globale_id)

values(id_raster_VS) <- NA
template_VS <- terra::rasterize(area_buffer_fix, id_raster_VS, field = 1, background = 0)

rm(globale_id_raster, id_raster_VS_globale_values, id_raster_VS_masked, df_extractie)
gc()

# 1. LAAD DE BRON-CROSSWALK IN
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

if(nrow(resultaten_gegroepeerd) > 0) {
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
    
    tabel_VS <- tabel_gefilterd[cel_id %in% studiegebied_globale_ids]
    tabel_VS_unique <- unique(tabel_VS, by = c("cel_id", "CODE"))
    
    tabel_cel_som <- tabel_VS_unique[, .(Oppervlakte = pmin(sum(BWK_FRAC, na.rm = TRUE), 1.0)), by = .(cel_id)]
    
    r_match_type <- id_raster_VS * NA
    r_opp_type   <- id_raster_VS * NA
    
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
}

lijst_matches      <- purrr::compact(lijst_matches)
lijst_oppervlaktes <- purrr::compact(lijst_oppervlaktes)

if (length(lijst_matches) == 0 || length(lijst_oppervlaktes) == 0) {
  message("⚠️ KRIJT: Geen enkele doel-BWK code voor '", soort, "' aangetroffen in de database. Failsafe ingeschakeld.")
  masker_bron <- template_VS * NA
  opp_bron    <- template_VS * NA
} else {
  masker_bron <- lijst_matches[[1]]
  opp_bron    <- lijst_oppervlaktes[[1]]
}

names(masker_bron) <- "Match"
names(opp_bron)    <- "Oppervlakte"

rm(tabel_vlaanderen, vertaal_df, lijst_matches, lijst_oppervlaktes, df_nieuw, resultaten_gegroepeerd)
gc()

drempel_m2 <- oppervlakte_ha * 10000

if(!all(is.na(suppressWarnings(terra::minmax(masker_bron))))) {
  masker_bron[masker_bron == 0] <- NA
}

# 1. Bereken de clusters op basis van de maximale binaire potentie
res_max <- cluster_filter_compleet(masker_bron, opp_bron, drempel_m2, afstand_m, werkelijk = FALSE)

final_max <- res_max$raster
cl_max    <- res_max$clusters

# 2. Evalueer de werkelijke oppervlakte binnen de cl_max clusters
if (!all(is.na(suppressWarnings(terra::minmax(cl_max))))) {
  zonal_werkelijk <- terra::zonal(opp_bron, cl_max, fun = "sum", na.rm = TRUE)
  colnames(zonal_werkelijk) <- c("ID", "Werk_m2")
  zonal_werkelijk$Werk_m2 <- zonal_werkelijk$Werk_m2 * 100
  
  ids_werkelijk_voldoet <- zonal_werkelijk$ID[zonal_werkelijk$Werk_m2 >= drempel_m2]
  
  if(length(ids_werkelijk_voldoet) > 0) {
    final_opp <- terra::mask(final_max, cl_max %in% ids_werkelijk_voldoet)
    cl_opp    <- terra::mask(cl_max, cl_max %in% ids_werkelijk_voldoet)
  } else {
    final_opp <- final_max * NA
    cl_opp    <- cl_max * NA
  }
} else {
  final_opp <- final_max * NA
  cl_opp    <- cl_max * NA
}

# FIX: Snelle minmax-check zonder RAM-inlaad crash
if (is.null(cl_opp) || all(is.na(suppressWarnings(terra::minmax(cl_opp))))) {
  if (!is.null(final_opp) && !all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
    cl_opp <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE)
  } else {
    cl_opp <- template_VS * NA
  }
}

# ==============================================================================
# SCHONE EXPORT BIOTOOP EN ANALYTISCH ID-RASTER (TURNHOUTS VENNEGEBIED)
# ==============================================================================
base_dir <- here::here("data/output/Voerstreek/Rasters_Soorten", scenario_naam)

folders <- list(
  potentie  = file.path(base_dir, "01_Maximale_Potentie"),
  werkelijk = file.path(base_dir, "02_Werkelijke_Oppervlaktes"),
  id_raster = file.path(base_dir, "00_ID_Rasters")
)
purrr::walk(folders, ~if (!dir.exists(.x)) dir.create(.x, showWarnings = FALSE, recursive = TRUE))

# 1. Bepaal Maximale Potentie Raster
potentie_export_rast <- if (exists("final_max") && !is.null(final_max) && !all(is.na(suppressWarnings(terra::minmax(final_max))))) {
  terra::ifel(!is.na(final_max) & final_max > 0, 1, NA)
} else {
  terra::rast(template_VS, vals = NA)
}

# 2. Bepaal Werkelijke Oppervlakte Raster
werkelijk_export_rast <- if (exists("final_opp") && !is.null(final_opp) && !all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  terra::ifel(!is.na(final_opp) & final_opp > 0, 1, NA)
} else {
  terra::rast(template_VS, vals = NA)
}

# 3. Bepaal Analytisch Metacluster ID-raster
if (exists("cl_opp") && !is.null(cl_opp) && !all(is.na(suppressWarnings(terra::minmax(cl_opp))))) {
  id_export_rast <- cl_opp
} else if (exists("final_opp") && !is.null(final_opp) && !all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  id_export_rast <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE)
} else {
  id_export_rast <- terra::rast(template_VS, vals = NA)
}

# Opschonen van de soortnaam voor bestandsnamen (zonder spaties)
bestands_soort_naam <- tolower(gsub(" ", "", soort))

# --- EXPORT LUS VOOR VISUELE BINAIR RASTERS ---
export_config <- list(
  list(rast = potentie_export_rast,  naam = "Maximale_Potentie",        folder = folders$potentie),
  list(rast = werkelijk_export_rast, naam = "Werkelijke_Oppervlaktes", folder = folders$werkelijk)
)

for(item in export_config) {
  suffix <- if (exists("zoek_sleutel") && grepl("_wv$", zoek_sleutel)) "_wv.tif" else ".tif"
  file_path <- file.path(item$folder, paste0("Habitat_", item$naam, "_", bestands_soort_naam, suffix))
  
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
  message(paste("    [OK] Geëxporteerd naar scenariomap:", basename(file_path)))
}

# --- EXPORT VOOR ANALYTISCH ID-RASTER ---
suffix <- if (exists("zoek_sleutel") && grepl("_wv$", zoek_sleutel)) "_wv.tif" else ".tif"
file_path_id <- file.path(folders$id_raster, paste0("ID_Netwerken_", bestands_soort_naam, suffix))

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

suppressWarnings(rm(potentie_export_rast, werkelijk_export_rast, id_export_rast, export_config, export_rast, masker_bron, opp_bron))
gc()

message(paste("🏁 SCENARIO EXPORT VOLLEDIG AFGEROND VOOR:", toupper(soort)))
