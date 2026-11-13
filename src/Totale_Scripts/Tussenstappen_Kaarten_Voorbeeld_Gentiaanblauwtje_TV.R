library(here)
library(knitr)
library(tidyverse)
library(sf)
library(terra)
library(readxl)
library(tidyterra)
library(data.table)
library(ggplot2)
library(gridExtra)

conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::select)
conflicted::conflicts_prefer(terra::intersect)
conflicted::conflicts_prefer(terra::any)

terra::terraOptions(memfrac = 0.6, memmin = 1, tempdir = tempdir(), verbose = FALSE)

# ==============================================================================
# 0. MAPSTRUCTUUR & HULPFUNCTIES
# ==============================================================================
folder_rapport_kaarten <- here::here("data/output/Rapport_Kaarten")
if (!dir.exists(folder_rapport_kaarten)) dir.create(folder_rapport_kaarten, recursive = TRUE)

# Hulpfunctie voor het genereren van schone rapportkaarten (PNG)
maak_tussenkaart <- function(raster_laag, titel, bestandsnaam, type_schaal = "continu", grens_poly = NULL) {
  if (is.null(raster_laag) || all(is.na(suppressWarnings(terra::minmax(raster_laag))))) {
    message(paste("⚠️ [OVERGESLAGEN]: Geen data om kaart te maken voor", bestandsnaam))
    return(invisible(NULL))
  }
  
  if (type_schaal == "cluster") {
    raster_plot <- terra::as.factor(raster_laag)
  } else {
    raster_plot <- raster_laag
  }
  
  p <- ggplot() +
    geom_spatraster(data = raster_plot)
  
  if (type_schaal == "continu") {
    p <- p + scale_fill_viridis_c(na.value = "transparent", name = "Bedekking / Fractie", option = "viridis")
  } else if (type_schaal == "cluster") {
    p <- p + scale_fill_viridis_d(
      na.value = "transparent", 
      na.translate = FALSE, 
      name = "Cluster ID", 
      option = "turbo"
    )
  } else if (type_schaal == "categorie") {
    p <- p + scale_fill_discrete(
      na.value = "transparent", 
      na.translate = FALSE, 
      name = "BWK Code"
    )
  }
  
  if (!is.null(grens_poly)) {
    p <- p + geom_spatvector(data = grens_poly, fill = NA, color = "#2c3e50", linewidth = 0.8)
  }
  
  p <- p + theme_minimal() +
    labs(
      title = titel,
      subtitle = "Turnhouts Vennengebied (Grens weergegeven als referentie)",
      x = NULL, y = NULL
    ) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(face = "bold", size = 12, color = "#1a1a1a"),
      plot.subtitle = element_text(size = 9, color = "#4a4a4a"),
      legend.position = "right"
    )
  
  ggsave(
    filename = file.path(folder_rapport_kaarten, bestandsnaam),
    plot = p, width = 8, height = 6, dpi = 300, device = "png"
  )
  message(paste("    [PNG GEËXPORTEERD]:", bestandsnaam))
}

# Hulpfunctie voor roterende clusterkleuren (zonder legende)
maak_cluster_roterend_kaart <- function(raster_laag, titel, bestandsnaam, grens_poly = NULL, n_kleuren = 8) {
  if (is.null(raster_laag) || all(is.na(suppressWarnings(terra::minmax(raster_laag))))) {
    message(paste("⚠️ [OVERGESLAGEN]: Geen data voor roterende clusterkaart", bestandsnaam))
    return(invisible(NULL))
  }
  
  r_roterend <- terra::deepcopy(raster_laag)
  
  terra::values(r_roterend) <- ifelse(
    !is.na(terra::values(raster_laag)), 
    (as.numeric(terra::values(raster_laag)) %% n_kleuren) + 1, 
    NA
  )
  
  r_plot <- terra::as.factor(r_roterend)
  
  p <- ggplot() +
    geom_spatraster(data = r_plot) +
    scale_fill_brewer(
      palette = "Set3", 
      na.value = "transparent", 
      guide = "none"
    )
  
  if (!is.null(grens_poly)) {
    p <- p + geom_spatvector(data = grens_poly, fill = NA, color = "#2c3e50", linewidth = 0.8)
  }
  
  p <- p + theme_minimal() +
    labs(
      title = titel,
      subtitle = "Turnhouts Vennengebied — Onderscheidende netwerkclusters (Roterende kleurschaal)",
      x = NULL, y = NULL
    ) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(face = "bold", size = 12, color = "#1a1a1a"),
      plot.subtitle = element_text(size = 9, color = "#4a4a4a")
    )
  
  ggsave(
    filename = file.path(folder_rapport_kaarten, bestandsnaam),
    plot = p, width = 8, height = 6, dpi = 300, device = "png"
  )
  message(paste("    [ROTATING CLUSTER PNG GEËXPORTEERD]:", bestandsnaam))
}

# Hulpfunctie voor uniforme GeoTIFF export op MasterGrid formaat
exporteer_geotiff <- function(raster_laag, bestandsnaam, master_grid, datatype = "FLT4S", na_flag = -9999) {
  if (!is.null(raster_laag) && !all(is.na(suppressWarnings(terra::minmax(raster_laag))))) {
    export_rast <- terra::deepcopy(raster_laag)
    crs(export_rast) <- crs(master_grid)
    export_rast <- terra::extend(export_rast, master_grid, fill = NA)
    
    file_path <- file.path(folder_rapport_kaarten, bestandsnaam)
    
    terra::writeRaster(
      export_rast, 
      filename = file_path, 
      overwrite = TRUE, 
      gdal = c("COMPRESS=LZW"), 
      datatype = datatype, 
      NAflag = na_flag
    )
    message(paste("    [GeoTIFF GEËXPORTEERD]:", bestandsnaam))
  }
}

# Clusterfunctie (exact zoals in het basisscript van gentiaanblauwtje: width = dist_m / 2)
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

# ==============================================================================
# 1. PARAMETERS & INPUTS INLADEN
# ==============================================================================
soort <- "gentiaanblauwtje"
soort_titel <- tools::toTitleCase(soort)
bestands_soort_naam <- tolower(gsub(" ", "", soort))

scenario_rds_path <- "data/input/Scenario_rds/TV_Scenario_BWK_2025.rds"
scenario_path <- here::here(scenario_rds_path)
scenario_naam <- "BWK_2025"

message(paste("Verwerken van soort:", toupper(soort), "binnen scenario:", scenario_naam))

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

area_shape  <- vect(here("data/input/Turnhouts_Vennegebied.shp"))
master_grid <- rast(here("data/input/Raster_Vlaanderen/Vlaanderen_MasterGrid_10m.tif"))[[1]]

area_shape_proj <- project(area_shape, crs(master_grid))
crs(area_shape_proj) <- crs(master_grid)
area_buffer_fix <- buffer(area_shape_proj, width = buffer_m)

# Vertaalraster opbouwen (snelle C++ extractie uit basisscript)
id_raster_TV <- crop(master_grid, area_buffer_fix, snap = "near")
globale_id_raster <- terra::init(master_grid, fun = "cell")

id_raster_TV_globale_values <- crop(globale_id_raster, area_buffer_fix, snap = "near")
id_raster_TV_masked <- mask(id_raster_TV_globale_values, area_buffer_fix)

df_extractie <- terra::as.data.frame(id_raster_TV_masked, cells = TRUE, na.rm = TRUE)
vertaal_df <- as.data.table(df_extractie)
setnames(vertaal_df, c(1, 2), c("lokale_id", "globale_id"))

vertaal_df <- vertaal_df[!is.na(globale_id)]
studiegebied_globale_ids <- unique(vertaal_df$globale_id)

values(id_raster_TV) <- NA
template_TV <- terra::rasterize(area_buffer_fix, id_raster_TV, field = 1, background = 0)

rm(globale_id_raster, id_raster_TV_globale_values, id_raster_TV_masked, df_extractie)
gc()

# ==============================================================================
# 2. CROSSWALK, TABELGENERATIE & RASTER OPBOUW
# ==============================================================================
df_nieuw <- read_csv(
  here("data/input/Excel_files/Resultaten_Totaal_Samengevoegd.csv"), 
  col_types = cols(.default = "c"), 
  show_col_types = FALSE
)

# Export overzichtstabel (CSV + PNG)
tabel_bwk_soort <- df_nieuw %>%
  mutate(Soort_clean = tolower(trimws(Soort))) %>%
  filter(Soort_clean == soort) %>%
  mutate(
    Code_geformatteerd = ifelse(tolower(trimws(Match)) == "bevat", paste0(Code, "%"), Code)
  ) %>%
  select(
    `Biotoop Type` = Type, 
    `BWK Code` = Code_geformatteerd
  ) %>%
  distinct() %>%
  arrange(`Biotoop Type`, `BWK Code`)

write_csv(tabel_bwk_soort, file.path(folder_rapport_kaarten, paste0("Tabel_BWK_Codes_", soort_titel, ".csv")))

tabel_grafisch <- tableGrob(
  tabel_bwk_soort, 
  rows = NULL,
  theme = ttheme_minimal(
    core = list(bg_params = list(fill = c("#f8f9fa", "#ffffff")), fg_params = list(fontsize = 10, hjust = 0, x = 0.05)),
    colhead = list(bg_params = list(fill = "#2c3e50"), fg_params = list(col = "white", fontface = "bold", fontsize = 11, hjust = 0, x = 0.05))
  )
)

ggsave(
  filename = file.path(folder_rapport_kaarten, paste0("Tabel_BWK_Codes_", soort_titel, ".png")),
  plot = tabel_grafisch, width = 4.5, height = max(1.8, nrow(tabel_bwk_soort) * 0.35 + 0.8), dpi = 300
)

# Scenario-RDS koppelen
resultaten_gegroepeerd <- df_nieuw %>%
  mutate(Soort_clean = tolower(trimws(Soort))) %>%
  filter(Soort_clean == soort) %>%
  group_by(Type) %>%
  nest(Data = c(Code, Match))

tabel_vlaanderen <- readRDS(scenario_path)
setDT(tabel_vlaanderen)
tabel_vlaanderen[, CODE := tolower(trimws(CODE))]

lijst_matches      <- list()
lijst_oppervlaktes <- list()
lijst_codes_raster <- list()

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
    
    tabel_TV <- tabel_gefilterd[cel_id %in% studiegebied_globale_ids]
    tabel_TV_unique <- unique(tabel_TV, by = c("cel_id", "CODE"))
    
    tabel_eerste_code <- tabel_TV_unique[, .(Meest_Dominante_Code = dplyr::first(CODE)), by = .(cel_id)]
    tabel_cel_som     <- tabel_TV_unique[, .(Oppervlakte = pmin(sum(BWK_FRAC, na.rm = TRUE), 1.0)), by = .(cel_id)]
    
    r_match_type <- id_raster_TV * NA
    r_opp_type   <- id_raster_TV * NA
    r_code_type  <- id_raster_TV * NA
    
    if(nrow(tabel_cel_som) > 0) {
      tabel_cel_som[, Oppervlakte := pmin(Oppervlakte, 1)]
      tabel_cel_som[, Match := ifelse(Oppervlakte >= 0.01, 1, 0)]
      
      setnames(tabel_cel_som, "cel_id", "globale_id")
      setnames(tabel_eerste_code, "cel_id", "globale_id")
      
      tabel_finaal_mapping <- merge(tabel_cel_som, vertaal_df, by = "globale_id", all.x = TRUE)
      tabel_finaal_mapping <- merge(tabel_finaal_mapping, tabel_eerste_code, by = "globale_id", all.x = TRUE)
      tabel_finaal_mapping <- tabel_finaal_mapping[!is.na(lokale_id)]
      
      if(nrow(tabel_finaal_mapping) > 0) {
        tabel_matches_clean <- tabel_finaal_mapping[Match == 1]
        if(nrow(tabel_matches_clean) > 0) {
          r_match_type[tabel_matches_clean$lokale_id] <- tabel_matches_clean$Match
        }
        r_opp_type[tabel_finaal_mapping$lokale_id] <- tabel_finaal_mapping$Oppervlakte
        
        code_factors <- as.numeric(as.factor(tabel_finaal_mapping$Meest_Dominante_Code))
        r_code_type[tabel_finaal_mapping$lokale_id] <- code_factors
        levels(r_code_type) <- data.frame(ID = unique(code_factors), Code = unique(tabel_finaal_mapping$Meest_Dominante_Code))
        
        lijst_matches[[h_type]]      <- r_match_type
        lijst_oppervlaktes[[h_type]] <- r_opp_type
        lijst_codes_raster[[h_type]] <- r_code_type
      }
    }
  }
}

lijst_matches      <- purrr::compact(lijst_matches)
lijst_oppervlaktes <- purrr::compact(lijst_oppervlaktes)
lijst_codes_raster <- purrr::compact(lijst_codes_raster)

if (length(lijst_matches) == 0 || length(lijst_oppervlaktes) == 0) {
  message("⚠️ Geen enkele doel-BWK code voor '", soort, "' aangetroffen in de database.")
  masker_bron <- template_TV * NA
  opp_bron    <- template_TV * NA
  codes_bron  <- template_TV * NA
} else {
  masker_bron <- lijst_matches[[1]]
  opp_bron    <- lijst_oppervlaktes[[1]]
  codes_bron  <- lijst_codes_raster[[1]]
}

names(masker_bron) <- "Match"
names(opp_bron)    <- "Oppervlakte"

rm(tabel_vlaanderen, vertaal_df, lijst_matches, lijst_oppervlaktes, lijst_codes_raster, df_nieuw, resultaten_gegroepeerd)
gc()

# ==============================================================================
# STAP 1A & 1B: RUWE BWK-CODES EN BIOTOOPFRACTIES VISUALISEREN
# ==============================================================================
message("-> Stap 1a & 1b: Ruwe BWK-verspreiding en biotoopfracties exporteren...")

maak_tussenkaart(
  raster_laag = codes_bron, 
  titel = paste0("Stap 1a: Ruimtelijke Verspreiding van Geselecteerde BWK-codes (", soort_titel, ")"), 
  bestandsnaam = paste0("Stap1a_Verspreiding_BWK_Codes_", soort_titel, ".png"),
  type_schaal = "categorie",
  grens_poly = area_shape_proj
)

exporteer_geotiff(
  raster_laag = codes_bron, 
  bestandsnaam = paste0("Stap1a_BWK_Codes_", bestands_soort_naam, ".tif"), 
  master_grid = master_grid, datatype = "INT2U", na_flag = 0
)

maak_tussenkaart(
  raster_laag = opp_bron, 
  titel = paste0("Stap 1b: Geschikte BWK Biotoopfracties (", soort_titel, ")"), 
  bestandsnaam = paste0("Stap1b_Werkelijk_BWK_Habitat_", soort_titel, ".png"),
  type_schaal = "continu",
  grens_poly = area_shape_proj
)

exporteer_geotiff(
  raster_laag = opp_bron, 
  bestandsnaam = paste0("Stap1b_Werkelijk_BWK_Habitat_", bestands_soort_naam, ".tif"), 
  master_grid = master_grid, datatype = "FLT4S", na_flag = -9999
)

# ==============================================================================
# STAP 2: CLUSTERING OP MAXIMALE POTENTIE
# ==============================================================================
message("-> Stap 2: Clustering op maximale binaire potentie...")

drempel_m2 <- oppervlakte_ha * 10000

if(!all(is.na(suppressWarnings(terra::minmax(masker_bron))))) {
  masker_bron[masker_bron == 0] <- NA
}

res_max <- cluster_filter_compleet(masker_bron, opp_bron, drempel_m2, afstand_m, werkelijk = FALSE)
final_max <- res_max$raster
cl_max    <- res_max$clusters

maak_cluster_roterend_kaart(
  raster_laag  = cl_max, 
  titel        = paste0("Stap 2: Maximale Potentie Clusters (", soort_titel, " ≥ ", oppervlakte_ha, " ha)"), 
  bestandsnaam = paste0("Stap2_Maximale_Potentie_", soort_titel, ".png"),
  grens_poly   = area_shape_proj,
  n_kleuren    = 8
)

exporteer_geotiff(
  raster_laag = final_max, 
  bestandsnaam = paste0("Stap2_Maximale_Potentie_", bestands_soort_naam, ".tif"), 
  master_grid = master_grid, datatype = "INT1U", na_flag = 255
)

# ==============================================================================
# STAP 3: FILTERING OP WERKELIJKE OPPERVLAKTE (DEFINITIEF NETWERK)
# ==============================================================================
message("-> Stap 3: Evalueren van werkelijke oppervlakte binnen clusters...")

if (!all(is.na(suppressWarnings(terra::minmax(cl_max))))) {
  zonal_werkelijk <- terra::zonal(opp_bron, cl_max, fun = "sum", na.rm = TRUE)
  colnames(zonal_werkelijk) <- c("ID", "Werk_m2")
  zonal_werkelijk$Werk_m2 <- zonal_werkelijk$Werk_m2 * 100
  
  ids_werkelijk_voldoet <- zonal_werkelijk$ID[zonal_werkelijk$Werk_m2 >= drempel_m2]
  
  if(length(ids_werkelijk_voldoet) > 0) {
    final_opp <- terra::mask(final_max, cl_max %in% ids_werkelijk_voldoet)
    cl_opp    <- terra::mask(cl_max, cl_max %in% ids_werkelijk_voldoet)
    # Bewaar ook de werkelijke fracties binnen de definitieve clusters
    r_opp_definitief <- terra::mask(opp_bron, cl_opp)
  } else {
    final_opp        <- final_max * NA
    cl_opp           <- cl_max * NA
    r_opp_definitief <- opp_bron * NA
  }
} else {
  final_opp        <- final_max * NA
  cl_opp           <- cl_max * NA
  r_opp_definitief <- opp_bron * NA
}

# Fallback check voor ID-raster
if (is.null(cl_opp) || all(is.na(suppressWarnings(terra::minmax(cl_opp))))) {
  if (!is.null(final_opp) && !all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
    cl_opp <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE)
  } else {
    cl_opp <- template_TV * NA
  }
}

# Kaart van de discrete netwerkclusters die voldoen aan de werkelijke oppervlakte
maak_cluster_roterend_kaart(
  raster_laag  = cl_opp, 
  titel        = paste0("Stap 3: Definitief Werkelijk Habitat (Discrete Clusters ", soort_titel, ")"), 
  bestandsnaam = paste0("Stap3_Definitief_Habitat_Werkelijk_", soort_titel, ".png"),
  grens_poly   = area_shape_proj,
  n_kleuren    = 8
)

# GeoTIFF exports voor Stap 3
exporteer_geotiff(
  raster_laag = r_opp_definitief, 
  bestandsnaam = paste0("Stap3_Definitief_Habitat_Werkelijk_", bestands_soort_naam, ".tif"), 
  master_grid = master_grid, datatype = "FLT4S", na_flag = -9999
)

exporteer_geotiff(
  raster_laag = cl_opp, 
  bestandsnaam = paste0("Stap3_Definitief_Habitat_ID_", bestands_soort_naam, ".tif"), 
  master_grid = master_grid, datatype = "INT4U", na_flag = 0
)

# ==============================================================================
# STAP 4: INLADEN EN EXPORTEREN VAN DE BESTAANDE ARPL (03_ARPL & 00_ID_Rasters)
# ==============================================================================
message("-> Stap 4: Bestaande ARPL inladen en verwerken...")

pad_arpl_bron <- here::here(
  "data/output/Turnhouts_Vennegebied/Rasters_Soorten", 
  scenario_naam, 
  "03_ARPL", 
  paste0("Habitat_ARPL_", bestands_soort_naam, ".tif")
)

pad_id_bron <- here::here(
  "data/output/Turnhouts_Vennegebied/Rasters_Soorten", 
  scenario_naam, 
  "00_ID_Rasters", 
  paste0("ID_Netwerken_", bestands_soort_naam, ".tif")
)

if (file.exists(pad_arpl_bron)) {
  arpl_raster <- terra::rast(pad_arpl_bron)
  
  arpl_cropped <- terra::crop(arpl_raster, area_shape_proj, snap = "near")
  arpl_masked  <- terra::mask(arpl_cropped, area_shape_proj)
  
  if (file.exists(pad_id_bron)) {
    id_raster  <- terra::rast(pad_id_bron)
    id_cropped <- terra::crop(id_raster, area_shape_proj, snap = "near")
    id_masked  <- terra::mask(id_cropped, area_shape_proj)
    arpl_cl    <- terra::mask(id_masked, arpl_masked)
  } else {
    arpl_cl    <- terra::mask(cl_opp, arpl_masked)
  }
  
  maak_tussenkaart(
    raster_laag = arpl_cl, 
    titel = paste0("Stap 4: Actueel Relevant Potentieel Leefgebied (ARPL ", soort_titel, ")"), 
    bestandsnaam = paste0("Stap4_ARPL_Habitat_", soort_titel, ".png"),
    type_schaal = "cluster",
    grens_poly = area_shape_proj
  )
  
  exporteer_geotiff(
    raster_laag = arpl_masked, 
    bestandsnaam = paste0("Stap4_ARPL_Habitat_Werkelijk_", bestands_soort_naam, ".tif"), 
    master_grid = master_grid, 
    datatype = "FLT4S", 
    na_flag = -9999
  )
  
  exporteer_geotiff(
    raster_laag = arpl_cl, 
    bestandsnaam = paste0("Stap4_ARPL_Analytisch_ID_", bestands_soort_naam, ".tif"), 
    master_grid = master_grid, 
    datatype = "INT4U", 
    na_flag = 0
  )
  
} else {
  message("⚠️ [INFO]: Geen bestaand ARPL-bestand gevonden op ", pad_arpl_bron)
}

message(paste("\n🏁 ALLE PNG'S EN GEOTIFFS SUCCESVOL GEËXPORTEERD VOOR:", toupper(soort)))
