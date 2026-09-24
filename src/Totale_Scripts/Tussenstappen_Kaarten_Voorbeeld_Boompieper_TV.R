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

terra::terraOptions(memfrac = 0.6, memmin = 1, verbose = FALSE)

# ==============================================================================
# 0. MAPSTRUCTUUR & HULPFUNCTIES
# ==============================================================================
folder_rapport_kaarten <- here::here("data/output/Rapport_Kaarten")
if (!dir.exists(folder_rapport_kaarten)) dir.create(folder_rapport_kaarten, recursive = TRUE)

# Hulpfunctie voor het genereren van schone rapportkaarten (PNG)
maak_tussenkaart <- function(raster_laag, titel, bestandsnaam, type_schaal = "continu", grens_poly = NULL) {
  
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

cluster_filter_compleet <- function(masker, opp_laag, drempel_m2, dist_m, werkelijk = FALSE) {
  if (terra::global(is.na(masker), "sum")[[1]] == terra::ncell(masker)) {
    return(list(raster = masker * NA, clusters = masker * NA))
  }
  
  tmp_buf <- tempfile(pattern = "cl_buf_", fileext = ".tif")
  on.exit(unlink(tmp_buf), add = TRUE)
  
  if (dist_m > 0) {
    r_binair <- terra::ifel(!is.na(masker) & masker > 0, 1, NA)
    
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
soort <- "boompieper"

scenario_rds_path <- "data/input/Scenario_rds/TV_Scenario_BWK_2025.rds"
scenario_path <- here::here(scenario_rds_path)
scenario_naam <- "BWK_2025"

message(paste("Verwerken van soort:", toupper(soort)))

df <- read_excel(here::here("data/input/Excel_files/Soorten_bwk_afstanden.xlsx"))
resultaat <- df %>%
  filter(tolower(trimws(Soort)) == soort) %>%
  select(Type, MinOpp_ha, AfstandBiotopen_m, Dispersiecap_m)

oppervlakte_ha <- resultaat$MinOpp_ha[1]
afstand_m      <- resultaat$AfstandBiotopen_m[1]
buffer_m       <- resultaat$Dispersiecap_m[1]
rm(df, resultaat)

area_shape  <- vect(here("data/input/Turnhouts_Vennegebied.shp"))
master_grid <- rast(here("data/input/Raster_Vlaanderen/Vlaanderen_MasterGrid_10m.tif"))[[1]]

area_shape_proj <- project(area_shape, crs(master_grid))
crs(area_shape_proj) <- crs(master_grid)
area_buffer_fix <- buffer(area_shape_proj, width = buffer_m)

id_raster_TV <- crop(master_grid, area_buffer_fix, snap = "near")
template_TV_mask <- terra::rasterize(area_buffer_fix, id_raster_TV, field = 1)

lokale_ids  <- terra::cells(template_TV_mask)
coords      <- terra::xyFromCell(template_TV_mask, lokale_ids)
globale_ids <- terra::cellFromXY(master_grid, coords)

vertaal_df <- data.table(lokale_id = lokale_ids, globale_id = globale_ids)[!is.na(globale_id)]
studiegebied_globale_ids <- vertaal_df$globale_id

template_TV <- terra::classify(template_TV_mask, cbind(NA, 0))
values(id_raster_TV) <- NA

rm(template_TV_mask, lokale_ids, coords, globale_ids)
gc()

# ==============================================================================
# 2. CROSSWALK EN TABELGENERATIE
# ==============================================================================
df_nieuw <- read_csv(
  here("data/input/Excel_files/Resultaten_Totaal_Samengevoegd.csv"), 
  col_types = cols(.default = "c"), 
  show_col_types = FALSE
)

# Export overzichtstabel (CSV + PNG)
tabel_bwk_boompieper <- df_nieuw %>%
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

write_csv(tabel_bwk_boompieper, file.path(folder_rapport_kaarten, "Tabel_BWK_Codes_Boompieper.csv"))

tabel_grafisch <- tableGrob(
  tabel_bwk_boompieper, 
  rows = NULL,
  theme = ttheme_minimal(
    core = list(bg_params = list(fill = c("#f8f9fa", "#ffffff")), fg_params = list(fontsize = 10, hjust = 0, x = 0.05)),
    colhead = list(bg_params = list(fill = "#2c3e50"), fg_params = list(col = "white", fontface = "bold", fontsize = 11, hjust = 0, x = 0.05))
  )
)

ggsave(
  filename = file.path(folder_rapport_kaarten, "Tabel_BWK_Codes_Boompieper.png"),
  plot = tabel_grafisch, width = 4.5, height = max(1.8, nrow(tabel_bwk_boompieper) * 0.35 + 0.8), dpi = 300
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

bwk_max        <- lijst_matches[["bwk"]]
bwk_opp        <- lijst_oppervlaktes[["bwk"]]
bwk_codes      <- lijst_codes_raster[["bwk"]]
bwk_bossen_max <- lijst_matches[["bwk_bossen"]]
bwk_bossen_opp <- lijst_oppervlaktes[["bwk_bossen"]]
heide_max      <- lijst_matches[["heide"]]
heide_opp      <- lijst_oppervlaktes[["heide"]]
naaldbos_max   <- lijst_matches[["naaldbos"]]
naaldbos_opp   <- lijst_oppervlaktes[["naaldbos"]]

rm(tabel_vlaanderen, vertaal_df, lijst_matches, lijst_oppervlaktes)
gc()

# ==============================================================================
# STAP 1A & 1B: OUTPUTS
# ==============================================================================
maak_tussenkaart(
  raster_laag = bwk_codes, 
  titel = "Stap 1a: Ruimtelijke Verspreiding van Geselecteerde BWK-codes", 
  bestandsnaam = "Stap1a_Verspreiding_BWK_Codes_Boompieper.png",
  type_schaal = "categorie",
  grens_poly = area_shape_proj
)

exporteer_geotiff(
  raster_laag = bwk_codes, 
  bestandsnaam = paste0("Stap1a_BWK_Codes_", soort, ".tif"), 
  master_grid = master_grid, datatype = "INT2U", na_flag = 0
)

maak_tussenkaart(
  raster_laag = bwk_opp, 
  titel = "Stap 1b: Geschikte BWK Biotoopfracties (Werkelijke Oppervlaktes)", 
  bestandsnaam = "Stap1b_Werkelijk_BWK_Habitat_Boompieper.png",
  type_schaal = "continu",
  grens_poly = area_shape_proj
)

exporteer_geotiff(
  raster_laag = bwk_opp, 
  bestandsnaam = paste0("Stap1b_Werkelijk_BWK_Habitat_", soort, ".tif"), 
  master_grid = master_grid, datatype = "FLT4S", na_flag = -9999
)

# ==============================================================================
# STAP 2: ONAFHANKELIJK CLUSTEREN EN 100M BOSRAND UITSNIJDEN
# ==============================================================================
message("-> Clustering open biotoop en bos (5 ha) + 100m bosrand uitsnijden...")

r_binair_open_opp <- terra::ifel(!is.na(bwk_opp) & bwk_opp > 0, 1, NA)
open_cluster_5ha_opp <- cluster_filter_compleet(
  masker = r_binair_open_opp, opp_laag = bwk_opp, drempel_m2 = 5 * 10000, dist_m = 50, werkelijk = TRUE
)$raster

r_binair_bos_opp <- terra::ifel(!is.na(bwk_bossen_opp) & bwk_bossen_opp > 0, 1, NA)
bos_cluster_5ha_opp <- cluster_filter_compleet(
  masker = r_binair_bos_opp, opp_laag = bwk_bossen_opp, drempel_m2 = 5 * 10000, dist_m = 50, werkelijk = TRUE
)$raster

dist_tot_rand <- terra::distance(!is.na(bos_cluster_5ha_opp))
bos_100m_5ha_opp <- terra::mask(bos_cluster_5ha_opp, dist_tot_rand <= 100, maskvalue = FALSE)

maak_tussenkaart(
  raster_laag = bos_100m_5ha_opp, 
  titel = "Stap 2: Geselecteerde Bosranden (100m Randzone van Bosclusters ≥ 5 ha)", 
  bestandsnaam = "Stap2_Bosranden_100m_Boompieper.png",
  type_schaal = "continu",
  grens_poly = area_shape_proj
)

exporteer_geotiff(
  raster_laag = bos_100m_5ha_opp, 
  bestandsnaam = paste0("Stap2_Bosranden_100m_", soort, ".tif"), 
  master_grid = master_grid, datatype = "FLT4S", na_flag = -9999
)

# ==============================================================================
# STAP 3: CONTACTZONE-ANALYSE VOOR GROOT LEEFGEBIED
# ==============================================================================
message("-> Contactzones (20m edge) bepalen...")

w_edge_20m <- matrix(1, nrow = 3, ncol = 3)

buf_open_max <- terra::focal(terra::ifel(!is.na(open_cluster_5ha_opp), 1, NA), w = w_edge_20m, fun = "max", na.rm = TRUE)
buf_bos_max  <- terra::focal(terra::ifel(!is.na(bos_100m_5ha_opp), 1, NA), w = w_edge_20m, fun = "max", na.rm = TRUE)

contactzone_groot_max <- terra::ifel(!is.na(buf_open_max) & !is.na(buf_bos_max), 1, NA)
bosrand_gekoppeld_opp <- terra::mask(bos_100m_5ha_opp, contactzone_groot_max)
bwk_cluster_opp <- terra::cover(open_cluster_5ha_opp, bosrand_gekoppeld_opp)

maak_tussenkaart(
  raster_laag = bwk_cluster_opp, 
  titel = "Stap 3: Habitat Groot Leefgebied (Gekoppelde Bosrand en Open Biotoop)", 
  bestandsnaam = "Stap3_Habitat_Groot_Leefgebied_Boompieper.png",
  type_schaal = "continu",
  grens_poly = area_shape_proj
)

exporteer_geotiff(
  raster_laag = bwk_cluster_opp, 
  bestandsnaam = paste0("Stap3_Habitat_Groot_Leefgebied_", soort, ".tif"), 
  master_grid = master_grid, datatype = "FLT4S", na_flag = -9999
)

# ==============================================================================
# STAP 4: KLEIN LEEFGEBIED (HEIDE + NAALDBOS, 2 HA)
# ==============================================================================
message("-> Klein Leefgebied bepalen...")

if(!is.null(naaldbos_max) && !is.null(heide_max) && !all(is.na(terra::minmax(naaldbos_max)))) {
  buf_naaldbos <- terra::focal(terra::ifel(naaldbos_max > 0, 1, NA), w = w_edge_20m, fun = "max", na.rm = TRUE)
  buf_heide    <- terra::focal(terra::ifel(heide_max > 0, 1, NA), w = w_edge_20m, fun = "max", na.rm = TRUE)
  contactzone_klein <- terra::ifel(!is.na(buf_naaldbos) & !is.na(buf_heide), 1, NA)
  heide_naaldbos_opp <- terra::mask(heide_opp, contactzone_klein)
} else {
  heide_naaldbos_opp <- id_raster_TV * NA
}

r_binair_klein <- terra::ifel(!is.na(heide_naaldbos_opp) & heide_naaldbos_opp > 0, 1, NA)
bwk_cluster_klein_opp_res <- cluster_filter_compleet(r_binair_klein, heide_naaldbos_opp, dist_m = 50, drempel_m2 = 2 * 10000, werkelijk = TRUE)
bwk_cluster_klein_opp <- bwk_cluster_klein_opp_res$raster

maak_tussenkaart(
  raster_laag = bwk_cluster_klein_opp, 
  titel = "Stap 4: Habitat Klein Leefgebied (Heide-Naaldbos Overgang zones ≥ 2 ha)", 
  bestandsnaam = "Stap4_Habitat_Klein_Leefgebied_Boompieper.png",
  type_schaal = "continu",
  grens_poly = area_shape_proj
)

exporteer_geotiff(
  raster_laag = bwk_cluster_klein_opp, 
  bestandsnaam = paste0("Stap4_Habitat_Klein_Leefgebied_", soort, ".tif"), 
  master_grid = master_grid, datatype = "FLT4S", na_flag = -9999
)

# ==============================================================================
# HULPFUNCTIE VOOR ROTERENDE CLUSTERKLEUREN (ZONDER LEGENDE)
# ==============================================================================
maak_cluster_roterend_kaart <- function(raster_laag, titel, bestandsnaam, grens_poly = NULL, n_kleuren = 8) {
  
  # Maak een kopie van het raster om te transformeren
  r_roterend <- terra::deepcopy(raster_laag)
  
  # Roteer de unieke ID's over een beperkt aantal kleurklassen (modulo)
  # Bijvoorbeeld: ID 1 t/m 8 krijgen klasse 1 t/m 8, ID 9 krijgt weer klasse 1, etc.
  terra::values(r_roterend) <- ifelse(
    !is.na(terra::values(raster_laag)), 
    (as.numeric(terra::values(raster_laag)) %% n_kleuren) + 1, 
    NA
  )
  
  # Zet om naar factor voor discrete kleuring
  r_plot <- terra::as.factor(r_roterend)
  
  p <- ggplot() +
    geom_spatraster(data = r_plot) +
    # Gebruik een harmonieus, goed contrasterend discreet kleurenpalet
    scale_fill_brewer(
      palette = "Set3", 
      na.value = "transparent", 
      guide = "none" # Verberg de legende met honderden nummers
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

# ==============================================================================
# STAP 5: DEFINITIEF HABITAT COMBINEREN EN CLUSTEREN
# ==============================================================================
message("-> Stap 5: Definitief habitat combineren en clusteren...")

# 1. Combineer Groot en Klein Leefgebied (Stap 3 & Stap 4)
bwk_gecombineerd_opp <- terra::cover(bwk_cluster_opp, bwk_cluster_klein_opp)

# 2. Maak binair masker voor clustering
r_binair_definitief <- terra::ifel(!is.na(bwk_gecombineerd_opp) & bwk_gecombineerd_opp > 0, 1, NA)

# 3. Voer clusterfiltering uit op het gecombineerde habitat
# Geef 'oppervlakte_ha' en 'buffer_m' mee die in Stap 1 zijn ingeladen
definitief_cluster_res <- cluster_filter_compleet(
  masker    = r_binair_definitief, 
  opp_laag  = bwk_gecombineerd_opp, 
  drempel_m2 = oppervlakte_ha * 10000, 
  dist_m    = buffer_m, 
  werkelijk = TRUE
)

# Ken de resultaten expliciet toe aan de globale omgeving
r_finaal  <- definitief_cluster_res$raster
cl_finaal <- definitief_cluster_res$clusters

# Exporteer GeoTIFFs voor Stap 5
exporteer_geotiff(
  raster_laag = r_finaal, 
  bestandsnaam = paste0("Stap5_Definitief_Habitat_Werkelijk_", soort, ".tif"), 
  master_grid = master_grid, datatype = "FLT4S", na_flag = -9999
)

exporteer_geotiff(
  raster_laag = cl_finaal, 
  bestandsnaam = paste0("Stap5_Definitief_Habitat_ID_", soort, ".tif"), 
  master_grid = master_grid, datatype = "INT4U", na_flag = 0
)

# Render de roterende clusterkaart
maak_cluster_roterend_kaart(
  raster_laag  = cl_finaal, 
  titel        = "Stap 5: Definitief Werkelijk Habitat (Discrete Clusters Boompieper)", 
  bestandsnaam = "Stap5_Definitief_Habitat_Werkelijk_Boompieper.png",
  grens_poly   = area_shape_proj,
  n_kleuren    = 8
)

# ==============================================================================
# STAP 6: INLADEN EN EXPORTEREN VAN DE BESTAANDE ARPL (03_ARPL & 00_ID_Rasters)
# ==============================================================================
message("-> Stap 6: Bestaande ARPL inladen en verwerken...")

pad_arpl_bron <- here::here(
  "data/output/Turnhouts_Vennegebied/Rasters_Soorten", 
  scenario_naam, 
  "03_ARPL", 
  paste0("Habitat_ARPL_", soort, ".tif")
)

pad_id_bron <- here::here(
  "data/output/Turnhouts_Vennegebied/Rasters_Soorten", 
  scenario_naam, 
  "00_ID_Rasters", 
  paste0("ID_Netwerken_", soort, ".tif")
)

if (file.exists(pad_arpl_bron)) {
  # Inladen van het bestaande ARPL raster
  arpl_raster <- terra::rast(pad_arpl_bron)
  
  # Afknippen en maskeren op studiegebied
  arpl_cropped <- terra::crop(arpl_raster, area_shape_proj, snap = "near")
  arpl_masked  <- terra::mask(arpl_cropped, area_shape_proj)
  
  # Koppeling maken met de unieke cluster ID's (uit 00_ID_Rasters of Stap 5)
  if (file.exists(pad_id_bron)) {
    id_raster <- terra::rast(pad_id_bron)
    id_cropped <- terra::crop(id_raster, area_shape_proj, snap = "near")
    id_masked  <- terra::mask(id_cropped, area_shape_proj)
    arpl_cl    <- terra::mask(id_masked, arpl_masked)
  } else {
    arpl_cl    <- terra::mask(cl_finaal, arpl_masked)
  }
  
  # --- STAP 6 OUTPUTS ---
  # 1. Visualisatie met unieke clusterkleuren voor actieve netwerken
  maak_tussenkaart(
    raster_laag = arpl_cl, 
    titel = "Stap 6: Actueel Relevant Potentieel Leefgebied (ARPL Discrete Clusters)", 
    bestandsnaam = "Stap6_ARPL_Habitat_Boompieper.png",
    type_schaal = "cluster",
    grens_poly = area_shape_proj
  )
  
  # 2. GeoTIFF export van het ARPL raster naar de rapportmap
  exporteer_geotiff(
    raster_laag = arpl_masked, 
    bestandsnaam = paste0("Stap6_ARPL_Habitat_Werkelijk_", soort, ".tif"), 
    master_grid = master_grid, 
    datatype = "FLT4S", 
    na_flag = -9999
  )
  
  # 3. GeoTIFF export van het ARPL cluster ID-raster naar de rapportmap
  exporteer_geotiff(
    raster_laag = arpl_cl, 
    bestandsnaam = paste0("Stap6_ARPL_Analytisch_ID_", soort, ".tif"), 
    master_grid = master_grid, 
    datatype = "INT4U", 
    na_flag = 0
  )
  
} else {
  message("⚠️ [INFO]: Geen bestaand ARPL-bestand gevonden op ", pad_arpl_bron)
}

message(paste("\n🏁 ALLE PNG'S EN GEOTIFFS INCLUSIEF ARPL (STAP 6) SUCCESVOL GEËXPORTEERD VOOR:", toupper(soort)))
