# ==============================================================================
# LEEFGEBIEDENGRID VOOR HEIVLINDER - HEEL VLAANDEREN (WERKELIJKE OPPERVLAKTE)
# Ecologische analyse open biotoop, zandsubstraat, beschutting & 500m afstandsmatrix (INBO)
# Author: Bert Van Hecke / Update 2026
# ==============================================================================

# 1. BIBLIOTHEKEN LADEN & CONFLICTEN BEHEREN -----------------------------------
library(here)
library(knitr)
library(tidyverse)
library(sf)
library(terra)
library(readxl)
library(tidyterra)
library(data.table)
library(arrow)
library(conflicted)

conflicted::conflicts_prefer(dplyr::filter)
conflicted::conflicts_prefer(dplyr::select)
conflicted::conflicts_prefer(dplyr::first)
conflicted::conflicts_prefer(terra::intersect)
conflicted::conflicts_prefer(terra::any)

# Instellingen voor terra (RAM-geheugen beheer)
terraOptions(
  memfrac = 0.8,       # Max. 80% van het RAM-geheugen gebruiken
  todisk = TRUE,
  tempdir = tempdir(), 
  verbose = FALSE
)

# 2. HELPER FUNCTIES -----------------------------------------------------------
calc_ha_exact <- function(r) {
  if (is.null(r)) return(0)
  if (all(is.na(terra::values(r, mat = FALSE)))) return(0)
  area_raster <- r * terra::cellSize(r, unit = "ha")
  val <- terra::global(area_raster, "sum", na.rm = TRUE)[[1]]
  return(as.numeric(val))
}

cluster_filter_compleet <- function(masker, opp_laag, drempel_m2, dist_m, werkelijk = FALSE) {
  if (terra::global(is.na(masker), "sum")[[1]] == terra::ncell(masker)) {
    return(list(raster = masker * NA, clusters = masker * NA))
  }
  
  # 1. Binaire kaart maken
  r_binair <- terra::ifel(!is.na(masker) & masker > 0, 1, NA)
  
  # 2. Netwerkvorming via tijdelijke schijfbestanden
  if (dist_m > 0) {
    stral_cellen <- ceiling((dist_m / 2) / 10)
    f_matrix     <- matrix(1, nrow = (2 * stral_cellen + 1), ncol = (2 * stral_cellen + 1))
    
    tmp_focal  <- tempfile(fileext = ".tif")
    r_buffered <- terra::focal(r_binair, w = f_matrix, fun = "max", na.rm = TRUE, 
                               filename = tmp_focal, overwrite = TRUE)
    r_buffered <- terra::ifel(r_buffered > 0, 1, NA)
    
    tmp_patch  <- tempfile(fileext = ".tif")
    cl_network <- terra::patches(r_buffered, directions = 4, zeroAsNA = TRUE,
                                 filename = tmp_patch, overwrite = TRUE)
    
    cl_biotoop_only <- terra::mask(cl_network, masker)
    
    unlink(c(tmp_focal, tmp_patch))
    rm(r_buffered, cl_network, r_binair)
  } else {
    tmp_patch       <- tempfile(fileext = ".tif")
    cl_biotoop_only <- terra::patches(masker, directions = 8, zeroAsNA = TRUE,
                                      filename = tmp_patch, overwrite = TRUE)
  }
  
  # 3. Slimme extractie: Enkel actieve cellen verwerken
  df_cl <- terra::as.data.frame(cl_biotoop_only, cells = TRUE)
  if (nrow(df_cl) == 0) return(list(raster = masker * NA, clusters = masker * NA))
  colnames(df_cl) <- c("cell", "ID")
  
  if (werkelijk) {
    opp_waarden <- terra::extract(opp_laag, df_cl$cell)[[1]]
    dt_calc <- data.table(ID = df_cl$ID, Waarde = opp_waarden)
  } else {
    dt_calc <- data.table(ID = df_cl$ID, Waarde = 1)
  }
  
  stats_dt <- dt_calc[!is.na(ID) & !is.na(Waarde), .(Area_m2 = sum(Waarde, na.rm = TRUE) * 100), by = ID]
  voldoet_ids <- stats_dt[Area_m2 >= drempel_m2, ID]
  
  if (length(voldoet_ids) == 0) return(list(raster = masker * NA, clusters = masker * NA))
  
  # 4. Resultaat filteren
  cl_finaal <- terra::ifel(cl_biotoop_only %in% voldoet_ids, cl_biotoop_only, NA)
  r_finaal  <- terra::mask(masker, cl_finaal)
  
  rm(df_cl, dt_calc, stats_dt, cl_biotoop_only)
  gc()
  
  return(list(raster = r_finaal, clusters = cl_finaal))
}

# 3. SOORT INFORMATIE OPHALEN --------------------------------------------------
df <- read_excel(here("data/input/Excel_files/Soorten_bwk_afstanden.xlsx"))
soort <- "heivlinder"

resultaat <- df %>% filter(tolower(trimws(Soort)) == soort) %>% select(Type, MinOpp_ha, AfstandBiotopen_m, Dispersiecap_m)

oppervlakte_ha <- resultaat$MinOpp_ha[1]
afstand_m      <- resultaat$AfstandBiotopen_m[1]
buffer_m       <- resultaat$Dispersiecap_m[1]

print(resultaat)
rm(df, resultaat)

# 4. GEBIED VOORBEREIDEN (MASTERGRID VLAANDEREN) -------------------------------
message("-> Mastergrid Vlaanderen inlezen...")
master_grid <- rast(here("data/input/Raster_Vlaanderen/Vlaanderen_MasterGrid_10m.tif"))[[1]]

template_Vlaanderen <- master_grid
values(template_Vlaanderen) <- NA

# 5. BIOTOOPFILTERING VIA PARQUET (BWK TYPETABELLEN) ---------------------------
message("-> BWK Parquet dataset koppelen via Arrow...")

df_nieuw <- read_csv(
  here("data/input/Excel_files/Resultaten_Totaal_Samengevoegd.csv"), 
  show_col_types = FALSE,
  progress = FALSE
)

resultaten_gegroepeerd <- df_nieuw %>%
  mutate(Soort_clean = tolower(trimws(Soort))) %>%
  filter(Soort_clean == soort) %>%
  group_by(Type) %>%
  nest(Data = c(Code, Match))

ds_vlaanderen <- open_dataset(here("data/input/Raster_Vlaanderen/BWK_TidyTabel_Smal_Vlaanderen_2025.parquet"))

lijst_oppervlaktes <- list()

for(i in 1:nrow(resultaten_gegroepeerd)) {
  h_type <- resultaten_gegroepeerd$Type[i]
  h_data <- resultaten_gegroepeerd$Data[[i]]
  
  exact_codes <- tolower(trimws(h_data$Code[h_data$Match == "exact"]))
  bevat_codes <- tolower(trimws(h_data$Code[h_data$Match == "bevat"]))
  
  if (length(bevat_codes) > 0) {
    bevat_codes_escaped  <- gsub("([\\.\\^\\$\\*\\+\\?\\(\\)\\[\\{\\\\\\|])", "\\\\\\1", bevat_codes)
    bevat_codes_anchored <- paste0("^", bevat_codes_escaped)
    regex_term           <- paste0(bevat_codes_anchored, collapse = "|")
    
    tabel_sub <- ds_vlaanderen %>%
      filter(tolower(CODE) %in% exact_codes | str_detect(tolower(CODE), regex_term)) %>%
      collect() %>%
      as.data.table()
  } else {
    tabel_sub <- ds_vlaanderen %>%
      filter(tolower(CODE) %in% exact_codes) %>%
      collect() %>%
      as.data.table()
  }
  
  if (nrow(tabel_sub) > 0) {
    tabel_sub[, CODE := tolower(trimws(CODE))]
    tabel_unique  <- unique(tabel_sub, by = c("cel_id", "CODE"))
    tabel_cel_som <- tabel_unique[, .(Oppervlakte = pmin(sum(BWK_FRAC, na.rm = TRUE), 1.0)), by = .(cel_id)]
    
    r_opp_type <- template_Vlaanderen
    r_opp_type[tabel_cel_som$cel_id] <- tabel_cel_som$Oppervlakte
    lijst_oppervlaktes[[h_type]]     <- r_opp_type
  } else {
    lijst_oppervlaktes[[h_type]]     <- template_Vlaanderen
  }
}

heide_opp          <- if("bwk_heide" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["bwk_heide"]] else template_Vlaanderen
duin_opp           <- if("bwk_duin" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["bwk_duin"]] else template_Vlaanderen
grasland_opp       <- if("bwk_grasland" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["bwk_grasland"]] else template_Vlaanderen
groeve_terril_opp  <- if("bwk_groeve_terril" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["bwk_groeve_terril"]] else template_Vlaanderen
zand_opp           <- if("zand" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["zand"]] else template_Vlaanderen
beschutting_opp    <- if("bwk_beschutting" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["bwk_beschutting"]] else template_Vlaanderen

rm(ds_vlaanderen, resultaten_gegroepeerd, df_nieuw, lijst_oppervlaktes)
gc()

# 6. BIOTOOP & SUBSTRAAT CLUSTERING (WERKELIJKE OPPERVLAKTE) ------------------
message("-> Open biotooplagen samenvoegen en clusteren (10 ha, 50m netwerk)...")

h_clean <- terra::ifel(is.na(heide_opp), 0, heide_opp)
d_clean <- terra::ifel(is.na(duin_opp), 0, duin_opp)
g_clean <- terra::ifel(is.na(grasland_opp), 0, grasland_opp)
t_clean <- terra::ifel(is.na(groeve_terril_opp), 0, groeve_terril_opp)

heivlinder_biotoop_opp_raw <- terra::clamp(h_clean + d_clean + g_clean + t_clean, upper = 1.0)
heivlinder_biotoop_opp     <- terra::ifel(heivlinder_biotoop_opp_raw > 0, heivlinder_biotoop_opp_raw, NA)

rm(h_clean, d_clean, g_clean, t_clean, heide_opp, duin_opp, grasland_opp, groeve_terril_opp, heivlinder_biotoop_opp_raw)
gc()

r_binair_biotoop_opp    <- terra::ifel(!is.na(heivlinder_biotoop_opp) & heivlinder_biotoop_opp > 0, 1, NA)
heivlinder_clusters_opp <- cluster_filter_compleet(
  masker     = r_binair_biotoop_opp,
  opp_laag   = heivlinder_biotoop_opp,
  drempel_m2 = 100000,   # 10 ha
  dist_m     = 50,
  werkelijk  = TRUE
)

rm(r_binair_biotoop_opp)
gc()

message("-> Zandsubstraat clusteren (0.5 ha, 100m netwerk)...")
r_binair_zand_opp <- terra::ifel(!is.na(zand_opp) & zand_opp > 0, 1, NA)
zand_clusters_opp <- cluster_filter_compleet(
  masker     = r_binair_zand_opp,
  opp_laag   = zand_opp,
  drempel_m2 = 5000,      # 0.5 ha
  dist_m     = 100,
  werkelijk  = TRUE
)

rm(r_binair_zand_opp, zand_opp)
gc()

# 7. GROENKAART & BESCHUTTING KOPPELEN -----------------------------------------
message("-> Groenkaart inlezen en synchroniseren voor beschutting...")

hooggroen_raw <- rast(here("data/input/ASCI Files/Groenkaart_2021.tif"))
if(!terra::compareGeom(hooggroen_raw, template_Vlaanderen, stopOnError = FALSE)) {
  groenkaart_10m_synchronized <- terra::resample(hooggroen_raw, template_Vlaanderen, method = "near")
} else {
  groenkaart_10m_synchronized <- hooggroen_raw
}

groen_0 <- terra::ifel(is.na(groenkaart_10m_synchronized), 0, groenkaart_10m_synchronized)
bwk_opp_0 <- terra::ifel(is.na(beschutting_opp), 0, beschutting_opp)

heivlinder_beschutting_opp_raw <- terra::clamp(bwk_opp_0 + groen_0, upper = 1.0)
heivlinder_beschutting_opp     <- terra::ifel(heivlinder_beschutting_opp_raw > 0, heivlinder_beschutting_opp_raw, NA)

rm(hooggroen_raw, groenkaart_10m_synchronized, groen_0, bwk_opp_0, beschutting_opp, heivlinder_beschutting_opp_raw)
gc()

message("-> Beschuttingslagen clusteren (0.04 ha = 400 m², 100m netwerk)...")
r_binair_beschutting_opp <- terra::ifel(!is.na(heivlinder_beschutting_opp) & heivlinder_beschutting_opp > 0, 1, NA)
beschutting_clusters_opp <- cluster_filter_compleet(
  masker     = r_binair_beschutting_opp,
  opp_laag   = heivlinder_beschutting_opp,
  drempel_m2 = 400,      # 0.04 ha
  dist_m     = 100,
  werkelijk  = TRUE
)

rm(r_binair_beschutting_opp, heivlinder_beschutting_opp)
gc()

# 8. BASISHABITAT OPBOUWEN EN 500M AFSTANDSFILTER (VOORTPLANTING X FOERAGEREN) --
message("-> Basis open biotoop en aangrenzende 50m bosranden combineren...")

biotoop_opp         <- heivlinder_clusters_opp$raster
beschutting_ruw_opp <- beschutting_clusters_opp$raster

r_binair_biotoop_opp <- terra::ifel(!is.na(biotoop_opp) & biotoop_opp > 0, 1, NA)

tmp_dist_bio         <- tempfile(fileext = ".tif")
dist_tot_biotoop_opp <- terra::distance(r_binair_biotoop_opp, filename = tmp_dist_bio, overwrite = TRUE)

r_bin_besch_ruw      <- terra::ifel(!is.na(beschutting_ruw_opp) & beschutting_ruw_opp > 0, 1, NA)
bosrand_masker_opp   <- terra::ifel(r_bin_besch_ruw == 1 & dist_tot_biotoop_opp > 0 & dist_tot_biotoop_opp <= 50, 1, NA)
bosrand_direct_opp   <- terra::mask(beschutting_ruw_opp, bosrand_masker_opp)

bio_clean  <- terra::ifel(is.na(biotoop_opp), 0, biotoop_opp)
rand_clean <- terra::ifel(is.na(bosrand_direct_opp), 0, bosrand_direct_opp)

foerageer_base_opp_raw <- terra::clamp(bio_clean + rand_clean, upper = 1.0)
foerageer_base_opp     <- terra::ifel(foerageer_base_opp_raw > 0, foerageer_base_opp_raw, NA)
voortplanting_base_opp <- zand_clusters_opp$raster

unlink(tmp_dist_bio)
rm(heivlinder_clusters_opp, beschutting_clusters_opp, zand_clusters_opp, biotoop_opp, beschutting_ruw_opp, 
   r_binair_biotoop_opp, r_bin_besch_ruw, bosrand_masker_opp, bosrand_direct_opp, bio_clean, rand_clean, foerageer_base_opp_raw)
gc()

message("-> Wederzijdse afstandsfilter (500m) tussen voortplanting en foerageerhabitat...")

if (!all(is.na(suppressWarnings(terra::minmax(voortplanting_base_opp))))) {
  r_binair_vpt_base_opp  <- terra::ifel(!is.na(voortplanting_base_opp) & voortplanting_base_opp > 0, 1, NA)
  r_binair_foer_base_opp <- terra::ifel(!is.na(foerageer_base_opp) & foerageer_base_opp > 0, 1, NA)
  
  tmp_dist_vpt                <- tempfile(fileext = ".tif")
  dist_naar_voortplanting_opp <- terra::distance(r_binair_vpt_base_opp, filename = tmp_dist_vpt, overwrite = TRUE)
  zone_voortplanting_500m_opp <- terra::ifel(dist_naar_voortplanting_opp <= 500, 1, NA)
  
  tmp_dist_foer              <- tempfile(fileext = ".tif")
  dist_naar_foerageer_opp    <- terra::distance(r_binair_foer_base_opp, filename = tmp_dist_foer, overwrite = TRUE)
  zone_foerageer_500m_opp    <- terra::ifel(dist_naar_foerageer_opp <= 500, 1, NA)
  
  foerageer_finaal_opp     <- terra::mask(foerageer_base_opp, zone_voortplanting_500m_opp)
  voortplanting_finaal_opp <- terra::mask(voortplanting_base_opp, zone_foerageer_500m_opp)
  
  unlink(c(tmp_dist_vpt, tmp_dist_foer))
  rm(r_binair_vpt_base_opp, r_binair_foer_base_opp, dist_naar_voortplanting_opp, zone_voortplanting_500m_opp, 
     dist_naar_foerageer_opp, zone_foerageer_500m_opp, foerageer_base_opp, voortplanting_base_opp)
} else {
  foerageer_finaal_opp     <- template_Vlaanderen
  voortplanting_finaal_opp <- template_Vlaanderen
}

f_clean <- terra::ifel(is.na(foerageer_finaal_opp), 0, foerageer_finaal_opp)
v_clean <- terra::ifel(is.na(voortplanting_finaal_opp), 0, voortplanting_finaal_opp)

final_opp_raw <- terra::clamp(f_clean + v_clean, upper = 1.0)
final_opp     <- terra::ifel(final_opp_raw > 0, final_opp_raw, NA)

if (!all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  tmp_patch_cl <- tempfile(fileext = ".tif")
  cl_opp       <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE, filename = tmp_patch_cl, overwrite = TRUE)
} else {
  cl_opp <- template_Vlaanderen
}

rm(foerageer_finaal_opp, voortplanting_finaal_opp, f_clean, v_clean, final_opp_raw)
gc()

# 9. EXPORT NAAR DISK (ENKEL 00_ID_Rasters EN 01_Werkelijke_Oppervlaktes) -----
message("-> Geformatteerde export starten...")

base_dir <- here::here("data/output/Vlaanderen/Rasters_Soorten")

folders <- list(
  id_raster = file.path(base_dir, "00_ID_Rasters"),
  werkelijk = file.path(base_dir, "01_Werkelijke_Oppervlaktes")
)
purrr::walk(folders, ~if (!dir.exists(.x)) dir.create(.x, showWarnings = FALSE, recursive = TRUE))

# 1. Export Werkelijke Oppervlakte Raster
file_path_opp <- file.path(folders$werkelijk, paste0("Habitat_Werkelijke_Oppervlaktes_", soort, ".tif"))
terra::writeRaster(
  final_opp, 
  filename = file_path_opp, 
  overwrite = TRUE, 
  gdal = c("COMPRESS=LZW"), 
  datatype = "FLT4S", 
  NAflag = -9999
)
message(paste("    [OK] Geëxporteerd:", basename(file_path_opp)))

# 2. Export Analytisch Netwerk ID Raster
file_path_id <- file.path(folders$id_raster, paste0("ID_Netwerken_", soort, ".tif"))
terra::writeRaster(
  cl_opp, 
  filename = file_path_id, 
  overwrite = TRUE, 
  gdal = c("COMPRESS=LZW"), 
  datatype = "INT4U", 
  NAflag = 0
)
message(paste("    [OK] Analytisch ID-raster geëxporteerd:", basename(file_path_id)))

message("🏁 LEEFGEBIEDENGRID HEIVLINDER SUCCESVOL GEËXPORTEERD VOOR HEEL VLAANDEREN!")
