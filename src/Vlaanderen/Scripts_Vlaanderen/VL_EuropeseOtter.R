# ==============================================================================
# LEEFGEBIEDENGRID VOOR EUROPESE OTTER - HEEL VLAANDEREN (WERKELIJKE OPPERVLAKTE)
# Ecologische analyse water-landkoppeling & 1000 ha landmatrix (INBO)
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
soort <- "europeseotter"

resultaat <- df %>%
  filter(tolower(trimws(Soort)) == soort) %>%
  select(Type, MinOpp_ha, AfstandBiotopen_m, Dispersiecap_m)

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

d_opp <- if("d" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["d"]] else template_Vlaanderen
h_opp <- if("h" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["h"]] else template_Vlaanderen
k_opp <- if("k" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["k"]] else template_Vlaanderen
m_opp <- if("m" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["m"]] else template_Vlaanderen
r_opp <- if("r" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["r"]] else template_Vlaanderen
s_opp <- if("s" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["s"]] else template_Vlaanderen
t_opp <- if("t" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["t"]] else template_Vlaanderen
v_opp <- if("v" %in% names(lijst_oppervlaktes)) lijst_oppervlaktes[["v"]] else template_Vlaanderen

rm(ds_vlaanderen, resultaten_gegroepeerd, df_nieuw, lijst_oppervlaktes)
gc()

# 6. INTEGRAL WATERNETWERK (HUET-ZONES + WATERVLAKKEN) ------------------------
message("-> Waterlopen (Huetzon) en watervlakken inlezen voor heel Vlaanderen...")

r_huetzon_raw <- terra::rast(here("data/input/Raster_Vlaanderen/vlaanderen_huetzon_10m.tif"))
if(!terra::compareGeom(r_huetzon_raw, template_Vlaanderen, stopOnError = FALSE)) {
  r_huetzon <- terra::resample(r_huetzon_raw, template_Vlaanderen, method = "near")
} else {
  r_huetzon <- r_huetzon_raw
}

r_waterlopen1 <- r_huetzon == 1

bl_path <- here("data/input/Raster_Vlaanderen/vlaanderen_watervlakken_2024_10m.tif")
r_bl_raw <- terra::rast(bl_path)
if(!terra::compareGeom(r_bl_raw, template_Vlaanderen, stopOnError = FALSE)) {
  r_bl <- terra::resample(r_bl_raw, template_Vlaanderen, method = "near")
} else {
  r_bl <- r_bl_raw
}

r_plassen <- r_bl == 1

otter_water_max <- (r_waterlopen1 == 1) | (r_plassen == 1)
otter_water_max <- terra::ifel(otter_water_max == 1, 1, NA)

rm(r_huetzon_raw, r_huetzon, r_waterlopen1, r_bl_raw, r_bl, r_plassen)
gc()

# 7. LANDBIOTOPEN SMELTEN EN ZEVEN OP 1000 HA ---------------------------------
message("-> Alle landbiotopen samenvoegen en zeven op 1000 ha...")

clean_opp <- function(r) {
  if (is.null(r)) return(template_Vlaanderen)
  r[is.na(r) | r == 0] <- NA
  return(r)
}

opp_stack <- c(
  clean_opp(d_opp), clean_opp(h_opp), clean_opp(k_opp), clean_opp(m_opp),
  clean_opp(r_opp), clean_opp(s_opp), clean_opp(t_opp), clean_opp(v_opp)
)

otter_biotopen_opp <- sum(opp_stack, na.rm = TRUE)
otter_biotopen_opp <- terra::clamp(otter_biotopen_opp, upper = 1.0)
r_binair_land_opp  <- terra::ifel(!is.na(otter_biotopen_opp) & otter_biotopen_opp > 0, 1, NA)

rm(opp_stack, d_opp, h_opp, k_opp, m_opp, r_opp, s_opp, t_opp, v_opp)
gc()

if (!all(is.na(suppressWarnings(terra::minmax(r_binair_land_opp))))) {
  res_land_1000ha <- cluster_filter_compleet(
    masker     = r_binair_land_opp,
    opp_laag   = otter_biotopen_opp,
    drempel_m2 = 1000 * 10000, # 1000 ha
    dist_m     = 50,
    werkelijk  = TRUE
  )
  otter_biotopen1000_opp <- res_land_1000ha$raster
  rm(res_land_1000ha)
} else {
  otter_biotopen1000_opp <- template_Vlaanderen
}

rm(otter_biotopen_opp, r_binair_land_opp)
gc()

# 8. WATERKOPPELING (< 1000M) EN HABITAT-FINALESERING -------------------------
message("-> Wederzijdse afstandsmetingen uitvoeren tussen water en land (max. 1000m)...")

r_binair_land1000_opp <- terra::ifel(!is.na(otter_biotopen1000_opp) & otter_biotopen1000_opp > 0, 1, NA)

if (!all(is.na(suppressWarnings(terra::minmax(otter_water_max)))) && !all(is.na(suppressWarnings(terra::minmax(r_binair_land1000_opp))))) {
  dist_tot_water_opp <- terra::distance(otter_water_max)
  dist_tot_land_opp  <- terra::distance(r_binair_land1000_opp)
  
  leefgebied3_opp_mask <- r_binair_land1000_opp == 1 & dist_tot_water_opp <= 1000
  leefgebied3_opp      <- terra::mask(otter_biotopen1000_opp, leefgebied3_opp_mask)
  
  leefgebied4_opp_mask <- otter_water_max == 1 & dist_tot_land_opp <= 1000
  r_water_opp_full     <- terra::ifel(otter_water_max == 1, 1.0, NA)
  leefgebied4_opp      <- terra::mask(r_water_opp_full, leefgebied4_opp_mask)
  
  l3_clean <- terra::ifel(is.na(leefgebied3_opp), 0, leefgebied3_opp)
  l4_clean <- terra::ifel(is.na(leefgebied4_opp), 0, leefgebied4_opp)
  
  otter_leefgebied_opp_raw <- terra::clamp(l3_clean + l4_clean, upper = 1.0)
  final_opp                <- terra::ifel(!is.na(otter_leefgebied_opp_raw) & otter_leefgebied_opp_raw > 0, otter_leefgebied_opp_raw, NA)
  
  rm(dist_tot_water_opp, dist_tot_land_opp, leefgebied3_opp_mask, leefgebied3_opp, 
     leefgebied4_opp_mask, r_water_opp_full, l3_clean, l4_clean, otter_leefgebied_opp_raw)
} else {
  final_opp <- template_Vlaanderen
}

if (!all(is.na(suppressWarnings(terra::minmax(final_opp))))) {
  tmp_patch_cl <- tempfile(fileext = ".tif")
  cl_opp       <- terra::patches(final_opp, directions = 8, zeroAsNA = TRUE, filename = tmp_patch_cl, overwrite = TRUE)
} else {
  cl_opp <- template_Vlaanderen
}

rm(r_binair_land1000_opp, otter_biotopen1000_opp, otter_water_max)
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

message("🏁 LEEFGEBIEDENGRID EUROPESE OTTER SUCCESVOL GEËXPORTEERD VOOR HEEL VLAANDEREN!")
