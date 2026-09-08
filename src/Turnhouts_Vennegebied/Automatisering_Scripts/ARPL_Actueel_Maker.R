# ==============================================================================
# BATCH SCRIPT: Herbereken & Update ARPL + Actuele Verspreidingskaarten (DEF)
# Author: Bert Van Hecke
# ==============================================================================

library(here)
library(terra)
library(sf)
library(tidyverse)
library(readxl)

# ------------------------------------------------------------------------------
# 1. CONFIGURATIE & GEBIEDSKEUZE
# ------------------------------------------------------------------------------
GEBIED_CODE <- "TV"   # "DM", "HB", "MH", "KH", "TV", "VS"

gebieden_info <- list(
  "TV" = list(naam = "Turnhouts_Vennegebied", col_excel = "Turnhouts_Vennegebied")
)

gebied_naam <- gebieden_info[[GEBIED_CODE]]$naam
excel_kolom <- gebieden_info[[GEBIED_CODE]]$col_excel

RASTER_BASE_DIR <- here("data/output", gebied_naam, "Rasters_Soorten")
TERRITORIA_DIR  <- here("data/input/Territoria_Soorten")
EXCEL_AFSTANDEN <- here("data/input/Excel_files/Soorten_bwk_afstanden.xlsx")
EXCEL_MAATWERK  <- here("data/input/Excel_files/Soortenlijst_Maatwerkgebieden_Gefilterd.xlsx")

target_soorten <- c("zwartkopmeeuw")

dir_arpl   <- file.path(RASTER_BASE_DIR, "03_ARPL")
dir_verspr <- file.path(RASTER_BASE_DIR, "04_Actuele_Verspreiding")
dir.create(dir_arpl, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_verspr, recursive = TRUE, showWarnings = FALSE)

df_afstanden <- read_excel(EXCEL_AFSTANDEN)
df_maatwerk  <- read_excel(EXCEL_MAATWERK)

clean_key <- function(s) gsub("[^a-zA-Z0-9]", "", tolower(as.character(s)))
df_maatwerk$key <- clean_key(df_maatwerk$`Nederlandse naam`)

overzicht_df <- data.frame(
  Soort = character(),
  ARPL_ha = numeric(),
  Verspreiding_ha = numeric(),
  Status = character(),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------------------
# 2. BATCH LUS OVER ALLE SOORTEN
# ------------------------------------------------------------------------------
cat("============================================================\n")
cat(paste0("START HERBEREKENING EN UPDATE (", gebied_naam, " - ", GEBIED_CODE, ")\n"))
cat("============================================================\n\n")

for (s_naam in target_soorten) {
  
  formatted_soort <- paste0(
    toupper(substring(unlist(strsplit(s_naam, "[ _-]")), 1, 1)),
    tolower(substring(unlist(strsplit(s_naam, "[ _-]")), 2))
  ) %>% paste(collapse = "_")
  
  s_key <- clean_key(s_naam)
  
  # A. Check relevantie
  m_row <- df_maatwerk %>% filter(key == s_key)
  is_relevant <- (nrow(m_row) > 0 && excel_kolom %in% colnames(m_row) && m_row[[excel_kolom]][1] == 1)
  
  if (!is_relevant) {
    message(paste0("⏭️ Overgeslagen (Niet relevant voor ", GEBIED_CODE, "): ", formatted_soort))
    next
  }
  
  message(paste0("🔵 Verwerken & Updaten: ", formatted_soort))
  
  arpl_out_path   <- file.path(dir_arpl, paste0("Habitat_ARPL_", formatted_soort, ".tif"))
  verspr_out_path <- file.path(dir_verspr, paste0("Habitat_Actuele_Verspreiding_", formatted_soort, ".tif"))
  
  # B. Parameters ophalen
  s_info <- df_afstanden %>% filter(clean_key(Soort) == s_key) %>% slice(1)
  min_opp_ha  <- if(nrow(s_info) > 0) s_info$MinOpp_ha[1] else 0
  dispersie_m <- if(nrow(s_info) > 0) s_info$Dispersiecap_m[1] else 5000
  
  # C. Werkelijk oppervlakteraster inlezen
  r_werk_path <- file.path(RASTER_BASE_DIR, "02_Werkelijke_Oppervlakte", paste0("Habitat_Werkelijke_Oppervlaktes_", formatted_soort, ".tif"))
  
  if (!file.exists(r_werk_path)) {
    message(paste0("   └─ ⚠️ Werkelijk oppervlakteraster ontbreekt in '02_Werkelijke_Oppervlakte'."))
    overzicht_df <- rbind(overzicht_df, data.frame(Soort = formatted_soort, ARPL_ha = 0, Verspreiding_ha = 0, Status = "Raster ontbreekt"))
    next
  }
  
  r_werkelijk <- rast(r_werk_path)
  
  if (all(is.na(values(r_werkelijk, mat = FALSE)))) {
    message("   └─ ℹ️ Geen werkelijk habitat aanwezig in het studiegebied.")
    writeRaster(r_werkelijk, arpl_out_path, overwrite = TRUE, gdal = c("COMPRESS=LZW"), datatype = "INT1U")
    writeRaster(r_werkelijk, verspr_out_path, overwrite = TRUE, gdal = c("COMPRESS=LZW"), datatype = "INT1U")
    overzicht_df <- rbind(overzicht_df, data.frame(Soort = formatted_soort, ARPL_ha = 0, Verspreiding_ha = 0, Status = "Geen habitat"))
    next
  }
  
  # D. Territoria inlezen uit de Python export
  csv_path <- file.path(TERRITORIA_DIR, paste0(GEBIED_CODE, "_Territoria_", formatted_soort, ".csv"))
  has_obs  <- FALSE
  
  if (file.exists(csv_path)) {
    df_pts <- read_csv(csv_path, show_col_types = FALSE)
    if (nrow(df_pts) > 0) {
      colnames(df_pts) <- tolower(colnames(df_pts))
      if ("x" %in% colnames(df_pts) && "y" %in% colnames(df_pts)) {
        sf_punten   <- st_as_sf(df_pts, coords = c("x", "y"), crs = 31370)
        vect_punten <- vect(sf_punten)
        
        # FIX 1: Projecteer de punten zuiver naar de exacte CRS van het raster
        vect_punten <- terra::project(vect_punten, crs(r_werkelijk))
        has_obs     <- TRUE
      }
    }
  }
  
  # E. BEREKENING (ULTRASNEL EN RAM-VEILIG VIA TERRA)
  if (has_obs) {
    
    cl_patches <- terra::patches(r_werkelijk, directions = 8, zeroAsNA = TRUE)
    
    stats_p <- terra::zonal(r_werkelijk, cl_patches, fun = "sum", na.rm = TRUE)
    colnames(stats_p) <- c("ID", "Werk_ha")
    stats_p$Werk_ha <- stats_p$Werk_ha * 0.01
    
    valid_ids <- stats_p$ID[stats_p$Werk_ha >= min_opp_ha]
    
    if (length(valid_ids) > 0) {
      
      # 1. ACTUELE VERSPREIDING (Punt valt direct op de cluster)
      ext_pts <- terra::extract(cl_patches, vect_punten, ID = FALSE)
      ids_verspr <- if(!is.null(ext_pts) && nrow(ext_pts) > 0) unique(na.omit(ext_pts[[1]])) else c()
      ids_verspr_ok <- intersect(valid_ids, ids_verspr)
      
      if (length(ids_verspr_ok) > 0) {
        m_verspr     <- cl_patches %in% ids_verspr_ok
        r_verspr_raw <- terra::mask(r_werkelijk, m_verspr)
        ha_verspr    <- sum(stats_p$Werk_ha[stats_p$ID %in% ids_verspr_ok], na.rm = TRUE)
        r_verspr_export <- terra::ifel(!is.na(r_verspr_raw) & r_verspr_raw > 0, 1, NA)
      } else {
        r_verspr_export <- r_werkelijk * NA
        ha_verspr <- 0
      }
      
      # 2. ARPL BEREKENING (Via terra::buffer op punten -> voorkomt RAM crash!)
      vect_smear <- terra::buffer(vect_punten, width = dispersie_m)
      ext_arpl   <- terra::extract(cl_patches, vect_smear, ID = FALSE)
      ids_arpl   <- if(!is.null(ext_arpl) && nrow(ext_arpl) > 0) unique(na.omit(ext_arpl[[1]])) else c()
      ids_arpl_ok <- intersect(valid_ids, ids_arpl)
      
      if (length(ids_arpl_ok) > 0) {
        m_arpl        <- cl_patches %in% ids_arpl_ok
        r_arpl_raw    <- terra::mask(r_werkelijk, m_arpl)
        ha_arpl       <- sum(stats_p$Werk_ha[stats_p$ID %in% ids_arpl_ok], na.rm = TRUE)
        r_arpl_export <- terra::ifel(!is.na(r_arpl_raw) & r_arpl_raw > 0, 1, NA)
      } else {
        r_arpl_export <- r_werkelijk * NA
        ha_arpl <- 0
      }
      
      rm(vect_smear)
      
    } else {
      r_arpl_export   <- r_werkelijk * NA
      r_verspr_export <- r_werkelijk * NA
      ha_arpl   <- 0
      ha_verspr <- 0
    }
    
    rm(cl_patches)
    
  } else {
    r_arpl_export   <- r_werkelijk * NA
    r_verspr_export <- r_werkelijk * NA
    ha_arpl   <- 0
    ha_verspr <- 0
  }
  
  # F. Overschrijf de bestaande rasters
  writeRaster(r_arpl_export, arpl_out_path, overwrite = TRUE, gdal = c("COMPRESS=LZW"), datatype = "INT1U")
  writeRaster(r_verspr_export, verspr_out_path, overwrite = TRUE, gdal = c("COMPRESS=LZW"), datatype = "INT1U")
  
  status_msg <- if(has_obs) "Bijgewerkt (Nieuwe territoria)" else "Geen territoria (broedcode >=4)"
  message(paste0("   ├─ ARPL: ", round(ha_arpl, 2), " ha"))
  message(paste0("   └─ Verspreiding: ", round(ha_verspr, 2), " ha"))
  
  overzicht_df <- rbind(overzicht_df, data.frame(
    Soort = formatted_soort, ARPL_ha = round(ha_arpl, 2), Verspreiding_ha = round(ha_verspr, 2), Status = status_msg
  ))
  
  gc(verbose = FALSE)
}

# ------------------------------------------------------------------------------
# 3. EINDOVERZICHT EXPORTEREN
# ------------------------------------------------------------------------------
cat("\n============================================================\n")
cat("EINDOVERZICHT AANMAKEN\n")
cat("============================================================\n")

print(overzicht_df)

csv_uitvoer <- file.path(RASTER_BASE_DIR, paste0("Overzicht_ARPL_Verspreiding_", GEBIED_CODE, ".csv"))
write_csv(overzicht_df, csv_uitvoer)
cat(paste0("\n📊 Bijgewerkt overzicht opgeslagen in: ", csv_uitvoer, "\n"))
