# ==============================================================================
# ARPL & ACTUEEL MAKER PER SOORT (INCL. CRASH-PREVENTIE & EXTRA ROBUUST)
# AUTEUR: Bert Van Hecke
# ==============================================================================

library(here)
library(dplyr)
library(readxl)
library(readr)
library(purrr)
library(sf)
library(terra)
library(data.table)

# Dwing terra tot strikt geheugenbeheer
terra::terraOptions(memfrac = 0.2, tempdir = tempdir(), verbose = FALSE)

# ------------------------------------------------------------------------------
# 0. INSTELLINGEN & REGIONALE CONFIGURATIE
# ------------------------------------------------------------------------------
DEFAULT_SCENARIO_SELECTIE <- list(
  Turnhouts_Vennegebied = "TV_Scenario_BWK_2025.rds"
)

if (!exists("SCENARIO_SELECTIE") || !is.list(SCENARIO_SELECTIE)) {
  actieve_scenarios <- DEFAULT_SCENARIO_SELECTIE
} else {
  actieve_scenarios <- SCENARIO_SELECTIE
}

gebieden_info <- list(
  De_Maten             = list(code = "DM"),
  Heesbossen           = list(code = "HB"),
  Kalmthoutse_Heide     = list(code = "KH"),
  Mechelse_Heide       = list(code = "MH"),
  Turnhouts_Vennegebied = list(code = "TV"),
  Voerstreek           = list(code = "VS")
)

broedvogels_lijst <- c(
  "blauwborst", "boomleeuwerik", "boompieper", "bruinekiekendief", "fluiter",
  "grauweklauwier", "grutto", "ijsvogel", "kwak", "kwartelkoning", "matkop",
  "middelstebontespecht", "nachtegaal", "nachtzwaluw", "paapje", "porseleinhoen", 
  "roerdomp", "tapuit", "watersnip", "wespendief", "wielewaal", "woudaap", 
  "wulp", "zomertortel", "zwartespecht", "zwartkopmeeuw"
)

master_grid_pad <- here("data/input/Raster_Vlaanderen/Vlaanderen_MasterGrid_10m.tif")
master_grid     <- terra::rast(master_grid_pad)[[1]]

df_afstanden <- read_excel(here("data/input/Excel_files/Soorten_bwk_afstanden.xlsx")) %>% 
  mutate(Soort_clean = tolower(gsub(" ", "", trimws(Soort))))

message("==================================================")
message(" STARTEN ARPL & ACTUEEL MAKER PER SOORT")
message("==================================================")

for (gb_naam in names(actieve_scenarios)) {
  rds_naam <- actieve_scenarios[[gb_naam]]
  info     <- gebieden_info[[gb_naam]]
  if (is.null(info)) info <- list(code = "TV")
  
  rds_pad <- here("data/input/Scenario_rds", rds_naam)
  huidig_scenario <- gsub("^.*_Scenario_|^Scenario_|_wv\\.rds$|\\.rds$", "", basename(rds_pad), ignore.case = TRUE)
  
  message("\n--------------------------------------------------")
  message("-> Verwerken gebied  : ", gb_naam, " (Code: ", info$code, ")")
  message("   Scenario Naam     : ", huidig_scenario)
  message("--------------------------------------------------")
  
  base_rasters_dir <- here("data/output", gb_naam, "Rasters_Soorten", huidig_scenario)
  
  map_id        <- file.path(base_rasters_dir, "00_ID_Rasters")
  map_werkelijk <- file.path(base_rasters_dir, "02_Werkelijke_Oppervlaktes")
  out_dir_arpl  <- file.path(base_rasters_dir, "03_ARPL")
  out_dir_act   <- file.path(base_rasters_dir, "04_Actuele_Verspreiding")
  
  if (!dir.exists(out_dir_arpl)) dir.create(out_dir_arpl, recursive = TRUE)
  if (!dir.exists(out_dir_act))  dir.create(out_dir_act,  recursive = TRUE)
  
  werkelijk_tifs <- list.files(map_werkelijk, pattern = "\\.tif$", full.names = TRUE)
  
  if (length(werkelijk_tifs) == 0) {
    warning("⚠️ Geen werkelijke oppervlakte TIFs gevonden in: ", map_werkelijk)
    next
  }
  
  message("   Aantal te verwerken soortrasters: ", length(werkelijk_tifs))
  
  for (f_werkelijk in werkelijk_tifs) {
    
    # Voorkom dat RStudio's WebView2 volloopt door plots/grafische apparaten te sluiten
    graphics.off()
    
    # Geheugen opruimen direct aan het begin van elke iteratie
    if (requireNamespace("terra", quietly = TRUE)) {
      terra::tmpFiles(current = TRUE, orphan = TRUE, old = TRUE, remove = TRUE)
    }
    gc(verbose = FALSE)
    
    is_wintervogel <- grepl("_wv\\.tif$", f_werkelijk, ignore.case = TRUE)
    
    soort_clean <- basename(f_werkelijk) %>% 
      tolower() %>% 
      gsub("^habitat_werkelijke_oppervlaktes_|_wv\\.tif$|\\.tif$", "", .) %>% 
      trimws()
    
    bestands_suffix <- if (is_wintervogel) paste0(soort_clean, "_wv") else soort_clean
    
    f_out_arpl <- file.path(out_dir_arpl, paste0("Habitat_ARPL_", bestands_suffix, ".tif"))
    f_out_act  <- file.path(out_dir_act,  paste0("Habitat_Actuele_Verspreiding_", bestands_suffix, ".tif"))
    
    # Snel overslaan als BEIDE al op schijf staan
    if (file.exists(f_out_arpl) && file.exists(f_out_act)) {
      message("   ⏩ Reeds verwerkt (overgeslagen): ", bestands_suffix)
      next
    }
    
    f_id <- file.path(map_id, paste0("ID_Netwerken_", bestands_suffix, ".tif"))
    if (!file.exists(f_id)) {
      f_id <- file.path(map_id, paste0("ID_Netwerken_", soort_clean, ".tif"))
    }
    
    if (!file.exists(f_id)) {
      message("   ⚠️ ID-raster ontbreekt voor: ", bestands_suffix, " -> Overgeslagen")
      next
    }
    
    cl_max <- terra::rast(f_id)
    
    # --------------------------------------------------------------------------
    # HULPFUNCTIE VOOR ROBUUST INLEZEN VOOR WAARNEMINGEN/TERRITORIA
    # --------------------------------------------------------------------------
    laad_waarnemingen_punten <- function(pad, is_komma) {
      if (!file.exists(pad)) return(NULL)
      tryCatch({
        raw_df <- if (is_komma) {
          readr::read_csv(pad, show_col_types = FALSE)
        } else {
          readr::read_delim(pad, delim = ";", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)
        }
        colnames(raw_df) <- tolower(colnames(raw_df))
        
        if (all(c("x", "y") %in% colnames(raw_df)) && nrow(raw_df) > 0) {
          clean_df <- raw_df[!is.na(raw_df$x) & !is.na(raw_df$y), ]
          if (nrow(clean_df) > 0) {
            return(terra::vect(as.matrix(clean_df[, c("x", "y")]), type = "points", crs = "EPSG:31370"))
          }
        }
        return(NULL)
      }, error = function(e) {
        message("   ⚠️ Fout bij inlezen waarnemingen/territoria file: ", basename(pad), " (", e$message, ")")
        return(NULL)
      })
    }
    
    # --------------------------------------------------------------------------
    # 1. ARPL BEREKENEN
    # --------------------------------------------------------------------------
    if (is_wintervogel) {
      message("   ❄️ Wintervogel gedetecteerd (", bestands_suffix, "): Werkelijke Oppervlakte 1-op-1 overgenomen als ARPL.")
      r_werkelijk <- terra::rast(f_werkelijk)
      arpl_export_rast <- terra::ifel(!is.na(r_werkelijk) & r_werkelijk > 0, 1, NA)
      rm(r_werkelijk)
    } else {
      row_dist <- df_afstanden %>% filter(Soort_clean == soort_clean)
      buffer_m <- if (nrow(row_dist) > 0) row_dist$Dispersiecap_m[1] else 500
      
      is_broedvogel <- soort_clean %in% broedvogels_lijst
      soort_format  <- paste0(toupper(substr(soort_clean, 1, 1)), substr(soort_clean, 2, nchar(soort_clean)))
      
      if (is_broedvogel) {
        waarnemingen_path <- here("data/input/Territoria_Soorten", paste0(info$code, "_Territoria_", soort_format, ".csv"))
        is_csv_komma <- TRUE
      } else {
        waarnemingen_path <- here("data/input/Waarnemingen_Soorten", gb_naam, paste0("Waarnemingen_", soort_format, ".csv"))
        is_csv_komma <- FALSE
      }
      
      v_pts <- laad_waarnemingen_punten(waarnemingen_path, is_csv_komma)
      
      # UNIEKE REGEL VOOR BUFFER >= 250.000m (bijv. Porseleinhoen)
      if (!is.na(buffer_m) && buffer_m >= 250000) {
        if (!is.null(v_pts) && length(v_pts) > 0) {
          message("   💡 Buffer >= 250km met waarnemingen gedetecteerd (", bestands_suffix, "): Werkelijke Oppervlakte 1-op-1 overgenomen als ARPL.")
          r_werkelijk <- terra::rast(f_werkelijk)
          arpl_export_rast <- terra::ifel(!is.na(r_werkelijk) & r_werkelijk > 0, 1, NA)
          rm(r_werkelijk)
        } else {
          message("   ⚠️ Buffer >= 250km maar GEEN waarnemingen voor (", bestands_suffix, "): ARPL wordt leeg ingesteld.")
          arpl_export_rast <- terra::rast(master_grid, vals = NA)
        }
        if (!is.null(v_pts)) rm(v_pts)
      } else {
        # STANDAARD BUFFER BEREKENING VOOR NORMALE AFSTANDEN
        v_blobs <- NULL
        
        if (!is.null(v_pts)) {
          tryCatch({
            v_blobs <- terra::buffer(v_pts, width = buffer_m)
            v_blobs <- terra::aggregate(v_blobs)
          }, error = function(e) {
            v_blobs <<- NULL
          })
          rm(v_pts)
        }
        
        if (!is.null(v_blobs) && !all(is.na(suppressWarnings(terra::minmax(cl_max))))) {
          ext_blobs <- terra::extract(cl_max, v_blobs, ID = FALSE)
          ids_arpl  <- if (!is.null(ext_blobs) && nrow(ext_blobs) > 0) unique(na.omit(ext_blobs[[1]])) else c()
          
          if (length(ids_arpl) > 0) {
            arpl_export_rast <- cl_max %in% ids_arpl
            arpl_export_rast <- terra::ifel(arpl_export_rast == 1, 1, NA)
          } else {
            arpl_export_rast <- terra::rast(master_grid, vals = NA)
          }
          rm(v_blobs)
        } else {
          arpl_export_rast <- terra::rast(master_grid, vals = NA)
        }
      }
    }
    
    # --------------------------------------------------------------------------
    # 2. ACTUELE VERSPREIDING BEREKENEN
    # --------------------------------------------------------------------------
    is_broedvogel <- soort_clean %in% broedvogels_lijst
    soort_format  <- paste0(toupper(substr(soort_clean, 1, 1)), substr(soort_clean, 2, nchar(soort_clean)))
    
    if (is_broedvogel) {
      waarnemingen_path <- here("data/input/Territoria_Soorten", paste0(info$code, "_Territoria_", soort_format, ".csv"))
      is_csv_komma <- TRUE
    } else {
      waarnemingen_path <- here("data/input/Waarnemingen_Soorten", gb_naam, paste0("Waarnemingen_", soort_format, ".csv"))
      is_csv_komma <- FALSE
    }
    
    v_points <- laad_waarnemingen_punten(waarnemingen_path, is_csv_komma)
    
    if (!is.null(v_points) && !all(is.na(suppressWarnings(terra::minmax(cl_max))))) {
      ext_points  <- terra::extract(cl_max, v_points, ID = FALSE)
      ids_actueel <- if (!is.null(ext_points) && nrow(ext_points) > 0) unique(na.omit(ext_points[[1]])) else c()
      
      if (length(ids_actueel) > 0) {
        actueel_export_rast <- cl_max %in% ids_actueel
        actueel_export_rast <- terra::ifel(actueel_export_rast == 1, 1, NA)
      } else {
        actueel_export_rast <- terra::rast(master_grid, vals = NA)
      }
      rm(v_points)
    } else {
      actueel_export_rast <- terra::rast(master_grid, vals = NA)
    }
    
    # --------------------------------------------------------------------------
    # 3. EXPORTEN DIRECT NAAR SCHIJF
    # --------------------------------------------------------------------------
    arpl_export_rast    <- terra::extend(arpl_export_rast, master_grid, fill = NA)
    actueel_export_rast <- terra::extend(actueel_export_rast, master_grid, fill = NA)
    
    terra::writeRaster(arpl_export_rast,    f_out_arpl, overwrite = TRUE, gdal = c("COMPRESS=LZW"), datatype = "INT1U", NAflag = 255)
    terra::writeRaster(actueel_export_rast, f_out_act,  overwrite = TRUE, gdal = c("COMPRESS=LZW"), datatype = "INT1U", NAflag = 255)
    
    message("   [OK] Geëxporteerd voor: ", bestands_suffix)
    
    # Grondige opruiming
    rm(cl_max, arpl_export_rast, actueel_export_rast)
    if (requireNamespace("terra", quietly = TRUE)) {
      terra::tmpFiles(current = TRUE, orphan = TRUE, old = TRUE, remove = TRUE)
    }
    gc(verbose = FALSE)
  }
}

message("\n==================================================")
message(" 🎉 FASE 2 VOLLEDIG AFGEROND INCLUSIEF WINTERVOGEL SNELKOPPELING & HIGH-BUFFER BYPASS!")
message("==================================================")
