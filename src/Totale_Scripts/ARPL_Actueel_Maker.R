# ==============================================================================
# ARPL & ACTUEEL MAKER PER SOORT (MET AUTOMATISCHE OVERSLAG VAN BESTAANDE OUTPUT)
# ==============================================================================

library(here)
library(dplyr)
library(readxl)
library(readr)
library(purrr)
library(sf)
library(terra)
library(data.table)

terra::terraOptions(memfrac = 0.2, tempdir = tempdir(), verbose = FALSE)

# ------------------------------------------------------------------------------
# INSTELLINGEN:
#   - FORCE_OVERWRITE: Zet op FALSE zodat bestaande TIFs NETJES WORDEN OVERSLAGEN.
#                      Zet op TRUE als je toch een herberekening wilt afdwingen.
#   - SPECIFIEKE_SOORT: Vul een soortnaam in (bijv. "bruinekiekendief") of zet op NULL voor ALLE soorten.
# ------------------------------------------------------------------------------
FORCE_OVERWRITE  <- FALSE  
SPECIFIEKE_SOORT <- NULL 

DEFAULT_SCENARIO_SELECTIE <- list(
  Turnhouts_Vennegebied = "TV_Scenario_BWK_2025.rds"
)

actieve_scenarios <- if (exists("SCENARIO_SELECTIE") && is.list(SCENARIO_SELECTIE)) SCENARIO_SELECTIE else DEFAULT_SCENARIO_SELECTIE

gebieden_info <- list(
  De_Maten              = list(code = "DM"),
  Heesbossen            = list(code = "HB"),
  Kalmthoutse_Heide      = list(code = "KH"),
  Mechelse_Heide        = list(code = "MH"),
  Turnhouts_Vennegebied = list(code = "TV"),
  Voerstreek            = list(code = "VS")
)

broedvogels_lijst <- c(
  "blauwborst", "boomleeuwerik", "boompieper", "bruinekiekendief", "fluiter",
  "grauweklauwier", "grutto", "ijsvogel", "kwak", "kwartelkoning", "matkop",
  "middelstebontespecht", "nachtegaal", "nachtzwaluw", "paapje", "porseleinhoen", 
  "roerdomp", "tapuit", "watersnip", "wespendief", "wielewaal", "woudaap", 
  "wulp", "zomertortel", "zwartespecht", "zwartkopmeeuw"
)

master_grid_pad <- here("data/input/Raster_Vlaanderen/Vlaanderen_MasterGrid_10m.tif")
master_grid      <- terra::rast(master_grid_pad)[[1]]

df_afstanden <- read_excel(here("data/input/Excel_files/Soorten_bwk_afstanden.xlsx")) %>% 
  mutate(Soort_clean = tolower(gsub(" ", "", trimws(Soort))))

# HULPFUNCTIE: Zoekt flexibel naar bestanden (negeert hoofdletters, spaties en underscores)
vind_waarnemingen_bestand <- function(map_pad, patroon_zoek) {
  if (!dir.exists(map_pad)) return(NULL)
  alle_bestanden <- list.files(map_pad, full.names = TRUE)
  
  patroon_clean   <- gsub("[_ ]", "", tolower(patroon_zoek))
  bestanden_clean <- gsub("[_ ]", "", tolower(basename(alle_bestanden)))
  
  gematcht <- alle_bestanden[bestanden_clean == patroon_clean]
  if (length(gematcht) > 0) return(gematcht[1])
  return(NULL)
}

laad_waarnemingen_punten <- function(pad, is_komma) {
  if (is.null(pad) || !file.exists(pad)) return(NULL)
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
    message("    ⚠️ Fout bij inlezen CSV: ", basename(pad), " (", e$message, ")")
    return(NULL)
  })
}

message("==================================================")
message(" STARTEN ARPL & ACTUEEL MAKER PER SOORT")
if (!is.null(SPECIFIEKE_SOORT)) {
  message(" 🎯 FILTER ACTIEF: Enkel verwerken voor -> ", SPECIFIEKE_SOORT)
}
if (!FORCE_OVERWRITE) {
  message(" 🔒 OVERSCHRIJVEN UITGESCHAKELD: Bestaande bestanden worden overgeslagen.")
}
message("==================================================")

for (gb_naam in names(actieve_scenarios)) {
  rds_naam <- actieve_scenarios[[gb_naam]]
  info     <- gebieden_info[[gb_naam]]
  if (is.null(info)) info <- list(code = "TV")
  
  rds_pad <- here("data/input/Scenario_rds", rds_naam)
  huidig_scenario <- gsub("^.*_Scenario_|^Scenario_|_wv\\.rds$|\\.rds$", "", basename(rds_pad), ignore.case = TRUE)
  
  base_rasters_dir <- here("data/output", gb_naam, "Rasters_Soorten", huidig_scenario)
  
  map_id        <- file.path(base_rasters_dir, "00_ID_Rasters")
  map_werkelijk <- file.path(base_rasters_dir, "02_Werkelijke_Oppervlaktes")
  out_dir_arpl  <- file.path(base_rasters_dir, "03_ARPL")
  out_dir_act   <- file.path(base_rasters_dir, "04_Actuele_Verspreiding")
  
  if (!dir.exists(out_dir_arpl)) dir.create(out_dir_arpl, recursive = TRUE)
  if (!dir.exists(out_dir_act))  dir.create(out_dir_act,  recursive = TRUE)
  
  werkelijk_tifs <- list.files(map_werkelijk, pattern = "\\.tif$", full.names = TRUE)
  if (length(werkelijk_tifs) == 0) next
  
  for (f_werkelijk in werkelijk_tifs) {
    graphics.off()
    if (requireNamespace("terra", quietly = TRUE)) {
      terra::tmpFiles(current = TRUE, orphan = TRUE, old = TRUE, remove = TRUE)
    }
    gc(verbose = FALSE)
    
    is_wintervogel <- grepl("_wv\\.tif$", f_werkelijk, ignore.case = TRUE)
    soort_clean <- basename(f_werkelijk) %>% 
      tolower() %>% 
      gsub("^habitat_werkelijke_oppervlaktes_|_wv\\.tif$|\\.tif$", "", .) %>% 
      trimws()
    
    # --- 1. FILTER OP SPECIFIEKE SOORT ---
    if (!is.null(SPECIFIEKE_SOORT) && tolower(trimws(SPECIFIEKE_SOORT)) != soort_clean) {
      next
    }
    
    bestands_suffix <- if (is_wintervogel) paste0(soort_clean, "_wv") else soort_clean
    
    f_out_arpl <- file.path(out_dir_arpl, paste0("Habitat_ARPL_", bestands_suffix, ".tif"))
    f_out_act  <- file.path(out_dir_act,  paste0("Habitat_Actuele_Verspreiding_", bestands_suffix, ".tif"))
    
    # --- 2. CONTROLE OF BESTANDEN AL BESTAAN ---
    if (!FORCE_OVERWRITE && file.exists(f_out_arpl) && file.exists(f_out_act)) {
      message("    ⏩ Reeds verwerkt (overgeslagen): ", bestands_suffix)
      next
    }
    
    f_id <- file.path(map_id, paste0("ID_Netwerken_", bestands_suffix, ".tif"))
    if (!file.exists(f_id)) f_id <- file.path(map_id, paste0("ID_Netwerken_", soort_clean, ".tif"))
    
    if (!file.exists(f_id)) {
      message("    ⚠️ ID-raster ontbreekt voor: ", bestands_suffix)
      next
    }
    
    cl_max <- terra::rast(f_id)
    id_col_naam <- names(cl_max)[1]
    
    # --- LOCATIE WAARNEMINGEN BEPALEN ---
    is_broedvogel <- soort_clean %in% broedvogels_lijst
    
    if (is_broedvogel) {
      map_zoek <- here("data/input/Territoria_Soorten")
      bestand_zoek <- paste0(info$code, "_Territoria_", soort_clean, ".csv")
      is_csv_komma <- TRUE
    } else {
      map_zoek <- here("data/input/Waarnemingen_Soorten", gb_naam)
      bestand_zoek <- paste0("Waarnemingen_", soort_clean, ".csv")
      is_csv_komma <- FALSE
    }
    
    waarnemingen_path <- vind_waarnemingen_bestand(map_zoek, bestand_zoek)
    v_pts <- laad_waarnemingen_punten(waarnemingen_path, is_csv_komma)
    
    # --- 1. ARPL BEREKENEN ---
    if (is_wintervogel) {
      r_werkelijk <- terra::rast(f_werkelijk)
      arpl_export_rast <- terra::ifel(!is.na(r_werkelijk) & r_werkelijk > 0, 1, NA)
    } else {
      row_dist <- df_afstanden %>% filter(Soort_clean == soort_clean)
      buffer_m <- if (nrow(row_dist) > 0) row_dist$Dispersiecap_m[1] else 500
      
      if (!is.na(buffer_m) && buffer_m >= 250000) {
        if (!is.null(v_pts) && length(v_pts) > 0) {
          r_werkelijk <- terra::rast(f_werkelijk)
          arpl_export_rast <- terra::ifel(!is.na(r_werkelijk) & r_werkelijk > 0, 1, NA)
        } else {
          arpl_export_rast <- terra::rast(master_grid, vals = NA)
        }
      } else {
        if (!is.null(v_pts) && !all(is.na(suppressWarnings(terra::minmax(cl_max))))) {
          v_blobs <- tryCatch(terra::aggregate(terra::buffer(v_pts, width = buffer_m)), error = function(e) NULL)
          
          if (!is.null(v_blobs)) {
            ext_blobs <- terra::extract(cl_max, v_blobs, ID = FALSE)
            ids_arpl  <- if (!is.null(ext_blobs) && nrow(ext_blobs) > 0) unique(na.omit(ext_blobs[[id_col_naam]])) else c()
            
            if (length(ids_arpl) > 0) {
              arpl_export_rast <- terra::ifel(cl_max %in% ids_arpl, 1, NA)
            } else {
              arpl_export_rast <- terra::rast(master_grid, vals = NA)
            }
          } else {
            arpl_export_rast <- terra::rast(master_grid, vals = NA)
          }
        } else {
          arpl_export_rast <- terra::rast(master_grid, vals = NA)
        }
      }
    }
    
    # --- 2. ACTUEEL BEREKENEN ---
    if (!is.null(v_pts) && !all(is.na(suppressWarnings(terra::minmax(cl_max))))) {
      ext_points  <- terra::extract(cl_max, v_pts, ID = FALSE)
      ids_actueel <- if (!is.null(ext_points) && nrow(ext_points) > 0) unique(na.omit(ext_points[[id_col_naam]])) else c()
      
      if (length(ids_actueel) > 0) {
        actueel_export_rast <- terra::ifel(cl_max %in% ids_actueel, 1, NA)
      } else {
        actueel_export_rast <- terra::rast(master_grid, vals = NA)
      }
    } else {
      actueel_export_rast <- terra::rast(master_grid, vals = NA)
    }
    
    # --- 3. EXPORT ---
    arpl_export_rast    <- terra::extend(arpl_export_rast, master_grid, fill = NA)
    actueel_export_rast <- terra::extend(actueel_export_rast, master_grid, fill = NA)
    
    terra::writeRaster(arpl_export_rast,    f_out_arpl, overwrite = TRUE, gdal = c("COMPRESS=LZW"), datatype = "INT1U", NAflag = 255)
    terra::writeRaster(actueel_export_rast, f_out_act,  overwrite = TRUE, gdal = c("COMPRESS=LZW"), datatype = "INT1U", NAflag = 255)
    
    message("    [OK] Herberekend en geëxporteerd voor: ", bestands_suffix)
  }
}
