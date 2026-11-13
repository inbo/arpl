# ==============================================================================
# ENKEL ARPL MAKER PER SOORT (RERUN STRUKTUUR)
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

FORCE_OVERWRITE  <- FALSE  
SPECIFIEKE_SOORT <- NULL 

# Gebiedenlijst en bijhorende regiocodes voor territoria/waarnemingen
gebieden_info <- list(
  Turnhouts_Vennegebied = list(code = "TV"),
  De_Maten              = list(code = "DM"),
  Mechelse_Heide        = list(code = "MH"),
  Heesbossen            = list(code = "HB"),
  Kalmthoutse_Heide      = list(code = "KH"),
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

# --- HULPFUNCTIES ---
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

# --- HOOFDLUS ---
base_rerun_dir <- here("data/output/Rerun")
beschikbare_gebieden <- list.dirs(base_rerun_dir, full.names = FALSE, recursive = FALSE)
beschikbare_gebieden <- intersect(beschikbare_gebieden, names(gebieden_info))

message("==================================================")
message(" STARTEN ARPL MAKER (RERUN)")
if (!is.null(SPECIFIEKE_SOORT)) message(" 🎯 FILTER ACTIEF: ", SPECIFIEKE_SOORT)
if (!FORCE_OVERWRITE) message(" 🔒 OVERSCHRIJVEN UITGESCHAKELD: Bestaande bestanden worden overgeslagen.")
message("==================================================")

for (gb_naam in beschikbare_gebieden) {
  info <- gebieden_info[[gb_naam]]
  
  message("\n--------------------------------------------------")
  message(" 📍 GEBIED: ", gb_naam)
  message("--------------------------------------------------")
  
  map_potentieel <- file.path(base_rerun_dir, gb_naam, "Potentieel_Leefgebied")
  map_id         <- file.path(base_rerun_dir, gb_naam, "ID_Rasters") # Indien aanwezig in rerun, anders fallback
  out_dir_arpl   <- file.path(base_rerun_dir, gb_naam, "ARPL")
  
  if (!dir.exists(out_dir_arpl)) dir.create(out_dir_arpl, recursive = TRUE)
  if (!dir.exists(map_potentieel)) {
    message("    ℹ️ Geen map 'Potentieel_Leefgebied' gevonden voor dit gebied.")
    next
  }
  
  potentieel_tifs <- list.files(map_potentieel, pattern = "\\.tif$", full.names = TRUE)
  if (length(potentieel_tifs) == 0) {
    message("    ℹ️ Geen tif-bestanden gevonden in Potentieel_Leefgebied.")
    next
  }
  
  for (f_pot in potentieel_tifs) {
    graphics.off()
    if (requireNamespace("terra", quietly = TRUE)) {
      terra::tmpFiles(current = TRUE, orphan = TRUE, old = TRUE, remove = TRUE)
    }
    gc(verbose = FALSE)
    
    b_naam <- basename(f_pot)
    is_wintervogel <- grepl("_wv", b_naam, ignore.case = TRUE)
    
    # Soortnaam extraheren uit bestandsnaam:
    # "Potentieel_Leefgebied_<soort>_<gebied>_buffer50km.tif" of generieke variant
    soort_clean <- tolower(b_naam) %>%
      gsub("^potentieel_leefgebied_", "", .) %>%
      gsub(paste0("_", tolower(gb_naam), ".*$"), "", .) %>%
      gsub("_buffer.*$|_wv\\.tif$|\\.tif$", "", .) %>%
      trimws()
    
    if (!is.null(SPECIFIEKE_SOORT) && tolower(trimws(SPECIFIEKE_SOORT)) != soort_clean) next
    
    bestands_suffix <- if (is_wintervogel) paste0(soort_clean, "_wv") else soort_clean
    f_out_arpl <- file.path(out_dir_arpl, paste0("ARPL_", bestands_suffix, ".tif"))
    
    if (!FORCE_OVERWRITE && file.exists(f_out_arpl)) {
      message("    ⏩ Reeds verwerkt (overgeslagen): ", bestands_suffix)
      next
    }
    
    # ID-Raster opzoeken (controleert eerst in Rerun/<Gebied>/ID_Rasters, daarna fallback naar reguliere output)
    f_id <- file.path(map_id, paste0("ID_Netwerken_", bestands_suffix, ".tif"))
    if (!file.exists(f_id)) f_id <- file.path(map_id, paste0("ID_Netwerken_", soort_clean, ".tif"))
    
    # Fallback naar de klassieke outputmap mocht het ID-raster niet in de Rerun-map staan
    if (!file.exists(f_id)) {
      f_id <- here("data/output", gb_naam, "Rasters_Soorten/BWK_2025/00_ID_Rasters", paste0("ID_Netwerken_", bestands_suffix, ".tif"))
      if (!file.exists(f_id)) {
        f_id <- here("data/output", gb_naam, "Rasters_Soorten/BWK_2025/00_ID_Rasters", paste0("ID_Netwerken_", soort_clean, ".tif"))
      }
    }
    
    if (!file.exists(f_id)) {
      message("    ⚠️ ID-raster ontbreekt voor: ", bestands_suffix)
      next
    }
    
    cl_max <- tryCatch(terra::rast(f_id), error = function(e) NULL)
    if (is.null(cl_max)) {
      message("    ❌ Fout bij openen ID-raster: ", basename(f_id))
      next
    }
    
    id_col_naam <- names(cl_max)[1]
    
    # Waarnemingen ophalen
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
    
    # ARPL Berekening
    if (is_wintervogel) {
      r_pot <- terra::rast(f_pot)
      arpl_export_rast <- terra::ifel(!is.na(r_pot) & r_pot > 0, 1, NA)
    } else {
      row_dist <- df_afstanden %>% filter(Soort_clean == soort_clean)
      buffer_m <- if (nrow(row_dist) > 0) row_dist$Dispersiecap_m[1] else 500
      
      if (!is.na(buffer_m) && buffer_m >= 250000) {
        if (!is.null(v_pts) && length(v_pts) > 0) {
          r_pot <- terra::rast(f_pot)
          arpl_export_rast <- terra::ifel(!is.na(r_pot) & r_pot > 0, 1, NA)
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
    
    # Harmonisatie en export
    arpl_export_rast <- terra::extend(arpl_export_rast, master_grid, fill = NA)
    terra::writeRaster(
      arpl_export_rast,
      f_out_arpl,
      overwrite = TRUE,
      gdal = c("COMPRESS=LZW"),
      datatype = "INT1U",
      NAflag = 255
    )
    
    message("    [OK] ARPL opgeslagen: ", basename(f_out_arpl))
  }
}

message("\n==================================================")
message(" VERWERKING ARPL VOLTOOID")
message("==================================================")