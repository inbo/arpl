library(terra)
library(sf)
library(here)
library(readxl)
library(dplyr)

# --- CONFIGURATIE ---
FORCE_OVERWRITE <- TRUE    # TRUE zodat alle bestanden met de nieuwe naam worden aangemaakt
buffer_dist     <- 10000   # 10000 = 10 km (zet op 50000 voor 50 km)
te_negeren_soorten <- c("blauwekiekendief") # Strip _wv hier niet nodig, gebeurt automatisch

map_rasters <- here("data/output/Vlaanderen/Rasters_Soorten/01_Potentieel_Leefgebied")
map_ids     <- here("data/output/Vlaanderen/Rasters_Soorten/00_ID_Rasters")
map_shapes  <- here("data/input")
base_rerun  <- here("data/output/Rerun")
pad_excel   <- here("data/input/Excel_files/Soortenlijst_Maatwerkgebieden_Gefilterd.xlsx")

gebieden_lijst <- c(
  "Turnhouts_Vennegebied",
  "De_Maten",
  "Mechelse_Heide",
  "Heesbossen",
  "Kalmthoutse_Heide",
  "Voerstreek"
)

# --- EXCEL INLEZEN ---
df_excel <- read_excel(pad_excel)
soort_kolom <- colnames(df_excel)[2]

df_excel <- df_excel %>%
  mutate(
    # 1. Verwijder eventuele _wv suffix uit de Excel-naam zodat we de zuivere basisnaam hebben
    Soort_basis = sub("_wv$", "", trimws(as.character(.data[[soort_kolom]])), ignore.case = TRUE),
    # 2. Maak schoon (geen spaties/underscores, kleine letters) voor robuuste matching
    Soort_clean = tolower(gsub("[_ ]", "", Soort_basis))
  ) %>%
  filter(!is.na(Soort_clean) & Soort_clean != "")

alle_pot_tifs <- list.files(map_rasters, pattern = "\\.tif$", full.names = TRUE)
alle_id_tifs  <- if (dir.exists(map_ids)) list.files(map_ids, pattern = "\\.tif$", full.names = TRUE) else character(0)

# Functie die gegarandeerd EERST broedvogel en DAARNA wintervogel (_wv) teruggeeft
vind_vlaanderen_tifs <- function(soort_clean_id, bestanden_vec, prefix) {
  if (length(bestanden_vec) == 0) return(character(0))
  
  # Strip prefix en .tif van de bestandsnamen op schijf
  namen_zonder_ext <- tools::file_path_sans_ext(basename(bestanden_vec))
  namen_zonder_pref <- sub(paste0("^", prefix), "", namen_zonder_ext, ignore.case = TRUE)
  
  # Check of het een wintervogel is (_wv op het einde) en bepaal de schone basisnaam
  is_wv <- grepl("_wv$", namen_zonder_pref, ignore.case = TRUE)
  basis_schijf <- tolower(gsub("[_ ]", "", sub("_wv$", "", namen_zonder_pref, ignore.case = TRUE)))
  
  # Zoek los van elkaar de broedvogel-hit en de wintervogel-hit
  bv_hit <- bestanden_vec[basis_schijf == soort_clean_id & !is_wv]
  wv_hit <- bestanden_vec[basis_schijf == soort_clean_id & is_wv]
  
  # Retourneer in strikte volgorde: eerst broedvogel (indien aanwezig), dan wintervogel (indien aanwezig)
  return(c(bv_hit, wv_hit))
}

# --- VERWERKING PER GEBIED ---
for (gebied_naam in gebieden_lijst) {
  kol_match <- colnames(df_excel)[gsub("[_ ]", "", tolower(colnames(df_excel))) == gsub("[_ ]", "", tolower(gebied_naam))][1]
  if (is.na(kol_match)) next
  
  shp_pad <- file.path(map_shapes, paste0(gebied_naam, ".shp"))
  if (!file.exists(shp_pad)) next
  
  actieve_soorten <- df_excel %>%
    filter(suppressWarnings(as.numeric(.data[[kol_match]])) == 1) %>%
    pull(Soort_clean) %>%
    unique() %>%
    setdiff(tolower(gsub("[_ ]", "", sub("_wv$", "", te_negeren_soorten, ignore.case = TRUE))))
  
  message("\n==================================================")
  message(sprintf(" 📍 GEBIED: %s (%d actieve soorten in Excel)", gebied_naam, length(actieve_soorten)))
  message("==================================================")
  
  if (length(actieve_soorten) == 0) next
  
  map_out_pot <- file.path(base_rerun, gebied_naam, "Potentieel_Leefgebied")
  map_out_id  <- file.path(base_rerun, gebied_naam, "ID_Rasters")
  dir.create(map_out_pot, showWarnings = FALSE, recursive = TRUE)
  dir.create(map_out_id,  showWarnings = FALSE, recursive = TRUE)
  
  gebied <- st_read(shp_pad, quiet = TRUE)
  
  for (soort_id in actieve_soorten) {
    tiff_paden <- vind_vlaanderen_tifs(soort_id, alle_pot_tifs, "Potentieel_Leefgebied_")
    
    if (length(tiff_paden) == 0) {
      warning(sprintf("    ⚠️ Geen broed- of wintervogelraster gevonden voor: %s", soort_id))
      next
    }
    
    # Loopt nu automatisch: 1) alleen BV, 2) alleen WV, of 3) eerst BV en dan WV
    for (tiff_pad in tiff_paden) {
      soort_suffix <- sub("^Potentieel_Leefgebied_", "", tools::file_path_sans_ext(basename(tiff_pad)), ignore.case = TRUE)
      type_label   <- if (grepl("_wv$", soort_suffix, ignore.case = TRUE)) "Wintervogel" else "Broedvogel"
      out_pad_pot  <- file.path(map_out_pot, paste0("Habitat_Werkelijke_Oppervlaktes_", soort_suffix, ".tif"))
      
      if (!FORCE_OVERWRITE && file.exists(out_pad_pot)) {
        message(sprintf("    ⏩ [%s] Reeds aanwezig (overgeslagen): %s", type_label, basename(out_pad_pot)))
        next
      }
      
      message(sprintf("    ✂️ [%s] Bijsnijden: %s -> %s", type_label, soort_suffix, basename(out_pad_pot)))
      
      r <- rast(tiff_pad)
      gebied_proj <- if (st_crs(gebied) != st_crs(r)) st_transform(gebied, st_crs(r)) else gebied
      buffer_vect <- vect(st_buffer(st_union(gebied_proj), dist = buffer_dist))
      
      r_cropped <- crop(r, buffer_vect)
      r_cropped <- ifel(!is.na(r_cropped) & r_cropped > 0, 1L, NA_integer_)
      writeRaster(r_cropped, out_pad_pot, overwrite = TRUE, gdal = c("COMPRESS=LZW"), datatype = "INT1U", NAflag = 255)
      rm(r, r_cropped)
      
      # Bijhorend ID-raster bijsnijden (werkt automatisch voor zowel BV als _wv)
      id_naam_doel <- sub("^Potentieel_Leefgebied_", "ID_Netwerken_", basename(tiff_pad), ignore.case = TRUE)
      id_hit <- alle_id_tifs[tolower(basename(alle_id_tifs)) == tolower(id_naam_doel)]
      
      if (length(id_hit) > 0 && file.exists(id_hit[1])) {
        r_id       <- rast(id_hit[1])
        r_id_crop  <- crop(r_id, buffer_vect)
        out_pad_id <- file.path(map_out_id, basename(id_hit[1]))
        writeRaster(r_id_crop, out_pad_id, overwrite = TRUE, gdal = c("COMPRESS=LZW"), datatype = "INT4U", NAflag = 0)
        rm(r_id, r_id_crop)
      }
      
      gc(verbose = FALSE)
    }
  }
}

message("\n🏁 Alle maatwerkgebieden zijn bijgewerkt!")
