# ==============================================================================
# BENCHMARK: TESTMETHODES VOOR NETWERKCLUSTERING (ARGUSVLINDER DATA)
# ==============================================================================

library(here)
library(terra)
library(microbenchmark)
library(dplyr)

# 1. Omgeving & setup inladen
source(here::here("src", "Functies", "00_global_setup.R"))
source(here::here("src", "Functies", "extract_bwk_laag.R"))

master_grid   <- terra::rast(PATHS$master_grid)[[1]]
master_grid <- master_grid * NA

ds_parquet    <- arrow::open_dataset(PATHS$bwk_parquet)
df_bwk_codes  <- data.table::fread(PATHS$bwk_codes_csv, data.table = FALSE)

# 2. Haal de echte fracties voor argusvlinder op
message("1. Data ophalen voor argusvlinder...")
r_opp_compleet <- extract_bwk_laag(
  soort_naam    = "argusvlinder",
  type_naam     = "bwk",
  ds_parquet    = ds_parquet,
  df_bwk_codes  = df_bwk_codes,
  master_grid = master_grid
)

# 3. Testgebied selecteren (bijv. 1% van de oppervlakte / een uitsnede met echte data)
# We zoeken de bounding box rond een gebied met veel data (bijv. 10 x 10 km = 10.000 ha)
message("2. Testuitsnede maken (1% steekproef)...")

# Zoek het zwaartepunt van actieve cellen
actieve_cellen <- which(!is.na(values(r_opp_compleet, mat = FALSE)))
if (length(actieve_cellen) == 0) stop("Geen data gevonden om te benchmarken!")

midden_cel <- actieve_cellen[round(length(actieve_cellen) / 2)]
xy_midden  <- xyFromCell(r_opp_compleet, midden_cel)

# Uitsnede van 10 km x 10 km (1.000 x 1.000 cellen = 100 km² ≈ 0,7% van Vlaanderen)
# Vergroot 'straal_test' naar 15000 of 20000 voor 2% of 4%
# ==============================================================================
# TEST OP 16 MILJOEN CELLEN (~12% van Vlaanderen, 40 km x 40 km)
# ==============================================================================

straal_test <- 20000  # 20 km naar elke zijde = 40x40 km uitsnede

test_extent <- ext(xy_midden[1] - straal_test, xy_midden[1] + straal_test,
                   xy_midden[2] - straal_test, xy_midden[2] + straal_test)

r_test <- crop(r_opp_compleet, test_extent)
r_bin  <- ifel(!is.na(r_test) & r_test > 0, 1, NA)

n_cellen_actief <- global(!is.na(r_bin), "sum", na.rm = TRUE)[[1]]
message(paste0("   Test-raster cellen: ", ncell(r_test), 
               " | Actieve biotoopcellen: ", n_cellen_actief, 
               " (", round(n_cellen_actief * 0.01, 1), " ha)"))

# Netwerkparameters
dist_m       <- 100
straal_m     <- dist_m / 2
stral_cellen <- ceiling(straal_m / 10)

# ==============================================================================
# FUNCTIES DEFINIËREN PER METHODE
# ==============================================================================

# METHODE 1: Huidige terra::focal()
methode_focal <- function(r_in) {
  f_matrix   <- matrix(1, nrow = (2 * stral_cellen + 1), ncol = (2 * stral_cellen + 1))
  r_buffered <- terra::focal(r_in, w = f_matrix, fun = "max", na.rm = TRUE)
  r_buffered <- terra::ifel(r_buffered > 0, 1, NA)
  cl         <- terra::patches(r_buffered, directions = 4, zeroAsNA = TRUE)
  cl_out     <- terra::mask(cl, r_in)
  return(cl_out)
}

# METHODE 2: terra::buffer() (Raster expansie)
methode_terra_buffer <- function(r_in) {
  r_buffered <- terra::buffer(r_in, width = straal_m, background = 0)
  r_buffered <- terra::ifel(r_buffered > 0, 1, NA)
  cl         <- terra::patches(r_buffered, directions = 4, zeroAsNA = TRUE)
  cl_out     <- terra::mask(cl, r_in)
  return(cl_out)
}

# METHODE 3: Vectorbuffer via GEOS (as.polygons + buffer)
methode_vector_geos <- function(r_in) {
  # 1. Zet raster om naar polygonen
  v_poly <- terra::as.polygons(r_in, dissolve = TRUE, na.rm = TRUE)
  if (nrow(v_poly) == 0) return(r_in * NA)
  
  # 2. Buffer met GEOS en combineer overlappende polygonen
  v_buf <- terra::buffer(v_poly, width = straal_m)
  v_dis <- terra::aggregate(v_buf)
  v_dis$ID <- 1:nrow(v_dis)
  
  # 3. Terug rasterizen
  cl_raster <- terra::rasterize(v_dis, r_in, field = "ID")
  cl_out    <- terra::mask(cl_raster, r_in)
  return(cl_out)
}

# METHODE 4: qgisprocess (GDAL Proximity)
heeft_qgis <- requireNamespace("qgisprocess", quietly = TRUE) && qgisprocess::qgis_has_algorithm("gdal:proximity")

methode_qgis <- function(r_in) {
  if (!heeft_qgis) return(NULL)
  
  f_in  <- tempfile(fileext = ".tif")
  f_out <- tempfile(fileext = ".tif")
  writeRaster(r_in, f_in, overwrite = TRUE)
  
  qgisprocess::qgis_run_algorithm(
    "gdal:proximity",
    INPUT = f_in,
    BAND = 1,
    VALUES = "1",
    UNITS = 0,
    MAX_DISTANCE = straal_m, # was MAX_DIST
    NODATA = -9999,
    OUTPUT = f_out,
    .quiet = TRUE
  )
  
  r_prox <- rast(f_out)
  r_buf  <- ifel(!is.na(r_prox) & r_prox >= 0 & r_prox <= straal_m, 1, NA)
  cl     <- terra::patches(r_buf, directions = 4, zeroAsNA = TRUE)
  cl_out <- terra::mask(cl, r_in)
  
  unlink(c(f_in, f_out))
  return(cl_out)
}

# ==============================================================================
# 4. BENCHMARK UITVOEREN (16M cellen)
# ==============================================================================
message("\n3. Benchmark starten op 16 miljoen cellen...")

cat("\n--- Test Methode 1: terra::focal ---\n")
t1 <- system.time({ res1 <- methode_focal(r_bin) })
print(t1)

cat("\n--- Test Methode 3: Vector GEOS (as.polygons + buffer) ---\n")
t3 <- system.time({ res3 <- methode_vector_geos(r_bin) })
print(t3)

cat("\n--- Test Methode 4: qgisprocess (GDAL Proximity) ---\n")
t4 <- system.time({ res4 <- methode_qgis(r_bin) })
print(t4)

# ==============================================================================
# 5. OVERZICHT RESULTATEN
# ==============================================================================
cat("\n================ RESUMÉ TIJDEN (16M cellen) ================\n")
cat("1. focal()         : ", round(t1["elapsed"], 2), "s\n")
cat("3. Vector (GEOS)   : ", round(t3["elapsed"], 2), "s\n")
cat("4. qgisprocess     : ", round(t4["elapsed"], 2), "s\n")
cat("=============================================================\n")
