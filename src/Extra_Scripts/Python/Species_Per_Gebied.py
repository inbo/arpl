import os
import re
import geopandas as gpd
import pandas as pd
from shapely.geometry import Point

# --- PATH SETUP ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SHAPEFILE_PATH = os.path.join(SCRIPT_DIR, "Turnhouts_Vennegebied.shp")
EXCEL_PATH = os.path.join(SCRIPT_DIR, "Soortenlijst_Maatwerkgebieden.xlsx")

INPUT_FILES = [
    os.path.join(
        SCRIPT_DIR, "INBODATAVR-481_dumptemp20260531_met_toestemming.csv"
    ),
    os.path.join(SCRIPT_DIR, "INBODATAVR-481_dumptemp20260531_sws.csv"),
    os.path.join(
        SCRIPT_DIR,
        "INBODATAVR_494_basis_wbe_dumptem20260705_met_toestemming.csv",
    ),
    os.path.join(SCRIPT_DIR, "INBODATAVR_494_basis_wbe_dumptem20260705_sws.csv"),
]

# --- CONFIGURATIE ---
MIN_YEAR = 2010
BUFFER_AFSTAND_METERS = 0

SPECIES_COLUMN = "naam_nl"
LATIN_COLUMN = "naam_lat"
VALIDATION_COLUMN = "status"
DATE_COLUMN = "datum"
X_COORD_COLUMN = "x"
Y_COORD_COLUMN = "y"

CSV_CRS = "EPSG:31370"
TARGET_CRS = "EPSG:31370"

# --- STEP 0: SOORTENLIJST INLADEN EN FILTEREN ---
if not os.path.exists(EXCEL_PATH):
  raise FileNotFoundError(f"⚠️ Kan het Excel-bestand niet vinden op {EXCEL_PATH}")

print("📋 Soortenlijst uit Excel inladen en filteren...")
df_excel = pd.read_excel(EXCEL_PATH)
df_excel.columns = df_excel.columns.str.strip()

mask_gebied = df_excel["TurnhoutsVennegebied"] == 1
mask_gedaan = df_excel["Gedaan"].astype(str).str.strip().str.lower() == "ja"

df_doelsoorten = df_excel[mask_gebied & mask_gedaan].copy()

doelsoorten_nl = set(
    df_doelsoorten["Nederlandse naam"].dropna().astype(str).str.strip()
)

print(f"🎯 Aantal geselecteerde doelsoorten uit Excel: {len(doelsoorten_nl)}")

soorten_teller = {
    soort: {"aantal": 0, "lat_naam": ""} for soort in doelsoorten_nl
}

# Lookup-map op kleine letters
doelsoorten_lookup = {soort.lower(): soort for soort in doelsoorten_nl}

# --- STEP 1: SHAPEFILE PREPAREREN ---
if not os.path.exists(SHAPEFILE_PATH):
  raise FileNotFoundError(
      f"⚠️ Kan de shapefile niet vinden op {SHAPEFILE_PATH}"
  )

print("🗺️ Shapefile inladen...")
grenzen = gpd.read_file(SHAPEFILE_PATH).to_crs(TARGET_CRS)

if BUFFER_AFSTAND_METERS > 0:
  filter_gebied = grenzen.buffer(BUFFER_AFSTAND_METERS)
else:
  filter_gebied = grenzen.geometry

filter_gebied_union = filter_gebied.union_all()

# RegEx patroon om ondersoort-toevoegingen te verwijderen
# Verwijdert alles vanaf ssp, subsp, s.s., s.l., var. (met of zonder punt, hoofd- of kleine letters)
REGEX_ONDERSOORT = (
    r"\s+\b(ssp|subsp|s\.s|s\.l|var)\b.*$"  # \b zorgt voor exacte woordgrenzen
)


# --- STEP 2: VERWERK CSV'S IN CHUNKS ---
for file_path in INPUT_FILES:
  filename = os.path.basename(file_path)
  if not os.path.exists(file_path):
    print(f"⚠️ Bestand niet gevonden: {filename}. Overslaan.")
    continue

  print(f"🔄 Verwerken van {filename}...")
  chunk_count = 0

  for chunk in pd.read_csv(
      file_path, chunksize=100000, low_memory=False, decimal=","
  ):
    chunk_count += 1

    clean_validation = (
        chunk[VALIDATION_COLUMN].astype(str).str.strip().str.lower()
    )
    chunk = chunk[clean_validation.str.startswith("goedgekeurd", na=False)]
    if chunk.empty:
      continue

    chunk[DATE_COLUMN] = pd.to_datetime(chunk[DATE_COLUMN], errors="coerce")
    chunk = chunk[chunk[DATE_COLUMN].dt.year >= MIN_YEAR]
    if chunk.empty:
      continue

    chunk_coords = chunk.dropna(subset=[X_COORD_COLUMN, Y_COORD_COLUMN])
    if chunk_coords.empty:
      continue

    geometry = [
        Point(xy)
        for xy in zip(chunk_coords[X_COORD_COLUMN], chunk_coords[Y_COORD_COLUMN])
    ]
    chunk_gpd = gpd.GeoDataFrame(chunk_coords, geometry=geometry, crs=CSV_CRS)
    chunk_gpd = chunk_gpd.to_crs(TARGET_CRS)

    binnen_gebied_mask = chunk_gpd.geometry.within(filter_gebied_union)
    gefilterde_chunk = chunk_gpd[binnen_gebied_mask].copy()

    if gefilterde_chunk.empty:
      continue

    # 1. Opschonen Latijnse naam (Genus + Soortnaam)
    gefilterde_chunk["wetenschappelijk_clean"] = (
        gefilterde_chunk[LATIN_COLUMN]
        .astype(str)
        .str.strip()
        .str.split()
        .str[:2]
        .str.join(" ")
    )

    # 2. Flexible Opschoning Nederlandse naam
    gefilterde_chunk["nl_clean"] = (
        gefilterde_chunk[SPECIES_COLUMN]
        .astype(str)
        .str.replace(
            REGEX_ONDERSOORT, "", regex=True, flags=re.IGNORECASE
        )  # Stript ssp crecca etc.
        .str.strip()
    )

    # Groepeer binnen de chunk
    chunk_tellingen = (
        gefilterde_chunk.groupby(["nl_clean", "wetenschappelijk_clean"])
        .size()
        .reset_index(name="aantal")
    )

    for _, row in chunk_tellingen.iterrows():
      nl_naam_raw = row["nl_clean"]
      nl_naam_key = nl_naam_raw.lower()
      lat_naam = row["wetenschappelijk_clean"]
      aantal = row["aantal"]

      # Check op opgeschoonde naam
      if nl_naam_key in doelsoorten_lookup:
        officiele_nl_naam = doelsoorten_lookup[nl_naam_key]
        soorten_teller[officiele_nl_naam]["aantal"] += aantal

        if not soorten_teller[officiele_nl_naam]["lat_naam"]:
          soorten_teller[officiele_nl_naam]["lat_naam"] = lat_naam

  print(f"✅ Klaar met {filename} ({chunk_count} chunks verwerkt).")

# --- STEP 3: RESULTATEN WEGSCHRIJVEN ---
data_voor_df = [
    {
        "Nederlandse_Naam": nl_naam,
        "Wetenschappelijke_Naam": info["lat_naam"],
        "Aantal_Waarnemingen": info["aantal"],
    }
    for nl_naam, info in soorten_teller.items()
]

df_resultaat = pd.DataFrame(data_voor_df)
df_resultaat = df_resultaat.sort_values(
    by=["Aantal_Waarnemingen", "Nederlandse_Naam"], ascending=[False, True]
)

output_file = os.path.join(SCRIPT_DIR, "soorten_in_gebied.csv")
df_resultaat.to_csv(output_file, index=False, sep=";")

print("\n🎉 Klaar!")
print(
    f"Overzicht voor alle {len(df_resultaat)} doelsoorten is succesvol"
    " opgesteld."
)
print(f"Het resultaat is opgeslagen in: {output_file}")