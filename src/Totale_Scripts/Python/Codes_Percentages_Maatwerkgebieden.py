from pathlib import Path
import geopandas as gpd
import pandas as pd
import polars as pl
import numpy as np

# =========================================================================
# CONFIGURATIE & INSTELLINGEN
# =========================================================================
MAP_SHAPEFILES = Path("Shapefiles")  # Map met je maatwerkgebieden (.shp)
PAD_BWK_SHP = "BwkHab.shp"  # Originele BWK vector-shapefile
BUFFER_METER = 0  # 2 km buffer (0 = exact maatwerkgebied)
MINIMUM_OPP_HA = 1.0  # Drempel: minstens 1.0 hectare!

# =========================================================================
# 1. SHAPEFILES INLADEN & BUFFEREN (GEFIXTE CRS VERWERKING)
# =========================================================================
print("-> Maatwerkgebieden & BWK Shapefile inladen...")
DOEL_CRS = "EPSG:31370"  # Belge Lambert 72

shp_bestanden = list(MAP_SHAPEFILES.glob("*.shp"))
lijst_gebieden = []

for shp_pad in shp_bestanden:
    gdf = gpd.read_file(shp_pad)

    # 1. Herprojecteer direct naar EPSG:31370
    if gdf.crs is None or gdf.crs.to_epsg() != 31370:
        gdf = gdf.to_crs(DOEL_CRS)
    else:
        gdf.crs = DOEL_CRS  # Forceer exact dezelfde CRS definitie

    gdf['geometry'] = gdf.geometry.make_valid()
    gdf['Maatwerkgebied'] = shp_pad.stem
    lijst_gebieden.append(gdf[['Maatwerkgebied', 'geometry']])

# Nu hebben alle GeoDataFrames exact dezelfde CRS-definitie en slaagt concat
gdf_gebieden = gpd.GeoDataFrame(pd.concat(lijst_gebieden, ignore_index=True), crs=DOEL_CRS)

# Inladen van BWK shapefile en herprojecteren
gdf_bwk = gpd.read_file(PAD_BWK_SHP)
if gdf_bwk.crs is None or gdf_bwk.crs.to_epsg() != 31370:
    gdf_bwk = gdf_bwk.to_crs(DOEL_CRS)
else:
    gdf_bwk.crs = DOEL_CRS

gdf_bwk['geometry'] = gdf_bwk.geometry.make_valid()

# Buffer van 2 km toepassen op de maatwerkgebieden
if BUFFER_METER > 0:
    print(f"-> Buffer van {BUFFER_METER}m toepassen op maatwerkgebieden...")
    gdf_gebieden['geometry'] = gdf_gebieden.geometry.buffer(BUFFER_METER)

# =========================================================================
# 2. RUIMTELIJKE OVERLAY (SNIJDEN VAN POLYGONEN)
# =========================================================================
print("-> Polygonen snijden (Intersection van Maatwerkgebieden x BWK)...")
# overlay snijdt de BWK polygonen exact af op de (gebufferde) maatwerkgebieden
gdf_intersect = gpd.overlay(gdf_gebieden, gdf_bwk, how='intersection')

# Bereken de werkelijke oppervlakte in m² van elk gesneden polygoonfragment
gdf_intersect['intersect_opp_m2'] = gdf_intersect.geometry.area

# Omzetten naar Pandas DataFrame voor snelle verwerking van EENH / HAB kolommen
df_intersect = pd.DataFrame(gdf_intersect.drop(columns='geometry'))

# =========================================================================
# 3. VERDEELLOGICA UIT JOUW R-SCRIPT TOEPASSEN (EENH1-4 & HAB1-4)
# =========================================================================
print("-> Bedekkingspercentages en ecologische fracties verwerken...")


# Gewichten-functies uit jouw R-script
def get_eenh_pct(n, idx):
    if n == 1 and idx == 1: return 1.00
    if n == 2: return 0.70 if idx == 1 else (0.30 if idx == 2 else 0.0)
    if n == 3: return 0.60 if idx == 1 else (0.20 if idx in [2, 3] else 0.0)
    if n >= 4: return 0.60 if idx == 1 else (0.20 if idx == 2 else (0.10 if idx in [3, 4] else 0.0))
    return 0.0


records = []

# Doorlopen per gesneden polygoonfragment
for _, row in df_intersect.iterrows():
    m_gebied = row['Maatwerkgebied']
    opp_poly_m2 = row['intersect_opp_m2']

    # --- A. BWK EENHEDEN (EENH1 t/m EENH4) ---
    eenh_cols = [c for c in ['EENH1', 'EENH2', 'EENH3', 'EENH4'] if c in row and pd.notna(row[c])]
    eenh_codes = [str(row[c]).strip().lower() for c in eenh_cols if str(row[c]).strip() not in ['', 'nan', 'none']]
    n_eenh = len(eenh_codes)

    for idx, code in enumerate(eenh_codes, 1):
        frac = get_eenh_pct(n_eenh, idx)
        opp_ha = (opp_poly_m2 * frac) / 10000.0
        if opp_ha > 0:
            records.append({'Maatwerkgebied': m_gebied, 'CODE': code, 'opp_ha': opp_ha})

    # --- B. NATURA 2000 HABITATS (HAB1 t/m HAB4 + PHAB1 t/m PHAB4) ---
    for i in range(1, 5):
        hab_col = f'HAB{i}'
        phab_col = f'PHAB{i}'

        if hab_col in row and pd.notna(row[hab_col]):
            code = str(row[hab_col]).strip().lower()
            if code not in ['', 'nan', 'none']:
                # Als PHAB aanwezig is, pak die waarde, anders 0
                phab_val = float(row[phab_col]) if phab_col in row and pd.notna(row[phab_col]) else 0.0
                frac = phab_val / 100.0
                opp_ha = (opp_poly_m2 * frac) / 10000.0
                if opp_ha > 0:
                    records.append({'Maatwerkgebied': m_gebied, 'CODE': code, 'opp_ha': opp_ha})

# =========================================================================
# 4. AGGREGATIE EN FILTEREN OP MINSTENS 1 HECTARE
# =========================================================================
print("-> Samenvatten en filteren op drempelwaarde (>= 1 hectare)...")
df_polars = pl.DataFrame(records)

df_finaal = (
    df_polars
    .group_by(["Maatwerkgebied", "CODE"])
    .agg(pl.col("opp_ha").sum())
    .with_columns(
        # Relative percentage berekenen t.o.v. het totale toegekende habitat-oppervlak in de gebufferde zone
        (pl.col("opp_ha") / pl.col("opp_ha").sum().over("Maatwerkgebied") * 100).alias("Percentage")
    )
    # FILTER: Minstens 1 hectare oppervlakte!
    .filter(pl.col("opp_ha") >= MINIMUM_OPP_HA)
    .sort(["Maatwerkgebied", "opp_ha"], descending=[False, True])
    .with_columns([
        pl.col("Percentage").round(2),
        pl.col("opp_ha").round(2)
    ])
)

# =========================================================================
# 5. EXPORT
# =========================================================================
print("\n--- RESULTATEN DE MATEN (DIRECT UIT VECTOR SHAPEFILE) ---")
print(df_finaal.filter(pl.col("Maatwerkgebied") == "De_Maten"))

output_naam = f"Habitat_Overzicht_Vector_Buffer_{BUFFER_METER}m_Min_{MINIMUM_OPP_HA}ha.csv"
df_finaal.write_csv(output_naam)
print(f"\n✓ Succesvol verwerkt en opgeslagen als: {output_naam}")
