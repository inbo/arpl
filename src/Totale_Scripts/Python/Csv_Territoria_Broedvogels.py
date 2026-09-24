import geopandas as gpd
import pandas as pd
from pathlib import Path
from tqdm import tqdm

# Paden instellen
DIR_AVIMAP = Path("./AvimapKarteringen")
DIR_SOORTEN = Path("./Soortenkaarten")
OUTPUT_CSV = "alle_puntwaarnemingen_gecombineerd.csv"


def verwerk_map(folder_path, bron_naam):
    """Verwerkt alle shapefiles in een map met een tqdm voortgangsbalk."""
    shp_files = list(folder_path.glob("*.shp"))
    if not shp_files:
        print(f"Geen .shp bestanden gevonden in: {folder_path}")
        return []

    dfs = []
    # tqdm lus om de voortgang te tonen per bestand
    pbar = tqdm(shp_files, desc=f"Verwerken {bron_naam}", unit="bestand")
    for file_path in pbar:
        # Toon huidige bestandsnaam in de statusbalk
        pbar.set_postfix_str(file_path.name[:25])

        try:
            gdf = gpd.read_file(file_path)
            if gdf.empty:
                continue

            # Metadata toevoegen
            gdf['bron_map'] = bron_naam
            gdf['bron_bestand'] = file_path.name

            # Coördinatentransformaties (alleen uitvoeren als het nodig is)
            # 1. Lambert 72 (EPSG:31370)
            if gdf.crs != "EPSG:31370":
                gdf_lambert = gdf.to_crs(epsg=31370)
            else:
                gdf_lambert = gdf

            gdf['x_lambert'] = gdf_lambert.geometry.x.round(2)
            gdf['y_lambert'] = gdf_lambert.geometry.y.round(2)

            # 2. WGS84 (EPSG:4326)
            if gdf.crs != "EPSG:4326":
                gdf_wgs = gdf.to_crs(epsg=4326)
            else:
                gdf_wgs = gdf

            gdf['longitude'] = gdf_wgs.geometry.x
            gdf['latitude'] = gdf_wgs.geometry.y

            # Omzetten naar gewone Pandas DataFrame (verwijder Geopandas geometry)
            df = pd.DataFrame(gdf.drop(columns=['geometry']))
            dfs.append(df)

        except Exception as e:
            print(f"\nFout bij verwerken van {file_path.name}: {e}")

    return dfs


# Main uitvoer
if __name__ == "__main__":
    alle_dfs = []

    # 1. Verwerk Soortenkaarten
    print("--- Start verwerking Soortenkaarten ---")
    dfs_soorten = verwerk_map(DIR_SOORTEN, "Soortenkaarten")
    alle_dfs.extend(dfs_soorten)

    # 2. Verwerk AvimapKarteringen
    print("\n--- Start verwerking AvimapKarteringen ---")
    dfs_avimap = verwerk_map(DIR_AVIMAP, "AvimapKarteringen")
    alle_dfs.extend(dfs_avimap)

    if not alle_dfs:
        print("Geen data om samen te voegen.")
    else:
        # 3. Samenvoegen
        print("\nDatasets samenvoegen in geheugen...")
        df_totaal = pd.concat(alle_dfs, ignore_index=True)
        print(f"Totaal aantal rijen ingelezen: {len(df_totaal):,}")

        # 4. Optionele Ontdubbeling (Haal '#' weg op de volgende regel om in te schakelen)
        # print("Dubbelen verwijderen...")
        # df_totaal = df_totaal.drop_duplicates(subset=['euring', 'jaar', 'x_lambert', 'y_lambert', 'aantal'], keep='first')

        # 5. Exporteren naar CSV
        print(f"Schrijven naar CSV ({OUTPUT_CSV})... Dit kan even duren.")
        df_totaal.to_csv(OUTPUT_CSV, index=False, encoding='utf-8-sig')
        print("Klaar! Samenvoeging succesvol afgerond.")
