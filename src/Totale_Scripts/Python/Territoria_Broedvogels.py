import os
import re
import pandas as pd
import geopandas as gpd

# ==========================================
# 1. CONFIGURATIE EN PADEN
# ==========================================

# Mappen en bestands-paden
SHP_SOORTEN_MAP = "./soortenkaarten"  # Map met de shapefiles per soort
OUTPUT_MAP = "./output_territoria"  # Map voor de gefilterde resultaten
MAATWERK_EXCEL = "./Soortenlijst_Maatwerkgebieden_Gefilterd.xlsx"
AFSTANDEN_EXCEL = "./Soorten_bwk_afstanden.xlsx"
SHAPEFILE_GEBIEDEN_MAP = "./Shapefiles"  # Map met de gebieden-shapefiles

# Coördinatensysteem voor analyse (Lambert 72 voor België in meters)
CRS_BELGIE = "EPSG:31370"

# Standaard kolomnaam voor broedcode (wordt opschoond)
BROEDCODE_COL_TARGET = "broedcode"

# Koppeling van gebiedsnamen in Excel naar de gewenste output-codes
GEBIED_MAPPING = {
    "De_Maten": "DM",
    "Heesbossen": "HB",
    "Mechelse_Heide": "MH",
    "Kalmthoutse_Heide": "KH",
    "Turnhouts_Vennegebied": "TV",
    "Voerstreek": "VS"
}

# De doellijst van 25 soorten
TARGET_SOORTEN = [
    "blauwborst", "boomleeuwerik", "boompieper", "bruine kiekendief", "fluiter",
    "grauwe klauwier", "grutto", "ijsvogel", "kwak", "kwartelkoning", "matkop",
    "middelste bonte specht", "nachtegaal", "nachtzwaluw", "paapje", "porseleinhoen",
    "roerdomp", "tapuit", "watersnip", "wespendief", "wielewaal", "woudaap",
    "wulp", "zomertortel", "zwarte specht", "zwartkopmeeuw"
]

# Standaard buffer indien een soort niet gevonden wordt in de afstanden-excel (in meter)
STANDAARD_BUFFER_M = 5000


# ==========================================
# 2. HULPFUNCTIES
# ==========================================

def opschonen_key(naam_str: str) -> str:
    """
    Zet 'Middelste Bonte Specht', 'Middelste_Bonte_Specht' of 'middelstebontespecht'
    om naar 'middelstebontespecht' voor een vlekkeloze matching.
    """
    return re.sub(r'[^a-zA-Z0-9]', '', str(naam_str)).lower()


def formatteer_output_naam(naam_str: str) -> str:
    """
    Zet 'middelste bonte specht' of 'Middelste_Bonte_Specht' om naar 'Middelste_Bonte_Specht'.
    """
    schoon = str(naam_str).replace("_", " ")
    woorden = [w.capitalize() for w in schoon.strip().split()]
    return "_".join(woorden)


def vind_shapefile_voor_soort(soortnaam: str, shp_map: str):
    """
    Zoekt de juiste .shp van de soort in de map 'soortenkaarten'.
    """
    target_key = opschonen_key(soortnaam)
    if not os.path.exists(shp_map):
        return None

    for f in os.listdir(shp_map):
        if f.lower().endswith(".shp") and target_key in opschonen_key(f):
            return os.path.join(shp_map, f)
    return None


def vind_shapefile_voor_gebied(gebiedsnaam: str, shp_map: str):
    """
    Zoekt de juiste .shp van een gebied in de map 'Shapefiles'.
    """
    target_key = opschonen_key(gebiedsnaam)
    if not os.path.exists(shp_map):
        return None

    for f in os.listdir(shp_map):
        if f.lower().endswith(".shp") and target_key in opschonen_key(f):
            return os.path.join(shp_map, f)
    return None


# ==========================================
# 3. HOOFDSCRIPT
# ==========================================

def main():
    os.makedirs(OUTPUT_MAP, exist_ok=True)

    # ------------------------------------------
    # A. EXCEL BESTANDEN INLEZEN EN MAPPEN
    # ------------------------------------------
    print("Excel bestanden inlezen...")

    # 1. Buffers inlezen uit Soorten_bwk_afstanden.xlsx
    df_afstanden = pd.read_excel(AFSTANDEN_EXCEL)
    buffer_dict = {}
    for _, row in df_afstanden.iterrows():
        key = opschonen_key(row["Soort"])
        try:
            buffer_dict[key] = float(row["Dispersiecap_m"])
        except (ValueError, TypeError):
            pass

    # 2. Maatwerkgebieden inlezen uit Soortenlijst_Maatwerkgebieden_Gefilterd.xlsx
    df_maatwerk = pd.read_excel(MAATWERK_EXCEL)
    df_maatwerk["key"] = df_maatwerk["Nederlandse naam"].apply(opschonen_key)

    # Voor het eindoorsicht
    gevonden_per_gebied = {code: set() for code in GEBIED_MAPPING.values()}
    relevante_soorten_per_gebied = {code: set() for code in GEBIED_MAPPING.values()}

    # ------------------------------------------
    # B. GEBIEDEN SHAPEFILES INLEZEN
    # ------------------------------------------
    print("Gebieden shapefiles inlezen uit 'Shapefiles'...")
    gebieden_gdf = {}
    for gebied_naam in GEBIED_MAPPING.keys():
        shp_pad = vind_shapefile_voor_gebied(gebied_naam, SHAPEFILE_GEBIEDEN_MAP)
        if shp_pad:
            gdf = gpd.read_file(shp_pad)
            # Zorg dat de gebieden altijd in Lambert 72 staan
            if gdf.crs is None:
                gdf.set_crs(CRS_BELGIE, inplace=True)
            else:
                gdf = gdf.to_crs(CRS_BELGIE)
            gebieden_gdf[gebied_naam] = gdf
            print(f"  └─ Gebied geladen: {gebied_naam} ({os.path.basename(shp_pad)})")
        else:
            print(f"  └─ ⚠️ WAARSCHUWING: Geen shapefile gevonden voor gebied '{gebied_naam}'")

        # ------------------------------------------
        # C. SOORTEN-SHAPEFILES VERWERKEN
        # ------------------------------------------
        print("\nSoortenkaarten verwerken uit map 'soortenkaarten'...")
        for soort_ruw in TARGET_SOORTEN:
            soort_key = opschonen_key(soort_ruw)
            formatted_soort = formatteer_output_naam(soort_ruw)

            # 1. Bepaal EERST voor welke gebieden deze soort relevant is volgens Maatwerk-Excel
            maatwerk_row = df_maatwerk[df_maatwerk["key"] == soort_key]
            if maatwerk_row.empty:
                print(f"\n⚠️ Soort '{formatted_soort}' niet gevonden in {MAATWERK_EXCEL}")
                continue

            # Registreer alle gebieden waar de soort een '1' heeft
            is_relevant_voor_enig_gebied = False
            for gebied_naam, gebied_code in GEBIED_MAPPING.items():
                if gebied_naam in maatwerk_row.columns:
                    if maatwerk_row[gebied_naam].values[0] == 1:
                        relevante_soorten_per_gebied[gebied_code].add(formatted_soort)
                        is_relevant_voor_enig_gebied = True

            if not is_relevant_voor_enig_gebied:
                # Soort heeft nergens een 1 staan, kan overgeslagen worden
                continue

            # 2. Zoek de shapefile van de soort
            shp_soort_pad = vind_shapefile_voor_soort(soort_ruw, SHP_SOORTEN_MAP)
            if not shp_soort_pad:
                print(f"\n❌ Geen shapefile gevonden voor soort: {formatted_soort}")
                continue

            print(f"\nVerwerken: {formatted_soort} ({os.path.basename(shp_soort_pad)})")

            # 3. Bepaal buffer voor deze soort
            buffer_m = buffer_dict.get(soort_key, STANDAARD_BUFFER_M)

            # 4. Lees de soorten-shapefile in met GeoPandas
            gdf_soort = gpd.read_file(shp_soort_pad)

            if gdf_soort.empty:
                print("  └─ ℹ️ Shapefile is leeg.")
                continue

            # Zorg voor Lambert 72 CRS
            if gdf_soort.crs is None:
                gdf_soort.set_crs(CRS_BELGIE, inplace=True)
            else:
                gdf_soort = gdf_soort.to_crs(CRS_BELGIE)

            # Opschonen van kolomnamen
            clean_cols = {col: str(col).strip().lower() for col in gdf_soort.columns}
            gdf_soort.rename(columns=clean_cols, inplace=True)

            # Zoek broedcode kolom
            broed_col_found = next((c for c in gdf_soort.columns if "broedcode" in c), None)

            if not broed_col_found:
                print(f"  └─ ⚠️ Kolom 'broedcode' niet gevonden in shapefile van {formatted_soort}.")
                continue

            # Filter op broedcode >= 4
            gdf_soort[broed_col_found] = pd.to_numeric(gdf_soort[broed_col_found], errors='coerce')
            gdf_broed = gdf_soort[gdf_soort[broed_col_found] >= 4].copy()

            if gdf_broed.empty:
                print("  └─ ℹ️ Geen waarnemingen met broedcode >= 4.")
                continue

            # 5. Voer de ruimtelijke snijdingen per relevant gebied uit
            for gebied_naam, gebied_code in GEBIED_MAPPING.items():
                if gebied_naam not in maatwerk_row.columns:
                    continue

                is_relevant = maatwerk_row[gebied_naam].values[0] == 1
                if not is_relevant or gebied_naam not in gebieden_gdf:
                    continue

                # Buffer toepassen op het gebied
                gebied_poly = gebieden_gdf[gebied_naam].geometry.union_all()
                gebufferde_poly = gebied_poly.buffer(buffer_m)

                # Spatial filtering
                punten_in_gebied = gdf_broed[gdf_broed.geometry.within(gebufferde_poly)].copy()

                if not punten_in_gebied.empty:
                    punten_in_gebied["x"] = punten_in_gebied.geometry.x
                    punten_in_gebied["y"] = punten_in_gebied.geometry.y

                    df_export = pd.DataFrame(punten_in_gebied.drop(columns=["geometry"]))

                    uitvoer_naam = f"{gebied_code}_Territoria_{formatted_soort}.csv"
                    uitvoer_pad = os.path.join(OUTPUT_MAP, uitvoer_naam)

                    df_export.to_csv(uitvoer_pad, index=False)
                    print(
                        f"  └─ Opbrengst {gebied_code} (buffer {int(buffer_m)}m): {len(df_export)} rijen opgeslagen -> {uitvoer_naam}")

                    gevonden_per_gebied[gebied_code].add(formatted_soort)

    # ------------------------------------------
    # D. EINDOVERZICHT PRINTEN
    # ------------------------------------------
    print("\n" + "=" * 60)
    print("EINDSTATISTIEKEN PER GEBIED (RELEVANTE SOORTEN ZONDER WAARNEMINGEN)")
    print("=" * 60)

    for gebied_naam, gebied_code in GEBIED_MAPPING.items():
        relevante = relevante_soorten_per_gebied[gebied_code]
        gevonden = gevonden_per_gebied[gebied_code]
        ontbrekend = sorted(list(relevante - gevonden))

        print(f"\n📍 Gebied: {gebied_naam} ({gebied_code})")
        print(f"   ├─ Totaal relevante soorten: {len(relevante)}")
        print(f"   ├─ Succesvol geëxporteerd:   {len(gevonden)}")
        print(f"   └─ NUL waarnemingen voor ({len(ontbrekend)}):")

        if ontbrekend:
            for s in ontbrekend:
                print(f"       • {s}")
        else:
            print("       (Geen! Alle relevante soorten hadden waarnemingen)")


if __name__ == "__main__":
    main()
