import os
import re
import glob
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point

# ==============================================================================
# CONFIGURATIE & PADEN
# ==============================================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Referentiebestanden
EXCEL_AFSTANDEN_PATH = os.path.join(SCRIPT_DIR, "soorten_bwk_afstanden.xlsx")
EXCEL_HTS_PATH = os.path.join(SCRIPT_DIR, "HTS.xlsx")
SHAPEFILES_DIR = os.path.join(SCRIPT_DIR, "Shapefiles")

# Mapconfiguratie per gebied
GEBIEDEN_CONFIG = [
    {
        "code": "TV",
        "naam": "Turnhouts Vennegebied",
        "excel_col": "Turnhouts_Vennegebied",
        "sheet_name": "Turnhouts_Vennegebied",
        "data_dir": os.path.join(SCRIPT_DIR, "Data_TV"),
        "shape_keywords": ["Turnhouts_Vennegebied", "Turnhout"],
    },
    {
        "code": "DM",
        "naam": "De Maten",
        "excel_col": "De_Maten",
        "sheet_name": "De_Maten",
        "data_dir": os.path.join(SCRIPT_DIR, "Data_DM"),
        "shape_keywords": ["De_Maten", "Maten"],
    },
    {
        "code": "HB",
        "naam": "Heesbossen",
        "excel_col": "Heesbossen",
        "sheet_name": "Heesbossen",
        "data_dir": os.path.join(SCRIPT_DIR, "Data_HB"),
        "shape_keywords": ["Heesbossen", "Heesbos"],
    },
    {
        "code": "KH",
        "naam": "Kalmthoutse Heide",
        "excel_col": "Kalmthoutse_Heide",
        "sheet_name": "Kalmthoutse_Heide",
        "data_dir": os.path.join(SCRIPT_DIR, "Data_KH"),
        "shape_keywords": ["Kalmthoutse_Heide", "Kalmthout"],
    },
    {
        "code": "MH",
        "naam": "Mechelse Heide",
        "excel_col": "Mechelse_Heide",
        "sheet_name": "Mechelse_Heide",
        "data_dir": os.path.join(SCRIPT_DIR, "Data_MH"),
        "shape_keywords": ["Mechelse_Heide", "Mechelse"],
    },
    {
        "code": "VS",
        "naam": "Voerstreek",
        "excel_col": "Voerstreek",
        "sheet_name": "Voerstreek",
        "data_dir": os.path.join(SCRIPT_DIR, "Data_VS"),
        "shape_keywords": ["Voerstreek", "Voer"],
    },
]

# Kolomnamen uit soorten_bwk_afstanden.xlsx
EXCEL_SOORT_COL = "Soort"
EXCEL_BUFFER_COL = "Dispersiecap_m"
EXCEL_GROEP_COL = "Groep"

# Trefwoorden voor goed bestudeerde groepen
GOED_BESTUDEERD_KEYWORDS = ["plant", "vogel", "vlinder"]

# Kolomnamen invoerdata
SPECIES_COLUMN = "naam_nl"
LATIN_COLUMN = "naam_lat"
VALIDATION_COLUMN = "status"
DATE_COLUMN = "datum"
X_COORD_COLUMN = "x"
Y_COORD_COLUMN = "y"

# CRS & Instellingen
CSV_CRS = "EPSG:31370"
TARGET_CRS = "EPSG:31370"
MIN_YEAR_GLOBAAL = 2010

REGEX_ONDERSOORT = r"\s+\b(ssp|subsp|s\.s|s\.l|var)\b.*$"


# ==============================================================================
# HULPFUNCTIES
# ==============================================================================
def normalize_name(name):
    """Verwijdert alle spaties, koppeltekens, leestekens en zet om naar kleine letters."""
    if pd.isna(name):
        return ""
    return re.sub(r"[^a-zA-Z0-9]", "", str(name)).lower()


def check_is_goed_bestudeerd(groep_tekst):
    """Bepaalt of de soortgroep behoort tot Planten, Vogels of Vlinders."""
    groep_norm = normalize_name(groep_tekst)
    return any(keyword in groep_norm for keyword in GOED_BESTUDEERD_KEYWORDS)


def find_shapefile_for_area(area_config):
    """Zoekt naar de passende shapefile in de map Shapefiles/."""
    if not os.path.exists(SHAPEFILES_DIR):
        return None

    shp_files = glob.glob(os.path.join(SHAPEFILES_DIR, "*.shp"))
    for kw in area_config["shape_keywords"]:
        for shp in shp_files:
            if kw.lower() in os.path.basename(shp).lower():
                return shp

    if len(shp_files) == 1:
        return shp_files[0]

    return None


def read_file_in_chunks(file_path, chunksize=100000):
    """Generator die data in chunks oplevert (werkt voor CSV én Excel)."""
    ext = os.path.splitext(file_path)[1].lower()

    if ext == ".csv":
        for chunk in pd.read_csv(file_path, chunksize=chunksize, low_memory=False, decimal=','):
            yield chunk
    elif ext in [".xlsx", ".xls"]:
        full_df = pd.read_excel(file_path)
        for i in range(0, len(full_df), chunksize):
            yield full_df.iloc[i:i + chunksize].copy()


def load_hts_for_area(excel_path, target_sheet):
    """Laadt de HTS-soorten uit het specifieke tabblad voor dit gebied."""
    hts_map = {}
    if not os.path.exists(excel_path):
        print(f"⚠️ Waarschuwing: {excel_path} niet gevonden.")
        return hts_map

    try:
        xls = pd.ExcelFile(excel_path)
        sheet_names = xls.sheet_names

        matched_sheet = None
        for s in sheet_names:
            if normalize_name(s) == normalize_name(target_sheet):
                matched_sheet = s
                break

        if not matched_sheet:
            print(f"⚠️ Waarschuwing: Tabblad '{target_sheet}' niet gevonden in HTS.xlsx. Beschikbaar: {sheet_names}")
            return hts_map

        df_hts = pd.read_excel(excel_path, sheet_name=matched_sheet)
        df_hts.columns = df_hts.columns.str.strip()

        kolom_kandidaten = ["naam_nl", "nederlandse naam", "soort", "naam", "soortnaam"]
        gevonden_col = None

        for col in df_hts.columns:
            if col.lower().strip() in kolom_kandidaten:
                gevonden_col = col
                break
        if not gevonden_col and not df_hts.empty:
            gevonden_col = df_hts.columns[0]

        if gevonden_col:
            for val in df_hts[gevonden_col].dropna():
                raw_val = str(val).strip()
                if raw_val:
                    hts_map[normalize_name(raw_val)] = raw_val

        print(f"🏷️ HTS-soorten geladen uit tabblad '{matched_sheet}': {len(hts_map)} soorten.")
    except Exception as e:
        print(f"⚠️ Fout bij inladen HTS-tabblad '{target_sheet}': {e}")

    return hts_map


# ==============================================================================
# HOOFDLOOP PER GEBIED
# ==============================================================================
for gebied in GEBIEDEN_CONFIG:
    code = gebied["code"]
    naam = gebied["naam"]
    excel_col = gebied["excel_col"]
    sheet_name = gebied["sheet_name"]
    data_dir = gebied["data_dir"]

    print("=" * 80)
    print(f"🌍 ANALYSE STARTEN VOOR GEBIED: {naam} ({code})")
    print("=" * 80)

    # 1. Controleer data directorij & bestanden
    if not os.path.exists(data_dir):
        print(f"⚠️ Map '{data_dir}' bestaat niet. Overslaan.\n")
        continue

    input_files = (
            glob.glob(os.path.join(data_dir, "*.csv")) +
            glob.glob(os.path.join(data_dir, "*.xlsx")) +
            glob.glob(os.path.join(data_dir, "*.xls"))
    )

    if not input_files:
        print(f"⚠️ Geen CSV/Excel-bestanden gevonden in {data_dir}. Overslaan.\n")
        continue

    # 2. Shapefile zoeken
    shapefile_path = find_shapefile_for_area(gebied)
    if not shapefile_path or not os.path.exists(shapefile_path):
        print(f"⚠️ Geen shapefile gevonden voor '{naam}' in '{SHAPEFILES_DIR}'. Overslaan.\n")
        continue

    print(f"🗺️ Shapefile geladen: {os.path.basename(shapefile_path)}")
    grenzen_base = gpd.read_file(shapefile_path).to_crs(TARGET_CRS)
    gebied_union_base = grenzen_base.geometry.union_all()

    # 3. Gebiedsspecifieke HTS-soorten laden
    hts_exclusief_map = load_hts_for_area(EXCEL_HTS_PATH, sheet_name)

    # 4. Hoofdsoortenlijst inladen & filteren op `excel_col == 1`
    if not os.path.exists(EXCEL_AFSTANDEN_PATH):
        print(f"⚠️ Fout: Afstanden-bestand '{EXCEL_AFSTANDEN_PATH}' niet gevonden!")
        continue

    df_excel = pd.read_excel(EXCEL_AFSTANDEN_PATH)
    df_excel.columns = df_excel.columns.str.strip()

    if excel_col not in df_excel.columns:
        print(f"⚠️ Kolom '{excel_col}' niet gevonden in {EXCEL_AFSTANDEN_PATH}. Overslaan.\n")
        continue

    df_excel[excel_col] = pd.to_numeric(df_excel[excel_col], errors='coerce').fillna(0)
    df_excel_gefilterd = df_excel[df_excel[excel_col] == 1].copy()

    print(f"🔍 Totaal doelsoorten gefilterd op {excel_col} (==1): {len(df_excel_gefilterd)}")

    # 5. Metadata per soort opbouwen
    soorten_meta = {}
    doelsoorten_lookup = {}

    for _, row in df_excel_gefilterd.iterrows():
        raw_naam = str(row[EXCEL_SOORT_COL]).strip()
        if pd.isna(raw_naam) or raw_naam == "nan":
            continue

        norm_key = normalize_name(raw_naam)
        is_hts_exclusief = norm_key in hts_exclusief_map
        schone_naam = hts_exclusief_map[norm_key] if is_hts_exclusief else raw_naam

        buffer_m = row.get(EXCEL_BUFFER_COL, 0)
        buffer_m = 0 if pd.isna(buffer_m) else float(buffer_m)

        groep_raw = str(row.get(EXCEL_GROEP_COL, "")).strip()
        is_goed_bestudeerd = check_is_goed_bestudeerd(groep_raw)

        if buffer_m > 0:
            gebufferd_gebied = grenzen_base.buffer(buffer_m).union_all()
        else:
            gebufferd_gebied = gebied_union_base

        soorten_meta[raw_naam] = {
            "norm_key": norm_key,
            "schone_nl_naam": schone_naam,
            "lat_naam": "",
            "buffer_m": buffer_m,
            "groep": groep_raw,
            "is_hts_exclusief": is_hts_exclusief,
            "is_goed_bestudeerd": is_goed_bestudeerd,
            "gebufferd_gebied": gebufferd_gebied,
            "waarnemingen_2010_plus": 0,
            "waarnemingen_2015_plus": 0
        }
        doelsoorten_lookup[norm_key] = raw_naam

    # 6. Verwerk data van alle invoerbestanden van dit gebied
    for file_path in input_files:
        filename = os.path.basename(file_path)
        print(f"  🔄 Verwerken: {filename}...")

        try:
            for chunk in read_file_in_chunks(file_path):
                chunk.columns = chunk.columns.str.strip()

                if VALIDATION_COLUMN not in chunk.columns or DATE_COLUMN not in chunk.columns:
                    continue

                # Status filter
                clean_validation = chunk[VALIDATION_COLUMN].astype(str).str.strip().str.lower()
                chunk = chunk[clean_validation.str.startswith("goedgekeurd", na=False)]
                if chunk.empty:
                    continue

                # Datum filter
                chunk[DATE_COLUMN] = pd.to_datetime(chunk[DATE_COLUMN], errors="coerce")
                chunk["jaar"] = chunk[DATE_COLUMN].dt.year
                chunk = chunk[chunk["jaar"] >= MIN_YEAR_GLOBAAL]
                if chunk.empty:
                    continue

                # Coördinaten filter
                chunk_coords = chunk.dropna(subset=[X_COORD_COLUMN, Y_COORD_COLUMN])
                if chunk_coords.empty:
                    continue

                # GeoDataFrame
                geometry = [Point(xy) for xy in zip(chunk_coords[X_COORD_COLUMN], chunk_coords[Y_COORD_COLUMN])]
                chunk_gpd = gpd.GeoDataFrame(chunk_coords, geometry=geometry, crs=CSV_CRS).to_crs(TARGET_CRS)

                # Namen opschonen & matchen
                chunk_gpd["wetenschappelijk_clean"] = chunk_gpd[LATIN_COLUMN].astype(str).str.strip().str.split().str[
                    :2].str.join(" ")
                chunk_gpd["nl_clean"] = chunk_gpd[SPECIES_COLUMN].astype(str).str.replace(
                    REGEX_ONDERSOORT, "", regex=True, flags=re.IGNORECASE
                ).str.strip()
                chunk_gpd["norm_key"] = chunk_gpd["nl_clean"].apply(normalize_name)

                matched_chunk = chunk_gpd[chunk_gpd["norm_key"].isin(doelsoorten_lookup.keys())]
                if matched_chunk.empty:
                    continue

                for norm_key, group in matched_chunk.groupby("norm_key"):
                    excel_key = doelsoorten_lookup[norm_key]
                    meta = soorten_meta[excel_key]

                    mask_binnen = group.geometry.within(meta["gebufferd_gebied"])
                    binnen_punten = group[mask_binnen]

                    if not binnen_punten.empty:
                        if not binnen_punten["nl_clean"].empty and " " in binnen_punten["nl_clean"].iloc[0]:
                            meta["schone_nl_naam"] = binnen_punten["nl_clean"].iloc[0]

                        if not meta["lat_naam"]:
                            meta["lat_naam"] = binnen_punten["wetenschappelijk_clean"].iloc[0]

                        n_2010 = len(binnen_punten)
                        n_2015 = len(binnen_punten[binnen_punten["jaar"] >= 2015])

                        meta["waarnemingen_2010_plus"] += n_2010
                        meta["waarnemingen_2015_plus"] += n_2015

        except Exception as e:
            print(f"⚠️ Fout bij verwerken bestand {filename}: {e}")

    # 7. Evaluatie & Resultaten exporteren per gebied
    resultaten = []
    for raw_naam, info in soorten_meta.items():
        is_hts = info["is_hts_exclusief"]
        is_goed = info["is_goed_bestudeerd"]
        n_2010 = info["waarnemingen_2010_plus"]
        n_2015 = info["waarnemingen_2015_plus"]

        behouden = False
        cat_label = ""
        reden = ""

        # REGEL 1: SBZ/SBP doelsoort -> ALTIJD BEHOUDEN
        if not is_hts:
            behouden = True
            cat_label = "SBZ / SBP Doelsoort"
            reden = f"Drempelvrij (SBZ/SBP): Behouden ongeacht aantal waarnemingen (Aantal 2015+: {n_2015})"

        # REGEL 2: HTS-exclusieve soort -> DREMPEL TOEPASSEN
        else:
            cat_label = "HTS Exclusief"
            if is_goed:
                # Vogels, Planten, Vlinders -> >= 10 sinds 2015
                if n_2015 >= 10:
                    behouden = True
                    reden = f"Goed onderzocht ({info['groep']}): >=10 waarnemingen sinds 2015 (Aantal: {n_2015})"
                else:
                    reden = f"Afgekeurd (Goed onderzocht): <10 waarnemingen sinds 2015 (Aantal: {n_2015})"
            else:
                # Overige groepen -> >= 10 sinds 2010
                if n_2010 >= 10:
                    behouden = True
                    reden = f"Matig onderzocht ({info['groep']}): >=10 waarnemingen sinds 2010 (Aantal: {n_2010})"
                else:
                    reden = f"Afgekeurd (Matig onderzocht): <10 waarnemingen sinds 2010 (Aantal: {n_2010})"

        resultaten.append({
            "Nederlandse_Naam": info["schone_nl_naam"],
            "Excel_Naam_Origineel": raw_naam,
            "Wetenschappelijke_Naam": info["lat_naam"],
            "Categorie": cat_label,
            "Soortgroep": info["groep"],
            "Dispersiebuffer_m": info["buffer_m"],
            "Onderzoeksgraad": "Goed (>=2015)" if is_goed else "Matig (>=2010)",
            "Waarnemingen_2015_2026": n_2015,
            "Waarnemingen_2010_2026": n_2010,
            "Behouden_In_Lijst": "JA" if behouden else "NEE",
            "Reden": reden
        })

    df_resultaat = pd.DataFrame(resultaten)
    df_resultaat = df_resultaat.sort_values(
        by=["Behouden_In_Lijst", "Categorie", "Nederlandse_Naam"],
        ascending=[False, True, True]
    )

    output_file = os.path.join(SCRIPT_DIR, f"soorten_geëvalueerd_met_buffers_{code}.csv")
    df_resultaat.to_csv(output_file, index=False, sep=";")

    print(f"🎉 Evaluatie voltooid voor {naam}! Resultaat opgeslagen in: {output_file}\n")

print("✨ ALLE ANALYSES ZIJN SUCCESVOL VOLTOOID!")
