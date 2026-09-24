import os
import re
import glob
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point

# ==============================================================================
# CONFIGURATIE & PADEN
# ==============================================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__)) if '__file__' in locals() else os.getcwd()

EXCEL_SPECIES_LIST = os.path.join(SCRIPT_DIR, "Soortenlijst_Maatwerkgebieden.xlsx")
EXCEL_BREEDING_SEASONS = os.path.join(SCRIPT_DIR, "Broedvogels_Datums.xlsx")
EXCEL_AFSTANDEN_PATH = os.path.join(SCRIPT_DIR, "soorten_bwk_afstanden.xlsx")
EXCEL_HTS_PATH = os.path.join(SCRIPT_DIR, "HTS.xlsx")
SHAPEFILES_DIR = os.path.join(SCRIPT_DIR, "Shapefiles")

GEBIEDEN_CONFIG = [
    {"code": "TV", "naam": "Turnhouts Vennegebied", "sheet_name": "Turnhouts_Vennegebied",
     "data_dir": os.path.join(SCRIPT_DIR, "Data_TV"), "shape_keywords": ["Turnhouts_Vennegebied", "Turnhout"]},
    {"code": "DM", "naam": "De Maten", "sheet_name": "De_Maten", "data_dir": os.path.join(SCRIPT_DIR, "Data_DM"),
     "shape_keywords": ["De_Maten", "Maten"]},
    {"code": "HB", "naam": "Heesbossen", "sheet_name": "Heesbossen", "data_dir": os.path.join(SCRIPT_DIR, "Data_HB"),
     "shape_keywords": ["Heesbossen", "Heesbos"]},
    {"code": "KH", "naam": "Kalmthoutse Heide", "sheet_name": "Kalmthoutse_Heide",
     "data_dir": os.path.join(SCRIPT_DIR, "Data_KH"), "shape_keywords": ["Kalmthoutse_Heide", "Kalmthout"]},
    {"code": "MH", "naam": "Mechelse Heide", "sheet_name": "Mechelse_Heide",
     "data_dir": os.path.join(SCRIPT_DIR, "Data_MH"), "shape_keywords": ["Mechelse_Heide", "Mechelse"]},
    {"code": "VS", "naam": "Voerstreek", "sheet_name": "Voerstreek", "data_dir": os.path.join(SCRIPT_DIR, "Data_VS"),
     "shape_keywords": ["Voerstreek", "Voer"]},
]

SPECIES_NL_COLUMNS = ["naam_nl", "nederlandse naam", "soort", "soortnaam"]
SPECIES_LAT_COLUMNS = ["naam_lat", "wetenschappelijke naam", "latijnse naam", "scientific_name"]
BEHAVIOR_COLUMNS = ["gedrag", "activity", "behavior"]

VALIDATION_COLUMN = "status"
DATE_COLUMN = "datum"
X_COORD_COLUMN = "x"
Y_COORD_COLUMN = "y"

CSV_CRS = "EPSG:31370"
TARGET_CRS = "EPSG:31370"

START_YEAR_EVAL = 2011
END_YEAR_EVAL = 2025
MIN_YEARS_HTS = 5
MIN_YEAR_GLOBAAL = 2011

# COMPLETE WINTERGASTENLIJST (Inclusief Bergeend, Blauwe kiekendief, Krakeend, etc.)
WINTER_SPECIES_RAW = {
    "bergeend", "tadorna tadorna",
    "wintertaling", "anas crecca",
    "smient", "mareca penelope",
    "tafeleend", "aythya ferina",
    "kuifeend", "aythya fuligula",
    "krakeend", "mareca strepera",
    "slobeend", "spatula clypeata",
    "pijlstaart", "anas acuta",
    "kolgans", "anser albifrons",
    "toendrarietgans", "anser fabalis",
    "kemphaan", "calidris pugnax",
    "regenwulp", "numenius phaeopus",
    "grote zilverreiger", "ardea alba", "egretta alba", "casmerodius albus",
    "blauwe kiekendief", "circus cyaneus"
}

WINTER_START = "10-15"
WINTER_END = "03-15"

WINTER_ALLOWED_BEHAVIOR = {
    'ter plaatse', 'foeragerend', 'jagend',
    'rustend', 'pleisterend', 'slaapplaats',
    '(kleur)ringdragend'
}

MAANDEN_DICT = {'jan': 1, 'feb': 2, 'mrt': 3, 'apr': 4, 'mei': 5, 'jun': 6, 'jul': 7, 'aug': 8, 'sep': 9, 'okt': 10,
                'nov': 11, 'dec': 12}


def normalize_name(name):
    if pd.isna(name): return ""
    clean_str = str(name).replace('\xa0', ' ')
    return re.sub(r"[^a-zA-Z0-9]", "", clean_str).lower()


WINTER_SPECIES_NORM = {normalize_name(s) for s in WINTER_SPECIES_RAW}


def parse_mixed_date_robust(series):
    """
    Parseert de datumkolom 100% correct:
    - DD/MM/YYYY of DD-MM-YYYY -> dayfirst=True (Europees)
    - YYYY-MM-DD -> ISO
    - Numerieke Excel datums -> origin='1899-12-30'
    """
    s_clean = series.astype(str).str.strip()

    numeric_mask = s_clean.str.match(r'^\d{5}(\.\d+)?$')
    iso_mask = s_clean.str.match(r'^\d{4}[-/]\d{1,2}[-/]\d{1,2}')
    euro_mask = s_clean.str.match(r'^\d{1,2}[-/]\d{1,2}[-/]\d{4}')

    res = pd.Series(index=series.index, dtype='datetime64[ns]')

    if numeric_mask.any():
        num_vals = pd.to_numeric(series[numeric_mask], errors='coerce')
        res.loc[numeric_mask] = pd.to_datetime(num_vals, unit='D', origin='1899-12-30', errors='coerce')

    if iso_mask.any():
        res.loc[iso_mask] = pd.to_datetime(series[iso_mask], format='mixed', dayfirst=False, errors='coerce')

    if euro_mask.any():
        res.loc[euro_mask] = pd.to_datetime(series[euro_mask], format='mixed', dayfirst=True, errors='coerce')

    unparsed_mask = res.isna() & series.notna() & (s_clean != '') & (s_clean != 'nan')
    if unparsed_mask.any():
        res.loc[unparsed_mask] = pd.to_datetime(series[unparsed_mask], format='mixed', dayfirst=True, errors='coerce')

    return res


def parse_dutch_date_to_mm_dd(date_str):
    parts = str(date_str).strip().lower().split()
    if len(parts) < 2: return None
    day_str = parts[0].zfill(2)
    month_text = parts[1][:3]
    if month_text in ['maa', 'maa']: month_text = 'mrt'
    month_code = MAANDEN_DICT.get(month_text)
    if not month_code:
        for k, v in MAANDEN_DICT.items():
            if k in month_text:
                month_code = v
                break
    if month_code and day_str.isdigit():
        return month_code, int(day_str)
    return None


def find_shapefile_for_area(area_config):
    if not os.path.exists(SHAPEFILES_DIR): return None
    shp_files = glob.glob(os.path.join(SHAPEFILES_DIR, "*.shp"))
    for kw in area_config["shape_keywords"]:
        for shp in shp_files:
            if kw.lower() in os.path.basename(shp).lower(): return shp
    return shp_files[0] if len(shp_files) == 1 else None


def read_file_in_chunks(file_path, chunksize=100000):
    ext = os.path.splitext(file_path)[1].lower()
    if ext == ".csv":
        with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
            first_line = f.readline()
        sep_char = ',' if first_line.count(',') > first_line.count(';') else ';'
        for chunk in pd.read_csv(file_path, chunksize=chunksize, sep=sep_char, low_memory=False):
            yield chunk
    elif ext in [".xlsx", ".xls"]:
        full_df = pd.read_excel(file_path)
        for i in range(0, len(full_df), chunksize):
            yield full_df.iloc[i:i + chunksize].copy()


def load_hts_for_area(excel_path, target_sheet):
    hts_set = set()
    if not os.path.exists(excel_path): return hts_set
    try:
        xls = pd.ExcelFile(excel_path)
        matched_sheet = next((s for s in xls.sheet_names if normalize_name(s) == normalize_name(target_sheet)), None)
        if matched_sheet:
            df_hts = pd.read_excel(excel_path, sheet_name=matched_sheet)
            df_hts.columns = df_hts.columns.astype(str).str.strip()
            gevonden_col = next(
                (c for c in df_hts.columns if c.lower() in SPECIES_NL_COLUMNS or c.lower() in SPECIES_LAT_COLUMNS),
                df_hts.columns[0] if not df_hts.empty else None)
            if gevonden_col:
                for val in df_hts[gevonden_col].dropna():
                    hts_set.add(normalize_name(val))
    except Exception as e:
        print(f"⚠️ Fout bij het lezen van HTS: {e}")
    return hts_set


def load_species_library_for_area(excel_path, gebied_config):
    df_species = pd.read_excel(excel_path)
    df_species.columns = df_species.columns.astype(str).str.strip()
    col_nl = next((c for c in df_species.columns if c.lower() in ["nederlandse naam", "naam_nl", "soort"]), None)
    col_lat = next((c for c in df_species.columns if
                    c.lower() in ["wetenschappelijke naam", "naam_lat", "latijnse naam", "scientific_name"]), None)

    area_col = next((c for c in df_species.columns if normalize_name(c) in [normalize_name(gebied_config["code"]),
                                                                            normalize_name(gebied_config["sheet_name"]),
                                                                            normalize_name(gebied_config["naam"])]),
                    None)
    df_filtered = df_species if not area_col else df_species[
        pd.to_numeric(df_species[area_col], errors='coerce') == 1].copy()

    species_library = {}
    for _, row in df_filtered.iterrows():
        raw_nl = str(row[col_nl]).strip() if pd.notna(row[col_nl]) else ""
        raw_lat = str(row[col_lat]).strip() if pd.notna(row[col_lat]) else ""
        if raw_lat and raw_nl:
            species_library[normalize_name(raw_lat)] = {
                "official_nl": raw_nl,
                "official_lat": raw_lat,
                "norm_nl": normalize_name(raw_nl),
                "synonyms_found": set([raw_nl])
            }
    return species_library


# ==============================================================================
# INLADEN METADATA
# ==============================================================================
print("=== LAAD HULP-METADATA IN ===")

raw_buffers = {}
if os.path.exists(EXCEL_AFSTANDEN_PATH):
    df_afstanden = pd.read_excel(EXCEL_AFSTANDEN_PATH)
    df_afstanden.columns = df_afstanden.columns.astype(str).str.strip()
    col_soort = df_afstanden.columns[1] if len(df_afstanden.columns) > 1 else df_afstanden.columns[0]
    col_buffer = df_afstanden.columns[5] if len(df_afstanden.columns) > 5 else df_afstanden.columns[-1]
    for _, row in df_afstanden.iterrows():
        raw_sname = row.get(col_soort)
        buf_val = row.get(col_buffer, 0)
        try:
            buf_val = 0.0 if pd.isna(buf_val) else float(buf_val)
        except ValueError:
            buf_val = 0.0
        if pd.notna(raw_sname):
            norm_s = normalize_name(raw_sname)
            if norm_s: raw_buffers[norm_s] = buf_val

raw_breeding = {}
if os.path.exists(EXCEL_BREEDING_SEASONS):
    df_breeding = pd.read_excel(EXCEL_BREEDING_SEASONS)
    df_breeding.columns = df_breeding.columns.astype(str).str.strip()
    col_br_s = next((c for c in df_breeding.columns if c.lower() in ["nederlandse naam", "soort", "naam"]),
                    df_breeding.columns[1])
    col_br_d = next(
        (c for c in df_breeding.columns if "datum" in c.lower() or "periode" in c.lower() or "broed" in c.lower()),
        df_breeding.columns[2])

    for _, row in df_breeding.iterrows():
        s_name = str(row[col_br_s]).strip()
        d_range = str(row[col_br_d]).strip()
        norm_s = normalize_name(s_name)
        if not norm_s: continue

        parts = re.split(r'[-–—]', d_range)
        if len(parts) == 2:
            res_start = parse_dutch_date_to_mm_dd(parts[0])
            res_end = parse_dutch_date_to_mm_dd(parts[1])
            if res_start and res_end:
                raw_breeding[norm_s] = {
                    "start_doy": pd.Timestamp(2023, res_start[0], res_start[1]).dayofyear,
                    "eind_doy": pd.Timestamp(2023, res_end[0], res_end[1]).dayofyear,
                    "raw_text": d_range
                }

raw_breeding["wulp"] = {"start_doy": pd.Timestamp(2023, 3, 10).dayofyear,
                        "eind_doy": pd.Timestamp(2023, 5, 31).dayofyear, "raw_text": "10 mrt - 31 mei"}
raw_breeding["boomleeuwerik"] = {"start_doy": pd.Timestamp(2023, 3, 5).dayofyear,
                                 "eind_doy": pd.Timestamp(2023, 6, 15).dayofyear, "raw_text": "5 mrt - 15 jun"}

print(
    f"✅ Broedvogel datumgrenzen verwerkt voor {len(raw_breeding)} soorten uit {os.path.basename(EXCEL_BREEDING_SEASONS)}.")
print("✅ Generieke hulp-metadata succesvol geladen.\n")

# ==============================================================================
# HOOFDLOOP PER GEBIED
# ==============================================================================
for gebied in GEBIEDEN_CONFIG:
    code, naam, sheet_name, data_dir = gebied["code"], gebied["naam"], gebied["sheet_name"], gebied["data_dir"]
    output_dir = os.path.join(SCRIPT_DIR, f"output_species_files_{code}")

    print("=" * 80)
    print(f"🌍 VERWERKEN GEBIED: {naam} ({code})")
    print("=" * 80)

    species_library = load_species_library_for_area(EXCEL_SPECIES_LIST, gebied)
    target_lat_keys = list(species_library.keys())

    if not os.path.exists(data_dir): continue
    input_files = glob.glob(os.path.join(data_dir, "*.csv")) + glob.glob(os.path.join(data_dir, "*.xlsx")) + glob.glob(
        os.path.join(data_dir, "*.xls"))
    if not input_files: continue

    shapefile_path = find_shapefile_for_area(gebied)
    if not shapefile_path or not os.path.exists(shapefile_path): continue

    grenzen_base = gpd.read_file(shapefile_path).to_crs(TARGET_CRS)
    gebied_union_base = grenzen_base.geometry.union_all()
    hts_exclusief_set = load_hts_for_area(EXCEL_HTS_PATH, sheet_name)

    species_meta = {}
    breeding_seasons = {}

    for lat_k, data in species_library.items():
        buffer_m = 0
        for norm_k in (lat_k, data["norm_nl"]):
            if norm_k in raw_buffers:
                buffer_m = raw_buffers[norm_k]
                break

        matched_br = None
        for norm_k in (lat_k, data["norm_nl"]):
            if norm_k in raw_breeding:
                matched_br = raw_breeding[norm_k]
                break

        if matched_br:
            breeding_seasons[lat_k] = matched_br

        is_hts = (lat_k in hts_exclusief_set) or (data["norm_nl"] in hts_exclusief_set)
        effective_buffer = buffer_m if buffer_m > 0 else 10.0
        species_meta[lat_k] = {
            "official_nl": data["official_nl"],
            "official_lat": data["official_lat"],
            "buffer_m": buffer_m,
            "is_hts": is_hts,
            "geometry": grenzen_base.buffer(effective_buffer).union_all()
        }

    os.makedirs(output_dir, exist_ok=True)
    species_collected_data = {lat_k: [] for lat_k in species_library.keys()}
    species_collected_wulp_wv = []
    template_columns = None

    for file_path in input_files:
        filename = os.path.basename(file_path)
        print(f"  📖 Lezen: {filename}...")

        try:
            for chunk in read_file_in_chunks(file_path):
                chunk.columns = chunk.columns.astype(str).str.strip()

                if template_columns is None and not chunk.empty:
                    template_columns = [c for c in chunk.columns if
                                        c not in ['_matched_lat_key', '_mm_dd', '_day_of_year', '_norm_data_lat',
                                                  '_norm_data_nl', 'jaar']]

                if VALIDATION_COLUMN not in chunk.columns or DATE_COLUMN not in chunk.columns: continue

                clean_validation = chunk[VALIDATION_COLUMN].astype(str).str.strip().str.lower()
                chunk = chunk[clean_validation.str.startswith("goedgekeurd", na=False)]
                if chunk.empty: continue

                # DATUM PARSING MET ROBUUSTE FIX
                chunk[DATE_COLUMN] = parse_mixed_date_robust(chunk[DATE_COLUMN])
                chunk = chunk.dropna(subset=[DATE_COLUMN])
                chunk['jaar'] = chunk[DATE_COLUMN].dt.year
                chunk = chunk[chunk['jaar'] >= MIN_YEAR_GLOBAAL]
                if chunk.empty: continue

                chunk[X_COORD_COLUMN] = pd.to_numeric(chunk[X_COORD_COLUMN].astype(str).str.replace(',', '.'),
                                                      errors='coerce')
                chunk[Y_COORD_COLUMN] = pd.to_numeric(chunk[Y_COORD_COLUMN].astype(str).str.replace(',', '.'),
                                                      errors='coerce')
                chunk = chunk.dropna(subset=[X_COORD_COLUMN, Y_COORD_COLUMN])
                if chunk.empty: continue

                geometry = [Point(xy) for xy in zip(chunk[X_COORD_COLUMN], chunk[Y_COORD_COLUMN])]
                chunk_gpd = gpd.GeoDataFrame(chunk, geometry=geometry, crs=CSV_CRS).to_crs(TARGET_CRS)

                data_col_nl = next((c for c in chunk_gpd.columns if c.lower() in SPECIES_NL_COLUMNS), None)
                data_col_lat = next((c for c in chunk_gpd.columns if c.lower() in SPECIES_LAT_COLUMNS), None)
                data_col_gedrag = next((c for c in chunk_gpd.columns if c.lower() in BEHAVIOR_COLUMNS), None)

                chunk_gpd['_matched_lat_key'] = None

                if data_col_lat:
                    chunk_gpd['_norm_data_lat'] = chunk_gpd[data_col_lat].apply(normalize_name)
                    for lat_k in target_lat_keys:
                        mask = chunk_gpd['_matched_lat_key'].isna() & (chunk_gpd['_norm_data_lat'] == lat_k)
                        chunk_gpd.loc[mask, '_matched_lat_key'] = lat_k
                        if data_col_nl and mask.any():
                            species_library[lat_k]["synonyms_found"].update(
                                chunk_gpd.loc[mask, data_col_nl].dropna().unique())

                if data_col_nl:
                    chunk_gpd['_norm_data_nl'] = chunk_gpd[data_col_nl].apply(normalize_name)
                    for lat_k, data in species_library.items():
                        mask = chunk_gpd['_matched_lat_key'].isna() & (chunk_gpd['_norm_data_nl'] == data["norm_nl"])
                        chunk_gpd.loc[mask, '_matched_lat_key'] = lat_k

                chunk_gpd = chunk_gpd[chunk_gpd['_matched_lat_key'].notna()]
                if chunk_gpd.empty: continue

                spatial_keep_mask = pd.Series(False, index=chunk_gpd.index)
                for lat_k, group_data in chunk_gpd.groupby('_matched_lat_key'):
                    target_poly = species_meta[lat_k]["geometry"]
                    spatial_keep_mask.loc[group_data.index[group_data.geometry.intersects(target_poly)]] = True

                chunk_gpd = chunk_gpd[spatial_keep_mask]
                if chunk_gpd.empty: continue

                chunk_gpd['_mm_dd'] = chunk_gpd[DATE_COLUMN].dt.strftime('%m-%d')
                chunk_gpd['_day_of_year'] = chunk_gpd[DATE_COLUMN].dt.dayofyear

                # A. WULP WINTERGAST LOGICA
                is_wulp_key = normalize_name("Numenius arquata")
                is_wulp_rows = chunk_gpd['_matched_lat_key'] == is_wulp_key
                is_in_winter_period = (chunk_gpd['_mm_dd'] >= WINTER_START) | (chunk_gpd['_mm_dd'] <= WINTER_END)

                if data_col_gedrag:
                    raw_g_wulp = chunk_gpd[data_col_gedrag]
                    clean_g_wulp = raw_g_wulp.astype(str).str.strip().str.lower()
                    g_ok_wulp = clean_g_wulp.isin(WINTER_ALLOWED_BEHAVIOR) | raw_g_wulp.isna() | clean_g_wulp.isin(
                        ['nan', 'none', '', 'onbekend'])
                    wulp_wv_mask = is_wulp_rows & is_in_winter_period & g_ok_wulp
                else:
                    wulp_wv_mask = is_wulp_rows & is_in_winter_period

                wulp_winter_chunk = chunk_gpd[wulp_wv_mask].copy()
                if not wulp_winter_chunk.empty:
                    species_collected_wulp_wv.append(pd.DataFrame(wulp_winter_chunk.drop(
                        columns=['geometry', '_matched_lat_key', '_mm_dd', '_day_of_year', 'jaar', '_norm_data_lat',
                                 '_norm_data_nl'], errors='ignore')))

                # B. ALGEMENE SEIZOENSFILTERING
                vogel_keep_mask = pd.Series(False, index=chunk_gpd.index)

                for lat_k, group_data in chunk_gpd.groupby('_matched_lat_key'):
                    norm_nl = species_library[lat_k]["norm_nl"]

                    PLANT_KEYWORDS = ["ranunculus", "sparganium", "waterranonkel", "egelskop", "gentiaandebout",
                                      "zonnedauw"]
                    is_plant = any(kw in lat_k for kw in PLANT_KEYWORDS) or any(kw in norm_nl for kw in PLANT_KEYWORDS)

                    is_wintergast = (lat_k in WINTER_SPECIES_NORM) or (norm_nl in WINTER_SPECIES_NORM)
                    has_breeding = lat_k in breeding_seasons

                    if is_plant:
                        vogel_keep_mask.loc[group_data.index] = True

                    elif lat_k == is_wulp_key:
                        start_doy = breeding_seasons[lat_k]["start_doy"] if lat_k in breeding_seasons else pd.Timestamp(
                            2023, 3, 10).dayofyear
                        eind_doy = breeding_seasons[lat_k]["eind_doy"] if lat_k in breeding_seasons else pd.Timestamp(
                            2023, 5, 31).dayofyear
                        doy_s = group_data['_day_of_year']
                        vogel_keep_mask.loc[group_data.index] = (doy_s >= start_doy) & (doy_s <= eind_doy)

                    elif is_wintergast:
                        m_date = is_in_winter_period.loc[group_data.index]
                        if data_col_gedrag:
                            raw_g = group_data[data_col_gedrag]
                            clean_g = raw_g.astype(str).str.strip().str.lower()
                            m_gedrag = clean_g.isin(WINTER_ALLOWED_BEHAVIOR) | raw_g.isna() | clean_g.isin(
                                ['nan', 'none', '', 'onbekend'])
                            vogel_keep_mask.loc[group_data.index] = m_date & m_gedrag
                        else:
                            vogel_keep_mask.loc[group_data.index] = m_date

                    elif has_breeding:
                        dates = breeding_seasons[lat_k]
                        s_doy, e_doy = dates["start_doy"], dates["eind_doy"]
                        doy_s = group_data['_day_of_year']
                        if s_doy <= e_doy:
                            vogel_keep_mask.loc[group_data.index] = (doy_s >= s_doy) & (doy_s <= e_doy)
                        else:
                            vogel_keep_mask.loc[group_data.index] = (doy_s >= s_doy) | (doy_s <= e_doy)

                    else:
                        vogel_keep_mask.loc[group_data.index] = True

                chunk_gpd = chunk_gpd[vogel_keep_mask]
                if chunk_gpd.empty: continue

                for lat_k, group_data in chunk_gpd.groupby('_matched_lat_key'):
                    clean_data = pd.DataFrame(group_data.drop(
                        columns=['geometry', '_matched_lat_key', '_mm_dd', '_day_of_year', '_norm_data_lat',
                                 '_norm_data_nl'], errors='ignore'))
                    species_collected_data[lat_k].append(clean_data)

        except Exception as e:
            print(f"⚠️ Fout bij verwerken van bestand {filename}: {e}")

    # EVALUATIE & OUTPUT
    print(f"\n📊 Evaluatie van waarnemingsjaren ({START_YEAR_EVAL}–{END_YEAR_EVAL}) voor {naam}...")
    overzicht_resultaten = []

    for lat_k, list_of_dfs in species_collected_data.items():
        meta = species_meta[lat_k]
        official_nl = meta["official_nl"]
        is_hts = meta["is_hts"]

        full_species_df = pd.concat(list_of_dfs, ignore_index=True) if list_of_dfs else pd.DataFrame()

        if not full_species_df.empty and "jaar" in full_species_df.columns:
            eval_df = full_species_df[
                (full_species_df["jaar"] >= START_YEAR_EVAL) & (full_species_df["jaar"] <= END_YEAR_EVAL)]
            unieke_jaren = sorted(eval_df["jaar"].dropna().unique().astype(int).tolist())
            aantal_jaren = len(unieke_jaren)
            totaal_waarnemingen = len(full_species_df)
        else:
            unieke_jaren, aantal_jaren, totaal_waarnemingen = [], 0, 0

        export_goedgekeurd = False
        reden = ""

        if is_hts:
            if aantal_jaren >= MIN_YEARS_HTS:
                export_goedgekeurd = True
                reden = f"HTS-Exclusief GOEDGEKEURD ({aantal_jaren}/{MIN_YEARS_HTS}+ unieke jaren)"
            else:
                export_goedgekeurd = False
                reden = f"HTS-Exclusief AFGEKEURD ({aantal_jaren}/{MIN_YEARS_HTS} unieke jaren)"
        else:
            export_goedgekeurd = True
            reden = f"Niet-HTS: Geëxporteerd ({totaal_waarnemingen} waarnemingen)" if totaal_waarnemingen > 0 else "Niet-HTS: Geëxporteerd (0 waarnemingen - leeg bestand)"

        if export_goedgekeurd:
            title_case_name = str(official_nl).title()
            safe_filename = "Waarnemingen_" + "".join([c for c in title_case_name if c.isalnum()]) + ".csv"
            out_path = os.path.join(output_dir, safe_filename)

            clean_export_df = full_species_df.drop(columns=['jaar'],
                                                   errors='ignore') if not full_species_df.empty else pd.DataFrame(
                columns=template_columns if template_columns else ["datum", "status", "x", "y", "naam_nl", "naam_lat"])
            clean_export_df.to_csv(out_path, index=False, sep=";", decimal=".")

        overzicht_resultaten.append({
            "Nederlandse_Naam": official_nl,
            "Wetenschappelijke_Naam": meta["official_lat"],
            "Gevonden_Synoniemen": ", ".join(sorted(species_library[lat_k]["synonyms_found"])),
            "Categorie": "HTS Exclusief" if is_hts else "SBZ / SBP Doelsoort",
            "Dispersiebuffer_m": meta["buffer_m"],
            "Aantal_Unieke_Jaren_2011_2025": aantal_jaren,
            "Unieke_Jaren_Lijst": ", ".join(map(str, unieke_jaren)),
            "Totaal_Waarnemingen_Sinds_2011": totaal_waarnemingen,
            "CSV_Gegenereerd": "JA" if export_goedgekeurd else "NEE",
            "Reden": reden
        })

    if species_collected_wulp_wv:
        full_wulp_wv = pd.concat(species_collected_wulp_wv, ignore_index=True)
        full_wulp_wv.to_csv(os.path.join(output_dir, "Waarnemingen_Wulp_wv.csv"), index=False, sep=";", decimal=".")

    df_overzicht = pd.DataFrame(overzicht_resultaten)
    df_overzicht["Unieke_Jaren_Lijst"] = df_overzicht["Unieke_Jaren_Lijst"].astype(object)
    df_overzicht["Reden"] = df_overzicht["Reden"].astype(object)
    df_overzicht = df_overzicht.sort_values(by=["CSV_Gegenereerd", "Categorie", "Nederlandse_Naam"],
                                            ascending=[False, True, True])

    overzicht_path = os.path.join(SCRIPT_DIR, f"soorten_jaren_overzicht_{code}.csv")
    df_overzicht.to_csv(overzicht_path, index=False, sep=";")

    print(f"📋 Overzichtstabel opgeslagen in: {overzicht_path}")
    print(f"🎉 Succes voor {naam}! Soort-CSV's opgeslagen in: {output_dir}\n")

print("✨ ALLE GEBIEDEN ZIJN SUCCESVOL VERWERKT!")
