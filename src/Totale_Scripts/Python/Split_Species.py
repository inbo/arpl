import os
import pandas as pd
import re

# --- PATH SETUP ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "output_species_files")
EXCEL_SPECIES_LIST = os.path.join(SCRIPT_DIR, "Soortenlijst_TV.xlsx")
EXCEL_BREEDING_SEASONS = os.path.join(SCRIPT_DIR, "Broedvogels_Datums.xlsx")

INPUT_FILES = [
    os.path.join(SCRIPT_DIR, "INBODATAVR-481_dumptemp20260531_met_toestemming.csv"),
    os.path.join(SCRIPT_DIR, "INBODATAVR-481_dumptemp20260531_sws.csv"),
    os.path.join(SCRIPT_DIR, "INBODATAVR_494_basis_wbe_dumptem20260705_met_toestemming.csv"),
    os.path.join(SCRIPT_DIR, "INBODATAVR_494_basis_wbe_dumptem20260705_sws.csv")
]

# --- COLUMN NAME CONFIGURATION ---
SPECIES_COLUMN = "naam_nl"
VALIDATION_COLUMN = "status"
SPECIES_GROUP_COLUMN = "soortgroep"
DATE_COLUMN = "datum"

# --- DATE FILTERS PER SPECIES GROUP ---
PERIOD_FILTERS = {
    # Eventuele algemene groepsfilters (Vogels hier weglaten)
}

# --- GENERIEKE JAARGRENS (Voor álle soorten) ---
MIN_YEAR = 2020

# --- SPECIFIEKE WINTERGAST-FILTER ---
WINTER_SPECIES = {
    "Bergeend", "Blauwe kiekendief", "Grote zilverreiger", "Kemphaan",
    "Krakeend", "Kuifeend", "Pijlstaart", "Regenwulp", "Slobeend",
    "Smient", "Tafeleend", "Wintertaling"
}
WINTER_START = "10-15"
WINTER_END = "03-15"


# --- HULPFUNCTIE: VERTAAL NEDERLANDSE DATUM NAAR MM-DD ---
def parse_dutch_date_to_mm_dd(date_str):
    """
    Zet een Nederlandse datumstring zoals '15 feb' of '1 apr' om naar 'MM-DD'.
    """
    months_mapping = {
        'jan': '01', 'feb': '02', 'mrt': '03', 'apr': '04',
        'mei': '05', 'jun': '06', 'jul': '07', 'aug': '08',
        'sep': '09', 'okt': '10', 'nov': '11', 'dec': '12'
    }

    # Maak de string schoon en splits op spatie (bijv. ['15', 'feb'])
    parts = date_str.strip().lower().split()
    if len(parts) < 2:
        return None

    day_str = parts[0].zfill(2)  # '5' wordt '05', '15' blijft '15'
    month_text = parts[1][:3]  # Pak de eerste 3 letters (bijv. 'maart' -> 'maa')

    # Correctie voor maart/mrt en mei
    if month_text == 'maa':
        month_text = 'mrt'
    elif month_text == 'mei':
        month_text = 'mei'

    month_code = months_mapping.get(month_text)

    if not month_code:
        # Als de maand niet herkend wordt, probeer een fallback
        for key, val in months_mapping.items():
            if key in month_text:
                month_code = val
                break

    if month_code and day_str.isdigit():
        return f"{month_code}-{day_str}"
    return None


os.makedirs(OUTPUT_DIR, exist_ok=True)
print(f"Output directory aangemaakt op: {OUTPUT_DIR}\n")

# --- LEES DE EXCEL SOORTENLIJST ---
if not os.path.exists(EXCEL_SPECIES_LIST):
    raise FileNotFoundError(
        f"⚠️ Kritieke fout: Kan de soortenlijst '{os.path.basename(EXCEL_SPECIES_LIST)}' niet vinden in {SCRIPT_DIR}.")

df_target_species = pd.read_excel(EXCEL_SPECIES_LIST)
target_species = df_target_species['Nederlandse naam'].dropna().astype(str).str.strip().unique().tolist()
print(f"📚 {len(target_species)} doelsoorten succesvol ingeladen uit Excel.\n")

target_species_sorted = sorted(target_species, key=len, reverse=True)

# --- LAAD EN PARSE DE BROEDVOGEL DATUMS ---
breeding_seasons = {}
if os.path.exists(EXCEL_BREEDING_SEASONS):
    print(f"📅 Broedvogel datums laden en vertalen uit {os.path.basename(EXCEL_BREEDING_SEASONS)}...")

    df_breeding = pd.read_excel(EXCEL_BREEDING_SEASONS)

    if 'Nederlandse Naam' in df_breeding.columns and 'Datumgrenzen' in df_breeding.columns:
        df_breeding_clean = df_breeding.dropna(subset=['Nederlandse Naam', 'Datumgrenzen'])

        for _, row in df_breeding_clean.iterrows():
            species_name = str(row['Nederlandse Naam']).strip()
            date_range_str = str(row['Datumgrenzen']).strip()

            # Verwacht formaat splitsen: "15 feb - 15 jun"
            if "-" in date_range_str:
                start_raw, end_raw = date_range_str.split("-")
                start_formatted = parse_dutch_date_to_mm_dd(start_raw)
                end_formatted = parse_dutch_date_to_mm_dd(end_raw)

                if start_formatted and end_formatted:
                    breeding_seasons[species_name] = {
                        "start": start_formatted,
                        "end": end_formatted
                    }
                else:
                    print(f"  ⚠️ Kon datums niet ontcijferen voor {species_name}: '{date_range_str}'")
            else:
                print(f"  ⚠️ Ongeldig formaat (missend '-' streepje) voor {species_name}: '{date_range_str}'")

        print(f"✅ {len(breeding_seasons)} broedvogel-seizoenen succesvol ingeladen en vertaald.\n")
    else:
        print("⚠️ Fout: Kolom 'Nederlandse Naam' of 'Datumgrenzen' niet gevonden in het broedvogelbestand!\n")
else:
    print(f"⚠️ Waarschuwing: Broedvogelbestand '{os.path.basename(EXCEL_BREEDING_SEASONS)}' niet gevonden.\n")

# --- PROCESS BOTH FILES ---
for file_path in INPUT_FILES:
    filename = os.path.basename(file_path)
    if not os.path.exists(file_path):
        print(f"⚠️ Warning: Could not find file {filename}. Skipping.")
        continue

    print(f"🔄 Processing {filename}...")
    chunk_count = 0

    for chunk in pd.read_csv(file_path, chunksize=100000, low_memory=False, decimal=','):
        chunk_count += 1

        # ----------------------------------------------------
        # FILTER 1: Alleen goedgekeurde waarnemingen
        # ----------------------------------------------------
        clean_validation = chunk[VALIDATION_COLUMN].astype(str).str.strip().str.lower()
        chunk = chunk[clean_validation.str.startswith("goedgekeurd", na=False)]

        if chunk.empty:
            continue

        # ----------------------------------------------------
        # FILTER 2: Datumfilters & Generieke Jaargrens (>= 2020)
        # ----------------------------------------------------
        chunk[DATE_COLUMN] = pd.to_datetime(chunk[DATE_COLUMN], errors='coerce')

        # A) Directe harde jaargrens toepassen voor ALLE waarnemingen
        chunk = chunk[chunk[DATE_COLUMN].dt.year >= MIN_YEAR]

        if chunk.empty:
            continue

        # Maak een masker om rijen bij te houden die we willen bewaren op basis van seizoenen
        keep_mask = pd.Series(True, index=chunk.index)

        # B) Dag/maand-filters toepassen (voor groepen anders dan Vogels)
        chunk['_mm_dd'] = chunk[DATE_COLUMN].dt.strftime('%m-%d')
        for group, dates in PERIOD_FILTERS.items():
            is_this_group = chunk[SPECIES_GROUP_COLUMN] == group
            is_in_period = (chunk['_mm_dd'] >= dates["start"]) & (chunk['_mm_dd'] <= dates["end"])
            keep_mask.loc[is_this_group] = is_in_period

        # Filter de chunk op basis van de algemene periodemaskers
        chunk = chunk[keep_mask]

        if chunk.empty:
            chunk = chunk.drop(columns=['_mm_dd'])
            continue

        # ----------------------------------------------------
        # CLEANING & ROBUUST MATCHEN MET EXCEL LIJST
        # ----------------------------------------------------
        chunk['_raw_species_str'] = chunk[SPECIES_COLUMN].astype(str).str.strip()
        chunk['_matched_target_species'] = None

        for species in target_species_sorted:
            escaped_species = re.escape(species)
            pattern = rf"^{escaped_species}(?:$|[^a-zA-Z0-9])"

            mask = chunk['_matched_target_species'].isna() & chunk['_raw_species_str'].str.contains(pattern, regex=True,
                                                                                                    case=False,
                                                                                                    na=False)
            chunk.loc[mask, '_matched_target_species'] = species

        chunk = chunk[chunk['_matched_target_species'].notna()]

        if chunk.empty:
            chunk = chunk.drop(columns=['_mm_dd'])
            continue

        # ----------------------------------------------------
        # FILTER 3: Wintergast- en Broedvogel-filters (Vogels)
        # ----------------------------------------------------

        # We splitsen de dataset hier op in een regulier deel en een speciaal "Wulp winter" deel
        is_wulp = chunk['_matched_target_species'] == "Wulp"
        is_in_winter_period = (chunk['_mm_dd'] >= WINTER_START) | (chunk['_mm_dd'] <= WINTER_END)

        # --- DEEL A: Speciaal Wulp Wintergast export ---
        wulp_winter_chunk = chunk[is_wulp & is_in_winter_period].copy()
        if not wulp_winter_chunk.empty:
            # Schrijf Wulp winter direct weg naar zijn eigen bestand om de hoofdloop niet te verstoren
            species_file_path = os.path.join(OUTPUT_DIR, "Waarnemingen_Wulp_wv.csv")
            clean_wulp_data = wulp_winter_chunk.drop(columns=['_raw_species_str', '_matched_target_species', '_mm_dd'])

            if not os.path.exists(species_file_path):
                clean_wulp_data.to_csv(species_file_path, index=False, sep=";", decimal=".")
            else:
                clean_wulp_data.to_csv(species_file_path, index=False, mode='a', header=False, sep=";", decimal=".")

        # --- DEEL B: Reguliere verwerking (inclusief Wulp als broedvogel) ---
        vogel_keep_mask = pd.Series(True, index=chunk.index)

        # 1. WINTERVOGEL CHECK (Wulp zit hier niet in, dus wordt hier niet gefilterd als wintervogel)
        is_winter_species = chunk['_matched_target_species'].isin(WINTER_SPECIES)
        vogel_keep_mask.loc[is_winter_species] = is_in_winter_period

        # 2. BROEDVOGEL CHECK (Inclusief Wulp!)
        if breeding_seasons:
            for species_name, dates in breeding_seasons.items():
                if species_name not in WINTER_SPECIES:  # Wulp mag hier doorheen, want hij zit niet in WINTER_SPECIES
                    is_this_species = chunk['_matched_target_species'] == species_name
                    is_in_breeding_season = (chunk['_mm_dd'] >= dates["start"]) & (chunk['_mm_dd'] <= dates["end"])
                    vogel_keep_mask.loc[is_this_species] = is_in_breeding_season

        # Pas de gecombineerde vogel-filter toe op de chunk
        chunk = chunk[vogel_keep_mask]

        # Nu kunnen we de tijdelijke datum-kolom veilig weggooien
        chunk = chunk.drop(columns=['_mm_dd'])

        if chunk.empty:
            continue

        # ----------------------------------------------------
        # GROUP BY TARGET SPECIES AND WRITE TO CSV
        # ----------------------------------------------------
        grouped = chunk.groupby('_matched_target_species')

        for official_species_name, group_data in grouped:
            title_case_name = str(official_species_name).title()
            safe_filename = "Waarnemingen_" + "".join([c for c in title_case_name if c.isalnum()])

            species_file_path = os.path.join(OUTPUT_DIR, f"{safe_filename}.csv")

            clean_group_data = group_data.drop(columns=['_raw_species_str', '_matched_target_species'])

            if not os.path.exists(species_file_path):
                clean_group_data.to_csv(species_file_path, index=False, sep=";", decimal=".")
            else:
                clean_group_data.to_csv(species_file_path, index=False, mode='a', header=False, sep=";", decimal=".")

    print(f"✅ Finished {filename} ({chunk_count} chunks processed).")

print(f"\n🎉 Succes! Alle waarnemingen die matchen met je Excel-lijst zijn opgeslagen in: {OUTPUT_DIR}")