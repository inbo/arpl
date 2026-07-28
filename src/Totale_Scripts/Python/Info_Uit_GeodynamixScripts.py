import os
import csv
import re
import pandas as pd  # Nodig voor het lezen van Excel

# --- CONFIGURATIE ---
folder_path = './GDX_Scripts'
legenda_file = './Soorten_bwk_afstanden.xlsx'  # Het nieuwe Excel bestand
output_folder = './Resultaten_Per_Groep'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# 1. LEES DE LEGENDA IN VANUIT EXCEL
legenda_data = {}

try:
    # Lees het Excel bestand in
    df_legenda = pd.read_excel(legenda_file)

    # Normaliseer kolomnamen: spaties weg, alles kleine letters
    # Dit zorgt dat 'MinOpp (ha)' gevonden wordt als we zoeken op 'minopp (ha)'
    df_legenda.columns = [str(col).strip().lower() for col in df_legenda.columns]

    # Zet de dataframe om naar een lijst van dictionaries
    records = df_legenda.to_dict('records')
    for row in records:
        # We gebruiken de kolom 'soort' als sleutel voor de koppeling
        s = str(row.get('soort', '')).strip()
        if s and s != 'nan':
            legenda_data[s] = row

    print(f"✅ Legenda ingeladen uit Excel: {len(legenda_data)} soorten herkend.")
except Exception as e:
    print(f"❌ Fout bij laden Excel-legenda: {e}")
    print("Check of de kolomnamen 'Soort' en 'Groep' aanwezig zijn in het Excel-bestand.")
    exit()

# 2. DATA VERZAMELEN
pattern = re.compile(r'assign\s*\(\s*([^,_]+)([^,]*)\s*,\s*bwk\s*\(\s*(.*?)\s*\)\s*\)', re.DOTALL)
data_per_groep = {}

bestanden = [f for f in os.listdir(folder_path) if f.lower().endswith('.txt')]
print(f"Bezig met het verwerken van {len(bestanden)} bestanden...")

for filename in bestanden:
    file_path = os.path.join(folder_path, filename)
    try:
        with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
            matches = pattern.findall(content)

            if not matches:
                continue

            # Master Soort Logica
            master_soort = matches[0][0].strip()

            # Haal info op uit de Excel-data (via de master_soort)
            info = legenda_data.get(master_soort, {})

            # Bepaal de groep (gebruik 'OVERIG' als soort niet in Excel staat)
            groep = str(info.get('groep', 'OVERIG')).strip()
            if groep == 'nan' or groep == '':
                groep = 'OVERIG'

            if groep not in data_per_groep:
                data_per_groep[groep] = []

            for match in matches:
                gevonden_naam = match[0].strip()
                type_suffix = match[1].strip().lstrip('_')
                codes_blok = match[2].strip()

                # Bepaal het Type veld (not_ logica)
                if gevonden_naam == master_soort:
                    final_type = type_suffix
                else:
                    final_type = f"not_{gevonden_naam}"
                    if type_suffix:
                        final_type += f"_{type_suffix}"

                individuele_codes = [c.strip().strip("'").strip('"') for c in codes_blok.split(',') if c.strip()]

                for code in individuele_codes:
                    # Wildcard Logica (%)
                    if code.endswith('%'):
                        match_type = 'bevat'
                        schone_code = code.rstrip('%')
                    else:
                        match_type = 'exact'
                        schone_code = code

                    # Rij samenstellen met data uit de Excel-legenda
                    data_per_groep[groep].append({
                        'Soort': master_soort,
                        'Type': final_type,
                        'Code': schone_code,
                        'Match': match_type,
                        # Koppeling naar de specifieke Excel-waarden
                        'MinOpp (ha)': info.get('minopp (ha)', ''),
                        'AfstandBiotopen (m)': info.get('afstandbiotopen (m)', ''),
                        'Dispersiecap (m)': info.get('dispersiecap (m)', '')
                    })
    except Exception as e:
        print(f"Fout in {filename}: {e}")

# 3. OPSLAAN PER GROEP (CSV)
for groep, rijen in data_per_groep.items():
    output_file = os.path.join(output_folder, f'biotopen_{groep}.csv')
    with open(output_file, 'w', newline='', encoding='utf-8') as csvfile:
        fieldnames = [
            'Soort', 'Type', 'Code', 'Match',
            'MinOpp (ha)', 'AfstandBiotopen (m)', 'Dispersiecap (m)'
        ]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rijen)

    status = "⚠️ CHECK DIT BESTAND!" if groep == "OVERIG" else "💾 Opgeslagen"
    print(f"{status} {groep}: {len(rijen)} regels.")

print(f"\nKlaar! De koppeling met '{legenda_file}' is voltooid.")