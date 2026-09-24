import pandas as pd
from pathlib import Path

# 1. Paden instellen
input_file = Path("P-30-001607_06_Veldkrekel_wbe_VL_DumpTem20260913.csv")
output_file = Path("Waarnemingen_Veldkrekel.csv")

print(f"=== START HERSTEL COÖRDINATEN & OPSCHONEN: {input_file.name} ===")

# CSV inlezen (probeer eerst ;, anders auto)
try:
    df = pd.read_csv(input_file, sep=';', engine='python')
except Exception:
    df = pd.read_csv(input_file, sep=None, engine='python')

df.columns = df.columns.str.strip().str.lower()
start_count = len(df)

# ------------------------------------------------------------------------------
# 1. HARD FILTER: Gedrag, Methode & Precisie
# ------------------------------------------------------------------------------
ongewenst_gedrag = ["dood"]
ongewenste_methode = ["binnenshuis", "uitgekweekt en losgelaten"]

mask_ongewenst = (
    df['gedrag'].astype(str).str.lower().isin(ongewenst_gedrag) |
    df['methode'].astype(str).str.lower().isin(ongewenste_methode) |
    (df['precisie'] > 100)
)

df_clean = df[~mask_ongewenst].copy()

# ------------------------------------------------------------------------------
# 2. HERSTEL EN REINIGING COÖRDINATEN (X en Y)
# ------------------------------------------------------------------------------
def fix_coördinaat(val, max_grens):
    if pd.isna(val):
        return None
    # Omzetten naar string en eventuele komma's vervangen door punten
    s = str(val).replace(',', '.').strip()
    try:
        f = float(s)
        # Als het getal veel te groot is (geen decimaal herkend), schalen
        while f > max_grens:
            f /= 10.0
        return f
    except ValueError:
        return None

# Lambert 72 X-coördinaten in Vlaanderen liggen tussen 20.000 en 260.000
# Lambert 72 Y-coördinaten in Vlaanderen liggen tussen 150.000 en 250.000
df_clean['x'] = df_clean['x'].apply(lambda v: fix_coördinaat(v, 300000.0))
df_clean['y'] = df_clean['y'].apply(lambda v: fix_coördinaat(v, 300000.0))

df_clean = df_clean.dropna(subset=['x', 'y'])

# ------------------------------------------------------------------------------
# 3. EXPORT MET PUNT ALS DECIMAALTEKEN
# ------------------------------------------------------------------------------
output_file.parent.mkdir(parents=True, exist_ok=True)

# float_format='%.2f' dwingt een punt af als decimaalscheiding
df_clean.to_csv(output_file, sep=';', index=False, float_format='%.2f')

print("--------------------------------------------------")
print(f" Resultaat : {len(df_clean)} van {start_count} waarnemingen over.")
print(f" Extensie X : Min = {df_clean['x'].min():.2f} | Max = {df_clean['x'].max():.2f}")
print(f" Extensie Y : Min = {df_clean['y'].min():.2f} | Max = {df_clean['y'].max():.2f}")
print("--------------------------------------------------")
print(f"🎉 Correct opgeslagen als: {output_file}")
