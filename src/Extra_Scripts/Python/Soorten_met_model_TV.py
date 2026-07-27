import os
import pandas as pd

# ==============================================================================
# CONFIGURATIE: Pas hier de bestandsnamen aan indien nodig
# ==============================================================================
hoofdbestand = "Soortenlijst_Maatwerkgebieden.xlsx"  # De Excel met alle tabbladen
referentiebestand = "Soorten_bwk_afstanden.xlsx"  # De Excel met model- en filterinfo
output_bestand = "Maatwerkgebieden_Ingevuld.xlsx"


def clean_species_name(name):
    """Helperfunctie om soortnamen te schonen (kleine letters, geen spaties)"""
    if pd.isna(name):
        return ""
    return str(name).replace(" ", "").lower().strip()


def main():
    print("-> Script opstarten...")

    # 1. Referentiebestand inlezen (Soorten_bwk_afstanden.xlsx)
    if not os.path.exists(referentiebestand):
        print(f"Fout: {referentiebestand} niet gevonden!")
        return

    df_ref = pd.read_excel(referentiebestand)

    # We bouwen twee 'sets' op voor een razendsnelle check achteraf
    model_soorten = set()
    automatisch_soorten = set()

    # Controleer of de kolom 'Soort' bestaat (index 1 / 2e kolom)
    if "Soort" in df_ref.columns:
        soort_kolom_ref = "Soort"
    else:
        # Fallback naar de 2e kolom op basis van positie
        soort_kolom_ref = df_ref.columns[1]

    # Pakt kolom G (7e kolom = index 6) voor de Automatisch-check
    kolom_g_naam = df_ref.columns[6]

    # Loop door het referentiebestand om de status per soort te bepalen
    for _, row in df_ref.iterrows():
        geschoonde_soort = clean_species_name(row[soort_kolom_ref])
        if geschoonde_soort:
            # Elke soort in deze lijst heeft in ieder casus een model
            model_soorten.add(geschoonde_soort)

            # Check of kolom G de waarde 'Simpel' bevat (ongevoelig voor hoofdletters/spaties)
            val_g = str(row[kolom_g_naam]).strip().lower()
            if val_g == "simpel":
                automatisch_soorten.add(geschoonde_soort)

    # 2. Hoofdbestand inlezen (alle tabbladen)
    if not os.path.exists(hoofdbestand):
        print(f"Fout: {hoofdbestand} niet gevonden!")
        return

    excel_sheets = pd.read_excel(hoofdbestand, sheet_name=None)

    if "Maatwerkgebieden" not in excel_sheets:
        print("Fout: Tabblad 'Maatwerkgebieden' niet gevonden!")
        return

    df_hoofd = excel_sheets["Maatwerkgebieden"]

    # Definieer de nieuwe kolomnamen op basis van posities (Python start bij 0)
    # Kolom 2 (index 1) = Soortnaam
    # Kolom 4 (index 3) = Model?
    # Kolom 5 (index 4) = Automatisch?
    soort_kolom = df_hoofd.columns[1]
    model_kolom = df_hoofd.columns[3]
    auto_kolom = df_hoofd.columns[4]

    print(f"-> Soorten uit hoofdtabel: '{soort_kolom}'")
    print(f"-> Modelstatus (Kolom 4): '{model_kolom}'")
    print(f"-> Automatisch status (Kolom 5): '{auto_kolom}'")

    # 3. STAP 1: 'Model?' kolom invullen (ja / nee)
    df_hoofd[model_kolom] = df_hoofd[soort_kolom].apply(
        lambda x: "ja" if clean_species_name(x) in model_soorten else "nee"
    )

    # 4. STAP 2: 'Automatisch?' kolom invullen (ja / nee)
    df_hoofd[auto_kolom] = df_hoofd[soort_kolom].apply(
        lambda x: "ja" if clean_species_name(x) in automatisch_soorten else "nee"
    )

    # 5. STAP 3: Gebiedskolommen invullen (1 / 0)
    # Door de extra kolom zijn de gebieden opgeschoven naar kolommen 6 tot 11 (index 5 tot 11)
    gebieds_kolommen = df_hoofd.columns[5:11]

    for gebied in gebieds_kolommen:
        print(f"-> Gebied verwerken: {gebied}")

        if gebied in excel_sheets:
            df_gebied = excel_sheets[gebied]
            gebied_soort_kolom = df_gebied.columns[1]  # Soortnaam in gebiedstabblad (2e kolom)

            soorten_in_gebied = set(
                df_gebied[gebied_soort_kolom].apply(clean_species_name)
            )

            df_hoofd[gebied] = df_hoofd[soort_kolom].apply(
                lambda x: 1 if clean_species_name(x) in soorten_in_gebied else 0
            )
        else:
            print(
                f"   [Waarschuwing] Geen tabblad gevonden voor '{gebied}'. Kolom krijgt overal 0."
            )
            df_hoofd[gebied] = 0

    # 6. Opslaan naar het nieuwe Excel-bestand
    print(f"-> Resultaat wegschrijven naar {output_bestand}...")
    # Denk eraan om dit bestand in Excel te sluiten voor je het script runt!
    with pd.ExcelWriter(output_bestand, engine="openpyxl") as writer:
        df_hoofd.to_excel(writer, sheet_name="Maatwerkgebieden", index=False)

        for sheet_name, df_sheet in excel_sheets.items():
            if sheet_name != "Maatwerkgebieden":
                df_sheet.to_excel(writer, sheet_name=sheet_name, index=False)

    print("-> Klaar! Alles is succesvol verwerkt en opgeschoven.")


if __name__ == "__main__":
    main()