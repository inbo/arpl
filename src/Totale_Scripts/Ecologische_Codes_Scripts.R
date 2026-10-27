library(readr)
library(dplyr)
library(here)

# 1. Lees het CSV-bestand in
# (Gebruik read_csv2 als de waarden gescheiden zijn door puntkomma's ';')
df <- read_csv(here("data/input/Excel_files/Resultaten_Totaal_Samengevoegd.csv"), show_col_types = FALSE)

# 2. Bepaal de uitvoermap en maak deze aan indien ze nog niet bestaat
output_dir <- here("data/references/Ecologisch_Verhaal_Soorten")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 3. Maak de opgemaakte code: voeg '%' toe als Match gelijk is aan 'bevat'
df <- df %>%
  mutate(
    FormattedCode = if_else(tolower(trimws(Match)) == "bevat", 
                            paste0(trimws(Code), "%"), 
                            trimws(Code))
  )

# 4. Groepeer per soort en genereer per soort een .txt-bestand
soorten <- unique(df$Soort)

for (s in soorten) {
  # Filter op de huidige soort
  df_soort <- df %>% filter(Soort == s)
  
  # Groepeer per Type en combineer de codes gescheiden door een komma
  lijnen <- df_soort %>%
    group_by(Type) %>%
    summarise(
      Regel = paste0(first(Type), ": ", paste(FormattedCode, collapse = ", ")),
      .groups = "drop"
    ) %>%
    pull(Regel)
  
  # Zorg dat de eerste letter altijd een hoofdletter is (bv. adder -> Adder)
  s_clean <- trimws(s)
  s_gekapitaliseerd <- paste0(toupper(substr(s_clean, 1, 1)), substr(s_clean, 2, nchar(s_clean)))
  
  # Bestandsnaam en volledig pad samenstellen
  bestandsnaam <- paste0(s_gekapitaliseerd, "_2026.txt")
  volledig_pad <- file.path(output_dir, bestandsnaam)
  
  # Schrijf weg naar bestand
  writeLines(lijnen, con = volledig_pad)
  cat("Bestand aangemaakt:", bestandsnaam, "\n")
}
