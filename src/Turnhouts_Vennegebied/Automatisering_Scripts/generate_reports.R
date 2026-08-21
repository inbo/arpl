library(here)
library(readxl)
library(dplyr)
library(rmarkdown)
library(stringr)

# 1. Mappenstructuur instellen en aanmaken
output_dir <- here("data/output/Turnhouts_Vennegebied/HTML")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  message("-> Map aangemaakt: ", output_dir)
}

# 2. Soortenlijst inlezen en filteren op Turnhouts Vennegebied
path_soortenlijst <- here("data/input/Excel_files/Soortenlijst_Maatwerkgebieden.xlsx")
df_soorten <- read_excel(path_soortenlijst)

# Filteren op soorten waar 'Turnhouts_Vennegebied' gelijk is aan 1
soorten_tv <- df_soorten %>%
  filter(Turnhouts_Vennegebied == 1) %>%
  pull(Soort) %>%
  tolower() %>%
  trimws() %>%
  unique()

cat("Aantal te verwerken soorten voor Turnhouts Vennegebied:", length(soorten_tv), "\n\n")

# 3. Lus uitvoeren over alle soorten
path_rmd <- here("Visualisatie_Habitatkaart.Rmd")

for (s in soorten_tv) {
  message(paste0("🚀 Genereren HTML voor: ", str_to_title(s), "..."))
  
  # Bestandsnaam opschonen voor de HTML-output
  bestand_naam <- paste0("Habitatkaart_", gsub(" ", "_", s), ".html")
  output_file_path <- file.path(output_dir, bestand_naam)
  
  tryCatch({
    rmarkdown::render(
      input = path_rmd,
      output_file = output_file_path,
      params = list(soort = s),
      quiet = TRUE
    )
    message(paste0("   [OK] Opgeslagen: ", bestand_naam))
  }, error = function(e) {
    message(paste0("   ⚠️ FOUT bij soort '", s, "': ", e$message))
  })
}

message("\n🏁 Alle rapporten succesvol gegenereerd in: ", output_dir)
