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

soorten_tv <- df_soorten %>%
  filter(Turnhouts_Vennegebied == 1) %>%
  pull('Nederlandse naam') %>%
  tolower() %>%
  trimws() %>%
  unique()

# Als 'wulp' in de lijst staat, splitsen we deze op in broedvogel en wintervogel
if ("wulp" %in% soorten_tv) {
  soorten_tv <- soorten_tv[soorten_tv != "wulp"] # Verwijder generieke 'wulp'
  soorten_tv <- c(soorten_tv, "wulp (broedvogel)", "wulp (wintervogel)")
}

cat("Aantal te verwerken soorten voor Turnhouts Vennegebied:", length(soorten_tv), "\n\n")

# 3. Lus uitvoeren over alle soorten
path_rmd <- here("src/Turnhouts_Vennegebied/Automatisering_Scripts/Visualisatie_Habitatkaart.Rmd")

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
