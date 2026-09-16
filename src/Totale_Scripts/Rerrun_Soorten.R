library(tidyverse)
library(rmarkdown)
library(here)

knitr::opts_knit$set(root.dir = here())

# 1. Definieer de lijst met uit te voeren combinaties
missing_tasks <- tribble(
  ~Gebied, ~Soort,
  "De_Maten", "blauwborst",
  "De_Maten", "bruinekiekendief",
  "De_Maten", "grotemodderkruiper",
  "De_Maten", "ijsvogel",
  "De_Maten", "kwak",
  "De_Maten", "matkop",
  "De_Maten", "nachtegaal",
  "De_Maten", "roerdomp",
  "De_Maten", "tapuit",
  "De_Maten", "watersnip",
  "De_Maten", "wielewaal",
  "De_Maten", "woudaap",
  "De_Maten", "zwartespecht",
  "Heesbossen", "grutto",
  "Heesbossen", "heivlinder",
  "Heesbossen", "kwartelkoning",
  "Heesbossen", "matkop",
  "Heesbossen", "nachtegaal",
  "Heesbossen", "paapje",
  "Heesbossen", "waterdrieblad",
  "Heesbossen", "watersnip",
  "Heesbossen", "wielewaal",
  "Heesbossen", "zomertortel",
  "Kalmthoutse_Heide", "adder",
  "Kalmthoutse_Heide", "blauwborst",
  "Kalmthoutse_Heide", "boomleeuwerik",
  "Kalmthoutse_Heide", "bruinekiekendief",
  "Kalmthoutse_Heide", "europeseotter",
  "Kalmthoutse_Heide", "fluiter",
  "Kalmthoutse_Heide", "grutto",
  "Kalmthoutse_Heide", "kamsalamander",
  "Kalmthoutse_Heide", "nachtzwaluw",
  "Kalmthoutse_Heide", "tapuit",
  "Kalmthoutse_Heide", "veldkrekel",
  "Kalmthoutse_Heide", "watersnip",
  "Kalmthoutse_Heide", "wespendief",
  "Kalmthoutse_Heide", "wulp",
  "Kalmthoutse_Heide", "zwartespecht",
  "Mechelse_Heide", "blauwborst",
  "Mechelse_Heide", "boomleeuwerik",
  "Mechelse_Heide", "bruineeikenpage",
  "Mechelse_Heide", "europeseotter",
  "Mechelse_Heide", "fluiter",
  "Mechelse_Heide", "grauweklauwier",
  "Mechelse_Heide", "kleinblaasjeskruid",
  "Mechelse_Heide", "knoflookpad",
  "Mechelse_Heide", "kwartelkoning",
  "Mechelse_Heide", "matkop",
  "Mechelse_Heide", "nachtegaal",
  "Mechelse_Heide", "nachtzwaluw",
  "Mechelse_Heide", "tapuit",
  "Mechelse_Heide", "vliegendhert",
  "Mechelse_Heide", "watersnip",
  "Mechelse_Heide", "wespendief",
  "Mechelse_Heide", "wielewaal",
  "Mechelse_Heide", "zwartespecht",
  "Turnhouts_Vennegebied", "bruinekiekendief",
  "Turnhouts_Vennegebied", "fluiter",
  "Turnhouts_Vennegebied", "grotepimpernel",
  "Turnhouts_Vennegebied", "heikikker",
  "Turnhouts_Vennegebied", "heivlinder",
  "Turnhouts_Vennegebied", "hoogveenglanslibel",
  "Turnhouts_Vennegebied", "kwartelkoning",
  "Turnhouts_Vennegebied", "matkop",
  "Turnhouts_Vennegebied", "paapje",
  "Turnhouts_Vennegebied", "porseleinhoen",
  "Turnhouts_Vennegebied", "spaanseruiter",
  "Turnhouts_Vennegebied", "tapuit",
  "Turnhouts_Vennegebied", "venglazenmaker",
  "Turnhouts_Vennegebied", "watersnip",
  "Turnhouts_Vennegebied", "wielewaal",
  "Turnhouts_Vennegebied", "zomertortel",
  "Turnhouts_Vennegebied", "zwartkopmeeuw",
  "Voerstreek", "geleanemoon",
  "Voerstreek", "kleinwarkruid",
  "Voerstreek", "kwartelkoning",
  "Voerstreek", "nachtegaal",
  "Voerstreek", "spaansevlag",
  "Voerstreek", "wielewaal",
  "Voerstreek", "zomertortel"
)

# Inventariseer alle Rmd-bestanden
all_rmd <- list.files(path = here(), pattern = "\\.[R|r]md$", recursive = TRUE, full.names = TRUE)

results_log <- data.frame(Gebied = character(), Soort = character(), Status = character(), File = character(), stringsAsFactors = FALSE)

# Tijdstip van gisterenochtend instellen (30 uur geleden)
start_cutoff <- Sys.time() - lubridate::hours(30)

# 2. Door de taken lussen
for (i in 1:nrow(missing_tasks)) {
  geb <- missing_tasks$Gebied[i]
  srt <- missing_tasks$Soort[i]
  
  message(sprintf("[%d/%d] Controle: %s - %s", i, nrow(missing_tasks), geb, srt))
  
  # EXACTE MATCHING (Voorkom wulp vs regenwulp)
  regex_geb <- paste0("(?i)", geb)
  regex_srt <- paste0("(?i)(^|[^a-z])", srt, "($|[^a-z])") 
  
  matches <- all_rmd[grepl(regex_geb, all_rmd, perl = TRUE) & grepl(regex_srt, all_rmd, perl = TRUE)]
  
  # Failsafe voor meerdere matches (bijv. exact wulp.Rmd kiezen boven wulp_wv.Rmd)
  if (length(matches) > 1) {
    exact_file_match <- matches[tolower(basename(matches)) == tolower(paste0(srt, ".rmd")) | 
                                  tolower(basename(matches)) == tolower(paste0("scenario_", srt, ".rmd"))]
    if (length(exact_file_match) > 0) {
      matches <- exact_file_match
    }
  }
  
  is_simpel_fallback <- FALSE
  
  # FALLBACK LOGICA: Geen specifiek script? Zoek Leefgebieden_Simpel
  if (length(matches) == 0) {
    simpel_matches <- all_rmd[grepl(regex_geb, all_rmd, perl = TRUE) & grepl("(?i)Leefgebieden_Simpel", all_rmd, perl = TRUE)]
    if (length(simpel_matches) > 0) {
      matches <- simpel_matches[1]
      is_simpel_fallback <- TRUE
      message(sprintf("  [INFO] Geen specifiek script voor '%s'. Fallback naar simpel script: %s", srt, basename(matches[1])))
    }
  }
  
  if (length(matches) >= 1) {
    rmd_file <- matches[1]
    
    # Bepaal het verwachte HTML-uitvoerbestand
    output_file_name <- if (is_simpel_fallback) paste0(geb, "_Leefgebied_", srt, ".html") else NULL
    expected_html <- if (is_simpel_fallback) file.path(dirname(rmd_file), output_file_name) else sub("\\.[R|r]md$", ".html", rmd_file)
    
    # CHECK: Sla over als de HTML bestaat EN is bijgewerkt na gisterenochtend
    if (file.exists(expected_html)) {
      mtime <- file.info(expected_html)$mtime
      if (!is.na(mtime) && mtime > start_cutoff) {
        message(sprintf("  [SKIP] Reeds succesvol gerund sinds gisterenochtend (%s): %s (%s)", format(mtime, "%H:%M"), basename(rmd_file), srt))
        results_log <- rbind(results_log, data.frame(Gebied = geb, Soort = srt, Status = "REEDS_RERENDERD", File = basename(rmd_file)))
        next
      }
    }
    
    # RENDEREN DYNAMISCH AFHANDELEN
    tryCatch({
      render_env <- new.env(parent = globalenv())
      render_env$soort <- srt
      render_env$soort_invoer <- srt
      render_env$geb <- geb
      
      # Gebruik 'soort_invoer' voor params bij Simpel-scripts
      render_params <- if (is_simpel_fallback) list(soort_invoer = srt) else NULL
      
      rmarkdown::render(
        input = rmd_file,
        output_file = output_file_name,
        params = render_params,
        envir = render_env,
        quiet = TRUE
      )
      
      message(sprintf("  [OK] Succesvol afgerond: %s (%s)", basename(rmd_file), srt))
      results_log <- rbind(results_log, data.frame(Gebied = geb, Soort = srt, Status = if (is_simpel_fallback) "OK (Simpel)" else "OK", File = basename(rmd_file)))
    }, error = function(e) {
      message(sprintf("  [FOUT] Fout bij verwerken van %s (%s): %s", basename(rmd_file), srt, e$message))
      results_log <- rbind(results_log, data.frame(Gebied = geb, Soort = srt, Status = paste("ERROR:", e$message), File = basename(rmd_file)))
    })
    
  } else {
    message(sprintf("  [WAARSCHUWING] Noch specifiek noch Simpel Rmd-bestand gevonden voor %s - %s", geb, srt))
    results_log <- rbind(results_log, data.frame(Gebied = geb, Soort = srt, Status = "NIET_GEVONDEN", File = "NA"))
  }
  
  gc()
}

print(results_log)
