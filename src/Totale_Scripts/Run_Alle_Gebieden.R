library(future)
library(furrr)

# 1. Definieer de hoofdmap waar de gebieden in staan
basis_pad <- "/src"

# 2. Zoek specifiek naar alle .R bestanden die binnen een Automatisering_Scripts map liggen
script_bestanden <- list.files(
  path = basis_pad, 
  pattern = "\\.R$", 
  full.names = TRUE, 
  recursive = TRUE
)

# Filteren zodat enkel de scripts uit de Automatisering_Scripts mappen meegenomen worden
script_bestanden <- script_bestanden[grep("/Automatisering_Scripts/", script_bestanden)]

# Controle vooraf: hoeveel scripts zijn er gevonden?
cat("Aantal gevonden automatisering scripts:", length(script_bestanden), "\n")

# 3. Stel exact 6 cores in voor de test
plan(multisession, workers = 6)

# 4. Functie om elk script afzonderlijk en veilig uit te voeren
run_script_veilig <- function(script_pad) {
  tijd_start <- Sys.time()
  
  # Haal het gebied uit het pad voor het overzicht
  gebied <- unlist(strsplit(script_pad, "/"))[3] 
  
  tryCatch({
    # Run in een schone, geïsoleerde omgeving om conflicten tussen variabelen te vermijden
    source(script_pad, local = new.env()) 
    
    tijd_eind <- Sys.time()
    duurtijd <- round(as.numeric(tijd_eind - tijd_start, units = "mins"), 2)
    
    return(data.frame(
      gebied = gebied, 
      script = basename(script_pad), 
      status = "Succes", 
      duurtijd_min = duurtijd
    ))
  }, error = function(e) {
    return(data.frame(
      gebied = gebied, 
      script = basename(script_pad), 
      status = paste("Fout:", e$message), 
      duurtijd_min = NA
    ))
  })
}

# 5. Voer de test uit op 6 cores en toon voortgang
cat("Starten van de parallelle test op 6 cores...\n")
start_totaal <- Sys.time()

# furrr regelt de verdeling van de scripts over de 6 cores
resultaten_df <- future_map_dfr(script_bestanden, run_script_veilig, .progress = TRUE)

eind_totaal <- Sys.time()
totale_tijd <- round(as.numeric(eind_totaal - start_totaal, units = "hours"), 2)

cat("\n--- TEST AFGEROND ---\n")
cat("Totale duurtijd op 6 cores:", totale_tijd, "uur\n")

# Bekijk samenvatting van fouten/succes per gebied
print(table(resultaten_df$gebied, resultaten_df$status))
