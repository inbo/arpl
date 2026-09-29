#!/bin/bash
#SBATCH --job-name=ARPL_Check
#SBATCH --cluster=wice
#SBATCH --partition=batch_icelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:05:00
#SBATCH --account=lp_inbo
#SBATCH --output=data/logs/check_gis_%j.out
#SBATCH --error=data/logs/check_gis_%j.err

module --force purge
module load cluster/wice/batch_icelake
module load R/4.5.1-gfbf-2025a
module load GDAL/3.11.1-foss-2025a
module load CMake/3.31.3-GCCcore-14.2.0
module load UDUNITS/2.2.28-GCCcore-14.2.0

Rscript -e '
lib_dir <- "/vsc-hard-mounts/leuven-data/392/vsc39293/Rlibs/rocky9/icelake/R-4.5.1"
.libPaths(c(lib_dir, .libPaths()))

pkgs <- c("tidyverse", "sf", "terra", "tidyterra", "leaflet", "readr", "exactextractr")

cat("=== CONTROLE OP BATCH ICELAKE NODE MET GIS MODULES ===\n")
for (p in pkgs) {
  res <- require(p, character.only = TRUE, quietly = TRUE)
  cat(sprintf("  %-15s : %s\n", p, ifelse(res, "[OK] GELADEN", "[FAIL] ONTBREKT")))
}
'
