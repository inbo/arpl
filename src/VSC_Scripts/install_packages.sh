#!/bin/bash
#SBATCH --job-name=ARPL_Install_Final
#SBATCH --cluster=wice
#SBATCH --partition=batch_icelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --account=lp_inbo
#SBATCH --output=data/logs/install_final_%j.out
#SBATCH --error=data/logs/install_final_%j.err

module --force purge
module load cluster/wice/batch_icelake
module load R/4.5.1-gfbf-2025a
module load GDAL/3.11.1-foss-2025a
module load CMake/3.31.3-GCCcore-14.2.0
module load UDUNITS/2.2.28-GCCcore-14.2.0

mkdir -p data/logs

Rscript -e '
lib_dir <- "/vsc-hard-mounts/leuven-data/392/vsc39293/Rlibs/rocky9/icelake/R-4.5.1"
dir.create(lib_dir, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(lib_dir, .libPaths()))

# Ruim achtergebleven lock-bestanden op
locks <- list.files(lib_dir, pattern = "^00LOCK", full.names = TRUE)
if (length(locks) > 0) unlink(locks, recursive = TRUE)

pkgs <- c("units", "vroom", "readr", "terra", "sf", "exactextractr", "leaflet", "tidyterra", "tidyverse")

cat("=== START ALLERLAATSTE INSTALLATIE-STAP ===\n")
install.packages(pkgs, lib = lib_dir, repos = "https://cloud.r-project.org", Ncpus = 8)
cat("=== COMPILATIE AFGEROND ===\n")
'
