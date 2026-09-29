#!/bin/bash
#SBATCH --job-name=ARPL_Maker
#SBATCH --cluster=wice
#SBATCH --partition=batch_icelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --account=lp_inbo
#SBATCH --output=data/logs/maker_%j.out
#SBATCH --error=data/logs/maker_%j.err

module --force purge
module load cluster/wice/batch_icelake
module load R/4.5.1-gfbf-2025a
module load GDAL/3.11.1-foss-2025a
module load CMake/3.31.3-GCCcore-14.2.0
module load UDUNITS/2.2.28-GCCcore-14.2.0

mkdir -p data/logs

Rscript src/Totale_Scripts/ARPL_Actueel_Maker.R
