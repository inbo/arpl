#!/bin/bash
#SBATCH --job-name=ARPL_Array
#SBATCH --cluster=wice
#SBATCH --partition=batch_icelake
#SBATCH --array=1-49
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00                # <--- AANGEPAST VAN 04:00:00 NAAR 10:00:00
#SBATCH --account=lp_inbo
#SBATCH --output=data/logs/array_%A_%a.out
#SBATCH --error=data/logs/array_%A_%a.err

module --force purge
module load cluster/wice/batch_icelake
module load R/4.5.1-gfbf-2025a
module load GDAL/3.11.1-foss-2025a
module load CMake/3.31.3-GCCcore-14.2.0
module load UDUNITS/2.2.28-GCCcore-14.2.0

mkdir -p data/logs

Rscript src/Totale_Scripts/02_run_enkele_taak.R $SLURM_ARRAY_TASK_ID