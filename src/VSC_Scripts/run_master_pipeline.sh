#!/bin/bash
#SBATCH --job-name=ARPL_Master_Pipeline
#SBATCH --cluster=wice
#SBATCH --partition=batch_icelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --account=lp_inbo
#SBATCH --output=data/logs/master_pipeline_%j.out
#SBATCH --error=data/logs/master_pipeline_%j.err

module --force purge
module load cluster/wice/batch_icelake
module load R/4.5.1-gfbf-2025a
module load GDAL/3.11.1-foss-2025a
module load CMake/3.31.3-GCCcore-14.2.0
module load UDUNITS/2.2.28-GCCcore-14.2.0

mkdir -p data/logs data/temp

echo "=================================================="
echo " STAP 1: GENEREREN VAN DYNAMISCHE TAKENLIJST"
echo "=================================================="
Rscript src/VSC_Scripts/00_Verzamel_Taken.R

TAKEN_AANTAL=$(Rscript -e "cat(length(readRDS('data/temp/globale_takenlijst.rds')))")

echo "Totaal aantal te verwerken taken: $TAKEN_AANTAL"

if [ "$TAKEN_AANTAL" -gt 0 ]; then
    echo "=================================================="
    echo " STAP 2: LAUNCHEN SLURM ARRAY ($TAKEN_AANTAL TAKEN PARALLEL)"
    echo "=================================================="
    
    ARRAY_JOB_ID=$(sbatch --parsable --cluster=wice --partition=batch_icelake \
      --job-name=ARPL_Array --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=128G \
      --time=08:00:00 --account=lp_inbo --array=1-$TAKEN_AANTAL \
      --output=data/logs/array_%A_%a.out --error=data/logs/array_%A_%a.err \
      --wrap="module load cluster/wice/batch_icelake R/4.5.1-gfbf-2025a GDAL/3.11.1-foss-2025a CMake/3.31.3-GCCcore-14.2.0 UDUNITS/2.2.28-GCCcore-14.2.0 && Rscript src/VSC_Scripts/01_Run_Enkele_Taak.R \$SLURM_ARRAY_TASK_ID")
      
    echo "Array Job ingediend met ID: $ARRAY_JOB_ID"
    echo "Wachten tot alle parallelle taken in $ARRAY_JOB_ID zijn afgerond..."
    
    sbatch --wait --cluster=wice --partition=batch_icelake \
      --job-name=ARPL_Wait --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=4G \
      --time=01:00:00 --account=lp_inbo --dependency=afterok:$ARRAY_JOB_ID \
      --output=data/logs/wait_%j.out --wrap="echo 'Alle parallelle taken zijn afgerond!'"
else
    echo "Alle rasters bestaan reeds op schijf. Slaan Stap 2 over."
fi

echo "=================================================="
echo " STAP 3: FINALE ARPL ACTUEEL MAKER RUNNEN"
echo "=================================================="
Rscript src/Totale_Scripts/ARPL_Actueel_Maker.R

echo "=================================================="
echo " 🎉 PIPELINE VOLLEDIG AFGEROND!"
echo "=================================================="