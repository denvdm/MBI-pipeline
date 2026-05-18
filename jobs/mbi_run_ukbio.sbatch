#!/bin/bash

#SBATCH --job-name=mbi_ukb
#SBATCH --account=<YOUR_ACCOUNT>
#SBATCH --mem-per-cpu=32G
#SBATCH --time=4:00:00
#SBATCH --output=logs/mbi_%j.out
#SBATCH --error=logs/mbi_%j.err

set -euo pipefail

mkdir -p logs

# Load R / modules as needed on cluster
module purge
module load R/4.2.0-foss-2021b

export PROJECT_ROOT="/path/to/mbi-pipeline"
export MAPS_DIR="/path/to/maps"
export FREESURFER_STATS_DIR="/path/to/freesurfer_stats"
export DEMOGRAPHICS_FILE="/path/to/demographics.tsv"

export OUTDIR="/path/to/runs/run_$(date +%Y%m%d_%H%M%S)"
export FEATURES_FILE="$OUTDIR/features/features_wide.tsv"
export FEATURE_MAP_FILE="$OUTDIR/features/feature_map.tsv"
export CFG="$PROJECT_ROOT/config_generic.yaml"
export TRAIN_IDS="/path/to/train_ids.txt"

export RUN_ADAPTER=1

cd $PROJECT_ROOT
bash "$PROJECT_ROOT/scripts/run_all.sh"


#ln -sfn "$(ls -1d run_* | sort | tail -n1)" current
