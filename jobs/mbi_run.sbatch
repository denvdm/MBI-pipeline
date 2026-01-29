#!/bin/bash
#SBATCH --job-name=mbi
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=logs/mbi_%j.out
#SBATCH --error=logs/mbi_%j.err

set -euo pipefail

mkdir -p logs

# Load R / modules as needed on cluster
module purge
module load R/4.2.0-foss-2021b

export PROJECT_ROOT="$PWD"
export MAPS_DIR="/cluster/projects/<PROJECT>/resources/maps"
export FREESURFER_STATS_DIR="/cluster/projects/<PROJECT>/ukb/fs_stats"
export DEMOGRAPHICS_FILE="/cluster/projects/<PROJECT>/ukb/demog.tsv"

export OUTDIR="/scratch/$USER/mbi/run_${SLURM_JOB_ID}"
export FEATURES_FILE="$OUTDIR/features/features_wide.tsv"
export FEATURE_MAP_FILE="$OUTDIR/features/feature_map.tsv"
export CFG="$PROJECT_ROOT/config_generic.yaml"

export RUN_ADAPTER=1

bash "$PROJECT_ROOT/scripts/run_all.sh"
