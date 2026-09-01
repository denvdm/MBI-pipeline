#!/bin/bash
#SBATCH --job-name=mbi_sens
#SBATCH --account=<YOUR_ACCOUNT>
#SBATCH --mem-per-cpu=32G
#SBATCH --time=2:00:00
#SBATCH --output=logs/mbi_sens_%j.out
#SBATCH --error=logs/mbi_sens_%j.err

# NOTE: --output/--error above are opened by SLURM before this script runs,
# relative to the submission directory. Run `mkdir -p logs` before sbatch.
set -euo pipefail
mkdir -p logs

module purge
module load R/4.2.0-foss-2021b

# -------------------------------------------------------------------
# Edit here: pick the alt template
# -------------------------------------------------------------------
export TEMPLATE_LABEL="${TEMPLATE_LABEL:-dellarosa2014}"

# Auto-discover the most recent canonical run unless caller overrides
RUNS_ROOT="/path/to/runs"
if [[ -z "${CANONICAL_OUTDIR:-}" ]]; then
  CANONICAL_OUTDIR="$(ls -dt "$RUNS_ROOT"/run_*/ 2>/dev/null | head -1)"
fi
[[ -n "${CANONICAL_OUTDIR:-}" ]] || { echo "[ERROR] No canonical run found under $RUNS_ROOT" >&2; exit 1; }
export CANONICAL_OUTDIR="${CANONICAL_OUTDIR%/}"

export PROJECT_ROOT="/path/to/mbi-pipeline"
export MAPS_DIR="/path/to/maps"

# Config the sensitivity run inherits its index parameters from. Leave unset to
# take the config the canonical run recorded in compute_mbi/run_info.txt.
# export CANONICAL_CFG="$PROJECT_ROOT/config_generic.yaml"

# Per-template PET file. Override ALT_PET_FILE explicitly for any template not
# listed below, or to use a different filename convention.
case "$TEMPLATE_LABEL" in
  dellarosa2014)
    export ALT_PET_FILE="${ALT_PET_FILE:-$MAPS_DIR/source-dellarosa2014_desc-fdg_space-MNI152_feature.nii.gz}"
    ;;
  cermep)
    export ALT_PET_FILE="${ALT_PET_FILE:-$MAPS_DIR/source-cermep_desc-fdg_space-MNI152_feature.nii.gz}"
    ;;
  *)
    : "${ALT_PET_FILE:?Set ALT_PET_FILE for template '$TEMPLATE_LABEL'}"
    ;;
esac

# Atlas + GM resampled to alt PET grid (made once per template by prep_petmap_alt.sh)
export ALT_ATLAS_FILE="$MAPS_DIR/aparc+aseg_in${TEMPLATE_LABEL}PETspace.nii.gz"
export ALT_GM_FILE="$MAPS_DIR/GMprob_in${TEMPLATE_LABEL}PETspace.nii.gz"
export ALT_LUT_FILE="$MAPS_DIR/lut_aparc_aseg_${TEMPLATE_LABEL}.csv"
# -------------------------------------------------------------------

echo "[INFO] TEMPLATE_LABEL:   $TEMPLATE_LABEL"
echo "[INFO] CANONICAL_OUTDIR: $CANONICAL_OUTDIR"
echo "[INFO] ALT_PET_FILE:     $ALT_PET_FILE"

cd "$PROJECT_ROOT"
bash "$PROJECT_ROOT/scripts/run_sens_all.sh"
