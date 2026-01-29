#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# User edit section
# -----------------------------
export PROJECT_ROOT="${PROJECT_ROOT:-$PWD}"

# Where shared resources live (PET maps, atlas, GM mask, LUT)
export MAPS_DIR="${MAPS_DIR:-/cluster/projects/<PROJECT>/resources/maps}"

# Where UKBB/FreeSurfer-exported tables live (only needed if running adapter)
export FREESURFER_STATS_DIR="${FREESURFER_STATS_DIR:-/cluster/projects/<PROJECT>/ukb/fs_stats}"

# Demographics file used by adapter (only needed if running adapter)
export DEMOGRAPHICS_FILE="${DEMOGRAPHICS_FILE:-/cluster/projects/<PROJECT>/ukb/demog.tsv}"

# Output directory for this run
export OUTDIR="${OUTDIR:-$PROJECT_ROOT/out/run_$(date +%Y%m%d_%H%M%S)}"

# Inputs consumed by normdev
export FEATURES_FILE="${FEATURES_FILE:-$OUTDIR/features/features_wide.tsv}"
export FEATURE_MAP_FILE="${FEATURE_MAP_FILE:-$OUTDIR/features/feature_map.tsv}"

# Config file (generic)
CFG="${CFG:-$PROJECT_ROOT/config_generic.yaml}"

# Toggle adapter step (1=yes, 0=no)
RUN_ADAPTER="${RUN_ADAPTER:-1}"

mkdir -p "$OUTDIR"/{logs,features}

echo "[INFO] OUTDIR: $OUTDIR"
echo "[INFO] CFG:    $CFG"
echo "[INFO] MAPS:   $MAPS_DIR"
echo "[INFO] RUN_ADAPTER: $RUN_ADAPTER"
echo "[INFO] FEATURES_FILE: $FEATURES_FILE"
echo "[INFO] FEATURE_MAP_FILE: $FEATURE_MAP_FILE"
echo

cd "$PROJECT_ROOT"
# -----------------------------
# Optional: build features_wide + feature_map
# -----------------------------
if [[ "$RUN_ADAPTER" == "1" ]]; then
  echo "[STEP] Adapter: FreeSurfer/UKBB -> features_wide.tsv + feature_map.tsv"
  Rscript "$PROJECT_ROOT/extras/prep_features.R" \
    -c "$CFG" -o "$OUTDIR" --overwrite \
    > "$OUTDIR/logs/00_adapter.log" 2>&1
  echo "[OK] Adapter done"
  echo
fi

# -----------------------------
# PET weights
# -----------------------------
echo "[STEP] PET weights extraction"
Rscript "$PROJECT_ROOT/scripts/01_extract_petvals.R" \
  -c "$CFG" -o "$OUTDIR" --overwrite \
  > "$OUTDIR/logs/01_pet.log" 2>&1
echo "[OK] PET done"
echo

# -----------------------------
# Normative deviations
# -----------------------------
echo "[STEP] Normative deviations"
Rscript "$PROJECT_ROOT/scripts/02_run_normdevs.R" \
  -c "$CFG" -o "$OUTDIR" --overwrite \
  > "$OUTDIR/logs/02_normdev.log" 2>&1
echo "[OK] Normdev done"
echo

# -----------------------------
# Compute MBI
# -----------------------------
echo "[STEP] Compute MBI"
Rscript "$PROJECT_ROOT/scripts/03_run_mbi.R" \
  -c "$CFG" -o "$OUTDIR" --overwrite \
  > "$OUTDIR/logs/03_mbi.log" 2>&1
echo "[OK] MBI done"
echo

# -----------------------------
# QC alignment
# -----------------------------
echo "[STEP] QC: relate MBI to PET weights"
Rscript "$PROJECT_ROOT/scripts/04_relate_mbi_pet.R" \
  -c "$CFG" -o "$OUTDIR" --overwrite \
  > "$OUTDIR/logs/04_qc.log" 2>&1
echo "[OK] QC done"
echo

echo "[DONE] All steps finished. Outputs in: $OUTDIR"


# ------------------------------------------------------------
# OPTIONAL: MBI calibration (EXTRAS)
# ------------------------------------------------------------
# This step performs a grid search over (gamma, beta) to
# calibrate the MBI weighting. It is NOT required for normal
# pipeline runs and is therefore not executed by default.
#
# To run manually after the pipeline completes:
#
#   Rscript extras/calibrate_mbi_grid.R \
#     --config ${CFG} \
#     --outdir ${OUTDIR}
#
# Results will be written to:
#   ${OUTDIR}/calibration/calibration_grid_results.tsv
# ------------------------------------------------------------

# ------------------------------------------------------------
# OPTIONAL (one-time): prepare PET map in fsaverage / PET space
# ------------------------------------------------------------
# This step prepares:
#   - aparc+aseg atlas in PET space
#   - GM probability map in PET space
#   - LUT CSV for region labels
#
# It is required ONCE per PET map / atlas combination and is
# NOT run automatically as part of the pipeline.
#
# To run manually (after editing paths):
#
#   bash pipeline/extras/prep_petmap.sh
#
# Outputs are written to MAPS_DIR and used by later stages.
# ------------------------------------------------------------
