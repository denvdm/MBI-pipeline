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
# Discovery/replication split
# -----------------------------
# Writes the assignment when split.mode is set and exits quietly when it is not,
# so this step is safe to leave in the driver unconditionally. It does not gate
# what stage 02 fits on: stage 02 recomputes the same assignment from the same
# two config values, so the file cannot fall out of step with the run. The file
# is for the phenotypic and genetic analyses downstream of this pipeline.
echo "[STEP] Discovery/replication split"
Rscript "$PROJECT_ROOT/extras/make_split.R" \
  -c "$CFG" -o "$OUTDIR" --overwrite \
  > "$OUTDIR/logs/01b_split.log" 2>&1
echo "[OK] Split step done"
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
# OPTIONAL: parameter sensitivity strip (EXTRAS)
# ------------------------------------------------------------
# Reports how the index behaves across emphasis (the only
# parameter that changes MBI) and beta. This is a sensitivity
# report, NOT a selection step: fix the parameters first and
# record them in extras/SELECTION_RULE.md, then run this to show
# the result does not hinge on them. Not executed by default.
#
# To run manually after the pipeline completes:
#
#   Rscript extras/sweep_emphasis.R \
#     --config ${CFG} \
#     --outdir ${OUTDIR}
#
# Results will be written to:
#   ${OUTDIR}/calibration/sweep_emphasis_results.tsv
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
#   bash extras/prep_petmap.sh
#
# Outputs are written to MAPS_DIR and used by later stages.
# ------------------------------------------------------------
