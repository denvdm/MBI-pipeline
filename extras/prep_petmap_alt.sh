#!/usr/bin/env bash
# Resample atlas + GM to alt PET grid. Run ONCE per alt template.
# Thin wrapper around extras/prep_petmap.sh; sets per-template output names.
#
# Required env vars:
#   PROJECT_ROOT      MBI-pipeline repo root
#   MAPS_DIR          where PET resources live (atlas/GM/LUT outputs land here)
#   FREESURFER_HOME   FreeSurfer install path
#   TEMPLATE_LABEL    short identifier, e.g. "vaishnavi2010"
#   ALT_PET_FILE      absolute path to alt PET .nii.gz
#
# Outputs (written under MAPS_DIR):
#   aparc+aseg_in${TEMPLATE_LABEL}PETspace.nii.gz
#   GMprob_in${TEMPLATE_LABEL}PETspace.nii.gz
#   lut_aparc_aseg_${TEMPLATE_LABEL}.csv
#
# Example:
#   export TEMPLATE_LABEL="vaishnavi2010"
#   export ALT_PET_FILE="$MAPS_DIR/source-vaishnavi2010_desc-cmrglc_space-MNI152_feature.nii.gz"
#   bash extras/prep_petmap_alt.sh
set -euo pipefail

: "${PROJECT_ROOT:?Need PROJECT_ROOT}"
: "${MAPS_DIR:?Need MAPS_DIR}"
: "${FREESURFER_HOME:?Need FREESURFER_HOME}"
: "${TEMPLATE_LABEL:?Need TEMPLATE_LABEL}"
: "${ALT_PET_FILE:?Need ALT_PET_FILE}"

[[ -f "$ALT_PET_FILE" ]] || { echo "[ERROR] ALT_PET_FILE not found: $ALT_PET_FILE" >&2; exit 1; }

# Override the env vars consumed by prep_petmap.sh
export PET_FILE="$ALT_PET_FILE"
export ATLAS_OUT="$MAPS_DIR/aparc+aseg_in${TEMPLATE_LABEL}PETspace.nii.gz"
export GM_OUT="$MAPS_DIR/GMprob_in${TEMPLATE_LABEL}PETspace.nii.gz"
export LUT_OUT="$MAPS_DIR/lut_aparc_aseg_${TEMPLATE_LABEL}.csv"
# MNI_FILE, GM_FILE, FSAVG_DIR, LUT_FILE inherit defaults from prep_petmap.sh

echo "[INFO] TEMPLATE_LABEL: $TEMPLATE_LABEL"
echo "[INFO] PET source:     $PET_FILE"
echo "[INFO] Atlas output:   $ATLAS_OUT"
echo "[INFO] GM output:      $GM_OUT"
echo "[INFO] LUT output:     $LUT_OUT"

bash "$PROJECT_ROOT/extras/prep_petmap.sh"

echo "[OK] Alt-template prep done for: $TEMPLATE_LABEL"
