#!/usr/bin/env bash
# Resample atlas + GM probability map to PET grid (3 mm) and write LUT.
# Run ONCE per atlas/PET-template combination. Outputs are consumed by
# scripts/01_extract_petvals.R.
#
# Environment overrides (set these in your shell or here at the top):
#   PROJECT_ROOT   repo root (where helpers/prep_fsaverage_to_pet.py lives)
#   PET_DIR        where raw PET .nii.gz files live
#   MAPS_DIR       where MRI templates (MNI152, GM prob map) live
#   OUT_DIR        where prepared outputs go (defaults to MAPS_DIR for
#                  backward compatibility with single-directory layouts)
#   FREESURFER_HOME path to FreeSurfer install (for fsaverage + LUT)
#   FSAVG_DIR      override fsaverage subject dir if not under FREESURFER_HOME
#   PYTHON_BIN     python executable (default: python3)
#
# Optional explicit overrides:
#   PET_FILE / MNI_FILE / GM_FILE / LUT_FILE  (inputs)
#   ATLAS_OUT / GM_OUT / LUT_OUT              (outputs)

set -euo pipefail

# ------------------------------------------------------------
# User edit section
# ------------------------------------------------------------
# Root of this repo (where helpers/prep_fsaverage_to_pet.py lives)
PROJECT_ROOT="${PROJECT_ROOT:-$PWD}" #/ess/p33/cluster/ukbio_users/dennisva/mbi-pipeline

# Where raw PET .nii.gz inputs live
PET_DIR="${PET_DIR:-/path/to/pet_raw}" #/ess/p33/cluster/ukbio_users/dennisva/mbi_ukb/neuromaps

# Where MRI template inputs (MNI152, GM probability) live
MAPS_DIR="${MAPS_DIR:-/path/to/templates}" #/ess/p33/cluster/ukbio_users/dennisva/mbi_ukb/neuromaps

# FreeSurfer fsaverage directory (set FREESURFER_HOME or set FSAVG_DIR directly)
FREESURFER_HOME="${FREESURFER_HOME:-/path/to/freesurfer}" #/cluster/software/EL9/amd/zen/easybuild/software/FreeSurfer/5.3.0-centos6_x86_64
FSAVG_DIR="${FSAVG_DIR:-${FREESURFER_HOME}/subjects/fsaverage}"

# Where prepared outputs go. Defaults to MAPS_DIR (legacy single-directory
# layout). Override if you keep inputs and outputs separate.
OUT_DIR="${OUT_DIR:-$MAPS_DIR}"

# Python executable
PY="${PYTHON_BIN:-python3}"

# (Optional) override any of these if you use different filenames
PET_FILE="${PET_FILE:-$PET_DIR/source-castrillon2023_desc-cmrglc_space-MNI152_res-3mm_feature.nii.gz}"
MNI_FILE="${MNI_FILE:-$MAPS_DIR/MNI152_T1_1mm_brain.nii.gz}"
GM_FILE="${GM_FILE:-$MAPS_DIR/avg152T1_gray.nii.gz}"
LUT_FILE="${LUT_FILE:-${FREESURFER_HOME}/FreeSurferColorLUT.txt}"

ATLAS_OUT="${ATLAS_OUT:-$OUT_DIR/aparc+aseg_inPETspace_3mm.nii.gz}"
GM_OUT="${GM_OUT:-$OUT_DIR/GMprob_inPETspace_3mm.nii.gz}"
LUT_OUT="${LUT_OUT:-$OUT_DIR/lut_aparc_aseg.csv}"
# ------------------------------------------------------------

PY_SCRIPT="$PROJECT_ROOT/helpers/prep_fsaverage_to_pet.py"

mkdir -p "$OUT_DIR"

echo "[INFO] PROJECT_ROOT: $PROJECT_ROOT"
echo "[INFO] PET_FILE:     $PET_FILE"
echo "[INFO] MNI_FILE:     $MNI_FILE"
echo "[INFO] GM_FILE:      $GM_FILE"
echo "[INFO] LUT_FILE:     $LUT_FILE"
echo "[INFO] FSAVG_DIR:    $FSAVG_DIR"
echo "[INFO] ATLAS_OUT:    $ATLAS_OUT"
echo "[INFO] GM_OUT:       $GM_OUT"
echo "[INFO] LUT_OUT:      $LUT_OUT"
echo

[[ -f "$PET_FILE"  ]] || { echo "[ERROR] PET_FILE not found"  >&2; exit 1; }
[[ -f "$MNI_FILE"  ]] || { echo "[ERROR] MNI_FILE not found"  >&2; exit 1; }
[[ -f "$GM_FILE"   ]] || { echo "[ERROR] GM_FILE not found"   >&2; exit 1; }
[[ -f "$LUT_FILE"  ]] || { echo "[ERROR] LUT_FILE not found"  >&2; exit 1; }
[[ -d "$FSAVG_DIR" ]] || { echo "[ERROR] FSAVG_DIR not found" >&2; exit 1; }
[[ -f "$PY_SCRIPT" ]] || { echo "[ERROR] PY_SCRIPT not found: $PY_SCRIPT" >&2; exit 1; }

"$PY" "$PY_SCRIPT" \
  --pet       "$PET_FILE" \
  --mni       "$MNI_FILE" \
  --gm        "$GM_FILE"  \
  --lut       "$LUT_FILE" \
  --fsavg     "$FSAVG_DIR" \
  --atlas-out "$ATLAS_OUT" \
  --gm-out    "$GM_OUT"   \
  --lut-out   "$LUT_OUT"

echo "[OK] Done"
