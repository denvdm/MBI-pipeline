#!/usr/bin/env bash
set -euo pipefail

# module load Python/3.10.4-GCCcore-11.3.0
# python -V
# python -m venv .venv_prep_pet
# source .venv_prep_pet/bin/activate
# python -m pip install -U pip setuptools wheel
# python -m pip install numpy nibabel scipy typing_extensions nilearn pandas antspyx

# ------------------------------------------------------------
# USER: EDIT HERE (keep it simple)
# ------------------------------------------------------------
# Root of this repo (where helpers/prep_fsaverage_to_pet.py lives)
PROJECT_ROOT="${PROJECT_ROOT:-$PWD}" #/ess/p33/cluster/ukbio_users/dennisva/mbi-pipeline

# Where your PET/map resources live
MAPS_DIR="${MAPS_DIR:-/path/to/neuromaps}" #/ess/p33/cluster/ukbio_users/dennisva/mbi_ukb/neuromaps

# FreeSurfer fsaverage directory (set FREESURFER_HOME or set FSAVG_DIR directly)
FREESURFER_HOME="${FREESURFER_HOME:-/path/to/freesurfer}" #/cluster/software/EL9/amd/zen/easybuild/software/FreeSurfer/5.3.0-centos6_x86_64
FSAVG_DIR="${FSAVG_DIR:-${FREESURFER_HOME}/subjects/fsaverage}" 

# Python executable
PY="${PYTHON_BIN:-python3}"

# (Optional) override any of these if you use different filenames
PET_FILE="${PET_FILE:-$MAPS_DIR/source-castrillon2023_desc-cmrglc_space-MNI152_res-3mm_feature.nii.gz}"
MNI_FILE="${MNI_FILE:-$MAPS_DIR/MNI152_T1_1mm_brain.nii.gz}"
GM_FILE="${GM_FILE:-$MAPS_DIR/avg152T1_gray.nii.gz}"
LUT_FILE="${LUT_FILE:-${FREESURFER_HOME}/FreeSurferColorLUT.txt}"

ATLAS_OUT="${ATLAS_OUT:-$MAPS_DIR/aparc+aseg_inPETspace_3mm.nii.gz}"
GM_OUT="${GM_OUT:-$MAPS_DIR/GMprob_inPETspace_3mm.nii.gz}"
LUT_OUT="${LUT_OUT:-$MAPS_DIR/lut_aparc_aseg.csv}"
# ------------------------------------------------------------

PY_SCRIPT="$PROJECT_ROOT/helpers/prep_fsaverage_to_pet.py"

# Small sanity checks
for f in "$PY_SCRIPT" "$PET_FILE" "$MNI_FILE" "$GM_FILE"; do
  [[ -f "$f" ]] || { echo "[ERROR] Missing file: $f" >&2; exit 1; }
done
[[ -d "$FSAVG_DIR" ]] || { echo "[ERROR] Missing FSAVG_DIR: $FSAVG_DIR" >&2; exit 1; }
[[ -f "$LUT_FILE" ]]  || { echo "[ERROR] Missing LUT_FILE: $LUT_FILE" >&2; exit 1; }

# main command:
"$PY" "$PY_SCRIPT" \
  --pet "$PET_FILE" \
  --mni "$MNI_FILE" \
  --fsavg-dir "$FSAVG_DIR" \
  --atlas-out "$ATLAS_OUT" \
  --gm "$GM_FILE" \
  --gm-out "$GM_OUT" \
  --lut "$LUT_FILE" \
  --lut-out "$LUT_OUT"

echo "[OK] Wrote:"
echo "  - $ATLAS_OUT"
echo "  - $GM_OUT"
echo "  - $LUT_OUT"
