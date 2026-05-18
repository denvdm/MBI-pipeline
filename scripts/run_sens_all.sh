#!/usr/bin/env bash
# Sensitivity run for an alternative volumetric PET template.
#
# Reuses the canonical run's normative-deviation outputs (PET-independent,
# slow). Re-runs only PET-dependent stages plus a comparison stage.
#
# Stages:
#   01  PET extraction with alt template
#   03  Compute MBI (alt PET weights x canonical normdev)
#   04  PET-MRI alignment QC (alt)
#   05  Compare alt vs canonical (Tier 0 weight-vector + Tier 1 subject-level)
#
# Required env vars (set by submit script or interactively):
#   PROJECT_ROOT       MBI-pipeline repo dir
#   MAPS_DIR           PET resources dir
#   CANONICAL_OUTDIR   existing canonical run with normdev/, extract_pet/, compute_mbi/
#   TEMPLATE_LABEL     short label, used in OUTDIR name (e.g. "dellarosa2014")
#   ALT_PET_FILE       alt PET .nii.gz (in MNI152 space)
#   ALT_ATLAS_FILE     atlas resampled to alt PET grid (from prep_petmap_alt.sh)
#   ALT_GM_FILE        GM-prob resampled to alt PET grid
#   ALT_LUT_FILE       LUT csv for the alt-grid atlas
#
# Optional env vars:
#   OUTDIR             default: $PROJECT_ROOT/runs/sens/${TEMPLATE_LABEL}_<timestamp>
set -euo pipefail

: "${PROJECT_ROOT:?Need PROJECT_ROOT}"
: "${MAPS_DIR:?Need MAPS_DIR}"
: "${CANONICAL_OUTDIR:?Need CANONICAL_OUTDIR}"
: "${ALT_PET_FILE:?Need ALT_PET_FILE}"
: "${ALT_ATLAS_FILE:?Need ALT_ATLAS_FILE}"
: "${ALT_GM_FILE:?Need ALT_GM_FILE}"
: "${ALT_LUT_FILE:?Need ALT_LUT_FILE}"
: "${TEMPLATE_LABEL:?Need TEMPLATE_LABEL}"

CANONICAL_OUTDIR="${CANONICAL_OUTDIR%/}"

NORMDEV_LONG="$CANONICAL_OUTDIR/normdev/deviations_z_long.txt"
FEAT_SIZES="$CANONICAL_OUTDIR/normdev/feats_avgSize.txt"
CANON_PET="$CANONICAL_OUTDIR/extract_pet/pet_fdg_featMedians.txt"
CANON_MBI="$CANONICAL_OUTDIR/compute_mbi/mbi_flagship_scores.txt"
for f in "$NORMDEV_LONG" "$FEAT_SIZES" "$CANON_PET" "$CANON_MBI"; do
  [[ -f "$f" ]] || { echo "[ERROR] Missing canonical input: $f" >&2; exit 1; }
done

for f in "$ALT_PET_FILE" "$ALT_ATLAS_FILE" "$ALT_GM_FILE" "$ALT_LUT_FILE"; do
  [[ -f "$f" ]] || { echo "[ERROR] Missing alt input: $f" >&2; exit 1; }
done

OUTDIR="${OUTDIR:-$(dirname "$CANONICAL_OUTDIR")/sens/${TEMPLATE_LABEL}_$(date +%Y%m%d_%H%M%S)}"
CFG="$OUTDIR/config_sens.yaml"

mkdir -p "$OUTDIR/logs"

echo "[INFO] PROJECT_ROOT:     $PROJECT_ROOT"
echo "[INFO] OUTDIR:           $OUTDIR"
echo "[INFO] CANONICAL_OUTDIR: $CANONICAL_OUTDIR"
echo "[INFO] TEMPLATE_LABEL:   $TEMPLATE_LABEL"
echo "[INFO] ALT_PET_FILE:     $ALT_PET_FILE"
echo "[INFO] CFG (generated):  $CFG"
echo

# Auto-generated sidecar config: override PET inputs, reuse canonical normdev
cat > "$CFG" <<EOF
# Auto-generated config for sensitivity run: $TEMPLATE_LABEL
# Generated: $(date)
# Canonical run used for normdev reuse: $CANONICAL_OUTDIR
paths:
  maps_dir: "$MAPS_DIR"

inputs:
  pet_file:   "$ALT_PET_FILE"
  atlas_file: "$ALT_ATLAS_FILE"
  gm_file:    "$ALT_GM_FILE"
  lut_file:   "$ALT_LUT_FILE"

  # Reuse canonical normative deviations (PET-independent, slow to compute)
  deviations_z_long: "$NORMDEV_LONG"
  feat_sizes:        "$FEAT_SIZES"
  # pet_weights left unset: stage 03 uses alt extract_pet/ output

mbi:
  beta: 0.25
  gamma: 2.5
  alpha: null
  winsor_q: 0.99
  coverage_thr: 0.80
  subvol_gate_pct: 0.20

qc:
  mods_use: ["CT", "SA", "SubVol"]
EOF
echo "[OK] Wrote sidecar config"
echo

cd "$PROJECT_ROOT"

echo "[STEP] 01 PET weights (alt: $TEMPLATE_LABEL)"
Rscript "$PROJECT_ROOT/scripts/01_extract_petvals.R" \
  -c "$CFG" -o "$OUTDIR" --overwrite \
  > "$OUTDIR/logs/01_pet.log" 2>&1
echo "[OK] 01 done"
echo

echo "[STEP] 03 Compute MBI"
Rscript "$PROJECT_ROOT/scripts/03_run_mbi.R" \
  -c "$CFG" -o "$OUTDIR" --overwrite \
  > "$OUTDIR/logs/03_mbi.log" 2>&1
echo "[OK] 03 done"
echo

echo "[STEP] 04 QC alignment"
Rscript "$PROJECT_ROOT/scripts/04_relate_mbi_pet.R" \
  -c "$CFG" -o "$OUTDIR" --overwrite \
  > "$OUTDIR/logs/04_qc.log" 2>&1
echo "[OK] 04 done"
echo

echo "[STEP] 05 Compare alt vs canonical"
Rscript "$PROJECT_ROOT/scripts/05_compare_pet_templates.R" \
  -c "$CFG" -o "$OUTDIR" \
  --canonical_outdir "$CANONICAL_OUTDIR" \
  --template_label "$TEMPLATE_LABEL" \
  --overwrite \
  2>&1 | tee "$OUTDIR/logs/05_compare.log"
echo "[OK] 05 done"
echo

echo "[DONE] Sensitivity run complete: $TEMPLATE_LABEL"
echo "       Outputs: $OUTDIR"
echo "       Headline correlation: $OUTDIR/compare/subject_level_correlations.tsv"
