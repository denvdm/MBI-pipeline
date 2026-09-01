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
#   CANONICAL_CFG      config the sensitivity run inherits its index parameters
#                      from. Defaults to the config the canonical run recorded
#                      in compute_mbi/run_info.txt, so the two runs cannot drift
#                      apart. A sensitivity run built at different index
#                      parameters from the run it is compared against is not a
#                      template comparison; it is a comparison of two indices.
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

# Which config the canonical run used. Recorded by write_run_info(), so the
# inheritance below points at the exact file that produced the run being
# compared against rather than at whatever config happens to be lying around.
CANON_RUNINFO="$CANONICAL_OUTDIR/compute_mbi/run_info.txt"
if [[ -z "${CANONICAL_CFG:-}" ]]; then
  if [[ -f "$CANON_RUNINFO" ]]; then
    CANONICAL_CFG="$(awk '$1 == "config" { $1 = ""; sub(/^ +/, ""); print; exit }' "$CANON_RUNINFO")"
  fi
fi
[[ -n "${CANONICAL_CFG:-}" && -f "$CANONICAL_CFG" ]] || {
  echo "[ERROR] Cannot locate the canonical run's config." >&2
  echo "        Looked in: $CANON_RUNINFO" >&2
  echo "        Set CANONICAL_CFG explicitly. Refusing to guess index parameters:" >&2
  echo "        a sensitivity run at different parameters is not a template comparison." >&2
  exit 1
}

mkdir -p "$OUTDIR/logs"

echo "[INFO] PROJECT_ROOT:     $PROJECT_ROOT"
echo "[INFO] OUTDIR:           $OUTDIR"
echo "[INFO] CANONICAL_OUTDIR: $CANONICAL_OUTDIR"
echo "[INFO] CANONICAL_CFG:    $CANONICAL_CFG"
echo "[INFO] TEMPLATE_LABEL:   $TEMPLATE_LABEL"
echo "[INFO] ALT_PET_FILE:     $ALT_PET_FILE"
echo "[INFO] CFG (generated):  $CFG"
echo

# Sidecar config: the canonical config with ONLY the input paths overridden.
# Every index parameter — emphasis, ratio, clamp_z, beta, coverage_thr,
# required_modalities, subvol_min_mm3 — is inherited rather than restated, so
# the two runs cannot drift apart as the canonical config evolves.
CANONICAL_CFG="$CANONICAL_CFG" CFG_OUT="$CFG" \
MAPS_DIR="$MAPS_DIR" ALT_PET_FILE="$ALT_PET_FILE" ALT_ATLAS_FILE="$ALT_ATLAS_FILE" \
ALT_GM_FILE="$ALT_GM_FILE" ALT_LUT_FILE="$ALT_LUT_FILE" \
NORMDEV_LONG="$NORMDEV_LONG" FEAT_SIZES="$FEAT_SIZES" \
TEMPLATE_LABEL="$TEMPLATE_LABEL" CANONICAL_OUTDIR="$CANONICAL_OUTDIR" \
Rscript -e '
  cfg <- yaml::read_yaml(Sys.getenv("CANONICAL_CFG"))
  if (is.null(cfg$paths))  cfg$paths  <- list()
  if (is.null(cfg$inputs)) cfg$inputs <- list()
  cfg$paths$maps_dir    <- Sys.getenv("MAPS_DIR")
  cfg$inputs$pet_file   <- Sys.getenv("ALT_PET_FILE")
  cfg$inputs$atlas_file <- Sys.getenv("ALT_ATLAS_FILE")
  cfg$inputs$gm_file    <- Sys.getenv("ALT_GM_FILE")
  cfg$inputs$lut_file   <- Sys.getenv("ALT_LUT_FILE")
  # Reuse the canonical normative deviations (PET-independent, slow to compute).
  cfg$inputs$deviations_z_long <- Sys.getenv("NORMDEV_LONG")
  cfg$inputs$feat_sizes        <- Sys.getenv("FEAT_SIZES")
  # pet_weights must NOT be inherited: stage 03 has to read the alt extraction.
  cfg$inputs$pet_weights <- NULL
  hdr <- c(sprintf("# Sensitivity run: %s", Sys.getenv("TEMPLATE_LABEL")),
           sprintf("# Generated: %s", format(Sys.time())),
           sprintf("# Inherited from: %s", Sys.getenv("CANONICAL_CFG")),
           sprintf("# Canonical run reused for normdev: %s", Sys.getenv("CANONICAL_OUTDIR")),
           "# Only inputs.* differ from the inherited config. Index parameters are",
           "# inherited so this run and the run it is compared against are the same index.")
  writeLines(c(hdr, yaml::as.yaml(cfg)), Sys.getenv("CFG_OUT"))
' || { echo "[ERROR] Failed to build sidecar config from $CANONICAL_CFG" >&2; exit 1; }

echo "[OK] Wrote sidecar config. Index parameters inherited from the canonical run:"
CFG_OUT="$CFG" Rscript -e \
  'cat(yaml::as.yaml(list(mbi = yaml::read_yaml(Sys.getenv("CFG_OUT"))$mbi)))'
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
