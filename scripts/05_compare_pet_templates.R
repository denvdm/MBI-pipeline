#!/usr/bin/env Rscript
# Stage 05: Compare alt PET-template MBI outputs vs the canonical run.
#
# Tier 0: weight-vector regional correlation (PET weights, raw + after the
#         pipeline's map-weight transform).
# Tier 1: subject-level construct correlations (MBI_raw, GBI, MBI).
#
# Workflow:
#   1. Run the full pipeline once with the canonical PET template (e.g.
#      castrillon2023), output to canonical OUTDIR.
#   2. Prepare the alt PET template (e.g. vaishnavi2010) via
#      extras/prep_petmap_alt.sh, then run the full pipeline pointing at the
#      alt atlas/GM/LUT outputs, into an alt OUTDIR.
#   3. Run this script with --canonical_outdir pointing at the canonical run.
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(optparse)
})

source("helpers/run_helpers.R")

extra <- list(
  make_option("--canonical_outdir", type = "character", default = NULL,
              help = "Path to the canonical (e.g. Castrillon) run OUTDIR. Must contain extract_pet/ and compute_mbi/.",
              metavar = "DIR"),
  make_option("--template_label", type = "character", default = "alt",
              help = "Short label for the alt template, used in headers [default: %default]")
)

opt <- parse_common_args(extra_options = extra)
cfg <- read_config(opt$config)

source_local("compute_mbi.R")        # defines transform_map_weights(), aliased below
source_local("compare_pet_funcs.R")

if (is.null(opt$canonical_outdir) || !nzchar(opt$canonical_outdir)) {
  stop("Need --canonical_outdir", call. = FALSE)
}

alt_outdir    <- opt$outdir
alt_weights   <- file.path(alt_outdir, "extract_pet", "pet_fdg_featMedians.txt")
alt_scores    <- file.path(alt_outdir, "compute_mbi", "mbi_flagship_scores.txt")

canon_weights <- file.path(opt$canonical_outdir, "extract_pet", "pet_fdg_featMedians.txt")
canon_scores  <- file.path(opt$canonical_outdir, "compute_mbi", "mbi_flagship_scores.txt")

for (f in c(alt_weights, alt_scores, canon_weights, canon_scores)) {
  assert_file_exists(f, basename(f))
}

out_stage <- stage_dir(alt_outdir, "compare")

w_canon <- fread(canon_weights, data.table = FALSE)
w_alt   <- fread(alt_weights,   data.table = FALSE)
s_canon <- fread(canon_scores,  data.table = FALSE)
s_alt   <- fread(alt_scores,    data.table = FALSE)

# Resolve once from the CANONICAL map so both templates are transformed on
# the same weight scale; otherwise each would self-scale and the Tier 0
# comparison would silently be between two differently-scaled vectors.
ratio    <- resolve_ratio(w_canon$weight, cfg$mbi$ratio %||% "auto")
clamp_z  <- cfg$mbi$clamp_z  %||% 3
emphasis <- cfg$mbi$emphasis %||% 1

# Tier 0: weight-vector regional correlation
W <- compare_weights_regional(w_canon, w_alt,
                              ratio = ratio, clamp_z = clamp_z, emphasis = emphasis,
                              drop_ventricles = TRUE)
write_tsv(W$summary,  file.path(out_stage, "weight_vector_correlations.tsv"))
write_tsv(W$regional, file.path(out_stage, "weight_vector_regional.tsv"))

# Tier 1: subject-level construct correlation
S <- compare_scores_subject(s_canon, s_alt)
write_tsv(S$summary, file.path(out_stage, "subject_level_correlations.tsv"))
fwrite(S$subject,
       file = file.path(out_stage, "subject_level_scores.tsv.gz"),
       sep = "\t", compress = "gzip")

# Run-level metadata for downstream provenance
meta <- data.frame(
  field = c("template_label", "canonical_outdir", "alt_outdir",
            "ratio", "clamp_z", "emphasis",
            "n_regions_compared", "n_subjects_compared", "timestamp"),
  value = c(opt$template_label, opt$canonical_outdir, alt_outdir,
            as.character(ratio), as.character(clamp_z), as.character(emphasis),
            as.character(nrow(W$regional)), as.character(nrow(S$subject)),
            format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
)
write_tsv(meta, file.path(out_stage, "compare_meta.tsv"))

write_run_info(out_stage, "05_compare_pet_templates", opt$config, params = list(
  template_label = opt$template_label, canonical_outdir = opt$canonical_outdir,
  ratio = ratio, clamp_z = clamp_z, emphasis = emphasis,
  n_regions_compared = nrow(W$regional), n_subjects_compared = nrow(S$subject)
))

cat("\n==== Tier 0: PET weight-vector regional correlation ====\n")
cat("Template:", opt$template_label, "vs canonical\n\n")
print(W$summary, row.names = FALSE)

cat("\n==== Tier 1: Subject-level construct correlation ====\n\n")
print(S$summary, row.names = FALSE)

cat("\nAll outputs written to:", out_stage, "\n")
