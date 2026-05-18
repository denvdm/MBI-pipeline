#!/usr/bin/env Rscript
# Stage 05: Compare alt PET-template MBI outputs vs the canonical run.
#
# Tier 0: weight-vector regional correlation (PET weights, raw + after the
#         pipeline's gamma + p99-winsor transform).
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

gamma    <- cfg$mbi$gamma    %||% 2.5
winsor_q <- cfg$mbi$winsor_q %||% 0.99

# Tier 0: weight-vector regional correlation
W <- compare_weights_regional(w_canon, w_alt,
                              gamma = gamma, winsor_q = winsor_q,
                              drop_ventricles = TRUE)
write.table(W$summary,
            file = file.path(out_stage, "weight_vector_correlations.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(W$regional,
            file = file.path(out_stage, "weight_vector_regional.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Tier 1: subject-level construct correlation
S <- compare_scores_subject(s_canon, s_alt)
write.table(S$summary,
            file = file.path(out_stage, "subject_level_correlations.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
fwrite(S$subject,
       file = file.path(out_stage, "subject_level_scores.tsv.gz"),
       sep = "\t", compress = "gzip")

# Run-level metadata for downstream provenance
meta <- data.frame(
  field = c("template_label", "canonical_outdir", "alt_outdir",
            "gamma", "winsor_q",
            "n_regions_compared", "n_subjects_compared", "timestamp"),
  value = c(opt$template_label, opt$canonical_outdir, alt_outdir,
            as.character(gamma), as.character(winsor_q),
            as.character(nrow(W$regional)), as.character(nrow(S$subject)),
            format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
)
write.table(meta, file = file.path(out_stage, "compare_meta.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n==== Tier 0: PET weight-vector regional correlation ====\n")
cat("Template:", opt$template_label, "vs canonical\n\n")
print(W$summary, row.names = FALSE)

cat("\n==== Tier 1: Subject-level construct correlation ====\n\n")
print(S$summary, row.names = FALSE)

cat("\nAll outputs written to:", out_stage, "\n")
