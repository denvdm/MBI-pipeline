#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# make_split.R  -  write the discovery/replication assignment.
#
# The assignment itself is computed in helpers/split_assign.R and depends only
# on the subject id and the label, so this script does not CREATE the split so
# much as write it down: running it twice on different cohorts gives the same
# label to every subject present in both. Nothing downstream has to read the
# file to stay consistent with it — stage 02 recomputes the same assignment from
# the same two config values — so the file exists to be consumed by the
# phenotypic and genetic analyses that live outside this repository, and to be
# archived with the run.
#
# WHAT THE SPLIT DOES AND DOES NOT DO HERE. It labels subjects. It does not
# decide who is scored: every subject who clears the inclusion rules receives
# GBI, MBI_raw and MBI, and the index is built identically for both halves. It
# governs the normative fit only if split.fit_on is set to "discovery", which is
# a sensitivity arm rather than the default, because the normative model
# contains no outcome and fitting it on the whole cohort is the more precise
# estimate.
#
# Run:
#   Rscript extras/make_split.R --config config.yaml --outdir <OUT>
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
})

source("helpers/run_helpers.R")
opt <- parse_common_args()
cfg <- read_config(opt$config)
source_local("split_assign.R")

sp <- read_split_cfg(cfg)

if (identical(sp$mode, "none")) {
  message("[split] split.mode is \"none\"; nothing to write. ",
          "Set split.mode: hash and split.label to produce an assignment.")
  quit(status = 0L)
}

features_file <- cfg$inputs$features_file %||% NULL
if (is.null(features_file)) stop("Missing config: inputs: features_file", call. = FALSE)
assert_file_exists(features_file, "features_file")

colmap  <- cfg$inputs$colmap %||% list()
col_id  <- colmap$id  %||% "id"
col_age <- colmap$age %||% "age"
col_sex <- colmap$sex %||% "sex"

hdr  <- names(fread(features_file, nrows = 0L))
want <- intersect(c(col_id, col_age, col_sex), hdr)
if (!col_id %in% want)
  stop("features_file has no id column '", col_id, "'", call. = FALSE)

dt <- fread(features_file, select = want, data.table = FALSE)
ids <- as.character(dt[[col_id]])

assign_df <- assign_split(ids, label = sp$label, discovery_frac = sp$discovery_frac)

# Balance across the covariates a stratified draw would have balanced on,
# reported rather than imposed.
covars <- list()
if (col_sex %in% names(dt)) covars[["sex"]] <- dt[[col_sex]]
if (col_age %in% names(dt)) {
  a    <- suppressWarnings(as.numeric(dt[[col_age]]))
  brks <- unique(stats::quantile(a, probs = seq(0, 1, 0.1), na.rm = TRUE))
  if (length(brks) > 2L)
    covars[["age_decile"]] <- cut(a, breaks = brks, include.lowest = TRUE, labels = FALSE)
}
bal <- split_balance(assign_df, covars)

out_stage <- stage_dir(opt$outdir, "split")
f_assign  <- file.path(out_stage, "split_assignment.txt")
f_train   <- file.path(out_stage, "train_ids.txt")
f_bal     <- file.path(out_stage, "split_balance.txt")

write_tsv(assign_df, f_assign)
write_tsv(assign_df[assign_df$split == "discovery", "id", drop = FALSE], f_train)
write_tsv(bal, f_bal)

n_disc <- sum(assign_df$split == "discovery")
cat(sprintf("[split] %d subjects | %d discovery (%.4f, requested %.4f) | %d replication\n",
            nrow(assign_df), n_disc, n_disc / nrow(assign_df), sp$discovery_frac,
            nrow(assign_df) - n_disc))
cat(sprintf("[split] label = \"%s\"  |  id checksum of discovery set = %s\n",
            sp$label, split_checksum(assign_df$id[assign_df$split == "discovery"])))
cat("\n[split] balance (a stratified draw would have forced these to the requested fraction):\n")
print(as.data.frame(bal), row.names = FALSE)

worst <- bal[bal$variable != "(overall)" & bal$n >= 30, , drop = FALSE]
if (nrow(worst)) {
  d <- max(abs(worst$frac_discovery - sp$discovery_frac))
  cat(sprintf("\n[split] largest deviation from the requested fraction, over levels with n >= 30: %.4f\n", d))
  if (d > 0.05)
    warning("[split] a covariate level deviates from the requested discovery ",
            "fraction by more than 0.05. At small N supply an external, ",
            "stratified assignment through TRAIN_IDS instead.", call. = FALSE)
}

cat("\nWrote:\n- ", f_assign, "\n- ", f_train, "\n- ", f_bal, "\n", sep = "")

write_run_info(out_stage, "make_split", opt$config, params = list(
  mode = sp$mode, label = sp$label,
  discovery_frac_requested = sp$discovery_frac,
  discovery_frac_realised  = round(n_disc / nrow(assign_df), 6),
  n_subjects = nrow(assign_df), n_discovery = n_disc,
  discovery_checksum = split_checksum(assign_df$id[assign_df$split == "discovery"])
))
