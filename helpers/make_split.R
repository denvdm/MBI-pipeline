#!/usr/bin/env Rscript

# -------------------------------------------------------------------
# make_split.R
# Create stratified 70/30 discovery/replication split.
# Stratified by sex and age decile for demographic balance.
# Reads features_wide.txt (needs id, sex, age columns).
# Writes split_ids.txt (id, split) and train_ids.txt (id, discovery only).
#
# Environment overrides:
#   FEATURES_FILE  path to features_wide.txt (required)
#   SPLIT_FILE     output path for split_ids.txt (default: ./split_ids.txt)
#   SPLIT_SEED     RNG seed (default: 42)
#   SPLIT_FRAC     discovery fraction (default: 0.7)
# -------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))

feats_file <- Sys.getenv("FEATURES_FILE")
out_file   <- Sys.getenv("SPLIT_FILE", "split_ids.txt")
seed       <- as.integer(Sys.getenv("SPLIT_SEED", "42"))
frac_disc  <- as.numeric(Sys.getenv("SPLIT_FRAC", "0.7"))

if (!nzchar(feats_file) || !file.exists(feats_file))
  stop("Set FEATURES_FILE to features_wide.txt path")

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

dt <- fread(feats_file, select = c("id", "sex", "age"))
dt <- dt[complete.cases(dt)]

# stratify by sex x age decile
# unique(): with integer ages, decile breaks can coincide and cut() errors on
# duplicated breaks. Deduplicating merges the affected strata, which is the
# sane degradation.
brks <- unique(quantile(dt$age, probs = seq(0, 1, 0.1), na.rm = TRUE))
if (length(brks) < 2L) stop("age has no spread; cannot stratify")
dt[, age_dec := cut(age, breaks = brks, include.lowest = TRUE, labels = FALSE)]
dt[, stratum := paste0(sex, "_", age_dec)]

set.seed(seed)
dt[, split := "replication"]
dt[, split := {
  n <- .N
  n_disc <- round(n * frac_disc)
  idx <- sample(.N, n_disc)
  s <- rep("replication", .N)
  s[idx] <- "discovery"
  s
}, by = stratum]

out <- dt[, .(id, split)]
fwrite(out, out_file, sep = "\t")

# Also write discovery-only IDs in the format expected by 02_run_normdevs.R's TRAIN_IDS hook.
train_file <- file.path(dirname(out_file), "train_ids.txt")
fwrite(dt[split == "discovery", .(id)], train_file, sep = "\t")

n_disc <- sum(out$split == "discovery")
n_repl <- sum(out$split == "replication")
message(sprintf("Split: %d discovery (%.1f%%), %d replication (%.1f%%). Seed=%d",
                n_disc, 100 * n_disc / nrow(out),
                n_repl, 100 * n_repl / nrow(out), seed))
message("Wrote: ", out_file)
message("Wrote: ", train_file)
