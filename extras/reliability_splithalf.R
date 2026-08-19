#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# reliability_splithalf.R  -  internal consistency of MBI and GBI.
#
# WHY THIS EXISTS. R^2(mbi_raw, gbi) is high, which invites the objection that
# residualising two near-collinear variables leaves only noise. That objection
# is correct when the two variables are separate NOISY MEASUREMENTS of the same
# thing. It does not apply here: mbi_raw and gbi are two exact weighted averages
# of the SAME regional deviations, so the residual is a linear contrast of those
# deviations, not a remainder.
#
# The question that DOES survive is whether that contrast direction has usable
# signal-to-noise. R^2 cannot answer it — it moves with `ratio`, which does not
# change MBI at all. Reliability can. This script builds the index twice from
# DISJOINT feature sets, left hemisphere and right hemisphere, and correlates
# the two. Independent noise gives zero; systematic biology does not.
#
# Report the GBI figure alongside MBI's as a benchmark. Expect MBI to be the
# less reliable of the two: a contrast discards the shared component that makes
# a global mean stable. `attenuation` in the output is the factor by which a
# true correlation is expected to shrink at the measured reliability, which is
# the number to quote when interpreting effect sizes.
#
# CAVEAT to state wherever this is reported: the two hemispheres share scanner,
# session and head-motion noise, so this is an UPPER BOUND on reliability.
# Test-retest ICC across a repeat-imaging visit is the clean measurement.
# Unpaired structures (e.g. Brain_Stem) fall in neither half and are excluded
# from both, so the full index is if anything slightly better than this
# estimate.
#
# Run:
#   Rscript extras/reliability_splithalf.R --config config.yaml --outdir <OUT>
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(optparse)
})

source("helpers/run_helpers.R")
opt <- parse_common_args()
cfg <- read_config(opt$config)
source_local("compute_mbi.R")

ratio          <- cfg$mbi$ratio          %||% "auto"
clamp_z        <- cfg$mbi$clamp_z        %||% 3
emphasis       <- cfg$mbi$emphasis       %||% 1
beta           <- cfg$mbi$beta           %||% 0
coverage_thr   <- cfg$mbi$coverage_thr   %||% 0.80
subvol_min_mm3 <- cfg$mbi$subvol_min_mm3 %||% 1000

# Same subject rule as stage 03. A reliability estimate computed on a different
# sample from the reported index is not a reliability estimate of the reported
# index.
required_modalities <- cfg$mbi$required_modalities %||% NULL

normdev_dir <- file.path(opt$outdir, "normdev")
pet_dir     <- file.path(opt$outdir, "extract_pet")

zs  <- fread(cfg$inputs$deviations_z_long %||% file.path(normdev_dir, "deviations_z_long.txt"),
             data.table = FALSE)
sz  <- fread(cfg$inputs$feat_sizes        %||% file.path(normdev_dir, "feats_avgSize.txt"),
             data.table = FALSE)
pet <- fread(cfg$inputs$pet_weights       %||% file.path(pet_dir, "pet_fdg_featMedians.txt"),
             data.table = FALSE)

keep_sub <- sz %>% filter(modality == "SubVol", is.finite(size),
                          size >= subvol_min_mm3) %>% pull(region)
zs <- zs %>% filter(modality != "SubVol" | region %in% keep_sub)

# Hemisphere assignment. "M" = midline / unpaired, excluded from both halves.
hemi <- function(r) {
  ifelse(grepl("^lh[._]|^Left|_lh$|^L_", r, ignore.case = TRUE), "L",
  ifelse(grepl("^rh[._]|^Right|_rh$|^R_", r, ignore.case = TRUE), "R", "M"))
}
zs$h <- hemi(zs$region)

feat <- zs %>% distinct(region, modality, h)
cat("[INFO] features per hemisphere:\n")
print(table(feat$h))
if (any(feat$h == "M"))
  cat("[INFO] unpaired (excluded from both halves): ",
      paste(sort(unique(zs$region[zs$h == "M"])), collapse = ", "), "\n", sep = "")

go <- function(sub) suppressMessages(compute_flagship_abg(
  zs = sub, fdg = pet, sizes = sz,
  ratio = ratio, clamp_z = clamp_z, emphasis = emphasis,
  beta = beta, coverage_thr = coverage_thr,
  required_modalities = required_modalities))

L <- go(zs %>% filter(h == "L"))
R <- go(zs %>% filter(h == "R"))

sb <- function(r) 2 * r / (1 + r)   # Spearman-Brown, correcting for the halving

m  <- inner_join(L$scores[, c("id","mbi")], R$scores[, c("id","mbi")], by = "id")
g  <- inner_join(L$scores[, c("id","gbi")], R$scores[, c("id","gbi")], by = "id")
r_m <- cor(m$mbi.x, m$mbi.y, use = "complete.obs")
r_g <- cor(g$gbi.x, g$gbi.y, use = "complete.obs")

out <- data.frame(
  index          = c("MBI", "GBI"),
  n_subjects     = c(sum(complete.cases(m)), sum(complete.cases(g))),
  n_features_L   = c(L$concentration$n_features, L$concentration$n_features),
  n_features_R   = c(R$concentration$n_features, R$concentration$n_features),
  r_split_half   = c(r_m, r_g),
  spearman_brown = c(sb(r_m), sb(r_g)),
  attenuation    = c(sqrt(sb(r_m)), sqrt(sb(r_g)))
)

stage <- file.path(opt$outdir, "reliability")
dir.create(stage, recursive = TRUE, showWarnings = FALSE)
f <- file.path(stage, "splithalf_reliability.tsv")
write_tsv(out, f)

cat("\n[OK] Wrote: ", f, "\n\n", sep = "")
print(out, row.names = FALSE)
cat("\nUPPER BOUND: hemispheres share scanner/session/motion noise.\n")
cat("`attenuation` is the factor by which a true correlation is expected to\n")
cat("shrink at this reliability (sqrt of Spearman-Brown).\n")
