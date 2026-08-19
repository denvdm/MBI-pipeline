#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# sweep_emphasis.R  -  parameter SENSITIVITY strip.
#
# This is NOT a selection step. Fix the parameters first and record them in
# extras/SELECTION_RULE.md, before any outcome is inspected; then run this to
# report how the index behaves as the two reportable parameters move, so a
# reader can see the result does not hinge on them.
#
#   emphasis : the only parameter that changes MBI
#   beta     : topography (0) versus burden (1)
#
# `ratio` is not swept: it is algebraically inert for MBI, since residualising
# on GBI removes the term it scales. extras/verify_invariance.R checks that.
#
# Per cell it reports:
#   r2_raw_gbi  shared variance of mbi_raw with gbi. A property of a discarded
#               intermediate - it moves while MBI does not. Reported for
#               completeness, NOT as a criterion.
#   n_eff       effective number of contributing features (inverse Herfindahl)
#   mbi_map_r   cor(MBI, per-subject r_tilt), the construct check
#   top_share   single largest feature influence share
#   share_lo    total influence share of the bottom CMRglc quintile
#   share_hi    total influence share of the top CMRglc quintile
#   tail_bal    share_lo / share_hi   (>1 means the low tail dominates)
#   share_mid   total share of the middle three quintiles
#
# Writes:
#   <outdir>/calibration/sweep_emphasis_results.tsv
#   <outdir>/calibration/sweep_emphasis_toploaders.tsv
#
# Run:
#   Rscript extras/sweep_emphasis.R --config config.yaml --outdir <OUT>
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(optparse)
})

# ------------------------------------------------------------
# USER: EDIT HERE
# ------------------------------------------------------------
emphases <- c(1, 1.5, 2, 3)     # config default is 1
betas    <- c(0, 1)             # config default is 0
top_n_detail <- 12
mods_use <- c("CT", "SA", "SubVol")   # modalities entering the r_tilt target
# ------------------------------------------------------------

source("helpers/run_helpers.R")
opt <- parse_common_args()
cfg <- read_config(opt$config)

source_local("compute_mbi.R")

stage_dir_local <- function(outdir, stage) {
  d <- file.path(outdir, stage); dir.create(d, recursive = TRUE, showWarnings = FALSE); d
}

subvol_min_mm3 <- cfg$mbi$subvol_min_mm3 %||% 1000
ratio          <- cfg$mbi$ratio          %||% "auto"
clamp_z        <- cfg$mbi$clamp_z        %||% 3
coverage_thr   <- cfg$mbi$coverage_thr   %||% 0.80

# Read the subject-inclusion rule from the config too. Omitting it would score
# every subject who clears coverage_thr on any modality, which is a LARGER and
# differently composed sample than stage 03 uses -- so the sensitivity strip
# would describe a different cohort from the index it is meant to characterise.
required_modalities <- cfg$mbi$required_modalities %||% NULL

normdev_dir <- file.path(opt$outdir, "normdev")
pet_dir     <- file.path(opt$outdir, "extract_pet")
cal_dir     <- stage_dir_local(opt$outdir, "calibration")

zs_file    <- cfg$inputs$deviations_z_long %||% file.path(normdev_dir, "deviations_z_long.txt")
sizes_file <- cfg$inputs$feat_sizes        %||% file.path(normdev_dir, "feats_avgSize.txt")
fdg_file   <- cfg$inputs$pet_weights       %||% file.path(pet_dir, "pet_fdg_featMedians.txt")

assert_file_exists(zs_file,    "deviations_z_long")
assert_file_exists(sizes_file, "feats_avgSize")
assert_file_exists(fdg_file,   "pet_fdg_featMedians")

zs      <- fread(zs_file,    data.table = FALSE) %>%
  mutate(region = as.character(region), modality = as.character(modality))
size_df <- fread(sizes_file, data.table = FALSE) %>%
  mutate(region = as.character(region), modality = as.character(modality))
fdg     <- fread(fdg_file,   data.table = FALSE) %>%
  mutate(region = as.character(region))

fdg_use <- if ("weight" %in% names(fdg)) {
  fdg %>% transmute(region, w_pet = as.numeric(weight))
} else if ("w_pet" %in% names(fdg)) {
  fdg %>% transmute(region, w_pet = as.numeric(w_pet))
} else stop("[ERROR] pet_fdg_featMedians.txt must contain 'weight' or 'w_pet'", call. = FALSE)

# --- SubVol gate: identical to stage 03 -----------------------------------
keep_sub <- size_df %>% filter(modality == "SubVol",
                               is.finite(size), size >= subvol_min_mm3) %>% pull(region)
zs <- zs %>% filter(modality != "SubVol" | region %in% keep_sub)

# --- r_tilt target: same orientation logic as MBI -------------------------
zo <- compute_oriented(zs)
r_tilt_df <- zo %>%
  filter(modality %in% mods_use) %>%
  select(id, region, z_oriented) %>%
  inner_join(fdg_use, by = "region") %>%
  group_by(id) %>%
  summarise(r_tilt = suppressWarnings(
    cor(z_oriented, w_pet, method = "spearman", use = "pairwise.complete.obs")),
    .groups = "drop")

# --- CMRglc quintiles, fixed ONCE so cells stay comparable ----------------
pet_region <- zs %>% distinct(region) %>%
  inner_join(fdg_use, by = "region") %>%
  filter(is.finite(w_pet), w_pet > 0)
qbrk <- quantile(pet_region$w_pet, probs = seq(0, 1, 0.2), na.rm = TRUE)
pet_region <- pet_region %>%
  mutate(cmr_quintile = as.integer(cut(w_pet, breaks = qbrk, include.lowest = TRUE)),
         cmr_pctile   = rank(w_pet) / length(w_pet))

cat(sprintf("[INFO] %d regions with a usable map weight | ratio=%s clamp_z=%.2f\n",
            nrow(pet_region), as.character(ratio), clamp_z))

grid <- expand.grid(emphasis = emphases, beta = betas,
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
cat(sprintf("[INFO] Sweeping %d cells\n", nrow(grid)))

detail_rows <- list()
score_cols  <- list()

eval_one <- function(E, B) {
  res <- suppressMessages(compute_flagship_abg(
    zs = zs, fdg = fdg_use %>% rename(weight = w_pet), sizes = size_df,
    ratio = ratio, clamp_z = clamp_z, emphasis = E,
    beta = B, coverage_thr = coverage_thr,
    required_modalities = required_modalities
  ))
  scores <- as.data.frame(res$scores)

  r_rg <- suppressWarnings(cor(scores$gbi, scores$mbi_raw, use = "pairwise.complete.obs"))
  tmp  <- inner_join(scores[, c("id", "mbi")], r_tilt_df, by = "id")
  mbi_map_r <- suppressWarnings(cor(tmp$mbi, tmp$r_tilt, use = "pairwise.complete.obs"))

  sh <- res$concentration$shares %>%
    left_join(pet_region, by = "region") %>%
    arrange(desc(share))

  qs <- sh %>% group_by(cmr_quintile) %>%
    summarise(s = sum(share, na.rm = TRUE), .groups = "drop")
  gq <- function(k) { v <- qs$s[qs$cmr_quintile == k]; if (length(v)) v else 0 }
  share_lo <- gq(1); share_hi <- gq(5)

  detail_rows[[length(detail_rows) + 1]] <<- sh %>%
    slice_head(n = top_n_detail) %>%
    mutate(emphasis = E, beta = B, rank = row_number()) %>%
    select(emphasis, beta, rank, modality, region, share,
           cmrglc = w_pet, cmr_pctile, cmr_quintile)

  # Keep the subject-level MBI for each cell so the cells can be compared as
  # indices, not only as summary statistics. Two settings can look similar on
  # N_eff and tail balance and still rank subjects differently; the correlation
  # between them is the number that says whether the choice matters.
  score_cols[[sprintf("emphasis_%g__beta_%g", E, B)]] <<-
    stats::setNames(scores$mbi, as.character(scores$id))

  data.frame(
    emphasis = E, beta = B,
    r2_raw_gbi = r_rg^2,
    n_eff      = res$concentration$n_eff,
    mbi_map_r  = mbi_map_r,
    top_share  = max(sh$share, na.rm = TRUE),
    share_lo   = share_lo,
    share_hi   = share_hi,
    tail_bal   = ifelse(share_hi > 0, share_lo / share_hi, NA_real_),
    share_mid  = gq(2) + gq(3) + gq(4),
    n_features = res$concentration$n_features,
    ratio_used     = res$diagnostics$ratio,
    realised_ratio = res$diagnostics$realised_ratio
  )
}

evals <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  cat(sprintf("  ... cell %d/%d (emphasis=%g, beta=%g)\n",
              i, nrow(grid), grid$emphasis[i], grid$beta[i]))
  eval_one(grid$emphasis[i], grid$beta[i])
}))

out_file <- file.path(cal_dir, "sweep_emphasis_results.tsv")
write_tsv(evals, out_file)
cat("[OK] Wrote: ", out_file, "\n", sep = "")

det <- bind_rows(detail_rows)
det_file <- file.path(cal_dir, "sweep_emphasis_toploaders.tsv")
write_tsv(det, det_file)
cat("[OK] Wrote: ", det_file, "\n", sep = "")

cat("\n[OK] Sensitivity strip (the reported row is the one in your config):\n")
print(evals %>% arrange(beta, emphasis) %>%
        mutate(across(where(is.numeric), ~round(.x, 4))), row.names = FALSE)

# --- Do the cells actually rank subjects differently? ----------------------
# The summary statistics above describe the WEIGHTS. This describes the INDEX.
# If two settings correlate at 0.99 the choice between them is presentational;
# if they correlate at 0.8 they are different measurements and the choice needs
# defending in the methods rather than a sensitivity row.
ids_all <- Reduce(intersect, lapply(score_cols, names))
if (length(ids_all) > 2 && length(score_cols) > 1) {
  S <- do.call(cbind, lapply(score_cols, function(v) v[ids_all]))
  colnames(S) <- names(score_cols)
  Cmat <- suppressWarnings(stats::cor(S, use = "pairwise.complete.obs"))

  cor_file <- file.path(cal_dir, "sweep_emphasis_score_correlations.tsv")
  write_tsv(data.frame(cell = rownames(Cmat), round(Cmat, 6),
                       check.names = FALSE), cor_file)
  cat("[OK] Wrote: ", cor_file, "\n", sep = "")

  cat(sprintf("\n[OK] Subject-level cor(MBI) between cells, n = %d subjects:\n",
              length(ids_all)))
  print(round(Cmat, 4))

  ref <- sprintf("emphasis_%g__beta_%g", cfg$mbi$emphasis %||% 1, cfg$mbi$beta %||% 0)
  if (ref %in% colnames(Cmat)) {
    off <- Cmat[, ref][setdiff(rownames(Cmat), ref)]
    cat(sprintf("\nAgainst the configured cell (%s): min r = %.4f (%s)\n",
                ref, min(off), names(off)[which.min(off)]))
    cat("A cell correlating above ~0.95 with the configured one is the same\n",
        "index differently presented. Below that, it is a different index.\n", sep = "")
  }
}

write_run_info(cal_dir, "sweep_emphasis", opt$config, params = list(
  emphases = emphases, betas = betas, ratio = as.character(ratio), clamp_z = clamp_z,
  coverage_thr = coverage_thr, subvol_min_mm3 = subvol_min_mm3,
  mods_use = mods_use, n_cells = nrow(grid)
))
