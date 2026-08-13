#!/usr/bin/env Rscript
# Calibrate (gamma, beta) for the MBI PET/size weighting on an existing run.
# Sweeps a grid, evaluates each (gamma, beta) via subject-level alignment
# between MBI and the per-subject r_tilt target (Spearman of regional
# deviations vs PET weights), writes a results TSV.
#
# Reads:
#   <outdir>/normdev/deviations_z_long.txt
#   <outdir>/normdev/feats_avgSize.txt
#   <outdir>/extract_pet/pet_fdg_featMedians.txt
# Writes:
#   <outdir>/calibration/calibration_grid_results.tsv

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(optparse)
})

# ------------------------------------------------------------
# USER: EDIT HERE (calibration grid)
# ------------------------------------------------------------
# `ratio` = maximum weight ratio between the highest- and lowest-demand region.
# It replaces gamma: under the robust_z / rank parameterisation, w^gamma is
# exactly (ratio^gamma)^(u/2), so exponentiating IS a change of ratio. Swept on
# a roughly geometric grid because the parameter is multiplicative.
ratios <- c(1.25, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10, 13, 16, 20, 26)
betas  <- seq(0.0, 0.5, by = 0.05)   # 11 values

# Constraints applied when reporting the recommended cell (see bottom of file).
# n_eff  : effective number of features carrying the index (concentration).
# r2 band: MBI_raw vs GBI shared variance. Too LOW and the two weightings barely
#          differ, so the residual is noise; too HIGH and MBI is a rescaling of
#          GBI. This bounds `ratio` from both sides and makes the choice
#          well-posed rather than a preference.
n_eff_min <- 30
r2_band   <- c(0.75, 0.95)

# Modalities used for the r_tilt target (alignment of deviations with PET weights)
#
# !! CHANGED from c("CT","SA") to all three modalities.
# !! The calibration target was previously computed on CORTEX ONLY, while the
# !! index it is tuned against is overwhelmingly subcortical (SubVol carries
# !! r = 0.46 with MBI; CT = 0.03, SA = 0.01). gamma and beta were therefore
# !! selected to maximise agreement between a subcortically-driven index and a
# !! cortex-only alignment target - which is part of why mbi_map_r ceilings
# !! around 0.55. The target must span the same features as the index.
mods_use <- c("CT", "SA", "SubVol")
# ------------------------------------------------------------

source("helpers/run_helpers.R")
opt <- parse_common_args()
cfg <- read_config(opt$config)

source_local("compute_mbi.R")

# Local helper for stage dir (avoids tight coupling with run_helpers internals)
stage_dir_local <- function(outdir, stage) {
  d <- file.path(outdir, stage)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

subvol_min_mm3 <- cfg$mbi$subvol_min_mm3 %||% 1000
map_method     <- cfg$mbi$map_method     %||% "robust_z"
clamp_z        <- cfg$mbi$clamp_z        %||% 2
coverage_thr   <- cfg$mbi$coverage_thr   %||% 0.80
gamma          <- cfg$mbi$gamma          %||% 2.5   # legacy path only
winsor_q       <- cfg$mbi$winsor_q       %||% 0.99  # legacy path only

if (!is.null(cfg$mbi$subvol_gate_pct))
  warning("config: 'subvol_gate_pct' is deprecated and ignored; ",
          "using absolute 'subvol_min_mm3' = ", subvol_min_mm3, call. = FALSE)

normdev_dir <- file.path(opt$outdir, "normdev")
pet_dir     <- file.path(opt$outdir, "extract_pet")
cal_dir     <- stage_dir_local(opt$outdir, "calibration")

zs_file    <- cfg$inputs$deviations_z_long %||% file.path(normdev_dir, "deviations_z_long.txt")
sizes_file <- cfg$inputs$feat_sizes        %||% file.path(normdev_dir, "feats_avgSize.txt")
fdg_file   <- cfg$inputs$pet_weights       %||% file.path(pet_dir, "pet_fdg_featMedians.txt")

assert_file_exists(zs_file,    "deviations_z_long")
assert_file_exists(sizes_file, "feats_avgSize")
assert_file_exists(fdg_file,   "pet_fdg_featMedians")

zs      <- fread(zs_file, data.table = FALSE)
size_df <- fread(sizes_file, data.table = FALSE)
fdg     <- fread(fdg_file, data.table = FALSE)

# Validate expected columns
req_zs <- c("id", "modality", "region", "z")
req_sz <- c("modality", "region", "size")
if (!all(req_zs %in% names(zs))) {
  stop("[ERROR] deviations_z_long.txt must contain: ", paste(req_zs, collapse = ", "), call. = FALSE)
}
if (!all(req_sz %in% names(size_df))) {
  stop("[ERROR] feats_avgSize.txt must contain: ", paste(req_sz, collapse = ", "), call. = FALSE)
}
if (!("region" %in% names(fdg))) {
  stop("[ERROR] pet_fdg_featMedians.txt must contain: region", call. = FALSE)
}

# Harmonize PET weight col name
fdg_use <- if ("weight" %in% names(fdg)) {
  fdg %>% transmute(region = as.character(region), w_pet = as.numeric(weight))
} else if ("w_pet" %in% names(fdg)) {
  fdg %>% transmute(region = as.character(region), w_pet = as.numeric(w_pet))
} else {
  stop("[ERROR] pet_fdg_featMedians.txt must contain 'weight' or 'w_pet'", call. = FALSE)
}

# -----------------------------
# SubVol gating (matches main pipeline): ABSOLUTE volume threshold, so this
# gate and align_feature_set() commute and the order never has to be defended.
# -----------------------------
keep_sub <- size_df %>%
  filter(modality == "SubVol") %>%
  filter(is.finite(size), size >= subvol_min_mm3) %>%
  pull(region)

zs <- zs %>% filter(modality != "SubVol" | region %in% keep_sub)

# -----------------------------
# Compute r_tilt using the SAME orientation logic as MBI
# -----------------------------
zo <- compute_oriented(zs)

mods_use <- unique(trimws(mods_use))
if (length(mods_use) == 0) stop("[ERROR] mods_use is empty", call. = FALSE)

r_tilt_df <- zo %>%
  filter(modality %in% mods_use) %>%
  select(id, region, z_oriented) %>%
  inner_join(fdg_use, by = "region") %>%
  group_by(id) %>%
  summarise(
    r_tilt = suppressWarnings(cor(z_oriented, w_pet, method = "spearman", use = "pairwise.complete.obs")),
    .groups = "drop"
  )

# -----------------------------
# Grid evaluation
# -----------------------------
grid <- expand.grid(ratio = ratios, beta = betas,
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
cat(sprintf("[INFO] Sweeping %d cells (method: %s, ratio: %.2f-%.2f, beta: %.2f-%.2f)\n",
            nrow(grid), map_method, min(ratios), max(ratios), min(betas), max(betas)))

eval_one <- function(R, b) {
  # suppressMessages: align_feature_set() and compute_flagship_abg() report per
  # call, which is useful once and unreadable across a full sweep.
  res <- suppressMessages(compute_flagship_abg(
    zs = zs, fdg = fdg_use %>% rename(weight = w_pet), sizes = size_df,
    beta = b, map_method = map_method, ratio = R, clamp_z = clamp_z,
    alpha = cfg$mbi$alpha %||% NULL,
    coverage_thr = coverage_thr,
    gamma = gamma, winsor_q = winsor_q
  ))
  scores <- as.data.frame(res$scores)

  corr_unw_pet <- suppressWarnings(cor(scores$gbi, scores$mbi_raw, use = "pairwise.complete.obs"))

  tmp <- inner_join(scores[, c("id", "mbi")], r_tilt_df, by = "id")
  mbi_map_r <- suppressWarnings(cor(tmp$mbi, tmp$r_tilt, use = "pairwise.complete.obs"))

  # Three numbers, three failure modes:
  #   mbi_map_r : does the index track the map alignment it claims to?
  #   n_eff     : or is it buying that alignment from a handful of regions?
  #   r2_raw_gbi: is there enough contrast left for the residual to be signal?
  data.frame(ratio = R, beta = b,
             corr_unw_pet   = corr_unw_pet,
             r2_raw_gbi     = corr_unw_pet^2,
             mbi_map_r      = mbi_map_r,
             n_eff          = res$concentration$n_eff,
             n_features     = res$concentration$n_features,
             top_share      = max(res$concentration$shares$share),
             realised_ratio = res$diagnostics$realised_ratio)
}

evals <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  if (i %% 20 == 0) cat(sprintf("  ... cell %d/%d\n", i, nrow(grid)))
  eval_one(grid$ratio[i], grid$beta[i])
}))
tab <- evals %>% arrange(desc(mbi_map_r))

out_file <- file.path(cal_dir, "calibration_grid_results.tsv")
write.table(tab, file = out_file, sep = "\t",
            row.names = FALSE, quote = FALSE)

cat("[OK] Wrote: ", out_file, "\n", sep = "")
cat("[OK] Best row by mbi_map_r, UNCONSTRAINED (do not adopt blindly):\n")
print(tab[1, , drop = FALSE])

feasible <- tab[is.finite(tab$n_eff) & tab$n_eff >= n_eff_min &
                is.finite(tab$r2_raw_gbi) &
                tab$r2_raw_gbi >= r2_band[1] & tab$r2_raw_gbi <= r2_band[2], ,
                drop = FALSE]

cat(sprintf("[OK] RECOMMENDED: best mbi_map_r subject to n_eff >= %d and R2(raw,GBI) in [%.2f, %.2f]:\n",
            n_eff_min, r2_band[1], r2_band[2]))
if (nrow(feasible)) {
  print(feasible[1, , drop = FALSE])
} else {
  cat("     NONE - no cell satisfies both constraints.\n")
  cat("     Widen the ratio grid, relax n_eff_min, or reconsider beta:\n")
  cat(sprintf("     best n_eff anywhere = %.2f; R2 range = [%.3f, %.3f]\n",
              max(tab$n_eff, na.rm = TRUE),
              min(tab$r2_raw_gbi, na.rm = TRUE), max(tab$r2_raw_gbi, na.rm = TRUE)))
}

cat("[OK] Concentration frontier (best mbi_map_r at each n_eff decile):\n")
tab$n_eff_bin <- cut(tab$n_eff, breaks = quantile(tab$n_eff, seq(0, 1, 0.1), na.rm = TRUE),
                     include.lowest = TRUE, labels = FALSE)
frontier <- do.call(rbind, lapply(split(tab, tab$n_eff_bin), function(d)
  d[which.max(d$mbi_map_r), c("ratio","beta","mbi_map_r","n_eff","r2_raw_gbi","top_share")]))
print(frontier, row.names = FALSE)
