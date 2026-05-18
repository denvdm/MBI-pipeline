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
gammas <- seq(1.0, 4.0, by = 0.25)   # 13 values
betas  <- seq(0.0, 0.5, by = 0.05)   # 11 values

# Modalities used for the r_tilt target (alignment of deviations with PET weights)
mods_use <- c("CT", "SA")

# SubVol size gate (passed through; the actual gate value comes from config)
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

subvol_gate_pct <- cfg$mbi$subvol_gate_pct %||% 0.20
winsor_q        <- cfg$mbi$winsor_q        %||% 0.99
coverage_thr    <- cfg$mbi$coverage_thr    %||% 0.80

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
# SubVol gating (matches main pipeline). Two-step filter so quantile()
# is computed within SubVol only.
# -----------------------------
keep_sub <- size_df %>%
  filter(modality == "SubVol") %>%
  filter(is.finite(size), size >= quantile(size, subvol_gate_pct, na.rm = TRUE)) %>%
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
grid <- expand.grid(gamma = gammas, beta = betas,
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
cat(sprintf("[INFO] Sweeping %d cells (gamma: %.2f-%.2f, beta: %.2f-%.2f)\n",
            nrow(grid), min(gammas), max(gammas), min(betas), max(betas)))

eval_one <- function(g, b) {
  res <- compute_flagship_abg(
    zs = zs, fdg = fdg_use %>% rename(weight = w_pet), sizes = size_df,
    beta = b, gamma = g,
    alpha = cfg$mbi$alpha %||% NULL,
    winsor_q     = winsor_q,
    coverage_thr = coverage_thr
  )
  scores <- as.data.frame(res$scores)

  corr_unw_pet <- suppressWarnings(cor(scores$gbi, scores$mbi_raw, use = "pairwise.complete.obs"))

  tmp <- inner_join(scores[, c("id", "mbi")], r_tilt_df, by = "id")
  mbi_map_r <- suppressWarnings(cor(tmp$mbi, tmp$r_tilt, use = "pairwise.complete.obs"))

  data.frame(gamma = g, beta = b, corr_unw_pet = corr_unw_pet, mbi_map_r = mbi_map_r)
}

evals <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  if (i %% 20 == 0) cat(sprintf("  ... cell %d/%d\n", i, nrow(grid)))
  eval_one(grid$gamma[i], grid$beta[i])
}))
tab <- evals %>% arrange(desc(mbi_map_r))

out_file <- file.path(cal_dir, "calibration_grid_results.tsv")
write.table(tab, file = out_file, sep = "\t",
            row.names = FALSE, quote = FALSE)

cat("[OK] Wrote: ", out_file, "\n", sep = "")
cat("[OK] Best row (by mbi_map_r):\n")
print(tab[1, , drop = FALSE])
cat("[OK] Current settings (gamma=2.5, beta=0.25):\n")
print(tab[tab$gamma == 2.5 & tab$beta == 0.25, , drop = FALSE])
