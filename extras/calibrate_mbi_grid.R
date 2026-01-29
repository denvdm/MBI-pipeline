#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(data.table)
  library(dplyr)
  library(tidyr)
})

# ------------------------------------------------------------
# USER: EDIT HERE (calibration grid)
# ------------------------------------------------------------
# Grid-search hyperparameters for the PET/size weighting.
# Keep modest unless you have serious compute.
gammas <- c(2.2, 2.3, 2.4, 2.5)
betas  <- c(0.30, 0.25, 0.20)

# Modalities used for r_tilt target (alignment of deviations with PET weights)
mods_use <- c("CT", "SA")

# Gate smallest SubVol regions (quantile cutoff)
subvol_gate_pct <- 0.20
# ------------------------------------------------------------


# -----------------------------
# CLI
# -----------------------------
option_list <- list(
  make_option(c("--config"), type = "character", default = NULL,
              help = "Path to YAML config used by the main pipeline (required)."),
  make_option(c("--outdir"), type = "character", default = NULL,
              help = "Pipeline output directory (required)."),
  make_option(c("--overwrite"), action = "store_true", default = FALSE,
              help = "Overwrite outputs if they exist.")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$config) || is.null(opt$outdir)) {
  cat("\n[ERROR] --config and --outdir are required.\n\n", file = stderr())
  quit(status = 1)
}

cfg <- yaml::read_yaml(opt$config)

# -----------------------------
# Resolve repository root robustly (so sourcing works from anywhere)
# -----------------------------
get_script_dir <- function() {
  arg <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", arg[grepl("^--file=", arg)])
  if (length(f) == 0) return(getwd())
  dirname(normalizePath(f))
}
repo_root_path <- normalizePath(file.path(get_script_dir(), ".."))

# Source helpers
source(file.path(repo_root_path, "helpers", "run_helpers.R"))
source(file.path(repo_root_path, "helpers", "compute_mbi.R"))

# -----------------------------
# Local utilities (keep explicit, do not rely on run_helpers internals)
# -----------------------------
assert_file_exists <- function(path, label = NULL) {
  if (!file.exists(path)) {
    msg <- paste0("[ERROR] Missing file",
                  if (!is.null(label)) paste0(" (", label, ")") else "",
                  ": ", path, "\n")
    stop(msg, call. = FALSE)
  }
}
stage_dir_local <- function(outdir, stage) {
  dir <- file.path(outdir, stage)
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  dir
}

# -----------------------------
# Resolve inputs from staged pipeline outputs (CURRENT naming)
# -----------------------------
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

# Validate expected columns (current pipeline)
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
# SubVol gating (match main pipeline default)
# -----------------------------
keep_sub <- size_df %>%
  filter(modality == "SubVol",
         is.finite(size),
         size >= quantile(size, subvol_gate_pct, na.rm = TRUE)) %>%
  pull(region)

zs <- zs %>% filter(modality != "SubVol" | region %in% keep_sub)

# -----------------------------
# Compute r_tilt using the SAME orientation logic as MBI
# -----------------------------
# compute_oriented() is defined in helpers/compute_mbi.R
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
# Grid evaluation (CURRENT score names: id, mbi_raw, gbi, mbi)
# -----------------------------
grid <- expand.grid(gamma = gammas, beta = betas, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

eval_one <- function(g, b) {
  res <- compute_flagship_abg(
    zs = zs, fdg = fdg_use %>% rename(weight = w_pet), sizes = size_df,
    beta = b, gamma = g,
    alpha = cfg$mbi$alpha %||% NULL,
    cap_frac = cfg$mbi$cap_frac %||% 0.05,
    coverage_thr = cfg$mbi$coverage_thr %||% 0.80
  )
  scores <- as.data.frame(res$scores)

  # corr between size-weighted (gbi) and PET-weighted (mbi_raw)
  corr_unw_pet <- suppressWarnings(cor(scores$gbi, scores$mbi_raw, use = "pairwise.complete.obs"))

  # alignment between residualized MBI (mbi) and r_tilt
  tmp <- inner_join(scores[, c("id", "mbi")], r_tilt_df, by = "id")
  mbi_map_r <- suppressWarnings(cor(tmp$mbi, tmp$r_tilt, use = "pairwise.complete.obs"))

  data.frame(gamma = g, beta = b, corr_unw_pet = corr_unw_pet, mbi_map_r = mbi_map_r)
}

evals <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) eval_one(grid$gamma[i], grid$beta[i])))
tab <- evals %>% arrange(desc(mbi_map_r))

out_file <- file.path(cal_dir, "calibration_grid_results.tsv")
if (file.exists(out_file) && !opt$overwrite) {
  stop("[ERROR] Output exists: ", out_file, " (use --overwrite)", call. = FALSE)
}
write.table(tab, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE)

cat("[OK] Wrote: ", out_file, "\n", sep = "")
cat("[OK] Best row (by mbi_map_r):\n")
print(tab[1, , drop = FALSE])
