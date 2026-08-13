#!/usr/bin/env Rscript

# -------------------------------------------------------------------
# run_toy_example.R
#
# Self-contained smoke test for the MBI pipeline.
#
# Computes the three flagship indices (GBI, MBI_raw, MBI) on the
# bundled 10-subject toy cohort using the repository's real index
# code (helpers/compute_mbi.R). It needs only base R + dplyr (+ yaml
# if a config is present): no PET NIfTI maps, no FreeSurfer stats, no
# large cohort.
#
# Two stages cannot be run literally on a 10-subject toy, so they are
# replaced by transparent surrogates; the index computation itself is
# the unmodified production code path:
#   Stage 01 (PET extraction) -> a small, deterministic, strictly
#       positive per-region PET weight vector (no NIfTI required).
#   Stage 02 (normative GAM)  -> a simple per-feature linear normative
#       fit (value ~ age + sex [+ icv for SA/SubVol]); the real GAM
#       stage requires >= 50 subjects, which a toy cannot satisfy.
#   Stage 03 (compute MBI)    -> helpers/compute_mbi.R, unchanged.
#
# Usage (from the repository root, or anywhere):
#   Rscript examples/toy/run_toy_example.R
#
# Outputs are written to examples/toy/output/.
# -------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---- locate repo root from this script's own path (robust to CWD) ----
this_file <- local({
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) normalizePath(f[1]) else normalizePath("examples/toy/run_toy_example.R")
})
toy_dir   <- dirname(this_file)
repo_root <- normalizePath(file.path(toy_dir, "..", ".."))
cat("[toy] repo root:", repo_root, "\n")

# ---- read the toy config (lenient plain-YAML read, not the strict ----
# ---- generic loader; every key has a built-in default below) ---------
cfg_path <- file.path(toy_dir, "config_toy.yaml")
cfg <- if (requireNamespace("yaml", quietly = TRUE) && file.exists(cfg_path)) {
  yaml::read_yaml(cfg_path)
} else {
  list()
}

# resolve a possibly-relative path against the repo root
abspath <- function(p, default) {
  p <- p %||% default
  if (!grepl("^(/|[A-Za-z]:)", p)) file.path(repo_root, p) else p
}

features_file <- abspath(cfg$paths$features_wide, "examples/toy/input/toy_features_wide.tsv")
fmap_file     <- abspath(cfg$paths$feature_map,   "examples/toy/input/toy_feature_map.tsv")
outdir        <- abspath(cfg$paths$outdir,        "examples/toy/output")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# toy generic column names -> internal names
cm       <- cfg$columns %||% list()
col_id   <- cm$id   %||% "eid"
col_age  <- cm$age  %||% "age"
col_sex  <- cm$sex  %||% "sex"
col_site <- cm$site %||% "site"
col_icv  <- cm$icv  %||% "icv"

# toy generic modality labels -> CT/SA/SubVol used by compute_mbi.R
mm <- cfg$modality_map %||% list(thickness = "CT", area = "SA", volume = "SubVol")

# index parameters (flagship defaults)
mb              <- cfg$mbi %||% list()
beta            <- mb$beta            %||% 0.25
map_method      <- mb$map_method      %||% "robust_z"
ratio           <- mb$ratio           %||% 4
clamp_z         <- mb$clamp_z         %||% 2
coverage_thr    <- mb$coverage_thr    %||% 0.80
subvol_min_mm3  <- mb$subvol_min_mm3  %||% 0   # toy sizes are synthetic; 0 = keep all
icv_modalities  <- cfg$icv_modalities %||% c("SA", "SubVol")

# ---- load toy inputs --------------------------------------------------
fw <- read.delim(features_file, stringsAsFactors = FALSE, check.names = FALSE)
fm <- read.delim(fmap_file,     stringsAsFactors = FALSE, check.names = FALSE)

if (!all(c("feature_id", "modality", "region") %in% names(fm))) {
  stop("feature_map must have columns: feature_id, modality, region", call. = FALSE)
}
fm$modality_internal <- unlist(mm[fm$modality], use.names = FALSE)
if (anyNA(fm$modality_internal)) {
  bad <- unique(fm$modality[is.na(fm$modality_internal)])
  stop("Unmapped modality label(s): ", paste(bad, collapse = ", "),
       ". Add them to modality_map in config_toy.yaml.", call. = FALSE)
}

need <- c(col_id, col_age, col_sex, col_site)
miss <- setdiff(need, names(fw))
if (length(miss)) stop("features_wide missing columns: ", paste(miss, collapse = ", "), call. = FALSE)

fw$.id  <- as.character(fw[[col_id]])
fw$.age <- as.numeric(fw[[col_age]])
fw$.sex <- factor(fw[[col_sex]])
fw$.icv <- if (col_icv %in% names(fw)) as.numeric(fw[[col_icv]]) else NA_real_

feat_cols <- fm$feature_id
miss_f <- setdiff(feat_cols, names(fw))
if (length(miss_f)) stop("features_wide missing feature column(s): ", paste(miss_f, collapse = ", "), call. = FALSE)

mod_counts <- table(fm$modality_internal)
cat(sprintf("[toy] %d subjects, %d features (%s)\n",
            nrow(fw), length(feat_cols),
            paste(sprintf("%s:%d", names(mod_counts), as.integer(mod_counts)), collapse = ", ")))

# ---- Stage 02 surrogate: simple normative deviations ------------------
# Real stage 02 fits a GAM per feature (mgcv, REML, k=10) on >= 50
# subjects. A toy cohort is far too small, so we use a transparent
# linear normative fit per feature: the residual of value ~ age + sex
# (+ icv for SA/SubVol), standardised by the residual SD. Sign
# orientation (higher = healthier) is applied inside compute_mbi.R.
z_one <- function(feat) {
  mod_int <- fm$modality_internal[match(feat, fm$feature_id)]
  y  <- as.numeric(fw[[feat]])
  df <- data.frame(y = y, age = fw$.age, sex = fw$.sex)
  use_icv <- (mod_int %in% icv_modalities) && any(is.finite(fw$.icv))
  if (use_icv) df$icv <- fw$.icv
  frm <- if (use_icv) y ~ age + sex + icv else y ~ age + sex
  fit <- stats::lm(frm, data = df, na.action = stats::na.exclude)
  mu  <- stats::predict(fit, newdata = df)
  s   <- stats::sd(y - mu, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) s <- 1
  (y - mu) / s
}

zs <- do.call(rbind, lapply(feat_cols, function(f) {
  data.frame(
    id       = fw$.id,
    modality = fm$modality_internal[match(f, fm$feature_id)],
    region   = fm$region[match(f, fm$feature_id)],
    z        = z_one(f),
    stringsAsFactors = FALSE
  )
}))

# feature sizes = mean raw value per (modality, region)
sizes <- data.frame(
  modality = fm$modality_internal,
  region   = fm$region,
  size     = vapply(feat_cols, function(f) mean(as.numeric(fw[[f]]), na.rm = TRUE), numeric(1)),
  stringsAsFactors = FALSE
)

# ---- Stage 01 surrogate: synthetic per-region PET weights -------------
# The toy ships no PET NIfTI, so we assign a fixed, deterministic,
# strictly positive metabolic weight per region. This makes MBI_raw
# genuinely PET-weighted (hence distinct from GBI). Illustrative only;
# the values carry no biological meaning.
pet_default <- list(R01 = 1.60, R02 = 1.10, R03 = 0.70, R04 = 0.90)
pet_cfg     <- cfg$pet_weights %||% pet_default
regions     <- sort(unique(zs$region))
weight      <- vapply(regions, function(r) {
  v <- pet_cfg[[r]]
  if (is.null(v)) 1.0 else as.numeric(v)
}, numeric(1))
fdg <- data.frame(region = regions, weight = weight, stringsAsFactors = FALSE)

# ---- SubVol size gate (identical logic to scripts/03_run_mbi.R) -------
# Absolute threshold, matching the main pipeline. The toy default is 0
# (keep everything): its sizes are synthetic means with no mm^3 meaning,
# so the production default of 1000 would gate on noise.
keep_sub <- sizes %>%
  dplyr::filter(modality == "SubVol") %>%
  dplyr::filter(is.finite(size), size >= subvol_min_mm3) %>%
  dplyr::pull(region)
zs <- zs %>% dplyr::filter(modality != "SubVol" | region %in% keep_sub)

# ---- Stage 03: REAL index computation ---------------------------------
source(file.path(repo_root, "helpers", "compute_mbi.R"))
res <- compute_flagship_abg(
  zs = zs, fdg = fdg, sizes = sizes,
  beta = beta, map_method = map_method, ratio = ratio, clamp_z = clamp_z,
  alpha = NULL, coverage_thr = coverage_thr
)
scores <- as.data.frame(res$scores)

# ---- write outputs ----------------------------------------------------
w_tsv <- function(x, f) write.table(x, file.path(outdir, f), sep = "\t",
                                    quote = FALSE, row.names = FALSE)
w_tsv(fdg,    "pet_weights_toy.txt")
w_tsv(zs,     "deviations_z_long.txt")
w_tsv(sizes,  "feats_avgSize.txt")
w_tsv(scores, "mbi_flagship_scores.txt")

# ---- QC summary + smoke-test assertions -------------------------------
r_mbi_gbi    <- suppressWarnings(stats::cor(scores$mbi,     scores$gbi))
r_mbiraw_gbi <- suppressWarnings(stats::cor(scores$mbi_raw, scores$gbi))
qc <- data.frame(
  metric = c("n_subjects", "n_regions_pet", "subvol_regions_kept",
             "mean_mbi", "sd_mbi", "cor_mbi_gbi", "cor_mbiraw_gbi"),
  value  = c(nrow(scores), nrow(fdg), length(keep_sub),
             round(mean(scores$mbi), 6), round(stats::sd(scores$mbi), 6),
             round(r_mbi_gbi, 6), round(r_mbiraw_gbi, 6))
)
w_tsv(qc, "toy_qc_summary.txt")

cat("\n[toy] scores (head):\n"); print(utils::head(scores))
cat("\n[toy] QC summary:\n");    print(qc, row.names = FALSE)

ok <- TRUE
chk <- function(cond, msg) {
  if (!isTRUE(cond)) { cat("  FAIL: ", msg, "\n", sep = ""); ok <<- FALSE }
  else               { cat("  ok  : ", msg, "\n", sep = "") }
}
cat("\n[toy] smoke-test checks:\n")
chk(nrow(scores) == nrow(fw),                              "one score row per subject")
chk(all(c("id", "mbi_raw", "gbi", "mbi") %in% names(scores)), "expected columns present")
chk(all(is.finite(scores$mbi)),                            "MBI finite for all subjects")
chk(abs(mean(scores$mbi)) < 1e-6,                          "MBI mean ~ 0 (standardised)")
chk(abs(stats::sd(scores$mbi) - 1) < 1e-6,                 "MBI sd ~ 1 (standardised)")
chk(is.finite(r_mbi_gbi) && abs(r_mbi_gbi) < 1e-6,         "MBI orthogonal to GBI (residualised)")
chk(all(res$mbi_details$weight_table$w_pet > 0),           "no region zero-weighted in MBI_raw")
chk(nrow(res$mbi_details$weight_table) ==
    nrow(res$gbi_details$weight_table),                    "MBI and GBI share one feature set")

cat("\n[toy] outputs written to: ", outdir, "\n", sep = "")
cat("[toy] R ", as.character(getRversion()), "\n", sep = "")
if (ok) {
  cat("[toy] SMOKE TEST PASSED\n")
} else {
  cat("[toy] SMOKE TEST FAILED\n")
  quit(status = 1L)
}
