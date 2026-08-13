#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})

source("helpers/run_helpers.R")
opt <- parse_common_args()
cfg <- read_config(opt$config)

source_local("compute_mbi.R")

# in/out
normdev_long <- cfg$inputs$deviations_z_long %||% file.path(opt$outdir, "normdev", "deviations_z_long.txt")
sizes_file   <- cfg$inputs$feat_sizes        %||% file.path(opt$outdir, "normdev", "feats_avgSize.txt")
fdg_file     <- cfg$inputs$pet_weights       %||% file.path(opt$outdir, "extract_pet", "pet_fdg_featMedians.txt")

assert_file_exists(normdev_long, "deviations_z_long")
assert_file_exists(sizes_file,   "feat_sizes (feats_avgSize.txt)")
assert_file_exists(fdg_file,     "pet_weights (pet_fdg_featMedians.txt)")

out_stage    <- stage_dir(opt$outdir, "compute_mbi")
out_scores   <- file.path(out_stage, "mbi_flagship_scores.txt")
out_features <- file.path(out_stage, "mbi_feature_manifest.txt")
out_shares   <- file.path(out_stage, "mbi_influence_shares.txt")

if (file.exists(out_scores) && !opt$overwrite) {
  stop("Output exists: ", out_scores, " (use --overwrite)", call. = FALSE)
}

zs      <- fread(normdev_long, data.table = FALSE) %>%
  mutate(region = as.character(region), modality = as.character(modality))
size_df <- fread(sizes_file, data.table = FALSE) %>%
  mutate(region = as.character(region), modality = as.character(modality))
fdg     <- fread(fdg_file, data.table = FALSE) %>%
  mutate(region = as.character(region))

# Parameters, overridable via config
beta         <- cfg$mbi$beta         %||% 0.25
map_method   <- cfg$mbi$map_method   %||% "robust_z"
ratio        <- cfg$mbi$ratio        %||% 4
clamp_z      <- cfg$mbi$clamp_z      %||% 2
alpha        <- cfg$mbi$alpha        %||% NULL
coverage_thr <- cfg$mbi$coverage_thr %||% 0.80

# legacy path only (map_method = "legacy"): reproduces the pre-rebuild build
gamma        <- cfg$mbi$gamma        %||% 2.5
winsor_q     <- cfg$mbi$winsor_q     %||% 0.99

if (!identical(map_method, "legacy") &&
    (!is.null(cfg$mbi$gamma) || !is.null(cfg$mbi$winsor_q)))
  message("note: gamma / winsor_q are ignored unless map_method = 'legacy'. ",
          "Contrast is set by 'ratio' (max weight ratio, high vs low demand).")

# ---------------------------------------------------------------------------
# SubVol reliability gate.
#
# ABSOLUTE volume threshold in mm^3, replacing the previous 20th-percentile
# rule. The quantile version was unstable: it was recomputed over whatever
# SubVol set survived upstream filtering, so removing a region for unrelated
# reasons (e.g. no valid map weight) silently shifted which OTHER structures
# were gated out. An absolute cut is invariant to that, which means this gate
# and the map-validity filter in align_feature_set() commute and the order of
# the two no longer needs defending.
#
# Rationale: very small aseg structures are segmented least reliably, so their
# deviations are dominated by measurement noise. Size is a proxy for
# reliability; the principled replacement is a per-feature test-retest ICC gate
# (UKB repeat imaging visit), applied across ALL modalities rather than SubVol
# alone. Config hook left below for that.
#
# !! SET FROM THE OBSERVED VOLUME DISTRIBUTION - the value below is a
# !! set to 1000 mm^3: just below amygdala (1590), the smallest structure
# !! surviving map-validity filtering, and close to the old effective cut (752).
# ---------------------------------------------------------------------------
subvol_min_mm3 <- cfg$mbi$subvol_min_mm3 %||% 1000

if (!is.null(cfg$mbi$subvol_gate_pct))
  warning("config: 'subvol_gate_pct' is deprecated and ignored; ",
          "using absolute 'subvol_min_mm3' = ", subvol_min_mm3, call. = FALSE)

# Reserved for the ICC-based gate (all modalities). NULL = not in use.
icc_file    <- cfg$inputs$feat_icc  %||% NULL
icc_min     <- cfg$mbi$icc_min      %||% NULL

sub_sizes <- size_df %>% dplyr::filter(modality == "SubVol")
keep_sub  <- sub_sizes %>% dplyr::filter(size >= subvol_min_mm3) %>% dplyr::pull(region)
drop_sub  <- setdiff(sub_sizes$region, keep_sub)

message(sprintf("SubVol gate: keeping %d/%d structures at >= %.0f mm^3%s",
                length(keep_sub), nrow(sub_sizes), subvol_min_mm3,
                if (length(drop_sub))
                  paste0(" (dropped: ", paste(sort(drop_sub), collapse = ", "), ")") else ""))

zs <- zs %>% dplyr::filter(modality != "SubVol" | region %in% keep_sub)

if (!is.null(icc_file) && !is.null(icc_min)) {
  assert_file_exists(icc_file, "feat_icc")
  icc <- fread(icc_file, data.table = FALSE) %>%
    mutate(region = as.character(region), modality = as.character(modality))
  keep_icc <- icc %>% dplyr::filter(is.finite(icc), icc >= icc_min) %>%
    dplyr::select(region, modality)
  n_before <- nrow(dplyr::distinct(zs, region, modality))
  zs <- zs %>% dplyr::semi_join(keep_icc, by = c("region","modality"))
  message(sprintf("ICC gate: %d -> %d features at ICC >= %.2f",
                  n_before, nrow(dplyr::distinct(zs, region, modality)), icc_min))
}

# ---------------------------------------------------------------------------
res <- compute_flagship_abg(
  zs = zs, fdg = fdg, sizes = size_df,
  beta = beta, map_method = map_method, ratio = ratio, clamp_z = clamp_z,
  alpha = alpha, coverage_thr = coverage_thr,
  gamma = gamma, winsor_q = winsor_q
)

mbi <- as.data.frame(res$scores)
write.table(mbi, file = out_scores, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Wrote: ", out_scores, "\n")

# Feature manifest: exactly what entered both indices, and what was dropped.
manifest <- rbind(
  data.frame(region = res$features$kept,    status = "kept",    stringsAsFactors = FALSE),
  data.frame(region = res$features$dropped, status = "dropped", stringsAsFactors = FALSE)
)
write.table(manifest, file = out_features, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Wrote: ", out_features, "\n")

# Per-feature influence shares + N_eff, for QC and for the methods table.
write.table(as.data.frame(res$concentration$shares), file = out_shares,
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Wrote: ", out_shares, "\n")
cat(sprintf("N_eff = %.2f over %d features (b = %.3f)\n",
            res$concentration$n_eff, res$concentration$n_features, res$concentration$b))
cat(sprintf("map_method = %s | ratio = %.2f (realised %.2f) | beta = %.2f | R2(MBI_raw, GBI) = %.3f\n",
            res$diagnostics$map_method, res$diagnostics$ratio,
            res$diagnostics$realised_ratio, res$diagnostics$beta,
            res$diagnostics$r2_raw_gbi))
