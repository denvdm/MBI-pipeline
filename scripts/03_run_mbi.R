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

out_stage  <- stage_dir(opt$outdir, "compute_mbi")
out_scores <- file.path(out_stage, "mbi_flagship_scores.txt")

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
beta            <- cfg$mbi$beta            %||% 0.25
gamma           <- cfg$mbi$gamma           %||% 2.5
alpha           <- cfg$mbi$alpha           %||% NULL
winsor_q        <- cfg$mbi$winsor_q        %||% 0.99
coverage_thr    <- cfg$mbi$coverage_thr    %||% 0.80
subvol_gate_pct <- cfg$mbi$subvol_gate_pct %||% 0.20

# SubVol gate: drop the smallest p<gate> regions within SubVol.
# Two-step filter ensures the quantile is computed within SubVol only
# (not over the full size_df, which would mix CT/SA/SubVol scales).
keep_sub <- size_df %>%
  dplyr::filter(modality == "SubVol") %>%
  dplyr::filter(size >= quantile(size, subvol_gate_pct, na.rm = TRUE)) %>%
  dplyr::pull(region)
zs <- zs %>% dplyr::filter(modality != "SubVol" | region %in% keep_sub)

res <- compute_flagship_abg(
  zs = zs, fdg = fdg, sizes = size_df,
  beta = beta, gamma = gamma, alpha = alpha,
  winsor_q = winsor_q, coverage_thr = coverage_thr
)

mbi <- as.data.frame(res$scores)
write.table(mbi, file = out_scores, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Wrote: ", out_scores, "\n")
