#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(broom)
  library(tidyr)
})

source("helpers/run_helpers.R")
opt <- parse_common_args()
cfg <- read_config(opt$config)

source_local("relate_pet_funcs.R")

normdev_long <- cfg$inputs$deviations_z_long %||% file.path(opt$outdir, "normdev", "deviations_z_long.txt")
sizes_file   <- cfg$inputs$feat_sizes        %||% file.path(opt$outdir, "normdev", "feats_avgSize.txt")
fdg_file     <- cfg$inputs$pet_weights       %||% file.path(opt$outdir, "extract_pet", "pet_fdg_featMedians.txt")

assert_file_exists(normdev_long, "deviations_z_long")
assert_file_exists(sizes_file,   "feat_sizes")
assert_file_exists(fdg_file,     "pet_weights")

out_stage <- stage_dir(opt$outdir, "qc_pet_align")

zs <- fread(normdev_long, data.table = FALSE)
size_df <- fread(sizes_file, data.table = FALSE)
fdg <- fread(fdg_file, data.table = FALSE)

mods_use <- cfg$qc$mods_use %||% c("CT","SA","SubVol")
subvol_gate_pct <- cfg$qc$subvol_gate_pct %||% 0.20

fdg_use <- prep_fdg(fdg)
zs_gated <- gate_subvol(zs, size_df, pct = subvol_gate_pct)
Zc <- collapse_zs(zs_gated, fdg_use, mods_use = mods_use)

mods_tag <- paste(mods_use, collapse = "+")

align_tab <- bind_rows(
  pc_pet_alignment(Zc, fdg_use, mods_tag),
  loo_pet_alignment(Zc, fdg_use)
)

write.table(align_tab, file = file.path(out_stage, "pet_alignment_summary.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

row_subj <- subject_pet_alignment(Zc, fdg_use)
write.table(row_subj, file = file.path(out_stage, "subject_alignment_summary.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

pc1_tbl <- pc1_loadings_table(Zc, fdg_use, align_to_pet = TRUE)
write.table(pc1_tbl, file = file.path(out_stage, "pc1_regional_z.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

meanZ <- mean_z_per_region(Zc)
write.table(meanZ, file = file.path(out_stage, "mean_regional_z.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

looB <- loo_beta_map(Zc)
write.table(looB, file = file.path(out_stage, "roi_beta_LOO.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("Wrote QC outputs to: ", out_stage, "\n")
