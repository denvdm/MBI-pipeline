#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(tidyr)
})

source("helpers/run_helpers.R")
opt <- parse_common_args()
cfg <- read_config(opt$config)

source_local("compute_mbi.R")        # for compute_oriented(): one orientation rule
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

# subvol_min_mm3 is sourced from cfg$mbi (single source of truth, consumed
# by both stage 03 and this stage so the QC sees the same gate that MBI did).
mods_use       <- cfg$qc$mods_use        %||% c("CT","SA","SubVol")
subvol_min_mm3 <- cfg$mbi$subvol_min_mm3 %||% 1000

if (!is.null(cfg$mbi$subvol_gate_pct))
  warning("config: 'subvol_gate_pct' is not read by this pipeline and is ignored; ",
          "using absolute 'subvol_min_mm3' = ", subvol_min_mm3, call. = FALSE)

fdg_use <- prep_fdg(fdg)

# Orient before anything else, with the same rule the index uses. Averaging a
# raw `z` across modalities inside collapse_zs() would mix sign conventions the
# moment a negatively-oriented modality (MD, RD, WMH) is in mods_use.
zs_oriented <- compute_oriented(zs)
n_dropped   <- dplyr::n_distinct(zs$modality) - dplyr::n_distinct(zs_oriented$modality)
if (n_dropped > 0)
  message(sprintf("[qc] %d modality/modalities have no defined orientation and were dropped: %s",
                  n_dropped,
                  paste(setdiff(unique(zs$modality), unique(zs_oriented$modality)),
                        collapse = ", ")))

zs_gated <- gate_subvol(zs_oriented, size_df, min_mm3 = subvol_min_mm3)
Zc <- collapse_zs(zs_gated, fdg_use, mods_use = mods_use)

mods_tag <- paste(mods_use, collapse = "+")

# The PC1 decomposition and the leave-one-out slopes each feed two outputs.
# Compute each once and pass it in, rather than letting every consumer rebuild
# it from Zc.
pc1   <- region_pc1(Zc)
betas <- region_loo_betas(Zc)

align_tab <- bind_rows(
  pc_pet_alignment(Zc, fdg_use, mods_tag, pc1 = pc1),
  loo_pet_alignment(Zc, fdg_use, betas = betas)
)
write_tsv(align_tab, file.path(out_stage, "pet_alignment_summary.txt"))

row_subj <- subject_pet_alignment(Zc, fdg_use)
write_tsv(row_subj, file.path(out_stage, "subject_alignment_summary.txt"))

pc1_tbl <- pc1_loadings_table(Zc, fdg_use, align_to_pet = TRUE, pc1 = pc1)
write_tsv(pc1_tbl, file.path(out_stage, "pc1_regional_z.txt"))

meanZ <- mean_z_per_region(Zc)
write_tsv(meanZ, file.path(out_stage, "mean_regional_z.txt"))

looB <- loo_beta_map(Zc, betas = betas)
write_tsv(looB, file.path(out_stage, "roi_beta_LOO.tsv"))

cat("Wrote QC outputs to: ", out_stage, "\n")

write_run_info(out_stage, "04_relate_mbi_pet", opt$config, params = list(
  mods_use = mods_use, subvol_min_mm3 = subvol_min_mm3,
  n_regions = length(pc1$region),
  n_subjects = row_subj$n_subjects,
  var_explained_pc1 = round(pc1$var_explained, 6)
))
