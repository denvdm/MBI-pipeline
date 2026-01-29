#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(RNifti)
  library(dplyr)
})

source("helpers/run_helpers.R")

opt <- parse_common_args()

cfg <- read_config(opt$config)

# Paths
maps_dir <- cfg$paths$maps_dir %||% NULL
if (is.null(maps_dir)) stop("Missing config: paths: maps_dir", call. = FALSE)

pet_file   <- cfg$inputs$pet_file   %||% file.path(maps_dir, "source-castrillon2023_desc-cmrglc_space-MNI152_res-3mm_feature.nii.gz")
atlas_file <- cfg$inputs$atlas_file %||% file.path(maps_dir, "aparc+aseg_inPETspace_3mm.nii.gz")
gm_file    <- cfg$inputs$gm_file    %||% file.path(maps_dir, "GMprob_inPETspace_3mm.nii.gz")
lut_file   <- cfg$inputs$lut_file   %||% file.path(maps_dir, "lut_aparc_aseg.csv")

assert_file_exists(pet_file, "pet_file")
assert_file_exists(atlas_file, "atlas_file")
assert_file_exists(gm_file, "gm_file")
assert_file_exists(lut_file, "lut_file")

source_local("pet_parcel_extract.R")

out_stage <- stage_dir(opt$outdir, "extract_pet")
out_parc  <- file.path(out_stage, "pet_fdg_parcellated.csv")
out_wts   <- file.path(out_stage, "pet_fdg_featMedians.txt")

if (file.exists(out_wts) && !opt$overwrite) {
  stop("Output exists: ", out_wts, " (use --overwrite)", call. = FALSE)
}

lut <- read.csv(lut_file)

res <- save_pet_parcel_csv(
  pet = pet_file, atlas = atlas_file, gm_prob = gm_file,
  output_csv = out_parc,
  gm_thr = 0.20, use_median = TRUE, clamp_neg = TRUE, winsor_q = NULL,
  min_vox = 50L, min_cov = 0.30, lut = lut
)

fdg_raw <- res$parcels |>
  dplyr::rename(index = label, fdg_median = pet_value, fdg_mean = pet_mean,
                nvox_used = n_used, nvox_total = n_total, region = ROI_Name)

colnames(fdg_raw) <- c("index","fdg_median","fdg_mean","p05","p95","nvox_used","nvox_total","coverage","region")
fdg_raw$region <- gsub("^ctx-","",fdg_raw$region)
fdg_raw$region <- gsub("-","_",fdg_raw$region)
fdg_raw$region <- make.names(fdg_raw$region, unique = TRUE)

fdg_raw$fdg_z <- scale(fdg_raw$fdg_median)
fdg_raw$fdg_w <- fdg_raw$fdg_z - min(fdg_raw$fdg_z, na.rm = TRUE) + 1e-6
fdg_raw$fdg_w <- fdg_raw$fdg_w / mean(fdg_raw$fdg_w, na.rm = TRUE)

fdg <- fdg_raw |> dplyr::transmute(region, weight = fdg_w)

write.table(fdg, file = out_wts, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Wrote:\n- ", out_parc, "\n- ", out_wts, "\n")
