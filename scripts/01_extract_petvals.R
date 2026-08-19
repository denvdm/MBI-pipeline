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
# Rename BY NAME. A positional colnames()<- here would silently relabel columns
# by position, so any change to the column order of pet_parcel_extract()'s
# output would scramble every downstream value without erroring.

# Normalise region names to match the adapter's canonical form:
#   - strip "ctx-" prefix from cortical
#   - convert "-" to "_"
#   - rename "Left_*" / "Right_*" subcortical to "lh_*" / "rh_*"
fdg_raw$region <- gsub("^ctx-",   "",    fdg_raw$region)
fdg_raw$region <- gsub("-",       "_",   fdg_raw$region)
fdg_raw$region <- gsub("^Left_",  "lh_", fdg_raw$region)
fdg_raw$region <- gsub("^Right_", "rh_", fdg_raw$region)
fdg_raw$region <- make.names(fdg_raw$region)
# Note the absence of unique = TRUE: that would rename a duplicate to
# "<name>.1", a region name matching nothing downstream, and the parcel would
# vanish from the index without a word. A duplicate here means the LUT or the
# name normalisation above is broken, which is an error, not something to
# rename around.
if (any(duplicated(fdg_raw$region)))
  stop("Duplicated region names after normalisation: ",
       paste(unique(fdg_raw$region[duplicated(fdg_raw$region)]), collapse = ", "),
       call. = FALSE)

# ---------------------------------------------------------------------------
# Emit the RAW parcel statistic. No z-scoring, no anchoring, no normalisation:
# the map -> weight transform lives entirely in
# helpers/compute_mbi.R::transform_map_weights(), and this stage must not
# duplicate any part of it. Rescaling here would make the weights depend on the
# map's distributional shape and let a single outlier parcel move every weight,
# invisibly to stage 03.
#
# QC columns travel with the map so downstream validity filtering can use the
# extraction criteria rather than inferring validity from the value.
# ---------------------------------------------------------------------------
# Within-parcel PET dispersion: (p95 - p05) / median, over the same GM-masked
# voxels as the parcel median. This asks whether a single number is a
# meaningful summary of the parcel at all - a large, heterogeneous parcel
# (cerebellar cortex spans ~52,000 mm^3 at 3 mm resolution; brainstem lumps
# midbrain, pons and medulla; VentralDC is an anatomical grab-bag) can have a
# perfectly valid median that nonetheless represents nothing in particular.
#
# It is outcome-blind, map-side, and generalises to any biological prior, which
# makes it a principled basis for excluding a region rather than an aesthetic
# one.
fdg <- fdg_raw |>
  dplyr::transmute(region,
                   weight = fdg_median,
                   coverage, nvox_used, nvox_total,
                   p05, p95,
                   dispersion = ifelse(is.finite(fdg_median) & fdg_median > 0,
                                       (p95 - p05) / fdg_median, NA_real_))

write_tsv(fdg, out_wts)
cat("Wrote:\n- ", out_parc, "\n- ", out_wts, "\n")

write_run_info(out_stage, "01_extract_petvals", opt$config, params = list(
  pet_file = basename(pet_file), atlas_file = basename(atlas_file),
  gm_file = basename(gm_file), lut_file = basename(lut_file),
  gm_thr = 0.20, use_median = TRUE, min_vox = 50L, min_cov = 0.30,
  n_regions = nrow(fdg),
  n_regions_valid = sum(is.finite(fdg$weight) & fdg$weight > 0)
))

# Surface the dispersion ranking in the stage-01 log, so the decision about
# which parcels to carry into stage 03 can be made before stage 03 is run.
cat("\nWithin-parcel PET dispersion (p95-p05)/median, highest first:\n")
rank_tbl <- fdg[order(-fdg$dispersion), c("region","weight","dispersion","coverage","nvox_used")]
print(utils::head(as.data.frame(rank_tbl), 20), row.names = FALSE)
cat("\n(High dispersion = the parcel median is a poor summary of the parcel.\n",
    " Decide exclusions from this ranking, not from region identity.)\n", sep = "")
