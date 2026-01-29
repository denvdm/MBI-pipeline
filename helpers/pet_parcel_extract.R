# pet_parcel_extract.R — PET→parcel extraction (GM mask + median) with robust QC
# Dependencies: RNifti (for NIfTI I/O). Uses base R for speed and minimal deps.
# Exposes:
#   pet_parcel_extract(pet, atlas, gm_prob=NULL, gm_thr=0.20, use_median=TRUE,
#                      clamp_neg=TRUE, winsor_q=NULL, min_vox=50L, min_cov=0.30, lut=NULL)
#   save_pet_parcel_csv(pet, atlas, gm_prob=NULL, output_csv="roi_values_with_names.csv", ...)
#
# Notes:
# - Expect atlas and (optional) GM map already resampled to PET grid (3mm) with NN for labels.
# - If dims differ, the function stops with an informative error.
# - Pass a lookup table (lut) with columns c("label","name") to attach ROI names.

pet_parcel_extract <- function(pet, atlas, gm_prob=NULL,
                               gm_thr=0.20, use_median=TRUE,
                               clamp_neg=TRUE, winsor_q=NULL,
                               min_vox=50L, min_cov=0.30, lut=NULL) {
  if (!requireNamespace("RNifti", quietly=TRUE)) stop("Install RNifti: install.packages('RNifti')")
  rd <- function(x) {
    if (inherits(x, "niftiImage")) return(x)
    if (is.character(x) && file.exists(x)) return(RNifti::readNifti(x))
    stop("Provide file paths (.nii/.nii.gz) or niftiImage objects for 'pet'/'atlas'/'gm_prob'.")
  }
  pet_img <- rd(pet); atlas_img <- rd(atlas); gm_img <- if (!is.null(gm_prob)) rd(gm_prob) else NULL

  # --- basic checks
  if (!all(dim(pet_img) == dim(atlas_img))) stop("PET and atlas dims differ; resample atlas to PET grid (NN) before running.")
  if (!is.null(gm_img) && !all(dim(pet_img) == dim(gm_img))) stop("PET and GM-prob dims differ; resample GM to PET grid.")

  pet_arr <- as.array(pet_img); atlas_arr <- as.array(atlas_img)
  if (!is.null(gm_img)) gm_arr <- as.array(gm_img)

  # --- mask
  if (!is.null(gm_img)) {
    gm_mask <- is.finite(pet_arr) & is.finite(atlas_arr) & (atlas_arr > 0) & (gm_arr >= gm_thr)
  } else {
    gm_mask <- is.finite(pet_arr) & is.finite(atlas_arr) & (atlas_arr > 0)
  }
  if (!any(gm_mask)) stop("No voxels passed the mask; check atlas/GM threshold.")

  # --- clamp/winsor
  if (isTRUE(clamp_neg)) pet_arr[pet_arr < 0] <- 0
  if (!is.null(winsor_q)) {
    x <- pet_arr[gm_mask]
    ql <- as.numeric(stats::quantile(x, probs=winsor_q, na.rm=TRUE, names=FALSE))
    qh <- as.numeric(stats::quantile(x, probs=1-winsor_q, na.rm=TRUE, names=FALSE))
    pet_arr[pet_arr < ql] <- ql; pet_arr[pet_arr > qh] <- qh
  }

  # --- voxel vectors (only GM∩atlas)
  lab <- as.integer(atlas_arr[gm_mask]); val <- as.numeric(pet_arr[gm_mask])

  # --- labels present and counts
  labs_all <- sort(unique(as.integer(atlas_arr[atlas_arr > 0])))
  labs_use <- sort(unique(lab))
  # total voxels per label in atlas (denominator for coverage)
  total_tab <- tabulate(match(as.integer(atlas_arr[atlas_arr > 0]), labs_all), nbins=length(labs_all))
  names(total_tab) <- labs_all
  # used voxels per label after mask
  used_tab  <- tabulate(match(lab, labs_all), nbins=length(labs_all))
  names(used_tab) <- labs_all

  # --- central tendency per label
  split_idx <- split(seq_along(val), lab)
  agg_fun <- if (isTRUE(use_median)) stats::median else base::mean
  pet_central <- vapply(split_idx, function(idx) agg_fun(val[idx], na.rm=TRUE), numeric(1))
  pet_mean    <- vapply(split_idx, function(idx) base::mean(val[idx], na.rm=TRUE), numeric(1))
  p05         <- vapply(split_idx, function(idx) stats::quantile(val[idx], 0.05, na.rm=TRUE, names=FALSE), numeric(1))
  p95         <- vapply(split_idx, function(idx) stats::quantile(val[idx], 0.95, na.rm=TRUE, names=FALSE), numeric(1))

  lbls <- as.integer(names(pet_central))
  n_used <- used_tab[match(lbls, as.integer(names(used_tab)))]
  n_total <- total_tab[match(lbls, as.integer(names(total_tab)))]
  covg <- n_used / pmax(1L, n_total)

  out <- data.frame(
    label = lbls,
    pet_value = as.numeric(pet_central),
    pet_mean = as.numeric(pet_mean),
    p05 = as.numeric(p05),
    p95 = as.numeric(p95),
    n_used = as.integer(n_used),
    n_total = as.integer(n_total),
    coverage = as.numeric(covg),
    stringsAsFactors = FALSE
  )

  # --- QC: too few voxels or low coverage → NA pet_value
  bad <- (out$n_used < min_vox) | (out$coverage < min_cov) | !is.finite(out$pet_value)
  out$pet_value[bad] <- NA_real_

  # --- optional ROI names
  if (!is.null(lut)) {
    if (!all(c("label","name") %in% names(lut))) stop("lut must have columns: label, name")
    out <- merge(out, lut[, c("label","name")], by="label", all.x=TRUE, sort=FALSE)
    names(out)[names(out)=="name"] <- "ROI_Name"
  }

  # --- sort by label
  out <- out[order(out$label), ]
  rownames(out) <- NULL

  list(parcels=out,
       settings=list(gm_thr=gm_thr, use_median=use_median, clamp_neg=clamp_neg,
                     winsor_q=winsor_q, min_vox=min_vox, min_cov=min_cov))
}

save_pet_parcel_csv <- function(pet, atlas, gm_prob=NULL, output_csv="roi_values_with_names.csv", ...) {
  res <- pet_parcel_extract(pet=pet, atlas=atlas, gm_prob=gm_prob, ...)
  utils::write.csv(res$parcels, file=output_csv, row.names=FALSE)
  message("Wrote: ", normalizePath(output_csv, winslash = "/"))
  invisible(res)
}
