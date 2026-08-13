#!/usr/bin/env Rscript
# Compare canonical and alt PET-template MBI outputs (Tier 0 + Tier 1).
#
# Tier 0: weight-vector correlation across regions, raw and after the
#         pipeline's gamma + p99-winsor transform.
# Tier 1: subject-level construct correlations for MBI_raw, GBI, and the
#         GBI-residualised MBI.
#
# Ventricle regions are dropped from Tier 0 by default to match the
# main pipeline (which excludes ventricles upstream at prep_features.R).
# Pass drop_ventricles = FALSE to keep them, e.g. for a template-level
# sanity check on PVE handling.
#
# REQUIRES compute_mbi.R to be sourced first (for transform_map_weights()).

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

VENT_PAT <- "(?i)ventricle|ventricular|lat_vent|Vent"

# Fisher z confidence intervals for a single correlation
fisher_ci <- function(r, n, conf = 0.95) {
  if (!is.finite(r) || n <= 3) return(c(lo = NA_real_, hi = NA_real_))
  r  <- max(min(r, 0.999999), -0.999999)
  z  <- atanh(r)
  se <- 1 / sqrt(n - 3)
  zc <- qnorm(1 - (1 - conf) / 2)
  c(lo = tanh(z - zc * se), hi = tanh(z + zc * se))
}

# The PET-weight transform is defined ONCE, in compute_mbi.R. This file used to
# carry its own copy with a comment saying it "mirrors" the index - two code
# paths that must agree, with nothing enforcing it. Aliasing means any change to
# the transform reaches the cross-template comparison automatically, so SF3 can
# never silently validate a construction the index no longer uses.
if (!exists("transform_map_weights", mode = "function")) {
  stop("compare_pet_funcs.R requires compute_mbi.R to be sourced first ",
       "(transform_map_weights() not found).", call. = FALSE)
}
apply_pet_transform <- function(w, method = "robust_z", ratio = 4, clamp_z = 2,
                                gamma = 2.5, winsor_q = 0.99) {
  transform_map_weights(w, method = method, ratio = ratio, clamp_z = clamp_z,
                        gamma = gamma, winsor_q = winsor_q, global_norm = TRUE)
}

# Tier 0: regional weight-vector comparison
compare_weights_regional <- function(w_canon_df, w_alt_df,
                                     method = "robust_z", ratio = 4, clamp_z = 2,
                                     gamma = 2.5, winsor_q = 0.99,
                                     drop_ventricles = TRUE) {
  stopifnot(all(c("region", "weight") %in% names(w_canon_df)))
  stopifnot(all(c("region", "weight") %in% names(w_alt_df)))

  d <- inner_join(
    dplyr::rename(w_canon_df, w_canon = weight),
    dplyr::rename(w_alt_df,   w_alt   = weight),
    by = "region"
  )

  # Lock the comparison to regions where BOTH templates have valid data.
  # Matches align_feature_set(): finite AND positive, not merely finite.
  d <- d[is.finite(d$w_canon) & is.finite(d$w_alt) &
         d$w_canon > 0 & d$w_alt > 0, , drop = FALSE]

  if (isTRUE(drop_ventricles)) {
    d <- d[!grepl(VENT_PAT, d$region, perl = TRUE), , drop = FALSE]
  }

  # Each template is transformed on ITS OWN robust statistics, which is the
  # point of the comparison: agreement should reflect shared topography, not a
  # shared arbitrary scale imposed from outside.
  d$w_canon_g <- apply_pet_transform(d$w_canon, method, ratio, clamp_z, gamma, winsor_q)
  d$w_alt_g   <- apply_pet_transform(d$w_alt,   method, ratio, clamp_z, gamma, winsor_q)

  metrics <- function(x, y, label) {
    ok <- is.finite(x) & is.finite(y)
    n  <- sum(ok)
    if (n < 4) {
      return(data.frame(
        stage = label, n_regions = n,
        r_Pearson = NA_real_, rP_lo = NA_real_, rP_hi = NA_real_,
        r_Spearman = NA_real_, rS_lo = NA_real_, rS_hi = NA_real_
      ))
    }
    rP  <- suppressWarnings(cor(x[ok], y[ok], method = "pearson"))
    rS  <- suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
    ciP <- fisher_ci(rP, n)
    ciS <- fisher_ci(rS, n)
    data.frame(
      stage = label, n_regions = n,
      r_Pearson  = rP, rP_lo = ciP["lo"], rP_hi = ciP["hi"],
      r_Spearman = rS, rS_lo = ciS["lo"], rS_hi = ciS["hi"]
    )
  }

  tag <- if (identical(method, "legacy"))
    sprintf("after_legacy_gamma%.2f_winsorq%.2f", gamma, winsor_q)
  else
    sprintf("after_%s_ratio%.2f", method, ratio)
  out <- bind_rows(
    metrics(d$w_canon,   d$w_alt,   "raw_normalised"),
    metrics(d$w_canon_g, d$w_alt_g, tag)
  )
  rownames(out) <- NULL
  list(summary = out, regional = d)
}

# Tier 1: subject-level construct comparison
compare_scores_subject <- function(scores_canon, scores_alt) {
  needed <- c("id", "mbi_raw", "gbi", "mbi")
  stopifnot(all(needed %in% names(scores_canon)))
  stopifnot(all(needed %in% names(scores_alt)))

  d <- inner_join(
    scores_canon %>%
      dplyr::select(all_of(needed)) %>%
      dplyr::rename(mbi_raw_canon = mbi_raw, gbi_canon = gbi, mbi_canon = mbi),
    scores_alt %>%
      dplyr::select(all_of(needed)) %>%
      dplyr::rename(mbi_raw_alt = mbi_raw, gbi_alt = gbi, mbi_alt = mbi),
    by = "id"
  )

  metrics <- function(x, y, label) {
    ok <- is.finite(x) & is.finite(y)
    n  <- sum(ok)
    if (n < 4) {
      return(data.frame(
        score = label, n_subjects = n,
        r_Pearson = NA_real_, rP_lo = NA_real_, rP_hi = NA_real_,
        r_Spearman = NA_real_, rS_lo = NA_real_, rS_hi = NA_real_
      ))
    }
    rP  <- suppressWarnings(cor(x[ok], y[ok], method = "pearson"))
    rS  <- suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
    ciP <- fisher_ci(rP, n)
    ciS <- fisher_ci(rS, n)
    data.frame(
      score = label, n_subjects = n,
      r_Pearson  = rP, rP_lo = ciP["lo"], rP_hi = ciP["hi"],
      r_Spearman = rS, rS_lo = ciS["lo"], rS_hi = ciS["hi"]
    )
  }

  out <- bind_rows(
    metrics(d$mbi_raw_canon, d$mbi_raw_alt, "MBI_raw"),
    metrics(d$gbi_canon,     d$gbi_alt,     "GBI"),
    metrics(d$mbi_canon,     d$mbi_alt,     "MBI_residualised")
  )
  rownames(out) <- NULL
  list(summary = out, subject = d)
}
