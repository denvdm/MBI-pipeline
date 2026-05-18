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

# Apply the pipeline's PET-weight transform: clamp/winsor/normalize/exponentiate.
# Mirrors compute_mbi.R::compute_mbi_flagship() so the comparison reflects what
# actually enters MBI rather than just the raw map.
apply_pet_transform <- function(w, gamma = 2.5, winsor_q = 0.99) {
  w[!is.finite(w) | w <= 0] <- 0
  # winsor positives at p99 BEFORE mean-norm (same order as compute_mbi.R)
  if (!is.null(winsor_q) && is.finite(winsor_q) && winsor_q > 0 && winsor_q < 1) {
    pos <- w[w > 0]
    if (length(pos) > 1) {
      qhi <- as.numeric(quantile(pos, winsor_q, na.rm = TRUE))
      w <- pmin(w, qhi)
    }
  }
  m <- mean(w[w > 0], na.rm = TRUE)
  if (is.finite(m) && m > 0) w <- w / m
  if (!is.null(gamma) && is.finite(gamma) && gamma != 1) w <- w ^ gamma
  w
}

# Tier 0: regional weight-vector comparison
compare_weights_regional <- function(w_canon_df, w_alt_df,
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
  d <- d[is.finite(d$w_canon) & is.finite(d$w_alt), , drop = FALSE]

  if (isTRUE(drop_ventricles)) {
    d <- d[!grepl(VENT_PAT, d$region, perl = TRUE), , drop = FALSE]
  }

  d$w_canon_g <- apply_pet_transform(d$w_canon, gamma, winsor_q)
  d$w_alt_g   <- apply_pet_transform(d$w_alt,   gamma, winsor_q)

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

  tag <- sprintf("after_gamma%.2f_winsorq%.2f", gamma, winsor_q)
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
