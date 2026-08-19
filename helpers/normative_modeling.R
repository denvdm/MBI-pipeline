#!/usr/bin/env Rscript

# Small, self-contained helpers for normative modeling.
#
# Design goal: keep this file free of hidden dependencies on caller environments.
# All runtime state is passed in explicitly from the runner script.
#
# ---------------------------------------------------------------------------
# LOCATION-SCALE NORMATIVE MODEL (mgcv gaulss)
#
# The default normative model fits BOTH the conditional mean and the
# conditional spread as functions of age, sex and scanner, so each subject's
# deviation is z_i = (y_i - mu_i) / sigma_i with a sigma of their own.
#
# The scale half is not optional decoration. Regional morphometry becomes more
# variable with age, at different rates in different regions. A location-only
# model — one scalar sigma per feature — therefore inflates |z| systematically
# at the ages where a feature is most variable. On simulated data with a
# linearly age-increasing SD:
#
#     sd(z) by age tertile, single sigma  :  0.70   0.95   1.27
#     sd(z) by age tertile, gaulss        :  1.00   0.99   1.01
#
# For a mean index that mostly rescales. For MBI it does more: MBI is a
# CONTRAST across regions, so where the variance-by-age slope differs between
# regions the inflation is spatially patterned — and spatially patterned
# inflation of z is arithmetically indistinguishable from genuine topography in
# the deviation profile. Of every upstream choice in this pipeline, this is the
# one most able to manufacture the thing MBI claims to measure. Treat
# normative$variance_model: false as a diagnostic comparison, not as an
# alternative specification; a run made that way records itself as
# location-only in the QC file.
#
# Implementation note: gaulss fits two linear predictors, one for the mean and
# one for the PRECISION 1/sigma. predict(type = "response") returns a
# two-column matrix whose second column is the precision, so sigma_i is its
# reciprocal — not the second column itself.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(mgcv)
  library(matrixStats)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ------------------------------
# Winsorizing
# ------------------------------
winsorize <- function(x, lower = 0.001, upper = 0.999, na.rm = TRUE) {
  if (!is.numeric(x)) return(x)
  q <- quantile(x, probs = c(lower, upper), na.rm = na.rm)
  x[x < q[1]] <- q[1]
  x[x > q[2]] <- q[2]
  x
}

# ------------------------------
# Formula builders
# ------------------------------
make_formula <- function(feat_name,
                         k_age,
                         with_icv,
                         site_RE,
                         sex_int,
                         has_multiple_scanners) {
  base <- sprintf("`%s` ~ s(age, k=%d, bs='cs') + sex", feat_name, k_age)
  if (isTRUE(sex_int)) {
    base <- sprintf("%s + s(age, by=sex, k=%d, bs='cs')", base, k_age)
  }
  if (isTRUE(with_icv)) {
    base <- sprintf("%s + icv", base)
  }
  if (isTRUE(site_RE) && isTRUE(has_multiple_scanners)) {
    base <- sprintf("%s + s(scanner, bs='re')", base)
  }
  as.formula(base)
}

# Scale (log-precision) linear predictor for gaulss. Deliberately leaner than
# the mean model: age and sex, plus the scanner random effect where there is
# more than one scanner. ICV is NOT included — it belongs in the mean model as
# a nuisance covariate, and adding it to the variance model would let head size
# rescale a subject's deviation, which is not what a normative z-score should
# do. k is smaller than for the mean because a variance curve estimated at this
# N does not support 10 basis functions and over-flexible scale models are the
# usual cause of gaulss convergence failures.
make_scale_formula <- function(k_age_scale, site_RE, has_multiple_scanners) {
  base <- sprintf("~ s(age, k=%d, bs='cs') + sex", k_age_scale)
  if (isTRUE(site_RE) && isTRUE(has_multiple_scanners)) {
    base <- sprintf("%s + s(scanner, bs='re')", base)
  }
  as.formula(base)
}

# ------------------------------
# Fit & predict for a single feature
# ------------------------------
# Returns a list with ok / z / mu / sigma. Under the location-scale model
# `sigma` in the QC table is the MEDIAN of the per-subject sigmas, and
# sigma_p10 / sigma_p90 record how much it moves across the sample.
fit_one_feature <- function(feat,
                            dt,
                            dt_raw,
                            k_age = 10,
                            use_site_RE = TRUE,
                            sex_interact = FALSE,
                            needs_icv = function(x) FALSE,
                            train_mask = NULL,
                            variance_model = TRUE,
                            k_age_scale = 5) {
  # Defensive checks
  if (!feat %in% names(dt)) {
    return(list(ok = FALSE, feat = feat, msg = "feature not found"))
  }
  y <- dt[[feat]]
  if (!is.numeric(y)) return(list(ok = FALSE, feat = feat, msg = "non-numeric"))
  if (!any(is.finite(y))) return(list(ok = FALSE, feat = feat, msg = "all NA/inf"))

  with_icv <- isTRUE(needs_icv(feat)) && ("icv" %in% names(dt)) && any(is.finite(dt$icv))
  has_multiple_scanners <- ("scanner" %in% names(dt)) && length(levels(dt$scanner)) > 1L
  frm <- make_formula(
    feat_name = feat,
    k_age = k_age,
    with_icv = with_icv,
    site_RE = use_site_RE,
    sex_int = sex_interact,
    has_multiple_scanners = has_multiple_scanners
  )

  # mgcv will drop rows with missing RHS; we also require finite y and age.
  idx <- which(is.finite(y) & is.finite(dt$age))

  # If train_mask provided, fit only on training subjects
  train_idx <- if (!is.null(train_mask)) intersect(idx, which(train_mask)) else idx
  if (length(train_idx) < 50L) {
    return(list(ok = FALSE, feat = feat, msg = "too few observations"))
  }

  dsub <- dt[train_idx, , drop = FALSE]
  dsub[[feat]] <- y[train_idx]

  # ---- fit -----------------------------------------------------------------
  # gaulss can fail to converge on a badly behaved feature. Rather than let
  # that kill the feature outright we fall back to the location-only model and
  # SAY SO in the QC table, so a mixed run is visible rather than silent.
  used_variance_model <- isTRUE(variance_model)
  fallback_msg <- NA_character_

  fit <- NULL
  if (used_variance_model) {
    frm_scale <- make_scale_formula(k_age_scale, use_site_RE, has_multiple_scanners)
    fit <- try(mgcv::gam(list(frm, frm_scale), data = dsub, method = "REML",
                         family = mgcv::gaulss()), silent = TRUE)
    if (inherits(fit, "try-error")) {
      fallback_msg <- paste0("gaulss failed, fell back to location-only: ",
                             trimws(as.character(fit)))
      used_variance_model <- FALSE
      fit <- NULL
    }
  }
  if (is.null(fit)) {
    fit <- try(mgcv::gam(frm, data = dsub, method = "REML", family = gaussian()),
               silent = TRUE)
    if (inherits(fit, "try-error")) {
      return(list(ok = FALSE, feat = feat, msg = as.character(fit)))
    }
  }

  # ---- EDF of the MEAN age smooth -----------------------------------------
  # Under gaulss the coefficient names of the SCALE predictor are prefixed
  # "s.1(", so "^s\\(age" continues to isolate the mean smooth. The scale
  # smooth is reported separately as edf_age_scale.
  nm <- names(fit$edf)
  edf_age <- tryCatch(sum(fit$edf[grepl("^s\\(age", nm)]), error = function(e) NA_real_)
  edf_age_scale <- if (used_variance_model) {
    tryCatch(sum(fit$edf[grepl("^s\\.1\\(age", nm)]), error = function(e) NA_real_)
  } else NA_real_

  # ---- predictions ---------------------------------------------------------
  if (used_variance_model) {
    pr <- try(predict(fit, newdata = dt, type = "response"), silent = TRUE)
    if (inherits(pr, "try-error") || is.null(dim(pr)) || ncol(pr) < 2L) {
      return(list(ok = FALSE, feat = feat, msg = "gaulss prediction failed",
                  edf_age = edf_age))
    }
    mu_full  <- as.numeric(pr[, 1])
    prec     <- as.numeric(pr[, 2])          # gaulss second component IS 1/sigma
    prec[!is.finite(prec) | prec <= 0] <- NA_real_
    sigma_i  <- 1 / prec
  } else {
    mu_full <- as.numeric(predict(fit, newdata = dt, type = "response"))
    s0 <- sqrt(summary(fit)$scale)
    if (!is.finite(s0) || s0 <= 0) {
      r_train <- resid(fit)
      s0 <- matrixStats::mad(r_train, constant = 1.4826)
      if (!is.finite(s0) || s0 <= 0) s0 <- sd(r_train, na.rm = TRUE)
    }
    if (!is.finite(s0) || s0 <= 0) {
      return(list(ok = FALSE, feat = feat, msg = "non-positive sigma",
                  edf_age = edf_age))
    }
    sigma_i <- rep(s0, length(mu_full))
  }

  if (!any(is.finite(sigma_i) & sigma_i > 0)) {
    return(list(ok = FALSE, feat = feat, msg = "no usable sigma",
                edf_age = edf_age))
  }

  # ---- variance explained --------------------------------------------------
  # Computed explicitly rather than read from summary(fit)$r.sq. summary.gam
  # returns NULL for r.sq under gaulss, and the low-edf arm of the
  # degenerate-fit gate below tests is.finite(r2) — so an NA here would silently
  # switch half the gate off. Keep this computed locally.
  y_tr  <- dsub[[feat]]
  mu_tr <- mu_full[train_idx]
  ok_tr <- is.finite(y_tr) & is.finite(mu_tr)
  r2 <- if (sum(ok_tr) > 2L && sd(y_tr[ok_tr]) > 0) {
    1 - sum((y_tr[ok_tr] - mu_tr[ok_tr])^2) / sum((y_tr[ok_tr] - mean(y_tr[ok_tr]))^2)
  } else NA_real_

  # ---- degenerate-fit gate -------------------------------------------------
  # Two failure modes, both of which produce a z-score that is not a deviation:
  # an age smooth that has absorbed the signal (edf_age > 8 against k_age = 10),
  # and a fit that explains essentially nothing while using almost no
  # flexibility (edf_age < 1.5 with r2 < 0.02). Check normative_qc.txt if a
  # feature you expected is missing from the index.
  if (!is.finite(edf_age) || edf_age > 8 ||
      (edf_age < 1.5 && is.finite(r2) && r2 < 0.02)) {
    message(sprintf("[NORMATIVE] Skipping %s (edf_age=%.2f, r2=%.3f)", feat, edf_age, r2))
    return(list(ok = FALSE, feat = feat, msg = "degenerate fit",
                edf_age = edf_age, r2 = r2))
  }

  # z-scores use raw (pre-winsorized) values for interpretability
  if (!feat %in% names(dt_raw)) {
    return(list(ok = FALSE, feat = feat, msg = "dt_raw missing feature",
                edf_age = edf_age, r2 = r2))
  }
  z <- (dt_raw[[feat]] - mu_full) / sigma_i

  sig_q <- stats::quantile(sigma_i, c(0.10, 0.50, 0.90), na.rm = TRUE)

  list(
    ok             = TRUE,
    feat           = feat,
    z              = z,
    mu             = mu_full,
    sigma          = as.numeric(sig_q[2]),
    sigma_p10      = as.numeric(sig_q[1]),
    sigma_p90      = as.numeric(sig_q[3]),
    sigma_ratio    = as.numeric(sig_q[3] / sig_q[1]),   # 1.0 == homoscedastic
    variance_model = used_variance_model,
    edf_age        = edf_age,
    edf_age_scale  = edf_age_scale,
    r2             = r2,
    n_train        = nrow(dsub),
    msg            = fallback_msg,
    formula        = paste(deparse(frm), collapse = " ")
  )
}
