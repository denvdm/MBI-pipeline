#!/usr/bin/env Rscript

# Small, self-contained helpers for normative modeling.
#
# Design goal: keep this file free of hidden dependencies on caller environments.
# All runtime state is passed in explicitly from the runner script.

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
# Formula builder
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

# ------------------------------
# Fit & predict for a single feature
# ------------------------------
fit_one_feature <- function(feat,
                            dt,
                            dt_raw,
                            k_age = 10,
                            use_site_RE = TRUE,
                            sex_interact = FALSE,
                            needs_icv = function(x) FALSE) {
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
  if (length(idx) < 50L) {
    return(list(ok = FALSE, feat = feat, msg = "too few observations"))
  }

  dsub <- dt[idx, , drop = FALSE]
  dsub[[feat]] <- y[idx]

  fit <- try(mgcv::gam(frm, data = dsub, method = "REML", family = gaussian), silent = TRUE)
  if (inherits(fit, "try-error")) {
    return(list(ok = FALSE, feat = feat, msg = as.character(fit)))
  }

  sm <- summary(fit)
  r2 <- sm$r.sq %||% NA_real_

  # EDF of age smooth (including sex-by-age smooth if present)
  edf_age <- tryCatch({
    sum(fit$edf[grepl("^s\\(age", names(fit$edf))])
  }, error = function(e) NA_real_)

  # Basic gate to skip degenerate fits
  if (!is.finite(edf_age) || edf_age > 8 || (edf_age < 1.5 && is.finite(r2) && r2 < 0.02)) {
    message(sprintf("[NORMATIVE] Skipping %s (edf_age=%.2f, r2=%.3f)", feat, edf_age, r2))
    return(list(ok = FALSE, feat = feat, msg = "degenerate fit", edf_age = edf_age, r2 = r2))
  }

  # Predictions on all rows (mgcv handles missing RHS)
  mu_full <- predict(fit, newdata = dt, type = "response")

  # Prefer model residual scale for stability
  sigma <- sqrt(summary(fit)$scale)
  if (!is.finite(sigma) || sigma <= 0) {
    r_train <- resid(fit)
    sigma <- matrixStats::mad(r_train, constant = 1.4826)
    if (!is.finite(sigma) || sigma <= 0) sigma <- sd(r_train, na.rm = TRUE)
  }
  if (!is.finite(sigma) || sigma <= 0) {
    return(list(ok = FALSE, feat = feat, msg = "non-positive sigma", edf_age = edf_age, r2 = r2))
  }

  # z-scores use raw (pre-winsorized) values for interpretability
  if (!feat %in% names(dt_raw)) {
    return(list(ok = FALSE, feat = feat, msg = "dt_raw missing feature", edf_age = edf_age, r2 = r2))
  }
  z <- (dt_raw[[feat]] - mu_full) / sigma

  list(
    ok = TRUE,
    feat = feat,
    z = z,
    mu = mu_full,
    sigma = sigma,
    edf_age = edf_age,
    r2 = r2,
    n_train = nrow(dsub),
    formula = deparse(formula(fit))
  )
}
