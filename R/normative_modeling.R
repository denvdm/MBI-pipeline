# Normative modelling utilities using GAMs and mixed-effects components.
# Fits feature-wise deviation models and outputs subject-level standardized residuals.

## Winsorizing
winsorize <- function(x, lower = 0.001, upper = 0.999, na.rm = TRUE) {
  if (!is.numeric(x)) return(x)
  q <- quantile(x, probs = c(lower, upper), na.rm = na.rm)
  x[x < q[1]] <- q[1]
  x[x > q[2]] <- q[2]
  x}

#============================== MODEL SPEC ====================================

make_formula <- function(feat_name, with_icv, site_RE, sex_int) {
  base <- sprintf("`%s` ~ s(age, k=%d, bs='cs') + sex", feat_name, k_age)
  if (sex_int) base <- sprintf("%s + s(age, by=sex, k=%d, bs='cs')", base, k_age)
  if (with_icv) base <- sprintf("%s + icv", base)
  if (site_RE && length(levels(dt$scanner)) > 1L) base <- sprintf("%s + s(scanner, bs='re')", base)
  as.formula(base)
}

#============================= FIT & PREDICT ==================================

fit_one_feature <- function(feat) {
  y <- dt[[feat]]
  if (!is.numeric(y)) return(list(ok=FALSE, feat=feat, msg="non-numeric"))
  if (!any(is.finite(y))) return(list(ok=FALSE, feat=feat, msg="all NA/inf"))

  with_icv <- needs_icv(feat) && any(is.finite(dt$icv))
  frm <- make_formula(feat, with_icv, use_site_RE, sex_interact)

  # mgcv handles missingness on RHS; we drop rows (subjects) with non-finite y
  idx <- which(is.finite(y) & is.finite(dt$age))
  if (length(idx) < 50L) return(list(ok=FALSE, feat=feat, msg="too few observations"))

  dsub <- dt[idx, , drop = FALSE]
  dsub[[feat]] <- y[idx]

  # Fit GAM
  fit <- try(gam(frm, data = dsub, method = "REML", family = family_gaussian), silent = TRUE)
  if (inherits(fit, "try-error")) return(list(ok=FALSE, feat=feat, msg=as.character(fit)))

  sm <- summary(fit)
  edf_age <- tryCatch(sum(fit$edf[grepl("^s\\(age", names(fit$edf))]), error=function(e) NA_real_)
  r2 <- sm$r.sq

  if (!is.finite(edf_age) || edf_age > 8 || (edf_age < 1.5 && r2 < 0.02)) {
    message(sprintf("[NORMATIVE] Skipping %s (edf_age=%.2f, r2=%.3f)", feat, edf_age, r2))
    return(list(ok=FALSE, feat=feat))  
  }

  # Predictions & residuals on *all* rows where covariates are finite (and age finite)
  mu_full <- predict(fit, newdata = dt, type = "response")

  # Use model residual scale (sigma) rather than sample SD of resid to be stable
  sigma <- sqrt(summary(fit)$scale)
  if (!is.finite(sigma) || sigma <= 0) {
    # fallback: robust scale of residuals in train
    r_train <- resid(fit)
    sigma <- matrixStats::mad(r_train, constant = 1.4826)
    if (!is.finite(sigma) || sigma <= 0) sigma <- sd(r_train, na.rm = TRUE)
  }

  # Deviations (z)
  z <- (dt_raw[[feat]] - mu_full) / sigma

  list(
    ok = TRUE,
    feat = feat,
    z = z,
    mu = mu_full,
    sigma = sigma,
    edf_age = tryCatch(sum(fit$edf[grepl("^s\\(age", names(fit$edf))]), error = function(e) NA_real_),
    n_train = nrow(dsub),
    formula = deparse(formula(fit))
  )
}
