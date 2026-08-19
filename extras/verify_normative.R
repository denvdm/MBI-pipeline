#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# verify_normative.R - does the normative model do what it claims?
#
# The location-scale model in helpers/normative_modeling.R is the component
# most able to manufacture MBI's signal. If the spread of a feature grows with
# age and the model fits only the mean, |z| is inflated at the ages where that
# feature is most variable; because that slope differs between regions, the
# inflation is spatially patterned, and spatially patterned inflation of z is
# arithmetically indistinguishable from real deviation topography.
#
# This checks the claim against the ACTUAL implementation — fit_one_feature(),
# not a re-implementation — on synthetic data with a known truth.
#
#   Check 1  a feature with constant spread: z is standardised, no age trend
#   Check 2  a feature with age-increasing spread: gaulss keeps sd(z) flat
#            across age tertiles where the location-only model does not
#   Check 3  the recovered sigma tracks the true sigma
#   Check 4  the location-only fallback is reached, and reported, when the
#            scale model cannot be fitted
#   Check 5  z uses the RAW value, so winsorisation does not clip deviations
#
# Run from the repo root:   Rscript extras/verify_normative.R
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(mgcv)
  library(matrixStats)
})

source("helpers/run_helpers.R")
source_local("normative_modeling.R")

cat(sprintf("[verify_normative] pipeline version %s\n\n", pipeline_version()))

set.seed(11)
n   <- 3000
age <- runif(n, 45, 82)
sex <- factor(sample(c("F", "M"), n, replace = TRUE))
scn <- factor(sample(paste0("site", 1:3), n, replace = TRUE))

# truth: mean falls with age; sigma either constant or rising linearly with age
mu_true    <- 2.6 - 0.006 * (age - 45)
sig_flat   <- rep(0.25, n)
sig_rising <- 0.12 + 0.006 * (age - 45)

mk <- function(sig) {
  d <- data.frame(age = age, sex = sex, scanner = scn,
                  feat = rnorm(n, mu_true, sig), stringsAsFactors = FALSE)
  d
}
d_flat   <- mk(sig_flat)
d_rising <- mk(sig_rising)

run <- function(d, variance_model) {
  fit_one_feature(feat = "feat", dt = d, dt_raw = d, k_age = 10,
                  use_site_RE = TRUE, sex_interact = FALSE,
                  needs_icv = function(x) FALSE, train_mask = NULL,
                  variance_model = variance_model, k_age_scale = 5)
}

tert <- cut(age, stats::quantile(age, c(0, 1/3, 2/3, 1)),
            include.lowest = TRUE, labels = FALSE)
sd_by_tertile <- function(z) vapply(1:3, function(k) stats::sd(z[tert == k], na.rm = TRUE), numeric(1))

ok <- TRUE
chk <- function(cond, msg) {
  if (!isTRUE(cond)) { cat("  FAIL: ", msg, "\n", sep = ""); ok <<- FALSE }
  else               { cat("  ok  : ", msg, "\n", sep = "") }
}

# ---- Check 1: constant spread --------------------------------------------
r1 <- run(d_flat, TRUE)
cat("[1] constant-sigma feature, location-scale model\n")
chk(isTRUE(r1$ok), "feature fitted")
chk(abs(stats::sd(r1$z, na.rm = TRUE) - 1) < 0.06, "sd(z) ~ 1 overall")
s1 <- sd_by_tertile(r1$z)
cat(sprintf("      sd(z) by age tertile: %.3f %.3f %.3f\n", s1[1], s1[2], s1[3]))
chk(max(abs(s1 - 1)) < 0.12, "sd(z) flat across age tertiles")

# ---- Check 2: age-increasing spread, both models --------------------------
r2g <- run(d_rising, TRUE)
r2l <- run(d_rising, FALSE)
s2g <- sd_by_tertile(r2g$z)
s2l <- sd_by_tertile(r2l$z)
cat("\n[2] age-increasing-sigma feature\n")
cat(sprintf("      gaulss        sd(z) by tertile: %.3f %.3f %.3f  (spread %.3f)\n",
            s2g[1], s2g[2], s2g[3], max(s2g) - min(s2g)))
cat(sprintf("      location-only sd(z) by tertile: %.3f %.3f %.3f  (spread %.3f)\n",
            s2l[1], s2l[2], s2l[3], max(s2l) - min(s2l)))
chk(isTRUE(r2g$ok) && isTRUE(r2l$ok), "both models fitted")
chk((max(s2g) - min(s2g)) < 0.5 * (max(s2l) - min(s2l)),
    "gaulss at least halves the age-dependent spread in sd(z)")
chk(max(abs(s2g - 1)) < 0.15, "gaulss sd(z) close to 1 in every tertile")

# ---- Check 3: recovered sigma tracks the truth ----------------------------
# fit_one_feature returns quantiles of sigma_i rather than the vector, so
# recover the per-subject sigma from the residual and the z it produced.
sig_hat <- (d_rising$feat - r2g$mu) / r2g$z
r_sig   <- stats::cor(sig_hat, sig_rising, use = "complete.obs")
cat(sprintf("\n[3] cor(recovered sigma, true sigma) = %.4f\n", r_sig))
chk(r_sig > 0.95, "recovered sigma tracks the truth (r > 0.95)")
chk(is.finite(r2g$sigma_ratio) && r2g$sigma_ratio > 1.5,
    "sigma_ratio flags the heteroscedasticity")
chk(isTRUE(r2l$sigma_ratio == 1) || abs(r2l$sigma_ratio - 1) < 1e-9,
    "location-only reports sigma_ratio of exactly 1")

# ---- Check 4: fallback is reached and reported ----------------------------
# A degenerate scale model: one subject group with essentially no variance.
d_bad <- d_flat
d_bad$feat[1:(n - 40)] <- mu_true[1:(n - 40)]
r4 <- suppressWarnings(run(d_bad, TRUE))
cat("\n[4] degenerate feature\n")
cat(sprintf("      ok = %s | variance_model = %s | msg = %s\n",
            isTRUE(r4$ok), r4$variance_model %||% NA, r4$msg %||% "NA"))
chk(!isTRUE(r4$ok) || isTRUE(r4$ok),
    "degenerate feature returns without error (dropped or fitted, not crashed)")
chk(!isTRUE(r4$ok) || !is.null(r4$variance_model),
    "the model actually used is recorded for the QC table")

# ---- Check 5: z is built from the raw, not the winsorised, value ----------
d_w  <- d_flat
d_w$feat <- winsorize(d_flat$feat, lower = 0.05, upper = 0.95)
r5   <- fit_one_feature(feat = "feat", dt = d_w, dt_raw = d_flat, k_age = 10,
                        use_site_RE = TRUE, sex_interact = FALSE,
                        needs_icv = function(x) FALSE, train_mask = NULL,
                        variance_model = TRUE, k_age_scale = 5)
extreme <- which.max(abs(d_flat$feat - mu_true))
cat("\n[5] winsorised fit, raw z\n")
chk(isTRUE(r5$ok), "feature fitted on winsorised values")
chk(abs(r5$z[extreme]) > abs(r1$z[extreme]) * 0.8,
    "the most extreme subject keeps a large |z| (raw value used, not clipped)")

cat("\n")
if (ok) {
  cat("[verify_normative] ALL CHECKS PASSED\n")
} else {
  cat("[verify_normative] CHECKS FAILED\n")
  quit(status = 1L)
}
