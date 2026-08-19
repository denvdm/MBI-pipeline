#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(tidyr)
})

# --- column pickers / guards --------------------------------------------------

pick_z_col <- function(zs) {
  if ("z_oriented" %in% names(zs)) {
    "z_oriented"
  } else if ("z" %in% names(zs)) {
    "z"
  } else {
    stop("Need zs$z_oriented or zs$z")
  }
}

prep_fdg <- function(fdg) {
  out <- if ("w_pet" %in% names(fdg)) {
    fdg %>% dplyr::select(region, w_pet)
  } else if ("weight" %in% names(fdg)) {
    fdg %>% dplyr::transmute(region, w_pet = weight)
  } else {
    stop("fdg must have column 'w_pet' or 'weight'")
  }
  # Match align_feature_set() in compute_mbi.R: a region without a finite,
  # positive map weight does not enter the index, so it must not enter the
  # QC that is meant to describe the index either.
  out %>%
    dplyr::mutate(region = as.character(region),
                  w_pet  = as.numeric(w_pet)) %>%
    dplyr::filter(is.finite(w_pet), w_pet > 0)
}

# NOTE: w_pet here is the RAW map value, not the winsorised/normalised/^gamma
# vector that enters the index. Spearman correlations are unaffected (monotone
# transform); reported Pearson values describe alignment with the biological
# map itself rather than with the transformed weights. Intentional, but state
# it wherever these numbers are reported.

# --- Fisher z CIs -------------------------------------------------------------

fisher_ci <- function(r, n, conf = 0.95) {
  if (!is.finite(r) || n <= 3) return(c(lo = NA_real_, hi = NA_real_))
  r  <- max(min(r, 0.999999), -0.999999)
  z  <- atanh(r)
  se <- 1 / sqrt(n - 3)
  zc <- qnorm(1 - (1 - conf) / 2)
  c(lo = tanh(z - zc * se), hi = tanh(z + zc * se))
}

fisher_mean_ci <- function(r_vec, n_vec, conf = 0.95) {
  ok <- is.finite(r_vec) & is.finite(n_vec) & (n_vec > 3)
  if (sum(ok) < 2) return(c(mean = NA, lo = NA, hi = NA))
  z <- atanh(pmin(pmax(r_vec[ok], -0.999999), 0.999999))
  w <- pmax(n_vec[ok] - 3, 1e-8)
  mz <- sum(w * z) / sum(w)
  se <- 1 / sqrt(sum(w))
  zc <- qnorm(1 - (1 - conf) / 2)
  c(
    mean = tanh(mz),
    lo   = tanh(mz - zc * se),
    hi   = tanh(mz + zc * se)
  )
}

# --- SubVol gate + collapse to (id,region) once ------------------------------

gate_subvol <- function(zs, size_df, min_mm3 = 1000) {
  # ABSOLUTE volume threshold, replacing the previous within-SubVol quantile.
  # The quantile version was recomputed over whatever SubVol set survived
  # upstream filtering, so dropping a region for unrelated reasons silently
  # shifted which OTHER structures were gated out. An absolute cut is invariant
  # to that, and keeps this QC stage gated identically to stage 03 regardless
  # of the order in which filters are applied.
  keep_sub <- size_df %>%
    dplyr::filter(modality == "SubVol") %>%
    dplyr::filter(is.finite(size), size >= min_mm3) %>%
    dplyr::pull(region)

  zs %>%
    dplyr::filter(modality != "SubVol" | region %in% keep_sub)
}

# Collapse to one deviation per (subject, region) by averaging across the
# modalities that describe that region.
#
# NOTE THAT THIS IS NOT WHAT THE INDEX DOES. The index treats every (modality,
# region) as its own feature with its own weight, so a cortical region enters
# twice, via CT and via SA. Here the two are averaged into one number, because
# the diagnostics below are region-level by construction: a PCA over regions, a
# leave-one-region-out slope, a per-subject correlation against a region-wise
# map. The alignment these functions report is therefore an alignment of a
# region-level summary of the deviations with the map, not of the index itself.
# Close enough to be informative, different enough to state when reporting.
collapse_zs <- function(zs, fdg_use, mods_use = c("CT", "SA", "SubVol")) {
  zcol <- pick_z_col(zs)

  zs %>%
    dplyr::filter(
      modality %in% mods_use,
      is.finite(.data[[zcol]])
    ) %>%
    dplyr::group_by(id, region) %>%
    dplyr::summarise(
      z = mean(.data[[zcol]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::semi_join(fdg_use, by = "region")
}

# --- PC1, computed once ------------------------------------------------------
# Both the alignment summary and the loadings table need the same PC1. Compute
# it once with region_pc1() and pass the result to each, rather than letting
# every consumer rebuild the wide matrix and decompose it again.
#
# CAVEAT to state wherever these loadings are reported: X[is.na(X)] <- 0 after
# scaling is mean-imputation. A region observed in fewer subjects is pulled
# toward the sample mean, which shrinks its loading. Acceptable for a QC
# diagnostic, misleading if read as a substantive regional ranking.
region_pc1 <- function(Z_collapsed) {
  M <- Z_collapsed %>%
    tidyr::pivot_wider(names_from = region, values_from = z) %>%
    dplyr::arrange(id) %>%
    dplyr::select(-id) %>%
    as.matrix()

  if (!ncol(M)) stop("No regions after collapse")

  vars <- apply(M, 2, sd, na.rm = TRUE)
  keep <- is.finite(vars) & vars > 0
  if (sum(keep) < 2) stop("Too few variable regions for PCA")

  X <- scale(M[, keep, drop = FALSE])
  X[is.na(X)] <- 0

  # Eigen-decomposition of the region x region cross-product rather than a full
  # SVD of the subject x region matrix. Identical to prcomp()'s rotation and
  # sdev to machine precision (verified), at a fraction of the cost: the small
  # matrix is n_regions square, the large one has a row per subject.
  ev <- eigen(crossprod(X), symmetric = TRUE)
  lam <- pmax(ev$values, 0)

  list(
    region       = colnames(X),
    loading      = as.numeric(ev$vectors[, 1]),
    var_explained = unname(lam[1] / sum(lam))
  )
}

pc_pet_alignment <- function(Z_collapsed, fdg_use, mods_tag, pc1 = NULL) {
  if (is.null(pc1)) pc1 <- region_pc1(Z_collapsed)

  pc1_tbl <- dplyr::tibble(region = pc1$region, loading = pc1$loading) %>%
    dplyr::inner_join(fdg_use, by = "region")

  rP  <- suppressWarnings(cor(pc1_tbl$loading, pc1_tbl$w_pet, use = "pairwise"))
  rS  <- suppressWarnings(cor(pc1_tbl$loading, pc1_tbl$w_pet, method = "spearman", use = "pairwise"))
  n   <- nrow(pc1_tbl)
  ciP <- fisher_ci(rP, n)
  ciS <- fisher_ci(rS, n)

  data.frame(
    metric      = paste0("PCA_loadings_vs_PET_", mods_tag),
    r_Pearson   = rP,
    rP_lo       = ciP["lo"],
    rP_hi       = ciP["hi"],
    r_Spearman  = rS,
    rS_lo       = ciS["lo"],
    rS_hi       = ciS["hi"],
    n_regions   = n
  )
}

# --- LOO contribution map, computed once -------------------------------------
# For each region: regress the leave-that-region-out mean deviation on the
# region's own deviation. The slope says how much a region's deviation tracks
# the rest of the brain.
#
# Closed form rather than one lm() per region. This is a simple bivariate
# regression, so the slope is Sxy/Sxx and the standard error follows from the
# residual sum of squares — algebraically identical to lm(), and verified equal
# to it to ~1e-15 on the slope, standard error, t and p. The previous version
# One lm() per region inside dplyr::do() costs a full model fit over every
# subject, and stage 04 needs these betas in two places.
region_loo_betas <- function(Z_collapsed) {
  sub_sum <- Z_collapsed %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(S = sum(z), P = dplyr::n(), .groups = "drop")

  Z_collapsed %>%
    dplyr::inner_join(sub_sum, by = "id") %>%
    dplyr::mutate(GI_minus_r = (S - z) / pmax(P - 1L, 1L)) %>%
    dplyr::group_by(region) %>%
    dplyr::summarise(
      n_obs = dplyr::n(),
      Sxx   = sum((z - mean(z))^2),
      Sxy   = sum((z - mean(z)) * (GI_minus_r - mean(GI_minus_r))),
      Syy   = sum((GI_minus_r - mean(GI_minus_r))^2),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      beta_LOO = ifelse(Sxx > 0, Sxy / Sxx, NA_real_),
      rss      = pmax(Syy - ifelse(Sxx > 0, Sxy^2 / Sxx, NA_real_), 0),
      se       = ifelse(n_obs > 2 & Sxx > 0, sqrt(rss / (n_obs - 2) / Sxx), NA_real_),
      t        = beta_LOO / se,
      p        = ifelse(is.finite(t), 2 * stats::pt(-abs(t), n_obs - 2), NA_real_)
    ) %>%
    dplyr::select(region, beta_LOO, se, t, p, n_obs)
}

# --- LOO contribution map × PET ----------------------------------------------

loo_pet_alignment <- function(Z_collapsed, fdg_use, betas = NULL) {
  if (is.null(betas)) betas <- region_loo_betas(Z_collapsed)

  beta_pet <- betas %>%
    dplyr::select(region, beta_LOO) %>%
    dplyr::inner_join(fdg_use, by = "region")

  rP  <- suppressWarnings(cor(beta_pet$beta_LOO, beta_pet$w_pet, use = "pairwise"))
  rS  <- suppressWarnings(cor(beta_pet$beta_LOO, beta_pet$w_pet, method = "spearman", use = "pairwise"))
  n   <- nrow(beta_pet)
  ciP <- fisher_ci(rP, n)
  ciS <- fisher_ci(rS, n)

  data.frame(
    metric      = "LOO_beta_vs_PET",
    r_Pearson   = rP,
    rP_lo       = ciP["lo"],
    rP_hi       = ciP["hi"],
    r_Spearman  = rS,
    rS_lo       = ciS["lo"],
    rS_hi       = ciS["hi"],
    n_regions   = n
  )
}

# --- Subject-level map×PET correlation summary --------------------------------

subject_pet_alignment <- function(Z_collapsed, fdg_use) {
  r_df <- Z_collapsed %>%
    dplyr::inner_join(fdg_use, by = "region") %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(
      r_s  = suppressWarnings(cor(z, w_pet, method = "spearman", use = "pairwise")),
      r_p  = suppressWarnings(cor(z, w_pet, method = "pearson",  use = "pairwise")),
      n_reg = dplyr::n(),
      .groups = "drop"
    )

  ci_s <- fisher_mean_ci(r_df$r_s, r_df$n_reg)
  ci_p <- fisher_mean_ci(r_df$r_p, r_df$n_reg)

  # The Fisher interval treats the per-subject correlations as independent
  # estimates of one quantity, so at cohort scale it collapses to a width of
  # essentially nothing and says only that the mean is precisely estimated. It
  # is NOT a statement about how much subjects differ from one another. The
  # quartiles and SD below are the numbers to quote for that; report the
  # interval only alongside them.
  qs <- stats::quantile(r_df$r_s[is.finite(r_df$r_s)], c(0.25, 0.5, 0.75),
                        na.rm = TRUE, names = FALSE)

  data.frame(
    n_subjects = sum(is.finite(r_df$r_s)),
    mean_s     = unname(ci_s["mean"]),
    lo_s       = unname(ci_s["lo"]),
    hi_s       = unname(ci_s["hi"]),
    mean_p     = unname(ci_p["mean"]),
    lo_p       = unname(ci_p["lo"]),
    hi_p       = unname(ci_p["hi"]),
    sd_s       = stats::sd(r_df$r_s, na.rm = TRUE),
    q25_s      = qs[1],
    median_s   = qs[2],
    q75_s      = qs[3]
  )
}

# --- PC1 loadings table (optionally sign-align to PET) ------------------------

pc1_loadings_table <- function(Z_collapsed, fdg_use = NULL, align_to_pet = TRUE,
                               pc1 = NULL) {
  if (is.null(pc1)) pc1 <- region_pc1(Z_collapsed)

  load <- pc1$loading
  reg  <- pc1$region

  # The sign of an eigenvector is arbitrary, so fix it: flip so the correlation
  # with the map is non-negative. Without this the loadings table can change
  # sign between runs for no substantive reason.
  if (!is.null(fdg_use)) {
    tmp <- dplyr::tibble(region = reg, loading = load) %>%
      dplyr::inner_join(fdg_use, by = "region")
    rP <- suppressWarnings(cor(tmp$loading, tmp$w_pet, use = "pairwise"))
    if (is.finite(rP) && align_to_pet && rP < 0) load <- -load
  }

  dplyr::tibble(
    region            = reg,
    pc1_loading       = load,
    var_explained_pc1 = pc1$var_explained
  )
}

# --- Mean z per region (collapsed across modalities & subjects) ---------------

mean_z_per_region <- function(Z_collapsed) {
  Z_collapsed %>%
    dplyr::group_by(region) %>%
    dplyr::summarise(
      mean_z = mean(z, na.rm = TRUE),
      .groups = "drop"
    )
}

# --- LOO beta map (region contribution slope) ---------------------------------

loo_beta_map <- function(Z_collapsed, betas = NULL) {
  if (is.null(betas)) betas <- region_loo_betas(Z_collapsed)
  betas <- dplyr::select(betas, region, beta_LOO, se, t, p)
  betas$FDR <- stats::p.adjust(betas$p, method = "BH")
  betas
}
