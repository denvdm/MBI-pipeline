# Convenience wrappers for linear and logistic models linking MBI metrics to outcomes.
# Standardizes covariates, handles missingness, and outputs tidy summaries.

## 0) Inputs: df must contain columns: eid, gbi, nmi, covariates, and outcomes (wide).
## Example placeholders:
# covars <- c("age","sex","site","ICV","fasting","batch")
# y_cont <- c("BMI","Glucose","HbA1c","Triglycerides","HDL","CRP")        # edit to your set
# y_diag <- c("T2D_dx","Dyslipidemia_dx","MDD_dx")                         # 0/1

## 1) PET-tilt (standardized residual of nmi ~ gbi)
make_tilt <- function(df, save_prov_path = NULL) {
  fit <- lm(nmi ~ gbi, data = df)
  a <- unname(coef(fit)["(Intercept)"])
  b <- unname(coef(fit)["gbi"])
  s <- sd(residuals(fit), na.rm = TRUE)
  df$tilt_fdg <- (df$nmi - (a + b * df$gbi)) / s
  if (!is.null(save_prov_path)) {
    saveRDS(list(
      a = a, b = b, s = s,
      r = cor(df$nmi, df$gbi, use = "complete.obs"),
      n = sum(complete.cases(df[, c("gbi","nmi")]))
    ), save_prov_path)
  }
  df
}

## 2) Clean helper: log1p skewed positives, leave others
smart_xform <- function(v) {
  if (all(v >= 0, na.rm = TRUE)) {
    x <- v[is.finite(v)]
    skewness <- if (length(x) < 20) 0 else e1071::skewness(x)
    if (skewness > 1) return(log1p(v))
  }
  v
}

## 3) Linear panel (ΔR2 + FDR)
# Linear panel with the same R² logic as logistic:
# R2_cov (covariates only), R2_gbi (GBI only), R2_base (cov+GBI), R2_full (cov+GBI+TILT),
# and deltas dR2_gbi (GBI beyond covariates) and dR2_tilt (TILT beyond GBI).
# Keeps tilt stats (beta/se/t/p) from the FULL model.
run_linear_panel <- function(df, y_vars, covars, use_adj_r2 = FALSE) {
  # pick R² type
  get_r2 <- function(m) {
    s <- summary(m)
    if (isTRUE(use_adj_r2)) s$adj.r.squared else s$r.squared
  }

  out <- lapply(y_vars, function(y){
    d <- df[, c(y, "gbi","nmi","tilt_fdg", covars)]
    d <- d[complete.cases(d), , drop = FALSE]
    # Recompute tilt on this exact analysis subset
    stopifnot(all(c("gbi","nmi") %in% names(d)))
    d$tilt_fdg <- as.numeric(scale(residuals(lm(nmi ~ gbi, data = d))))


    # optional mild transform for skewed-positive outcomes (same as before)
    if (all(d[[y]] >= 0, na.rm = TRUE) && requireNamespace("e1071", quietly = TRUE)) {
      x <- d[[y]][is.finite(d[[y]])]
      if (length(x) >= 20 && e1071::skewness(x) > 1) d[[y]] <- log1p(d[[y]])
    }

    f_cov  <- stats::as.formula(paste(y, "~", paste(covars, collapse = "+")))
    f_gbi  <- stats::as.formula(paste(y, "~ gbi"))
    f_base <- stats::as.formula(paste(y, "~", paste(covars, collapse = "+"), "+ gbi"))
    f_full <- stats::as.formula(paste(y, "~", paste(covars, collapse = "+"), "+ gbi + tilt_fdg"))

    m_cov  <- stats::lm(f_cov,  data = d)
    m_gbi  <- stats::lm(f_gbi,  data = d)
    m_base <- stats::lm(f_base, data = d)
    m_full <- stats::lm(f_full, data = d)

    R2_cov  <- get_r2(m_cov)
    R2_gbi  <- get_r2(m_gbi)
    R2_base <- get_r2(m_base)
    R2_full <- get_r2(m_full)
    dR2_gbi  <- R2_base - R2_cov          # GBI beyond covariates
    dR2_tilt <- R2_full - R2_base         # TILT beyond GBI (increment of interest)

    bs <- coef(summary(m_full))
    beta_TILT <- bs["tilt_fdg","Estimate"]; se_TILT <- bs["tilt_fdg","Std. Error"]
    t_TILT    <- bs["tilt_fdg","t value"];  p_TILT  <- bs["tilt_fdg","Pr(>|t|)"]
    beta_GBI  <- bs["gbi","Estimate"];      se_GBI  <- bs["gbi","Std. Error"]
    t_GBI     <- bs["gbi","t value"];       p_GBI   <- bs["gbi","Pr(>|t|)"]

    # effect-size extras for tilt (optional but handy)
    dfres <- stats::df.residual(m_full)
    partial_R2_TILT <- t_TILT^2 / (t_TILT^2 + dfres)
    f2_TILT         <- dR2_tilt / (1 - R2_full)
    q               <- stats::quantile(d$tilt_fdg, c(.10, .90), na.rm = TRUE)
    delta_10to90_TILT <- beta_TILT * as.numeric(q[2] - q[1])

    data.frame(
      outcome = y, n = nrow(d),
      R2_cov = R2_cov, R2_gbi = R2_gbi, R2_base = R2_base, R2_full = R2_full,
      dR2_gbi = dR2_gbi, dR2_tilt = dR2_tilt,
      beta_GBI = beta_GBI, se_GBI = se_GBI, t_GBI = t_GBI, p_GBI = p_GBI,
      beta_TILT = beta_TILT, se_TILT = se_TILT, t_TILT = t_TILT, p_TILT = p_TILT,
      partial_R2_TILT = partial_R2_TILT, f2_TILT = f2_TILT, delta_10to90_TILT = delta_10to90_TILT,
      stringsAsFactors = FALSE
    )
  })

  tab <- dplyr::bind_rows(out)
  if (nrow(tab)) {
    tab$FDR_TILT <- p.adjust(tab$p_TILT, method = "BH")
    tab$FDR_GBI  <- p.adjust(tab$p_GBI,  method = "BH")
  }
  tab[order(tab$FDR_TILT, tab$p_TILT), ]
}

## 4) Logistic panel (OR + ΔAUC optional)
# Bias-reduced GLM to avoid separation; returns R2 for cov, gbi, base (cov+gbi), full (+tilt),
# plus deltas (gbi beyond cov; tilt beyond gbi), and AUCs in the same structure.
run_logistic_panel <- function(df, y_vars, covars, auc = TRUE, r2_type = c("McFadden","Tjur")) {
  r2_type <- match.arg(r2_type)
  if (!requireNamespace("brglm2", quietly=TRUE)) stop("install.packages('brglm2')")
  if (!requireNamespace("pscl",   quietly=TRUE)) stop("install.packages('pscl')")
  if (isTRUE(auc) && !requireNamespace("pROC", quietly=TRUE)) stop("install.packages('pROC')")

  get_r2 <- function(m, y) {
    if (r2_type == "McFadden") as.numeric(pscl::pR2(m)["McFadden"]) else {
      p <- as.numeric(fitted(m)); y <- as.numeric(y)
      mean(p[y == 1L], na.rm=TRUE) - mean(p[y == 0L], na.rm=TRUE)  # Tjur's R2
    }
  }

  res_list <- lapply(stats::setNames(y_vars, y_vars), function(yname) {
    d <- df[, c(yname, "gbi","nmi", "tilt_fdg", covars)]
    d <- d[stats::complete.cases(d), , drop = FALSE]

    # Recompute tilt on this exact analysis subset
    stopifnot(all(c("gbi","nmi") %in% names(d)))
    d$tilt_fdg <- as.numeric(scale(residuals(stats::lm(nmi ~ gbi, data = d))))

    if (length(unique(d[[yname]])) < 2L) return(NULL)
    yv <- d[[yname]]

    f_cov  <- stats::as.formula(paste(yname, "~", paste(covars, collapse = "+")))
    f_gbi  <- stats::as.formula(paste(yname, "~ gbi"))
    f_base <- stats::as.formula(paste(yname, "~", paste(covars, collapse = "+"), "+ gbi"))
    f_full <- stats::as.formula(paste(yname, "~", paste(covars, collapse = "+"), "+ gbi + tilt_fdg"))

    m_cov  <- stats::glm(f_cov,  data = d, family = binomial(), method = brglm2::brglmFit)
    m_gbi  <- stats::glm(f_gbi,  data = d, family = binomial(), method = brglm2::brglmFit)
    m_base <- stats::glm(f_base, data = d, family = binomial(), method = brglm2::brglmFit)
    m_full <- stats::glm(f_full, data = d, family = binomial(), method = brglm2::brglmFit)

    # Tilt effect from full model (Wald)
    b  <- stats::coef(m_full)["tilt_fdg"]; se <- sqrt(diag(stats::vcov(m_full)))["tilt_fdg"]
    OR <- exp(b); z <- b / se; p <- 2 * stats::pnorm(-abs(z))
    CI_lo <- exp(b - 1.96 * se); CI_hi <- exp(b + 1.96 * se)
    
    # GBI effect from full model (for symmetry)
    bg  <- stats::coef(m_full)["gbi"]; seg <- sqrt(diag(stats::vcov(m_full)))["gbi"]
    OR_GBI  <- exp(bg); z_GBI <- bg / seg; p_GBI <- 2 * stats::pnorm(-abs(z_GBI))
    CIg_lo  <- exp(bg - 1.96 * seg); CIg_hi <- exp(bg + 1.96 * seg)

    # R2s
    R2_cov  <- get_r2(m_cov,  yv)
    R2_gbi  <- get_r2(m_gbi,  yv)
    R2_base <- get_r2(m_base, yv)
    R2_full <- get_r2(m_full, yv)
    dR2_gbi  <- R2_base - R2_cov           # GBI beyond covariates
    dR2_tilt <- R2_full - R2_base          # TILT beyond GBI

    # AUCs
    if (isTRUE(auc)) {
      p_cov  <- as.numeric(fitted(m_cov))
      p_base <- as.numeric(fitted(m_base))
      p_full <- as.numeric(fitted(m_full))
      roc_cov  <- pROC::roc(yv, p_cov,  quiet = TRUE)
      roc_base <- pROC::roc(yv, p_base, quiet = TRUE)
      roc_full <- pROC::roc(yv, p_full, quiet = TRUE)
      AUC_cov  <- as.numeric(roc_cov$auc)
      AUC_base <- as.numeric(roc_base$auc)
      AUC_full <- as.numeric(roc_full$auc)
      dAUC_gbi  <- AUC_base - AUC_cov
      dAUC_tilt <- AUC_full - AUC_base
    } else {
      AUC_cov <- AUC_base <- AUC_full <- dAUC_gbi <- dAUC_tilt <- NA_real_
    }

    dplyr::tibble(
      outcome = yname, n = nrow(d),
      OR_TILT = OR, CI_lo = CI_lo, CI_hi = CI_hi, z_TILT = z, p_TILT = p,
      OR_GBI  = OR_GBI, CIg_lo = CIg_lo, CIg_hi = CIg_hi, z_GBI = z_GBI, p_GBI = p_GBI,
      R2_cov = R2_cov, R2_gbi = R2_gbi, R2_base = R2_base, R2_full = R2_full,
      dR2_gbi = dR2_gbi, dR2_tilt = dR2_tilt,
      AUC_cov = AUC_cov, AUC_base = AUC_base, AUC_full = AUC_full,
      dAUC_gbi = dAUC_gbi, dAUC_tilt = dAUC_tilt
    )
  })

  tab <- dplyr::bind_rows(res_list)
  if (nrow(tab)) tab$FDR_TILT <- p.adjust(tab$p_TILT, method = "BH")
  dplyr::arrange(tab, FDR_TILT, p_TILT)
}


# ---------- Regional effects for arbitrary outcomes (map-ready) ----------
# Computes regionwise linear effects: z_region ~ outcome + covariates
# - Handles multiple outcomes at once (continuous or binary).
# - Uses flagship collapse: mean across CT/SA/SubVol (with optional SubVol p20 gate).
# - Returns one long table: outcome, region, beta, se, t, p, FDR, n_used.
#
# Dependencies: dplyr, tidyr, broom

# Helpers (local copies to avoid cross-file deps)
.pick_z_col <- function(zs) {
  if ("z_oriented" %in% names(zs)) "z_oriented" else
  if ("z" %in% names(zs)) "z" else stop("Need zs$z_oriented or zs$z")
}
.gate_subvol <- function(zs, size_df, pct = 0.20) {
  if (is.null(size_df)) return(zs)
  keep_sub <- size_df |>
    dplyr::filter(modality == "SubVol") |>
    dplyr::filter(size >= stats::quantile(size, pct, na.rm = TRUE)) |>
    dplyr::pull(region)
  zs |>
    dplyr::filter(modality != "SubVol" | region %in% keep_sub)
}
.collapse_zs <- function(zs, mods_use = c("CT","SA","SubVol")) {
  zcol <- .pick_z_col(zs)
  zs |>
    dplyr::filter(modality %in% mods_use, is.finite(.data[[zcol]])) |>
    dplyr::group_by(eid, region) |>
    dplyr::summarise(z = mean(.data[[zcol]], na.rm = TRUE), .groups = "drop")
}

# Main function
region_betas_for_outcomes <- function(zs, dt, y_vars,
                                      covars = c("age","sex","site","ICV","fasting","batch"),
                                      mods_use = c("CT","SA","SubVol"),
                                      size_df = NULL, subvol_gate_pct = 0.20,
                                      scale_y = TRUE) {
  stopifnot(all(c("eid", "region", "modality") %in% names(zs)))
  stopifnot(all(y_vars %in% names(dt)))

  # Flagship gating + collapse to (eid, region)
  zs2 <- .gate_subvol(zs, size_df, pct = subvol_gate_pct)
  Zc  <- .collapse_zs(zs2, mods_use = mods_use)

  out_list <- lapply(y_vars, function(y) {
    # Build analysis frame for this outcome
    req <- c("eid", y, covars)
    d   <- dt[, intersect(req, names(dt)), drop = FALSE]
    if (!("eid" %in% names(d))) stop("dt must contain 'eid'")
    d <- d[stats::complete.cases(d), , drop = FALSE]

    # Construct y_std (z-scale numeric; for 2-level factors/logicals, cast to 0/1 then z if requested)
    yy <- d[[y]]
    if (is.numeric(yy)) {
      y_std <- if (scale_y) as.numeric(scale(yy)) else yy
    } else if (is.logical(yy) || (is.factor(yy) && nlevels(yy) == 2)) {
      y01 <- if (is.logical(yy)) as.integer(yy) else as.integer(yy) - 1L
      y_std <- if (scale_y) as.numeric(scale(y01)) else y01
    } else {
      stop(sprintf("Outcome '%s' must be numeric or binary (logical/2-level factor).", y))
    }
    d$y_std <- y_std

    # Join with collapsed zs and fit regionwise models
    form <- stats::as.formula(paste("z ~ y_std",
                                    if (length(covars)) paste("+", paste(covars, collapse = "+")) else ""))

    dat <- dplyr::inner_join(Zc, d, by = "eid")
    dat <- dat[stats::complete.cases(dat[, c("z", "y_std", covars[covars %in% names(dat)])]), , drop = FALSE]

    res <- dat |>
      dplyr::group_by(region) |>
      do(broom::tidy(stats::lm(form, data = .))) |>
      dplyr::ungroup() |>
      dplyr::filter(term == "y_std") |>
      dplyr::transmute(outcome = y, region,
                       beta = estimate, se = std.error, t = statistic, p = p.value,
                       n_used = dplyr::n())  # n rows of full dat per region group; OK as an approximation

    # FDR within outcome
    res$FDR <- p.adjust(res$p, method = "BH")
    res
  })

  dplyr::bind_rows(out_list)
}

