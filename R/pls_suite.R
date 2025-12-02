# PLS utilities for multivariate associations between MBI metrics and outcome sets.
# Includes fit, permutation testing, component extraction, and summary tools.

# ========================== pls_suite.R =======================================
# Minimal deps: pls (and optionally splines if you enable spline age in regression)

# Fit PLS2 on metabolite panel (X) and Y = c(NMI_unw, PET_tilt). Return model + PET_tilt VIP.
run_pls_panel <- function(data, metab_cols, tilt_col,
                          unweighted_col = "nmi_unw",
                          ncomp = 3, scaleX = TRUE) {
  if (!requireNamespace("pls", quietly = TRUE)) stop("Please install.packages('pls')")
  stopifnot(all(c(unweighted_col, tilt_col) %in% names(data)))
  stopifnot(all(metab_cols %in% names(data)))
  rows <- stats::complete.cases(data[, c(metab_cols, unweighted_col, tilt_col)])
  if (!any(rows)) stop("No complete rows for X+Y")
  X <- as.matrix(data[rows, metab_cols, drop = FALSE])
  if (scaleX) X <- scale(X)
  Y <- cbind(NMI_unw = data[[unweighted_col]][rows],
             PET_tilt = data[[tilt_col]][rows])
  colnames(Y) <- c("NMI_unw","PET_tilt")

  # PLS2 with CV
  fit <- pls::plsr(Y ~ X, ncomp = ncomp, validation = "CV")

  # VIP for PET_tilt (response index = 2), standard formulation
  T  <- pls::scores(fit)[, 1:ncomp, drop = FALSE]                    # n x A
  Qy <- as.numeric(fit$Yloadings[2, 1:ncomp, drop = FALSE])          # length A
  W  <- fit$loading.weights[, 1:ncomp, drop = FALSE]                 # p x A
  SSY <- colSums((T * rep(Qy, each = nrow(T)))^2)                    # Y-var per comp
  if (!all(is.finite(SSY)) || sum(SSY) == 0) SSY <- rep(1, length(SSY))
  Wsq <- W^2
  Wsq_norm <- sweep(Wsq, 2, colSums(Wsq), "/")
  vip_pet <- sqrt(ncol(X) * as.numeric(Wsq_norm %*% (SSY / sum(SSY))))
  names(vip_pet) <- colnames(X)

  list(model = fit,
       rows_used = which(rows),
       vip_pet = sort(vip_pet, decreasing = TRUE))
}

# Extract cumulative R² / Q² for PET_tilt across 1..A components
# Extract cumulative R² / Q² for PET_tilt across 1..A components (robust to 2D/3D)
pls_r2q2 <- function(pls_fit, response = "PET_tilt") {
  read_r2 <- function(r2obj) {
    arr <- r2obj$val
    dn  <- dimnames(arr)
    # Find the "R2" row name robustly
    stat_names <- if (length(dn) >= 1) dn[[1]] else rownames(arr)
    r_idx <- which(tolower(stat_names) == "r2")
    if (length(r_idx) == 0L) r_idx <- 1L  # fallback

    if (length(dim(arr)) == 2L) {
      # shape: [stat x comps]
      as.numeric(arr[r_idx, , drop = TRUE])
    } else {
      # shape: [stat x comps x responses]
      # pick response by name; if not found, take the first
      resp_names <- dn[[3]]
      ridx <- if (!is.null(resp_names)) which(resp_names %in% response) else integer(0)
      if (length(ridx) == 0L) ridx <- 1L
      as.numeric(arr[r_idx, , ridx, drop = TRUE])
    }
  }
  R2Y <- read_r2(pls::R2(pls_fit, estimate = "train", intercept = TRUE))
  Q2Y <- read_r2(pls::R2(pls_fit, estimate = "CV",    intercept = TRUE))
  list(R2Y = R2Y, Q2Y = Q2Y)
}

# Choose number of components A by Q² peak (fallback to 1 if all <=0)
pls_select_A <- function(Q2Y) {
  if (all(!is.finite(Q2Y))) return(1L)
  if (all(Q2Y <= 0, na.rm = TRUE)) return(1L)
  which.max(Q2Y)
}

# Compute ΔR² of PET_tilt explained by PLS scores beyond covariates (regression)
pls_deltaR2_from_scores <- function(data, pls_obj, A, tilt_col,
                                    covars = NULL, use_spline_age = FALSE) {
  stopifnot(tilt_col %in% names(data))
  T_scores <- as.data.frame(pls::scores(pls_obj)[, 1:A, drop = FALSE])
  names(T_scores) <- paste0("PLS", seq_len(A))
  df <- cbind(data, T_scores)
  # model formulas
  rhs_cov <- if (length(covars)) paste(covars, collapse = " + ") else "1"
  if (use_spline_age && "age" %in% covars) {
    if (!requireNamespace("splines", quietly = TRUE)) stop("Please install/load 'splines'")
    rhs_cov <- sub("\\bage\\b", "splines::ns(age,3)", rhs_cov)
  }
  f0 <- stats::as.formula(paste(tilt_col, "~", rhs_cov))
  f1 <- stats::as.formula(paste(tilt_col, "~", rhs_cov, "+", paste(names(T_scores), collapse = " + ")))
  # fit on complete cases
  mm0 <- stats::model.frame(f0, data = df)
  mm1 <- stats::model.frame(f1, data = df)
  idx <- stats::complete.cases(mm1)   # ensures both models use same rows
  m0 <- stats::lm(f0, data = df[idx, , drop = FALSE])
  m1 <- stats::lm(f1, data = df[idx, , drop = FALSE])
  R2_0 <- summary(m0)$r.squared
  R2_1 <- summary(m1)$r.squared
  list(R2_0 = R2_0, R2_1 = R2_1, dR2 = R2_1 - R2_0, n = nrow(mm1[idx, , drop = FALSE]))
}

# Aggregate VIP to family pathways (optional)
pls_family_vip <- function(vip_pet, annotation_df,
                           metabolite_col = "metabolite",
                           family_col = "family") {
  stopifnot(all(c(metabolite_col, family_col) %in% names(annotation_df)))
  tmp <- data.frame(metabolite = names(vip_pet), VIP = as.numeric(vip_pet))
  names(tmp)[1] <- metabolite_col
  out <- merge(tmp, annotation_df[, c(metabolite_col, family_col)], by = metabolite_col, all.x = TRUE)
  agg <- stats::aggregate(list(VIP2_sum = out$VIP^2, mean_VIP = out$VIP),
                          by = list(family = out[[family_col]]),
                          FUN = function(x) c(sum = sum(x), mean = mean(x)))
  # unpack the list-cols
  agg$VIP2_sum  <- vapply(agg$VIP2_sum,  function(z) z["sum"],  numeric(1))
  agg$mean_VIP  <- vapply(agg$mean_VIP,  function(z) z["mean"], numeric(1))
  agg$n         <- as.integer(table(out[[family_col]])[as.character(agg$family)])
  agg[order(-agg$VIP2_sum), c("family","n","VIP2_sum","mean_VIP")]
}

# One-call wrapper: fit, pick A, compute R2/Q2, ΔR² (covariate-adjusted), family summary
pls_quick <- function(data, metab_cols, tilt_col,
                      unweighted_col = "nmi_unw",
                      covars = NULL, use_spline_age = FALSE,
                      max_comp = 3, scaleX = TRUE,
                      family_annotation = NULL) {
  fit <- run_pls_panel(data, metab_cols, tilt_col,
                       unweighted_col = unweighted_col,
                       ncomp = max_comp, scaleX = scaleX)
  r2q2 <- pls_r2q2(fit$model, "PET_tilt")
  A    <- pls_select_A(r2q2$Q2Y)
  dR2  <- pls_deltaR2_from_scores(data, fit$model, A, tilt_col, covars, use_spline_age)
  fam  <- if (!is.null(family_annotation))
            pls_family_vip(fit$vip_pet, family_annotation) else NULL
  summary <- list(
    A_star = A,
    R2Y_cum = r2q2$R2Y[A],
    Q2Y_cum = r2q2$Q2Y[A],
    dR2_metab = dR2$dR2,
    n_reg = dR2$n
  )
  # concise console message
  msg <- sprintf("PLS: A=%d | R2Y=%.3f | Q2Y=%.3f | ΔR2_metab=%.4f (n=%d)",
                 A, summary$R2Y_cum %||% NA_real_, summary$Q2Y_cum %||% NA_real_,
                 summary$dR2_metab %||% NA_real_, summary$n_reg %||% NA_integer_)
  message(msg)
  list(model = fit$model,
       rows_used = fit$rows_used,
       vip_pet = fit$vip_pet,
       R2Y = r2q2$R2Y, Q2Y = r2q2$Q2Y,
       A_star = A,
       dR2 = dR2,
       family_summary = fam,
       summary = summary)
}

# small infix for pretty printing (base-only fallback)
`%||%` <- function(a, b) if (!is.null(a) && length(a) && all(is.finite(a))) a else b
# ============================================================================
