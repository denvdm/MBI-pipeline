# Functions for relating brain-wide deviation patterns to PET metabolism maps.
# Provides correlation, alignment, and similarity measures between PET and MRI-derived metrics.

pacman::p_load(dplyr,data.table,broom,tidyr)

# --- column pickers / guards --------------------------------------------------
pick_z_col <- function(zs) {
  if ("z_oriented" %in% names(zs)) "z_oriented" else
    if ("z" %in% names(zs)) "z" else stop("Need zs$z_oriented or zs$z")
}

prep_fdg <- function(fdg) {
  if ("w_pet" %in% names(fdg)) fdg %>% select(region, w_pet)
  else if ("weight" %in% names(fdg)) fdg %>% transmute(region, w_pet = weight)
  else stop("fdg must have column 'w_pet' or 'weight'")
}

# --- Fisher z CIs -------------------------------------------------------------
fisher_ci <- function(r, n, conf = 0.95) {
  if (!is.finite(r) || n <= 3) return(c(lo = NA_real_, hi = NA_real_))
  r <- max(min(r, 0.999999), -0.999999); z <- atanh(r); se <- 1/sqrt(n - 3); zc <- qnorm(1 - (1-conf)/2)
  c(lo = tanh(z - zc*se), hi = tanh(z + zc*se))
}

fisher_mean_ci <- function(r_vec, n_vec, conf = 0.95) {
  ok <- is.finite(r_vec) & is.finite(n_vec) & (n_vec > 3); if (sum(ok) < 2) return(c(mean=NA, lo=NA, hi=NA))
  z <- atanh(pmin(pmax(r_vec[ok], -0.999999), 0.999999)); w <- pmax(n_vec[ok]-3, 1e-8)
  mz <- sum(w*z)/sum(w); se <- 1/sqrt(sum(w)); zc <- qnorm(1 - (1-conf)/2)
  c(mean = tanh(mz), lo = tanh(mz - zc*se), hi = tanh(mz + zc*se))
}

# --- SubVol gate + collapse to (eid,region) once ------------------------------
gate_subvol <- function(zs, size_df, pct = 0.20) {
  keep_sub <- size_df %>% filter(modality=="SubVol") %>% filter(size >= quantile(size, pct, na.rm=TRUE)) %>% pull(region)
  zs %>% filter(modality!="SubVol" | region %in% keep_sub)
}

collapse_zs <- function(zs, fdg_use, mods_use = c("CT","SA","SubVol")) {
  zcol <- pick_z_col(zs)
  zs %>% filter(modality %in% mods_use, is.finite(.data[[zcol]])) %>%
    group_by(eid, region) %>% summarise(z = mean(.data[[zcol]], na.rm=TRUE), .groups="drop") %>%
    semi_join(fdg_use, by="region")
}

# --- PC1 loadings × PET -------------------------------------------------------
pc_pet_alignment <- function(Z_collapsed, fdg_use, mods_tag) {
  M <- Z_collapsed %>% pivot_wider(names_from=region, values_from=z) %>% arrange(eid) %>% select(-eid) %>% as.matrix()
  if (!ncol(M)) stop("No regions after collapse")
  vars <- apply(M, 2, sd, na.rm=TRUE); keep <- is.finite(vars) & vars > 0
  if (sum(keep) < 2) stop("Too few variable regions for PCA")
  X <- scale(M[, keep, drop=FALSE]); X[is.na(X)] <- 0
  p <- prcomp(X, center=FALSE, scale.=FALSE)
  pc1_tbl <- tibble(region = colnames(X), loading = as.numeric(p$rotation[,1])) %>% inner_join(fdg_use, by="region")
  rP <- suppressWarnings(cor(pc1_tbl$loading, pc1_tbl$w_pet, use="pairwise"))
  rS <- suppressWarnings(cor(pc1_tbl$loading, pc1_tbl$w_pet, method="spearman", use="pairwise"))
  n  <- nrow(pc1_tbl); ciP <- fisher_ci(rP, n); ciS <- fisher_ci(rS, n)
  data.frame(metric = paste0("PCA_loadings_vs_PET_", mods_tag),
             r_Pearson=rP, rP_lo=ciP["lo"], rP_hi=ciP["hi"],
             r_Spearman=rS, rS_lo=ciS["lo"], rS_hi=ciS["hi"], n_regions=n)
}

# --- LOO contribution map × PET ----------------------------------------------
loo_pet_alignment <- function(Z_collapsed, fdg_use) {
  sub_sum <- Z_collapsed %>% group_by(eid) %>% summarise(S = sum(z), P = n(), .groups="drop")
  Z2 <- Z_collapsed %>% inner_join(sub_sum, by="eid") %>% mutate(GI_minus_r = (S - z)/pmax(P - 1L, 1L))
  betas <- Z2 %>% group_by(region) %>% do(tidy(lm(GI_minus_r ~ z, data = .))) %>%
    ungroup() %>% filter(term=="z") %>% transmute(region, beta_LOO = estimate)
  beta_pet <- betas %>% inner_join(fdg_use, by="region")
  rP <- suppressWarnings(cor(beta_pet$beta_LOO, beta_pet$w_pet, use="pairwise"))
  rS <- suppressWarnings(cor(beta_pet$beta_LOO, beta_pet$w_pet, method="spearman", use="pairwise"))
  n  <- nrow(beta_pet); ciP <- fisher_ci(rP, n); ciS <- fisher_ci(rS, n)
  data.frame(metric="LOO_beta_vs_PET", r_Pearson=rP, rP_lo=ciP["lo"], rP_hi=ciP["hi"],
             r_Spearman=rS, rS_lo=ciS["lo"], rS_hi=ciS["hi"], n_regions=n)
}

# --- Subject-level map×PET correlation summary --------------------------------
subject_pet_alignment <- function(Z_collapsed, fdg_use) {
  r_df <- Z_collapsed %>% inner_join(fdg_use, by="region") %>%
    group_by(eid) %>%
    summarise(r_s = suppressWarnings(cor(z, w_pet, method="spearman", use="pairwise")),
              r_p = suppressWarnings(cor(z, w_pet, method="pearson",  use="pairwise")),
              n_reg = n(), .groups="drop")
  ci_s <- fisher_mean_ci(r_df$r_s, r_df$n_reg); ci_p <- fisher_mean_ci(r_df$r_p, r_df$n_reg)
  data.frame(n_subjects = sum(is.finite(r_df$r_s)),
             mean_s = unname(ci_s["mean"]), lo_s = unname(ci_s["lo"]), hi_s = unname(ci_s["hi"]),
             mean_p = unname(ci_p["mean"]), lo_p = unname(ci_p["lo"]), hi_p = unname(ci_p["hi"]))
}

# --- Writers ------------------------------------------------------------------
write_alignment_rows <- function(rows_df, out_path) {
  dir.create(dirname(out_path), showWarnings=FALSE, recursive=TRUE)
  have <- file.exists(out_path)
  write.table(rows_df, file=out_path, sep="\t", row.names=FALSE, col.names=!have, append=have)
}

write_subject_summary <- function(row_df, out_path) {
  dir.create(dirname(out_path), showWarnings=FALSE, recursive=TRUE)
  write.table(row_df, file=out_path, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
}



# --- PC1 loadings table (optionally sign-align to PET) ------------------------
pc1_loadings_table <- function(Z_collapsed, fdg_use = NULL, align_to_pet = TRUE) {
  M <- Z_collapsed |>
    tidyr::pivot_wider(names_from = region, values_from = z) |>
    dplyr::arrange(eid) |>
    dplyr::select(-eid) |>
    as.matrix()
  if (!ncol(M)) stop("No regions after collapse")
  vars <- apply(M, 2, sd, na.rm = TRUE)
  keep <- is.finite(vars) & vars > 0
  if (sum(keep) < 2) stop("Too few variable regions for PCA")
  X <- scale(M[, keep, drop = FALSE]); X[is.na(X)] <- 0
  p <- prcomp(X, center = FALSE, scale. = FALSE)

  load <- as.numeric(p$rotation[, 1])
  reg  <- colnames(X)

  # optional: flip sign so correlation with PET is non-negative
  if (!is.null(fdg_use)) {
    tmp <- dplyr::tibble(region = reg, loading = load) |>
      dplyr::inner_join(fdg_use, by = "region")
    rP <- suppressWarnings(cor(tmp$loading, tmp$w_pet, use = "pairwise"))
    if (is.finite(rP) && align_to_pet && rP < 0) load <- -load
  }

  var_expl_pc1 <- unname(p$sdev[1]^2 / sum(p$sdev^2))
  dplyr::tibble(region = reg, pc1_loading = load, var_explained_pc1 = var_expl_pc1)
}

# --- Mean z per region (collapsed across modalities & subjects) ---------------
mean_z_per_region <- function(Z_collapsed) {
  Z_collapsed |>
    dplyr::group_by(region) |>
    dplyr::summarise(mean_z = mean(z, na.rm = TRUE), .groups = "drop")
}

# --- LOO beta map (region contribution slope) ---------------------------------
loo_beta_map <- function(Z_collapsed) {
  sub_sum <- Z_collapsed |>
    dplyr::group_by(eid) |>
    dplyr::summarise(S = sum(z), P = dplyr::n(), .groups = "drop")
  Z2 <- Z_collapsed |>
    dplyr::inner_join(sub_sum, by = "eid") |>
    dplyr::mutate(GI_minus_r = (S - z) / pmax(P - 1L, 1L))

  betas <- Z2 |>
    dplyr::group_by(region) |>
    do(broom::tidy(lm(GI_minus_r ~ z, data = .))) |>
    dplyr::ungroup() |>
    dplyr::filter(term == "z") |>
    dplyr::transmute(region, beta_LOO = estimate, se = std.error, t = statistic, p = p.value)
  betas$FDR <- p.adjust(betas$p, method = "BH")
  betas
}

