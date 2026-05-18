#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
})

# ---- Orientation (sign convention: higher z_oriented = healthier) --------
# Modalities not listed get orientation 0 and are dropped in compute_oriented().
get_orientation <- function(zs) {
  dplyr::case_when(
    zs$modality %in% c("CT","SA","SubVol","FA") ~  +1,
    zs$modality %in% c("MD","RD","WMH")         ~  -1,
    TRUE                                         ~   0
  )
}

# ---- Sign-aware oriented z ------------------------------------------------
# Rows from modalities with no defined orientation (s == 0) are dropped here
# so they cannot dilute Dm with zero-weighted contributions downstream.
compute_oriented <- function(zs) {
  zs %>%
    mutate(s = get_orientation(.), z_oriented = s * z) %>%
    filter(s != 0)
}

# ---- Prepare PET weights --------------------------------------------------
# If pet_w is NULL -> w_pet = 1 (uniform; used for GBI).
prepare_pet_weights <- function(zs, pet_w = NULL) {
  if (is.null(pet_w)) {
    zs %>% distinct(region) %>% mutate(w_pet = 1.0)
  } else {
    pet_w %>%
      transmute(region, w_raw = pmax(0, as.numeric(weight))) %>%
      mutate(w_pet = ifelse(w_raw > 0, w_raw / mean(w_raw[w_raw > 0], na.rm = TRUE), 0))
  }
}

# ---- Prepare size weights -------------------------------------------------
prepare_size_weights <- function(zs, size_df = NULL) {
  if (is.null(size_df)) {
    zs %>% distinct(region, modality) %>% mutate(size = 1.0)
  } else {
    size_df %>% mutate(size = as.numeric(size))
  }
}

# ---- Main computation -----------------------------------------------------
# Flagship settings:
#  - PET weights mean-normalised, p99-winsorised, then exponentiated by gamma = 2.5
#  - Size weights with beta = 0.25
#  - Equal alpha per available modality per subject (unless alpha supplied)
#  - coverage_thr = 0.80
#
# The p99 winsorisation of w_pet (before gamma) limits dominance by the most
# metabolically active region(s) without imposing the hard plateau a per-region
# cap on w_bar would produce. Empirically near-identical to no winsor under
# typical FDG distributions but provides a principled safeguard.
compute_mbi_flagship <- function(zs, pet_w = NULL, size_df = NULL, alpha = NULL, tune = NULL) {

  stopifnot(all(c("id","modality","region","z") %in% names(zs)))

  cfg <- list(
    use_global_norm   = TRUE,
    gamma             = 2.5,
    use_size          = TRUE,
    beta              = 0.25,
    per_modality_norm = FALSE,
    winsor_q          = 0.99,
    coverage_thr      = 0.80,
    debug             = FALSE
  )
  if (!is.null(tune) && length(tune)) {
    for (k in intersect(names(tune), names(cfg))) cfg[[k]] <- tune[[k]]
  }

  #--- orient z-scores (drops modalities with s == 0)
  zo <- compute_oriented(zs)

  #--- PET weights: prep, winsor at p99 of positives, normalize, exponentiate
  w_pet0 <- prepare_pet_weights(zs, pet_w) %>%
    mutate(w_pet = ifelse(is.finite(w_pet) & w_pet > 0, w_pet, 0))

  if (!is.null(cfg$winsor_q) && is.finite(cfg$winsor_q) &&
      cfg$winsor_q > 0 && cfg$winsor_q < 1) {
    pos <- w_pet0$w_pet[w_pet0$w_pet > 0]
    if (length(pos) > 1) {
      qhi <- as.numeric(quantile(pos, cfg$winsor_q, na.rm = TRUE))
      w_pet0$w_pet <- pmin(w_pet0$w_pet, qhi)
    }
  }

  if (cfg$use_global_norm) {
    m <- mean(w_pet0$w_pet[w_pet0$w_pet > 0], na.rm = TRUE)
    if (is.finite(m) && m > 0) w_pet0$w_pet <- w_pet0$w_pet / m
  }

  if (!is.null(cfg$gamma) && is.finite(cfg$gamma) && cfg$gamma != 1)
    w_pet0$w_pet <- w_pet0$w_pet ^ cfg$gamma

  #--- size weights (sublinear scaling)
  w_size0 <- prepare_size_weights(zs, size_df)
  if (cfg$use_size && !is.null(cfg$beta) && is.finite(cfg$beta) && cfg$beta != 1)
    w_size0 <- w_size0 %>% mutate(size = pmax(1e-12, size) ^ cfg$beta)

  #--- merge + combined weights
  zw <- zo %>%
    inner_join(w_pet0,  by = "region") %>%
    inner_join(w_size0, by = c("region","modality")) %>%
    mutate(w_comb = pmax(0, w_pet) * pmax(1e-12, size))

  #--- optional per-modality normalisation (cosmetic for w_bar display; does
  #--- NOT affect the subject-level score because Dm is a weighted mean and
  #--- a within-modality constant rescaling cancels in num/den).
  if (isTRUE(cfg$per_modality_norm)) {
    zw <- zw %>% group_by(modality) %>%
      mutate(w_bar = w_comb / mean(w_comb[w_comb > 0], na.rm = TRUE)) %>%
      ungroup()
  } else {
    zw <- zw %>% mutate(w_bar = w_comb)
  }

  #--- modality-level weighted averages
  Dm <- zw %>%
    group_by(id, modality) %>%
    summarise(
      num = sum(w_bar * z_oriented, na.rm = TRUE),
      den = sum(w_bar[!is.na(z_oriented)], na.rm = TRUE),
      Dm  = ifelse(den > 0, num / den, NA_real_),
      .groups = "drop"
    )

  #--- coverage rule (>= coverage_thr valid parcels)
  cover <- zw %>%
    group_by(id, modality) %>%
    summarise(cover = mean(!is.na(z_oriented), na.rm = TRUE), .groups = "drop")

  Dm <- Dm %>%
    left_join(cover, by = c("id","modality")) %>%
    mutate(Dm = ifelse(cover >= cfg$coverage_thr, Dm, NA_real_))

  #--- alpha weights: equal per available modality unless supplied
  alpha_tbl <- if (is.null(alpha)) {
    Dm %>%
      filter(!is.na(Dm)) %>%
      group_by(id) %>% mutate(alpha = 1 / n()) %>% ungroup() %>%
      select(id, modality, alpha)
  } else {
    Dm %>%
      mutate(alpha = alpha[modality]) %>%
      group_by(id) %>%
      mutate(alpha = alpha / sum(alpha, na.rm = TRUE)) %>%
      ungroup() %>%
      select(id, modality, alpha)
  }

  Dm <- Dm %>% left_join(alpha_tbl, by = c("id","modality"))

  #--- subject-level aggregation + z-standardisation
  D <- Dm %>%
    filter(!is.na(Dm), !is.na(alpha)) %>%
    group_by(id) %>%
    summarise(D = sum(alpha * Dm), .groups = "drop")

  D_mean <- mean(D$D, na.rm = TRUE)
  D_sd   <- sd(D$D,   na.rm = TRUE)

  MBI <- D %>%
    mutate(MBI = (D - D_mean) / D_sd) %>%
    arrange(id)

  #--- QC outputs
  qc_modality <- Dm %>%
    filter(!is.na(Dm), !is.na(alpha)) %>%
    group_by(modality) %>%
    summarise(mean_Dm = mean(Dm), sd_Dm = sd(Dm), .groups = "drop")

  weight_table <- zw %>%
    distinct(region, modality, w_pet, size, w_bar) %>%
    arrange(modality, region)

  diag <- NULL
  if (isTRUE(cfg$debug)) {
    join_stats <- zw %>%
      group_by(modality) %>%
      summarise(
        n        = n(),
        nz       = sum(w_bar > 0, na.rm = TRUE),
        frac_zero = 1 - nz / n,
        .groups  = "drop"
      )
    cover_stats <- cover %>%
      group_by(modality) %>%
      summarise(
        med_cover  = median(cover, na.rm = TRUE),
        p10        = quantile(cover, 0.10, na.rm = TRUE),
        p90        = quantile(cover, 0.90, na.rm = TRUE),
        pct_ge_thr = mean(cover >= cfg$coverage_thr, na.rm = TRUE),
        .groups    = "drop"
      )
    bad_subjects <- cover %>%
      group_by(id) %>%
      summarise(
        min_cover = min(cover, na.rm = TRUE),
        mods      = paste(modality[order(cover)], collapse = ","),
        .groups   = "drop"
      ) %>%
      arrange(min_cover) %>% head(10)
    diag <- list(join_stats = join_stats, cover_stats = cover_stats, bad_subjects = bad_subjects)
  }

  list(MBI = MBI, modality_scores = Dm, weight_table = weight_table,
       qc_modality = qc_modality, diag = diag)
}


# ---- Wrapper: yields GBI(alpha,beta), MBI_raw(alpha,beta,gamma) and MBI --
compute_flagship_abg <- function(
  zs, fdg, sizes,
  beta = 0.25, gamma = 2.5,
  alpha = NULL,
  winsor_q     = 0.99,
  coverage_thr = 0.80
) {
  # MBI_raw: PET-weighted
  mbi_out <- compute_mbi_flagship(
    zs      = zs,
    pet_w   = fdg,
    size_df = sizes,
    alpha   = alpha,
    tune    = list(beta = beta, gamma = gamma,
                   winsor_q = winsor_q, coverage_thr = coverage_thr)
  )
  mbi <- dplyr::rename(mbi_out$MBI, mbi_raw = MBI)

  # GBI: size-weighted only (PET weights = 1; gamma & winsor_q irrelevant when all w_pet = 1)
  gbi_out <- compute_mbi_flagship(
    zs      = zs,
    pet_w   = NULL,
    size_df = sizes,
    alpha   = alpha,
    tune    = list(beta = beta, coverage_thr = coverage_thr)
  )
  gbi <- dplyr::rename(gbi_out$MBI, gbi = MBI)

  # MBI: residual of MBI_raw ~ GBI (z-standardised)
  dt <- dplyr::inner_join(mbi, gbi, by = "id")
  dt$mbi <- as.numeric(scale(residuals(stats::lm(mbi_raw ~ gbi, data = dt))))

  list(
    scores      = dt[, c("id","mbi_raw","gbi","mbi")],
    mbi_details = mbi_out,
    gbi_details = gbi_out
  )
}
