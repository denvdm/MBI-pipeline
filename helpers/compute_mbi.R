#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
})

get_orientation <- function(zs) {
  s <- dplyr::case_when(
    zs$modality %in% c("CT","SA","SubVol","FA") ~  +1,
    zs$modality %in% c("MD","RD","WMH") ~  -1,
    zs$modality %in% c("AD","T1intensity","T2intensity","MyelinProxy","Gradient") ~ 0,
    TRUE ~ 0
  )
  vent_pat <- "(ventricle|ventricular|lat_vent)"
  is_vent  <- grepl(vent_pat, tolower(zs$region))
  s[is_vent] <- -1
  s
}

# ---- Prepare PET weights ------------------------------------------------------
# If pet_w is NULL → w_pet = 1 (simulated).
prepare_pet_weights <- function(zs, pet_w = NULL) {
  if (is.null(pet_w)) {zs %>% distinct(region) %>% mutate(w_pet = 1.0)
  } else {pet_w %>% transmute(region, w_raw = pmax(0, as.numeric(weight))) %>%
      mutate(w_pet = ifelse(w_raw > 0, w_raw / mean(w_raw[w_raw > 0], na.rm = TRUE), 0))}
}

# ---- Prepare size weights -----------------------------------------------------
# If NULL, default size = 1.
prepare_size_weights <- function(zs, size_df = NULL) {
  if (is.null(size_df)) {zs %>% distinct(region, modality) %>% mutate(size = 1.0)
  } else {size_df %>% mutate(size = as.numeric(size))}
}

# ---- Compute sign-aware oriented z (no truncation) ----------------------------
compute_oriented <- function(zs) {zs %>% mutate(s = get_orientation(.), z_oriented = s * z)}


# Defaults (flagship settings):
# - PET weights globally normalised, gamma = 2.5
# - Size weights on, beta = 0.25
# - No PET-based alpha (use_alpha_pet = FALSE)
# - No per-modality normalisation
# - cap_frac = 0.05, coverage_thr = 0.80
# - debug = FALSE unless explicitly turned on

compute_mbi_flagship <- function(zs, pet_w = NULL, size_df = NULL, alpha = NULL, tune = NULL) {

  stopifnot(all(c("id","modality","region","z") %in% names(zs)))

  cfg <- list(
    use_global_norm   = TRUE,
    gamma             = 2.5,
    use_size          = TRUE,
    beta              = 0.25,
    per_modality_norm = FALSE,
    cap_frac          = 0.05,
    coverage_thr      = 0.80,
    use_alpha_pet     = FALSE,
    top_pet_quantile  = NULL,
    debug             = FALSE
  )

  if (!is.null(tune) && length(tune)) {
    if (!is.null(tune$gamma))        cfg$gamma        <- tune$gamma
    if (!is.null(tune$beta))         cfg$beta         <- tune$beta
    if (!is.null(tune$cap_frac))     cfg$cap_frac     <- tune$cap_frac
    if (!is.null(tune$coverage_thr)) cfg$coverage_thr <- tune$coverage_thr
    if (!is.null(tune$debug))        cfg$debug        <- tune$debug
  }

  #--- orient z-scores (higher=better)
  zo <- compute_oriented(zs)

  #--- PET weights: prep + clamp + normalize + exponentiate
  w_pet0 <- prepare_pet_weights(zs, pet_w)
  if (!is.null(cfg$top_pet_quantile)) {
    qv <- quantile(w_pet0$w_pet[w_pet0$w_pet > 0], probs = cfg$top_pet_quantile, na.rm = TRUE)
    w_pet0 <- w_pet0 %>% mutate(w_pet = ifelse(w_pet >= qv, w_pet, 0))
  }
  w_pet0 <- w_pet0 %>%
    mutate(w_pet = ifelse(is.finite(w_pet) & w_pet > 0, w_pet, 0))

  if (cfg$use_global_norm) {
    m <- mean(w_pet0$w_pet[w_pet0$w_pet > 0], na.rm = TRUE)
    if (is.finite(m) && m > 0) w_pet0$w_pet <- w_pet0$w_pet / m
  }

  if (!is.null(cfg$gamma) && is.finite(cfg$gamma) && cfg$gamma != 1)
    w_pet0$w_pet <- w_pet0$w_pet ^ cfg$gamma

  #--- size weights (optional sublinear scaling)
  w_size0 <- prepare_size_weights(zs, size_df)
  if (cfg$use_size && !is.null(cfg$beta) && is.finite(cfg$beta) && cfg$beta != 1)
    w_size0 <- w_size0 %>% mutate(size = pmax(1e-12, size) ^ cfg$beta)

  #--- merge all + compute combined weights
  zw <- zo %>%
    inner_join(w_pet0,  by = "region") %>%
    inner_join(w_size0, by = c("region","modality")) %>%
    mutate(w_comb = pmax(0, w_pet) * pmax(1e-12, size))

  #--- optional per-modality normalization
  if (isTRUE(cfg$per_modality_norm)) {
    zw <- zw %>% group_by(modality) %>%
      mutate(w_bar = w_comb / mean(w_comb[w_comb > 0], na.rm = TRUE)) %>%
      ungroup()
  } else {
    zw <- zw %>% mutate(w_bar = w_comb)
  }

  #--- optional capping (avoid single-region dominance)
  if (!is.null(cfg$cap_frac) && is.finite(cfg$cap_frac) && cfg$cap_frac > 0) {
    zw <- zw %>%
      group_by(modality) %>%
      mutate(w_bar = pmin(w_bar, cfg$cap_frac * sum(w_bar, na.rm = TRUE))) %>%
      ungroup()
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

  #--- coverage rule (≥ coverage_thr valid parcels)
  cover <- zw %>%
    group_by(id, modality) %>%
    summarise(cover = mean(!is.na(z_oriented), na.rm = TRUE), .groups = "drop")

  Dm <- Dm %>%
    left_join(cover, by = c("id","modality")) %>%
    mutate(Dm = ifelse(cover >= cfg$coverage_thr, Dm, NA_real_))

  #--- alpha weights (currently: equal per available modality unless alpha supplied)
  if (isTRUE(cfg$use_alpha_pet) && is.null(alpha)) {
    pet_load <- zw %>%
      group_by(id, modality) %>%
      summarise(load = sum(w_bar, na.rm = TRUE), .groups = "drop")

    alpha_tbl <- Dm %>%
      select(id, modality) %>%
      left_join(pet_load, by = c("id","modality")) %>%
      group_by(id) %>%
      mutate(alpha = ifelse(is.finite(load) & load > 0, load / sum(load, na.rm = TRUE), NA_real_)) %>%
      ungroup() %>%
      select(id, modality, alpha)

  } else if (is.null(alpha)) {
    alpha_tbl <- Dm %>%
      filter(!is.na(Dm)) %>%
      group_by(id) %>%
      mutate(alpha = 1 / n()) %>%
      ungroup() %>%
      select(id, modality, alpha)

  } else {
    alpha_tbl <- Dm %>%
      mutate(alpha = alpha[modality]) %>%
      group_by(id) %>%
      mutate(alpha = alpha / sum(alpha, na.rm = TRUE)) %>%
      ungroup() %>%
      select(id, modality, alpha)
  }

  Dm <- Dm %>% left_join(alpha_tbl, by = c("id","modality"))

  #--- subject-level aggregation + z-standardization
  D <- Dm %>%
    filter(!is.na(Dm), !is.na(alpha)) %>%
    group_by(id) %>%
    summarise(D = sum(alpha * Dm), .groups = "drop")

  D_mean <- mean(D$D, na.rm = TRUE)
  D_sd   <- sd(D$D, na.rm = TRUE)

  MBI <- D %>%
    mutate(MBI = (D - D_mean) / D_sd) %>%
    arrange(id)

  #--- QC outputs (unchanged)
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
        med_cover = median(cover, na.rm = TRUE),
        p10       = quantile(cover, 0.10, na.rm = TRUE),
        p90       = quantile(cover, 0.90, na.rm = TRUE),
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
      arrange(min_cover) %>%
      head(10)

    diag <- list(join_stats = join_stats, cover_stats = cover_stats, bad_subjects = bad_subjects)
  }

  list(MBI = MBI, modality_scores = Dm, weight_table = weight_table, qc_modality = qc_modality, diag = diag)
}



# -- NEW: wrapper that yields GBI(α,β), MBI(α,β,γ), and MDA -------------------
compute_flagship_abg <- function(
  zs, fdg, sizes,
  beta = 0.25, gamma = 2.5,
  alpha = NULL,                  # named vector c(CT=..., SA=..., SubVol=...) or NULL -> equal per modality
  cap_frac = 0.05,
  coverage_thr = 0.80
) {
  # MBI(α,β,γ): PET-weighted
  mbi_out <- compute_mbi_flagship(
    zs      = zs,
    pet_w   = fdg,
    size_df = sizes,
    alpha   = alpha,
    tune    = list(beta = beta, gamma = gamma, cap_frac = cap_frac, coverage_thr = coverage_thr)
  )
  mbi <- dplyr::rename(mbi_out$MBI, mbi_raw = MBI)  # id, mbi_raw

  # GBI(α,β): size-weighted only (PET weights = 1)
  gbi_out <- compute_mbi_flagship(
    zs      = zs,
    pet_w   = NULL,
    size_df = sizes,
    alpha   = alpha,
    tune    = list(beta = beta, cap_frac = cap_frac, coverage_thr = coverage_thr)
  )
  gbi <- dplyr::rename(gbi_out$MBI, gbi = MBI) # id, gbi

  # MBI: residual of MBI_raw ~ GBI (z-standardised)
  dt <- dplyr::inner_join(mbi, gbi, by = "id")
  dt$mbi <- as.numeric(scale(residuals(stats::lm(mbi_raw ~ gbi, data = dt))))

  list(
    scores      = dt[, c("id","mbi_raw","gbi","mbi")],
    mbi_details = mbi_out,
    gbi_details = gbi_out
  )
}
