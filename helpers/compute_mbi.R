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

# ---- SINGLE definition of the biological-map weight transform ------------
# Turns a raw biological map statistic (one value per region) into index
# weights. This is the ONLY place the map becomes a weight: stage 01 emits the
# raw parcel statistic and nothing else, and helpers/compare_pet_funcs.R
# aliases this function rather than re-implementing it.
#
# Parameterisation (methods "robust_z" and "rank"):
#
#   w_r = ratio ^ (u_r / 2),   u_r in [-1, 1],   then mean-normalised
#
# so the weight vector spans [ratio^-0.5, ratio^+0.5] and `ratio` is the
# maximum possible ratio between the highest- and lowest-weighted region.
# One parameter, one sentence: the highest-demand region counts at most
# `ratio` times the lowest.
#
# Why this replaces the previous z-score -> subtract-min -> ^gamma chain:
#   * the old scale was pinned to the single coldest parcel, which landed at
#     ~1e-6 and, after ^gamma, ~1e-15 - present in GBI, absent from MBI_raw;
#   * it depended on the map's distributional shape, so gamma = 2.5 meant
#     something different for FDG than for any other prior;
#   * one outlier parcel moved every weight.
# Both methods below are bounded, cannot zero a region, and are defined by the
# map's own robust statistics, so the same `ratio` means the same thing for
# any biological prior.
#
# gamma is NOT a separate parameter here: w^gamma = (ratio^gamma)^(u/2), i.e.
# exponentiating is exactly a change of `ratio`. Likewise winsorisation is
# unnecessary - neither method has outliers to trim. Both survive only under
# method = "legacy", which reproduces the pre-rebuild construction.
#
# Applied to a constant vector (the GBI case, all weights 1) every method
# returns all ones, so GBI is invariant to ratio, gamma and winsor_q.
transform_map_weights <- function(x,
                                  method      = c("robust_z", "rank", "legacy"),
                                  ratio       = 4,
                                  clamp_z     = 2,
                                  gamma       = 2.5,
                                  winsor_q    = 0.99,
                                  global_norm = TRUE) {
  method <- match.arg(method)
  x  <- as.numeric(x)
  ok <- is.finite(x) & x > 0

  # ---- legacy path: the pre-rebuild construction, kept reproducible -------
  if (identical(method, "legacy")) {
    v <- x; v[!ok] <- 0
    if (!is.null(winsor_q) && is.finite(winsor_q) && winsor_q > 0 && winsor_q < 1) {
      pos <- v[v > 0]
      if (length(pos) > 1)
        v <- pmin(v, as.numeric(quantile(pos, winsor_q, na.rm = TRUE)))
    }
    if (isTRUE(global_norm)) {
      m <- mean(v[v > 0], na.rm = TRUE)
      if (is.finite(m) && m > 0) v <- v / m
    }
    if (!is.null(gamma) && is.finite(gamma) && gamma != 1) v <- v ^ gamma
    return(v)
  }

  if (!is.finite(ratio) || ratio < 1)
    stop("transform_map_weights(): ratio must be finite and >= 1", call. = FALSE)

  w  <- rep(0, length(x))
  xv <- x[ok]
  n  <- length(xv)
  if (n == 0L) return(w)

  u <- rep(0, n)
  if (n > 1L) {
    if (identical(method, "rank")) {
      # Topography only: depends on the map's ordering, nothing else.
      r <- (rank(xv, ties.method = "average") - 0.5) / n
      u <- 2 * r - 1
    } else {
      # robust_z: keeps relative magnitude within the bulk of the map, on a
      # median/MAD scale so it is comparable across maps, with a hard clamp at
      # +/- clamp_z robust SDs so no outlier parcel can stretch the range.
      med <- stats::median(xv, na.rm = TRUE)
      s   <- stats::mad(xv, center = med, constant = 1.4826, na.rm = TRUE)
      if (!is.finite(s) || s <= 0) s <- stats::sd(xv, na.rm = TRUE)
      if (is.finite(s) && s > 0)
        u <- pmax(-1, pmin(1, (xv - med) / (clamp_z * s)))
      # else: no spread -> u stays 0 -> uniform weights (the GBI case)
    }
  }

  wv <- ratio ^ (u / 2)
  wv <- wv / mean(wv)
  w[ok] <- wv
  w
}

# Realised (as opposed to maximum) weight ratio actually present in a vector.
# `ratio` is an upper bound: it is only attained if the clamp engages on both
# sides. Report this so the configured and the achieved contrast never get
# confused for one another.
realised_weight_ratio <- function(w) {
  pos <- w[is.finite(w) & w > 0]
  if (!length(pos)) return(NA_real_)
  max(pos) / min(pos)
}

# ---- Feature-set alignment ------------------------------------------------
# Both indices MUST be built on identical features. Previously, regions with a
# missing biological weight were silently zeroed in MBI_raw but carried at full
# weight in GBI, so their deviations loaded on MBI purely by construction.
#
# A region is retained only if it has (a) a finite, positive map weight,
# (b) deviation data, and (c) a size entry. Everything else is dropped from
# BOTH indices, and reported.
align_feature_set <- function(zs, sizes, map) {

  zs    <- zs    %>% dplyr::mutate(region = as.character(region))
  sizes <- sizes %>% dplyr::mutate(region = as.character(region))
  map   <- map   %>% dplyr::mutate(region = as.character(region))

  if (any(duplicated(map$region)))
    stop("align_feature_set(): map has duplicated region rows: ",
         paste(unique(map$region[duplicated(map$region)]), collapse = ", "),
         ". A duplicated region silently double-counts after the joins.",
         call. = FALSE)

  map_ok <- map %>%
    dplyr::filter(is.finite(as.numeric(weight)), as.numeric(weight) > 0) %>%
    dplyr::pull(region) %>% unique()

  in_zs    <- unique(zs$region)
  in_sizes <- unique(sizes$region)

  keep <- Reduce(intersect, list(map_ok, in_zs, in_sizes))

  dropped_nomap  <- setdiff(in_zs, map_ok)
  dropped_nosize <- setdiff(in_zs, in_sizes)
  dropped        <- setdiff(in_zs, keep)

  if (length(dropped_nomap))
    message("align_feature_set(): ", length(dropped_nomap),
            " region(s) dropped for no valid map weight: ",
            paste(sort(dropped_nomap), collapse = ", "))
  if (length(dropped_nosize))
    message("align_feature_set(): ", length(dropped_nosize),
            " region(s) dropped for no size entry: ",
            paste(sort(dropped_nosize), collapse = ", "))

  list(
    zs      = zs    %>% dplyr::filter(region %in% keep),
    sizes   = sizes %>% dplyr::filter(region %in% keep),
    map     = map   %>% dplyr::filter(region %in% keep),
    kept    = sort(keep),
    dropped = sort(dropped)
  )
}

# ---- Prepare PET weights --------------------------------------------------
# If pet_w is NULL -> w_pet = 1 (uniform; used for GBI).
prepare_pet_weights <- function(zs, pet_w = NULL) {
  if (is.null(pet_w)) {
    zs %>% distinct(region) %>% mutate(w_pet = 1.0)
  } else {
    # Pass the RAW map statistic straight through. All normalisation and
    # contrast shaping happens in transform_map_weights(), so there is exactly
    # one place where a map becomes a weight. (This used to mean-normalise
    # first, which was redundant under the new transform and load-bearing
    # under the old one - precisely the split-across-two-files problem.)
    pet_w %>% transmute(region, w_pet = as.numeric(weight))
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
#  - Map weights via transform_map_weights(method = "robust_z", ratio, clamp_z)
#  - Size weights with beta
#  - Equal alpha per available modality per subject (unless alpha supplied)
#  - coverage_thr = 0.80
#
# Note on where influence comes from. Writing both indices as normalised
# weighted means gives contrast weights e_r = W_r - b*U_r, and since
# W_r ~ w_r*s_r while U_r ~ s_r, this reduces to
#
#     e_r  ~  s_r * (w_r - c)
#
# i.e. a region's influence on MBI is its SIZE weight times how far its MAP
# weight sits from the mean. Both `beta` and `ratio` are therefore levers on
# concentration, and index_concentration() below measures the result.
compute_mbi_flagship <- function(zs, pet_w = NULL, size_df = NULL, alpha = NULL, tune = NULL) {

  stopifnot(all(c("id","modality","region","z") %in% names(zs)))

  cfg <- list(
    map_method        = "robust_z",
    ratio             = 4,
    clamp_z           = 2,
    use_size          = TRUE,
    beta              = 0.25,
    per_modality_norm = FALSE,
    coverage_thr      = 0.80,
    # legacy path only (map_method = "legacy"); ignored otherwise
    gamma             = 2.5,
    winsor_q          = 0.99,
    use_global_norm   = TRUE,
    debug             = FALSE
  )
  if (!is.null(tune) && length(tune)) {
    for (k in intersect(names(tune), names(cfg))) cfg[[k]] <- tune[[k]]
  }

  #--- orient z-scores (drops modalities with s == 0)
  zo <- compute_oriented(zs)

  #--- map weights via the single shared transform
  w_pet0 <- prepare_pet_weights(zs, pet_w)
  w_pet0$w_pet <- transform_map_weights(
    w_pet0$w_pet,
    method      = cfg$map_method,
    ratio       = cfg$ratio,
    clamp_z     = cfg$clamp_z,
    gamma       = cfg$gamma,
    winsor_q    = cfg$winsor_q,
    global_norm = cfg$use_global_norm
  )

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
    # Two fixes to the supplied-alpha branch (latent: flagship runs use NULL):
    # 1. cfg$mbi$alpha arrives from YAML as a *list*; alpha[modality] on a list
    #    yields a list-column and downstream arithmetic corrupts. Coerce and
    #    require names.
    # 2. Renormalise over modalities the subject actually has (filter NA Dm
    #    FIRST), matching the equal-alpha branch. Previously the sum included
    #    missing modalities, so subjects with incomplete data had alpha summing
    #    to < 1 and systematically shrunken scores.
    a <- unlist(alpha)
    if (is.null(names(a)) || !all(nzchar(names(a))))
      stop("alpha must be a NAMED numeric vector/list, e.g. c(CT=1, SA=1, SubVol=1)",
           call. = FALSE)
    miss <- setdiff(unique(Dm$modality[!is.na(Dm$Dm)]), names(a))
    if (length(miss))
      stop("alpha has no entry for modality: ", paste(miss, collapse = ", "),
           call. = FALSE)
    Dm %>%
      filter(!is.na(Dm)) %>%
      mutate(alpha = as.numeric(a[modality])) %>%
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


# ---- Concentration diagnostic --------------------------------------------
# Effective number of contributing features (inverse Herfindahl on the
# contrast weights e = W - b*U, where W and U are the normalised effective
# weights of MBI_raw and GBI). N_eff = n_features would mean perfectly even
# influence; small values mean a handful of regions drive the index.
#
# Computed from weights only - no subject data - so it can be used to screen
# candidate designs before running the pipeline.
index_concentration <- function(weight_table_mbi, weight_table_gbi) {
  key <- function(x) paste(x$modality, x$region, sep = "|")
  stopifnot(identical(sort(key(weight_table_mbi)), sort(key(weight_table_gbi))))

  m <- weight_table_mbi %>% arrange(modality, region)
  g <- weight_table_gbi %>% arrange(modality, region)

  nmod <- dplyr::n_distinct(m$modality)
  Wn <- m %>% group_by(modality) %>% mutate(W = (w_bar / sum(w_bar)) / nmod) %>% ungroup()
  Un <- g %>% group_by(modality) %>% mutate(U = (w_bar / sum(w_bar)) / nmod) %>% ungroup()

  W <- Wn$W / sqrt(sum(Wn$W^2))
  U <- Un$U / sqrt(sum(Un$U^2))
  # b here is the cosine between the two weight vectors, i.e. the regression
  # slope under an iid-z approximation. The empirical slope (score-level
  # cor(mbi_raw, gbi)) also reflects the z covariance structure, so the two
  # differ slightly; for design screening the approximation is sufficient and
  # requires no subject data. (Empirical check on the flagship build: both
  # give N_eff ~ 7.75.)
  b <- sum(W * U)
  e <- W - b * U
  p <- e^2 / sum(e^2)

  list(
    n_features = nrow(m),
    b          = b,
    n_eff      = 1 / sum(p^2),
    shares     = Wn %>% transmute(modality, region, share = p) %>% arrange(desc(share))
  )
}


# ---- Wrapper: yields GBI(alpha,beta), MBI_raw(alpha,beta,gamma) and MBI --
compute_flagship_abg <- function(
  zs, fdg, sizes,
  beta         = 0.25,
  map_method   = "robust_z",
  ratio        = 4,
  clamp_z      = 2,
  alpha        = NULL,
  coverage_thr = 0.80,
  # legacy path only (map_method = "legacy")
  gamma        = 2.5,
  winsor_q     = 0.99
) {

  # ---- Enforce ONE feature set for both indices --------------------------
  al    <- align_feature_set(zs, sizes, fdg)
  zs    <- al$zs
  sizes <- al$sizes
  fdg   <- al$map

  # MBI_raw: PET-weighted
  mbi_out <- compute_mbi_flagship(
    zs      = zs,
    pet_w   = fdg,
    size_df = sizes,
    alpha   = alpha,
    tune    = list(beta = beta, map_method = map_method,
                   ratio = ratio, clamp_z = clamp_z,
                   gamma = gamma, winsor_q = winsor_q,
                   coverage_thr = coverage_thr)
  )
  mbi <- dplyr::rename(mbi_out$MBI, mbi_raw = MBI)

  # GBI: size-weighted only (PET weights = 1; gamma & winsor_q irrelevant when all w_pet = 1)
  # GBI: uniform map weights. Every transform method returns all ones on a
  # constant vector, so GBI is invariant to map_method, ratio, gamma, winsor_q.
  gbi_out <- compute_mbi_flagship(
    zs      = zs,
    pet_w   = NULL,
    size_df = sizes,
    alpha   = alpha,
    tune    = list(beta = beta, map_method = map_method,
                   coverage_thr = coverage_thr)
  )
  gbi <- dplyr::rename(gbi_out$MBI, gbi = MBI)

  # ---- Guard: refuse to residualise mismatched feature sets --------------
  fkey <- function(x) sort(paste(x$modality, x$region, sep = "|"))
  if (!identical(fkey(mbi_out$weight_table), fkey(gbi_out$weight_table)))
    stop("MBI and GBI feature sets diverged - refusing to residualise. ",
         "Check align_feature_set() and the map/size/deviation inputs.",
         call. = FALSE)

  # MBI: residual of MBI_raw ~ GBI (z-standardised)
  # na.exclude, NOT the default na.omit: with na.omit, residuals() returns a
  # vector shorter than nrow(dt) whenever any score is NA, and the column
  # assignment below would either error or (worse) silently recycle and
  # misalign subjects. na.exclude pads residuals back to full length with NA.
  dt <- dplyr::inner_join(mbi, gbi, by = "id")
  fit <- stats::lm(mbi_raw ~ gbi, data = dt, na.action = stats::na.exclude)
  dt$mbi <- as.numeric(scale(stats::residuals(fit)))

  conc      <- index_concentration(mbi_out$weight_table, gbi_out$weight_table)
  w_ratio   <- realised_weight_ratio(mbi_out$weight_table$w_pet)
  r2_raw_g  <- suppressWarnings(stats::cor(dt$mbi_raw, dt$gbi,
                                           use = "pairwise.complete.obs")^2)

  message(sprintf(
    "compute_flagship_abg(): %d features | method=%s ratio=%.2f (realised %.2f) beta=%.2f | b=%.3f | R2(raw,GBI)=%.3f | N_eff=%.2f",
    conc$n_features, map_method, ratio, w_ratio, beta, conc$b, r2_raw_g, conc$n_eff))

  list(
    scores        = dt[, c("id","mbi_raw","gbi","mbi")],
    mbi_details   = mbi_out,
    gbi_details   = gbi_out,
    features      = list(kept = al$kept, dropped = al$dropped, n_kept = length(al$kept)),
    concentration = conc,
    diagnostics   = list(map_method = map_method, ratio = ratio,
                         realised_ratio = w_ratio, beta = beta,
                         r2_raw_gbi = r2_raw_g)
  )
}
