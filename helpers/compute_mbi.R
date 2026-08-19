#!/usr/bin/env Rscript
#
# MBI construction.
#
# WHAT MBI IS
#
#   mbi_raw = sum_r W_r z_r      map-weighted  mean regional deviation
#   gbi     = sum_r U_r z_r      unweighted    mean regional deviation
#   MBI     = standardised residual of mbi_raw on gbi
#
# Because both are weighted averages of the SAME regional deviations, the
# residual is exactly a linear contrast of those deviations:
#
#   MBI  =  sum_r (W_r - b * U_r) * z_r ,  standardised
#
# with coefficients that sum to ~0. It is a difference between two weighted
# averages — a profile score, like a laterality index — not a remainder. State
# it this way: it makes clear that a high R^2(mbi_raw, gbi) is a property of how
# the two averages were scaled, not of the index that is reported.
#
# PARAMETERS
#
#   emphasis   the only parameter that changes MBI. Report emphasis in {1,2,3}
#              as a sensitivity strip.
#   clamp_z    safety rail against a single outlier parcel. Whether it engages
#              at all depends on the map; check realised_weight_ratio().
#   ratio      sets the weight scale. ALGEBRAICALLY INERT for MBI — see
#              transform_map_weights() below. Defaults to "auto", i.e. read off
#              the supplied map, so it is a derived quantity rather than a
#              setting. A number is accepted for testing.
#   beta       size weighting exponent. 0 asks a topography question (every
#              standardised deviation is one observation); 1 asks a burden
#              question (volume weighting). Report as a sensitivity.
#   coverage_thr  subject inclusion.

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
# so they cannot dilute the index with zero-weighted contributions downstream.
compute_oriented <- function(zs) {
  zs %>%
    mutate(s = get_orientation(.), z_oriented = s * z) %>%
    filter(s != 0)
}

# ---- Resolving `ratio` ----------------------------------------------------
# `ratio` is the nominal max:min of the weight vector. It is inert for MBI, so
# its only job is to put mbi_raw on an interpretable scale — and the
# interpretable scale is the one the biology already has. "auto" therefore reads
# it off the map itself: the weights span the same multiplicative range as the
# supplied values.
#
# This is what makes the config honest. A hard-coded 1.73 looks like a chosen
# number, has to be re-derived by hand whenever the map changes, and goes stale
# silently if nobody does. Derived, it is not a setting at all.
#
# A literal number is still accepted, which is what extras/verify_invariance.R
# uses to sweep ratio over three orders of magnitude and show MBI does not move.
resolve_ratio <- function(x, ratio = "auto") {
  if (is.character(ratio)) {
    if (!tolower(ratio) %in% c("auto", "map", "observed"))
      stop("transform_map_weights(): ratio must be a number >= 1 or \"auto\"; got \"",
           ratio, "\"", call. = FALSE)
    v <- as.numeric(x)
    v <- v[is.finite(v) & v > 0]
    if (length(v) < 2L) return(1)
    r <- max(v) / min(v)
    return(if (is.finite(r) && r >= 1) r else 1)
  }
  if (!is.finite(ratio) || ratio < 1)
    stop("transform_map_weights(): ratio must be finite and >= 1, or \"auto\"",
         call. = FALSE)
  ratio
}


# ---- SINGLE definition of the biological-map weight transform ------------
# Turns a raw biological map statistic (one value per region) into index
# weights. This is the ONLY place a map becomes a weight: stage 01 emits the
# raw parcel statistic and nothing else, and helpers/compare_pet_funcs.R
# aliases this function rather than re-implementing it.
#
#   u_r = clamp( (x_r - median(x)) / (clamp_z * MAD(x)) , -1, +1 )
#   u_r <- sign(u_r) * |u_r|^emphasis
#   w_r = 1 + kappa * u_r ,   kappa = (ratio - 1) / (ratio + 1)
#   w   <- w / mean(w)
#
# The transform is ADDITIVE in u, and must be: influence on MBI goes as
# e_r ~ (w_r - mean(w)), so under this form a region one robust SD BELOW the
# median of the map pulls exactly as hard as one a SD above, in the opposite
# direction. A multiplicative form such as ratio^(u/2) is bounded below by zero,
# so after mean-normalisation a high-value region can sit arbitrarily far above
# the mean weight while a low-value region sits at most 1.0 below it: low-map
# regions could then only fail to contribute, never pull. That asymmetry is a
# property of the functional form rather than of the data, and it suppresses
# exactly the low-metabolism territory the index is meant to be sensitive to.
# Do not replace this with a multiplicative weighting.
#
# `emphasis` applies a convex, sign-symmetric emphasis to the tails. Since
# |u| <= 1, emphasis > 1 shrinks typical-demand regions toward the mean weight
# while leaving both extremes at +/-1. The hypothesis it encodes: metabolic
# demand relates non-linearly to structural vulnerability, with regions at
# either extreme of the demand distribution disproportionately affected and
# typical-demand regions buffered. emphasis = 1 recovers the plain linear form.
# As emphasis grows the index converges on a high- versus low-demand
# extreme-group contrast.
#
# `ratio` IS INERT FOR MBI, and this is provable rather than incidental:
# mbi_raw is proportional to (sum_r z_r) + kappa * (sum_r u_r z_r) over a
# constant denominator, and gbi is proportional to (sum_r z_r). Residualising
# on gbi removes the first term, leaving a direction with no kappa in it; and
# the influence shares e_r ~ kappa * u_r have kappa cancel in
# p = e^2 / sum(e^2). `ratio` therefore cannot be tuned to change a result,
# only to relabel the weight scale. extras/verify_invariance.R checks this
# against this implementation over three orders of magnitude in ratio.
#
# Applied to a constant vector (the GBI case) the transform returns all ones,
# so GBI is invariant to every parameter here.
transform_map_weights <- function(x,
                                  ratio    = "auto",
                                  clamp_z  = 3,
                                  emphasis = 1) {
  x     <- as.numeric(x)
  ratio <- resolve_ratio(x, ratio)
  ok    <- is.finite(x) & x > 0
  if (!is.finite(clamp_z) || clamp_z <= 0)
    stop("transform_map_weights(): clamp_z must be finite and > 0", call. = FALSE)
  if (!is.finite(emphasis) || emphasis <= 0)
    stop("transform_map_weights(): emphasis must be finite and > 0", call. = FALSE)

  w  <- rep(0, length(x))
  xv <- x[ok]
  n  <- length(xv)
  if (n == 0L) return(w)

  u <- rep(0, n)
  if (n > 1L) {
    med <- stats::median(xv, na.rm = TRUE)
    s   <- stats::mad(xv, center = med, constant = 1.4826, na.rm = TRUE)
    if (!is.finite(s) || s <= 0) s <- stats::sd(xv, na.rm = TRUE)
    if (is.finite(s) && s > 0)
      u <- pmax(-1, pmin(1, (xv - med) / (clamp_z * s)))
    # else: no spread -> u stays 0 -> uniform weights (the GBI case)
  }

  if (emphasis != 1) u <- sign(u) * abs(u)^emphasis

  kappa <- (ratio - 1) / (ratio + 1)   # < 1, so every weight stays positive
  wv    <- 1 + kappa * u
  wv    <- wv / mean(wv)

  w[ok] <- wv
  w
}

# Realised (as opposed to maximum) weight ratio actually present in a vector.
# `ratio` is an upper bound, attained only if the clamp engages on both sides.
# Report this alongside the configured value so the two are never confused.
realised_weight_ratio <- function(w) {
  pos <- w[is.finite(w) & w > 0]
  if (!length(pos)) return(NA_real_)
  max(pos) / min(pos)
}

# ---- Feature-set alignment ------------------------------------------------
# Both indices MUST be built on identical features. A region carried in gbi but
# not in mbi_raw — or the reverse — loads on the residual purely by
# construction, because the contrast e_r = W_r - b*U_r then has a term with
# only one of the two weights in it.
#
# A region is retained only if it has (a) a finite, positive map weight,
# (b) deviation data, and (c) a usable size entry. Everything else is dropped
# from BOTH indices, and reported.
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
  # a region needs a usable size for EVERY feature it contributes, because
  # under global pooling `size` is a single cross-modality scale rather than
  # something that cancels within a bucket
  bad_size <- unique(sizes$region[!is.finite(as.numeric(sizes$size)) |
                                   as.numeric(sizes$size) <= 0])
  in_sizes <- setdiff(unique(sizes$region), bad_size)

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

# ---- Prepare map weights --------------------------------------------------
# If pet_w is NULL -> w_pet = 1 (uniform; this is the GBI case).
prepare_pet_weights <- function(zs, pet_w = NULL) {
  if (is.null(pet_w)) {
    zs %>% distinct(region) %>% mutate(w_pet = 1.0)
  } else {
    # Pass the RAW map statistic straight through. All normalisation and
    # contrast shaping happens in transform_map_weights(), so there is exactly
    # one place where a map becomes a weight.
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
# GLOBAL POOLING:  D = sum_all(w * z) / sum_all(w), over every feature at once.
# There is no per-modality nesting and no per-modality mixing weight. Each
# feature contributes in proportion to its own weight, full stop.
#
# This is a constraint, not a default to be relaxed. Renormalising the map
# weights inside each modality would stop the map speaking ACROSS modalities,
# leaving a mixing weight rather than the map to decide how much subcortex
# counts. It also introduces leverage nobody chose: within a modality,
# sum(e^2) scales as spread^2 / n_m, and that 1/n_m term is invisible to any
# mixing weight, so a modality with few regions hands each of its regions
# several times the per-region leverage of a finely parcellated one.
#
# Interaction with beta: because pooling is global, `size` does not cancel
# within a bucket and must be on ONE scale across all features —
# 02_run_normdevs.R emits tissue volume in mm^3 for every feature. At beta = 0
# the size term is inert, but it must stay coherent for beta > 0 to mean
# anything.
#
# Where influence comes from. Writing both indices as normalised weighted means
# gives contrast weights e_r = W_r - b*U_r, and since W_r ~ w_r*s_r while
# U_r ~ s_r, this reduces to
#
#     e_r  ~  s_r * (w_r - c)
#
# i.e. a region's influence on MBI is its SIZE weight times how far its MAP
# weight sits from the mean. index_concentration() below measures the result.
compute_mbi_flagship <- function(zs, pet_w = NULL, size_df = NULL, tune = NULL) {

  stopifnot(all(c("id","modality","region","z") %in% names(zs)))

  cfg <- list(
    ratio        = "auto",
    clamp_z      = 3,
    emphasis     = 1,
    beta         = 0,
    coverage_thr = 0.80,
    required_modalities = NULL,
    quiet        = FALSE,   # the GBI pass reuses the same subject rule; report once
    debug        = FALSE
  )
  if (!is.null(tune) && length(tune)) {
    unknown <- setdiff(names(tune), names(cfg))
    if (length(unknown))
      stop("compute_mbi_flagship(): unknown tune parameter(s): ",
           paste(unknown, collapse = ", "),
           ". Accepted: ", paste(names(cfg), collapse = ", "), ".",
           call. = FALSE)
    # Assign through a temporary, because `cfg[[k]] <- NULL` DELETES the entry
    # rather than setting it. Without this, tune = list(beta = NULL) silently
    # removes beta and the next is.finite(cfg$beta) evaluates on logical(0),
    # which errors with "argument is of length zero" far from the cause.
    # required_modalities is the one parameter for which NULL is meaningful,
    # so it is set explicitly and everything else rejects NULL outright.
    for (k in names(tune)) {
      v <- tune[[k]]
      if (is.null(v)) {
        if (!identical(k, "required_modalities"))
          stop("compute_mbi_flagship(): tune$", k, " is NULL. Only ",
               "required_modalities may be NULL; omit a parameter to use its ",
               "default.", call. = FALSE)
        cfg["required_modalities"] <- list(NULL)
      } else {
        cfg[[k]] <- v
      }
    }
  }

  #--- orient z-scores (drops modalities with s == 0)
  zo <- compute_oriented(zs)

  #--- map weights via the single shared transform
  w_pet0 <- prepare_pet_weights(zs, pet_w)
  w_pet0$w_pet <- transform_map_weights(
    w_pet0$w_pet,
    ratio    = cfg$ratio,
    clamp_z  = cfg$clamp_z,
    emphasis = cfg$emphasis
  )

  #--- size weights
  w_size0 <- prepare_size_weights(zs, size_df)
  if (is.finite(cfg$beta) && cfg$beta != 1)
    w_size0 <- w_size0 %>% mutate(size = pmax(1e-12, size) ^ cfg$beta)

  #--- merge + combined weights
  zw <- zo %>%
    inner_join(w_pet0,  by = "region") %>%
    inner_join(w_size0, by = c("region","modality")) %>%
    mutate(w_bar = pmax(0, w_pet) * pmax(1e-12, size))

  #--- coverage rule, computed per (id, modality) so subject inclusion is
  #--- driven by per-modality data completeness rather than by the pooled total
  cover <- zw %>%
    group_by(id, modality) %>%
    summarise(cover = mean(!is.na(z_oriented), na.rm = TRUE), .groups = "drop")

  zwg <- zw %>%
    left_join(cover, by = c("id", "modality")) %>%
    filter(cover >= cfg$coverage_thr)

  # ---- REQUIRED-MODALITY GUARDRAIL ----------------------------------------
  # The coverage rule above is per (subject, modality): it decides whether a
  # subject contributes a modality, not whether they have enough modalities to
  # get a score at all. Without a further rule a subject with no subcortical
  # data still receives an index, built from CT and SA only.
  #
  # For a mean index an incomplete modality set would merely be a noisier
  # estimate of the same quantity. MBI is a CONTRAST, so a different feature set
  # is a different linear functional — a genuinely different index,
  # standardised into the same column as everyone else's. On simulated data with
  # a shared latent spatial pattern, an index built without SubVol correlates
  # r = 0.92 with the full index and one built without CT r = 0.77: close, but
  # not interchangeable.
  #
  # Set mbi.required_modalities to enforce a complete modality set. NULL scores
  # everyone who clears the coverage rule and warns if subjects differ in how
  # many modalities they contribute. A requested modality that is absent from
  # the data entirely is an ERROR, not a silently relaxed constraint.
  mod_count <- zwg %>%
    distinct(id, modality) %>%
    count(id, name = "n_mod")
  mods_all <- sort(unique(zw$modality))

  req <- cfg$required_modalities
  if (!is.null(req) && length(req)) {
    req <- as.character(req)
    missing_mods <- setdiff(req, mods_all)
    if (length(missing_mods))
      stop("compute_mbi_flagship(): required_modalities names modality/modalities ",
           "absent from the data: ", paste(missing_mods, collapse = ", "),
           ". Present: ", paste(mods_all, collapse = ", "),
           ". Refusing to silently relax the requirement.", call. = FALSE)

    keep_ids <- zwg %>%
      distinct(id, modality) %>%
      filter(modality %in% req) %>%
      count(id, name = "n_req") %>%
      filter(n_req == length(req)) %>%
      pull(id)

    n_before <- dplyr::n_distinct(zwg$id)
    n_after  <- length(keep_ids)
    if (n_after < n_before && !isTRUE(cfg$quiet)) {
      short <- zwg %>%
        distinct(id, modality) %>%
        filter(!id %in% keep_ids, modality %in% req) %>%
        count(modality, name = "n_subj_with")
      message(sprintf(
        "required_modalities: %d of %d subjects dropped for an incomplete modality set (required: %s)",
        n_before - n_after, n_before, paste(req, collapse = ", ")))
      if (nrow(short)) print(as.data.frame(short), row.names = FALSE)
    } else if (!isTRUE(cfg$quiet)) {
      message(sprintf("required_modalities: all %d subjects have the full set (%s)",
                      n_before, paste(req, collapse = ", ")))
    }
    zwg <- zwg %>% filter(id %in% keep_ids)
  } else {
    tab <- mod_count %>% count(n_mod, name = "n_subjects") %>% arrange(desc(n_mod))
    if (nrow(tab) > 1L && !isTRUE(cfg$quiet))
      warning("compute_mbi_flagship(): mbi.required_modalities is not set and subjects ",
              "differ in how many modalities they contribute, so MBI is not the same ",
              "linear contrast for every subject. Modality counts: ",
              paste(sprintf("%d mod -> %d subj", tab$n_mod, tab$n_subjects),
                    collapse = "; "), call. = FALSE)
  }

  modality_composition <- mod_count %>%
    left_join(zwg %>% distinct(id, modality) %>%
                group_by(id) %>%
                summarise(mods = paste(sort(modality), collapse = "+"), .groups = "drop"),
              by = "id") %>%
    filter(!is.na(mods)) %>%
    count(mods, name = "n_subjects") %>%
    arrange(desc(n_subjects))

  D <- zwg %>%
    group_by(id) %>%
    summarise(
      num = sum(w_bar * z_oriented, na.rm = TRUE),
      den = sum(w_bar[!is.na(z_oriented)], na.rm = TRUE),
      D   = ifelse(den > 0, num / den, NA_real_),
      .groups = "drop"
    ) %>%
    filter(!is.na(D)) %>%
    select(id, D)

  # per-modality scores, QC only — they are not combined into D
  Dm <- zwg %>%
    group_by(id, modality) %>%
    summarise(
      num = sum(w_bar * z_oriented, na.rm = TRUE),
      den = sum(w_bar[!is.na(z_oriented)], na.rm = TRUE),
      Dm  = ifelse(den > 0, num / den, NA_real_),
      .groups = "drop"
    ) %>%
    left_join(cover, by = c("id", "modality"))

  D_mean <- mean(D$D, na.rm = TRUE)
  D_sd   <- sd(D$D,   na.rm = TRUE)

  MBI <- D %>%
    mutate(MBI = (D - D_mean) / D_sd) %>%
    arrange(id)

  #--- QC outputs
  qc_modality <- Dm %>%
    filter(!is.na(Dm)) %>%
    group_by(modality) %>%
    summarise(mean_Dm = mean(Dm), sd_Dm = sd(Dm), .groups = "drop")

  weight_table <- zw %>%
    distinct(region, modality, w_pet, size, w_bar) %>%
    arrange(modality, region)

  diag <- NULL
  if (isTRUE(cfg$debug)) {
    join_stats <- zw %>%
      group_by(modality) %>%
      summarise(n = n(), nz = sum(w_bar > 0, na.rm = TRUE),
                frac_zero = 1 - nz / n, .groups = "drop")
    cover_stats <- cover %>%
      group_by(modality) %>%
      summarise(med_cover  = median(cover, na.rm = TRUE),
                p10        = quantile(cover, 0.10, na.rm = TRUE),
                p90        = quantile(cover, 0.90, na.rm = TRUE),
                pct_ge_thr = mean(cover >= cfg$coverage_thr, na.rm = TRUE),
                .groups    = "drop")
    bad_subjects <- cover %>%
      group_by(id) %>%
      summarise(min_cover = min(cover, na.rm = TRUE),
                mods = paste(modality[order(cover)], collapse = ","),
                .groups = "drop") %>%
      arrange(min_cover) %>% head(10)
    diag <- list(join_stats = join_stats, cover_stats = cover_stats,
                 bad_subjects = bad_subjects)
  }

  list(MBI = MBI, modality_scores = Dm, weight_table = weight_table,
       qc_modality = qc_modality, modality_composition = modality_composition,
       diag = diag)
}


# ---- Concentration diagnostic --------------------------------------------
# Effective number of contributing features (inverse Herfindahl on the contrast
# weights e = W - b*U, where W and U are the normalised effective weights of
# mbi_raw and gbi). N_eff = n_features would mean perfectly even influence;
# small values mean a handful of regions drive the index, which is the single
# most useful thing to know about a candidate weighting. Report it, and report
# the top influence shares alongside it.
#
# Computed from weights only — no subject data — so candidate designs can be
# screened before running the pipeline.
#
# WHICH SLOPE. e depends on b, and because W and U are nearly parallel, e is a
# near-cancellation: a small change in b is a large relative change in the
# contrast direction. Two slopes are available and they are not the same:
#
#   b_cosine     sum(W * U), the slope under an iid-z approximation. Needs no
#                subject data, so it is the one available at design time.
#   b_empirical  the realised slope of lm(mbi_raw ~ gbi), which also carries
#                the covariance structure of z. This is the slope the reported
#                index is actually built with.
#
# Pass b_empirical when you have it and n_eff will be computed at the index's
# own slope, with n_eff_cosine retained alongside for comparability. On a
# synthetic replica of a 143-feature build the two differ by roughly 2%, so
# this matters for reporting a figure, not for choosing a design.
index_concentration <- function(weight_table_mbi, weight_table_gbi,
                                b_empirical = NULL) {
  key <- function(x) paste(x$modality, x$region, sep = "|")
  stopifnot(identical(sort(key(weight_table_mbi)), sort(key(weight_table_gbi))))

  m <- weight_table_mbi %>% arrange(modality, region)
  g <- weight_table_gbi %>% arrange(modality, region)

  Wn <- m %>% mutate(W = w_bar / sum(w_bar))
  Un <- g %>% mutate(U = w_bar / sum(w_bar))

  W <- Wn$W / sqrt(sum(Wn$W^2))
  U <- Un$U / sqrt(sum(Un$U^2))

  b_cos  <- sum(W * U)
  shares <- function(b) { e <- W - b * U; e^2 / sum(e^2) }
  p_cos  <- shares(b_cos)

  use_emp <- !is.null(b_empirical) && is.finite(b_empirical)
  b_used  <- if (use_emp) b_empirical else b_cos
  p_used  <- if (use_emp) shares(b_used) else p_cos

  list(
    n_features    = nrow(m),
    b             = b_used,
    b_cosine      = b_cos,
    b_empirical   = if (use_emp) b_empirical else NA_real_,
    b_source      = if (use_emp) "empirical (lm slope)" else "cosine approximation",
    n_eff         = 1 / sum(p_used^2),
    n_eff_cosine  = 1 / sum(p_cos^2),
    shares        = Wn %>% transmute(modality, region, share = p_used) %>%
                      arrange(desc(share))
  )
}


# ---- Wrapper: yields gbi, mbi_raw and MBI --------------------------------
compute_flagship_abg <- function(
  zs, fdg, sizes,
  ratio        = "auto",
  clamp_z      = 3,
  emphasis     = 1,
  beta         = 0,
  coverage_thr = 0.80,
  required_modalities = NULL
) {

  # ---- Enforce ONE feature set for both indices --------------------------
  al    <- align_feature_set(zs, sizes, fdg)
  zs    <- al$zs
  sizes <- al$sizes
  fdg   <- al$map

  # Resolve ratio ONCE, from the aligned map, and pass the number down. Both
  # passes then share one weight scale, and the value that gets reported is the
  # value that was used.
  ratio_cfg  <- ratio
  ratio      <- resolve_ratio(fdg$weight, ratio)

  tune <- list(ratio = ratio, clamp_z = clamp_z, emphasis = emphasis,
               beta = beta, coverage_thr = coverage_thr,
               required_modalities = required_modalities)

  # mbi_raw: map-weighted
  mbi_out <- compute_mbi_flagship(zs = zs, pet_w = fdg, size_df = sizes, tune = tune)
  mbi     <- dplyr::rename(mbi_out$MBI, mbi_raw = MBI)

  # gbi: uniform map weights. The transform returns all ones on a constant
  # vector, so gbi is invariant to ratio, clamp_z and emphasis.
  gbi_out <- compute_mbi_flagship(zs = zs, pet_w = NULL, size_df = sizes,
                                  tune = utils::modifyList(tune, list(quiet = TRUE)))
  gbi     <- dplyr::rename(gbi_out$MBI, gbi = MBI)

  # ---- Guard: the two passes must score the SAME subjects ----------------
  # The subject rules (coverage, required_modalities) depend only on which
  # deviations exist, never on the weights, so the two sets are identical by
  # construction. Assert it anyway: if it ever stopped being true, the
  # residualisation below would silently be fitted on an inner join that is not
  # the analysis sample.
  if (!identical(sort(as.character(mbi$id)), sort(as.character(gbi$id))))
    stop("MBI and GBI scored different subject sets - refusing to residualise. ",
         "Check the coverage and required_modalities rules.", call. = FALSE)

  # ---- Guard: refuse to residualise mismatched feature sets --------------
  fkey <- function(x) sort(paste(x$modality, x$region, sep = "|"))
  if (!identical(fkey(mbi_out$weight_table), fkey(gbi_out$weight_table)))
    stop("MBI and GBI feature sets diverged - refusing to residualise. ",
         "Check align_feature_set() and the map/size/deviation inputs.",
         call. = FALSE)

  # MBI: residual of mbi_raw ~ gbi (z-standardised).
  #
  # THE RESULTING COLUMN IS SAMPLE-DEPENDENT. Both the slope b and the
  # standardisation are estimated on the subjects scored here, so `mbi` is
  # orthogonal to `gbi` and has SD 1 IN THIS SAMPLE AND NO OTHER. Take a subset
  # afterwards — a diagnostic group, the complete cases for some outcome — and
  # the column you are holding is no longer orthogonal to GBI and no longer
  # standardised. Analyses on a subset should re-residualise within that subset
  # rather than reuse this column.
  #
  # na.exclude, NOT the default na.omit: with na.omit, residuals() returns a
  # vector shorter than nrow(dt) whenever any score is NA, and the column
  # assignment below would either error or (worse) silently recycle and
  # misalign subjects. na.exclude pads residuals back to full length with NA.
  dt  <- dplyr::inner_join(mbi, gbi, by = "id")
  fit <- stats::lm(mbi_raw ~ gbi, data = dt, na.action = stats::na.exclude)
  dt$mbi <- as.numeric(scale(stats::residuals(fit)))
  b_emp  <- unname(stats::coef(fit)[2])

  # Concentration at the index's own realised slope, not the design-time
  # approximation; both are returned.
  conc     <- index_concentration(mbi_out$weight_table, gbi_out$weight_table,
                                  b_empirical = b_emp)
  w_ratio  <- realised_weight_ratio(mbi_out$weight_table$w_pet)
  r2_raw_g <- suppressWarnings(stats::cor(dt$mbi_raw, dt$gbi,
                                          use = "pairwise.complete.obs")^2)

  message(sprintf(
    "compute_flagship_abg(): %d features | ratio=%.3f %s (realised %.2f) clamp_z=%.2f emph=%.2f beta=%.2f | b=%.4f | R2(raw,GBI)=%.4f | N_eff=%.2f (at cosine b=%.4f: %.2f)",
    conc$n_features, ratio,
    if (is.character(ratio_cfg)) "[auto]" else "[config]",
    w_ratio, clamp_z, emphasis, beta, conc$b, r2_raw_g,
    conc$n_eff, conc$b_cosine, conc$n_eff_cosine))

  list(
    scores        = dt[, c("id","mbi_raw","gbi","mbi")],
    mbi_details   = mbi_out,
    gbi_details   = gbi_out,
    features      = list(kept = al$kept, dropped = al$dropped, n_kept = length(al$kept)),
    concentration = conc,
    diagnostics   = list(ratio = ratio, realised_ratio = w_ratio, clamp_z = clamp_z,
                         emphasis = emphasis, beta = beta, r2_raw_gbi = r2_raw_g,
                         ratio_source = if (is.character(ratio_cfg))
                           "derived from the map (auto)" else "set in config")
  )
}
