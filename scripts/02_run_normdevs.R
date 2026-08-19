#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(future.apply)
  library(tidyr)
})

source("helpers/run_helpers.R")
opt <- parse_common_args()
cfg <- read_config(opt$config)

source_local("normative_modeling.R")

# -----------------------------
# generic feature table + feature map
# -----------------------------
features_file    <- cfg$inputs$features_file    %||% NULL
feature_map_file <- cfg$inputs$feature_map_file %||% NULL

if (is.null(features_file))    stop("Missing config: inputs: features_file", call. = FALSE)
if (is.null(feature_map_file)) stop("Missing config: inputs: feature_map_file", call. = FALSE)

assert_file_exists(features_file, "features_file")
assert_file_exists(feature_map_file, "feature_map_file")

# Column mapping (allow callers to avoid UKB-specific names)
colmap <- cfg$inputs$colmap %||% list()
col_id      <- colmap$id      %||% "id"
col_age     <- colmap$age     %||% "age"
col_sex     <- colmap$sex     %||% "sex"
col_scanner <- colmap$scanner %||% "scanner"
col_icv     <- colmap$icv     %||% NULL  # optional

# Normative parameters
k_age <- cfg$normative$k_age %||% 10
use_site_RE <- isTRUE(cfg$normative$use_site_RE %||% TRUE)
sex_interact <- isTRUE(cfg$normative$sex_interact %||% FALSE)

# Location-SCALE normative model (mgcv gaulss). Each subject is z-scored
# against their own sigma_i = sigma(age, sex, scanner) rather than one scalar
# sigma per feature. See helpers/normative_modeling.R for why this matters more
# for a contrast index than for a mean index; false gives a location-only fit
# and is a diagnostic comparison, not an alternative specification.
variance_model <- isTRUE(cfg$normative$variance_model %||% TRUE)
k_age_scale    <- cfg$normative$k_age_scale %||% 5

cat(sprintf("[normative] model = %s (k_age = %d%s)\n",
            if (variance_model) "location-scale, mgcv gaulss" else "location-only, gaussian",
            k_age,
            if (variance_model) sprintf(", k_age_scale = %d", k_age_scale) else ""))

# -----------------------------
# Runtime
# -----------------------------
n_cores <- cfg$runtime$n_cores %||% max(1L, parallel::detectCores(logical = TRUE) - 1L)
future::plan(future::multisession, workers = n_cores)

# -----------------------------
# Outputs
# -----------------------------
out_stage <- stage_dir(opt$outdir, "normdev")
out_qc    <- file.path(out_stage, "normative_qc.txt")
out_z     <- file.path(out_stage, "deviations_z.txt")
out_mu    <- file.path(out_stage, "predicted_means.txt")
out_long  <- file.path(out_stage, "deviations_z_long.txt")
out_raw   <- file.path(out_stage, "feats_rawVals.txt")
out_sizes <- file.path(out_stage, "feats_avgSize.txt")

if (file.exists(out_long) && !opt$overwrite) {
  stop("Output exists: ", out_long, " (use --overwrite)", call. = FALSE)
}

# -----------------------------
# Load inputs
# -----------------------------
dt_in <- fread(features_file, data.table = FALSE)
fm <- fread(feature_map_file, data.table = FALSE)

req_fm_cols <- c("feature", "modality", "region")
miss_fm <- setdiff(req_fm_cols, names(fm))
if (length(miss_fm) > 0) stop("feature_map_file missing columns: ", paste(miss_fm, collapse = ", "), call. = FALSE)

fm$feature <- as.character(fm$feature)
fm$modality <- as.character(fm$modality)
fm$region <- as.character(fm$region)

# Ensure unique mapping per feature
if (any(duplicated(fm$feature))) {
  dups <- unique(fm$feature[duplicated(fm$feature)])
  stop("feature_map_file has duplicated 'feature' entries (must be unique). Examples: ", paste(head(dups, 5), collapse = ", "), call. = FALSE)
}

# ...and a unique (modality, region) per feature. Nothing downstream can detect
# a violation: two features sharing a (modality, region) make the cortical
# volume computation below produce list-columns, and if they sit in different
# modalities the region simply contributes twice to both indices, showing up in
# the influence shares with no explanation. It has to be caught here.
fm_key <- paste(fm$modality, fm$region, sep = "|")
if (any(duplicated(fm_key))) {
  dupk <- unique(fm_key[duplicated(fm_key)])
  stop("feature_map_file maps more than one feature to the same ",
       "(modality, region): ", paste(head(dupk, 5), collapse = "; "),
       if (length(dupk) > 5) sprintf(" ... and %d more", length(dupk) - 5) else "",
       ". Each (modality, region) must correspond to exactly one feature.",
       call. = FALSE)
}

# Validate required covariate columns
req_cols <- c(col_id, col_age, col_sex, col_scanner)
miss_cov <- setdiff(req_cols, names(dt_in))
if (length(miss_cov) > 0) {
  stop("features_file missing required columns: ", paste(miss_cov, collapse = ", "), call. = FALSE)
}

# Validate feature columns
feature_cols <- fm$feature
miss_feat <- setdiff(feature_cols, names(dt_in))
if (length(miss_feat) > 0) {
  stop("features_file is missing feature columns referenced in feature_map_file. Examples: ",
       paste(head(miss_feat, 8), collapse = ", "), call. = FALSE)
}

# Standardize column names internally
dt <- dt_in %>%
  dplyr::rename(
    id = !!col_id,
    age = !!col_age,
    sex = !!col_sex,
    scanner = !!col_scanner
  )
dt$id <- as.character(dt$id)

# Handle ICV if provided
if (!is.null(col_icv)) {
  if (!col_icv %in% names(dt)) stop("colmap: icv column not found in features_file: ", col_icv, call. = FALSE)
  dt <- dt %>% dplyr::rename(icv_raw = !!col_icv)
}

# Minimal completeness
keep_cols <- c("id", "age", "sex", "scanner", feature_cols)
if ("icv_raw" %in% names(dt)) keep_cols <- c(keep_cols, "icv_raw")
dt <- dt[, unique(keep_cols)]

dt$sex <- as.factor(dt$sex)
dt$scanner <- as.factor(dt$scanner)

# -----------------------------
# Optional: discovery/replication split
# If TRAIN_IDS is set, normative models are fit on those IDs only.
# Z-scores are still computed for ALL subjects.
# -----------------------------
train_ids_file <- Sys.getenv("TRAIN_IDS", "")
train_mask <- rep(TRUE, nrow(dt))

if (nzchar(train_ids_file) && file.exists(train_ids_file)) {
  train_ids <- fread(train_ids_file, select = "id")$id
  train_mask <- dt$id %in% as.character(train_ids)
  cat(sprintf("[split] Fitting normative models on %d/%d subjects (discovery)\n",
              sum(train_mask), nrow(dt)))

  # A scanner with no discovery subjects still receives predictions, because
  # subsetting a factor keeps its levels: mgcv sees the level, fits its random
  # effect on nothing, shrinks it to zero and says nothing. Those subjects are
  # then z-scored against a scanner effect that was never estimated. Warn
  # loudly rather than let it pass.
  sc_tab <- table(dt$scanner, ifelse(train_mask, "train", "held_out"))
  empty  <- rownames(sc_tab)[sc_tab[, "train"] == 0 &
                             (if ("held_out" %in% colnames(sc_tab)) sc_tab[, "held_out"] > 0 else FALSE)]
  if (length(empty))
    warning("[split] scanner level(s) with ZERO discovery subjects but held-out ",
            "subjects present: ", paste(empty, collapse = ", "),
            ". Their scanner random effect is not estimated; those deviations ",
            "are not comparable to the rest. Drop them or rebuild the split.",
            call. = FALSE)

  thin <- rownames(sc_tab)[sc_tab[, "train"] > 0 & sc_tab[, "train"] < 30]
  if (length(thin))
    message(sprintf("[split] scanner level(s) with fewer than 30 discovery subjects: %s",
                    paste(thin, collapse = ", ")))
}

# -----------------------------
# Optional: ICV residualization
#
# If icv_raw is present, compute residualized icv (centred residual) and use it
# as nuisance covariate for features that need ICV.
# -----------------------------
if ("icv_raw" %in% names(dt) && any(is.finite(dt$icv_raw))) {
  # Guard the scanner random effect: s(scanner, bs = "re") errors out on a
  # single-level factor, i.e. on any single-site dataset. fit_one_feature()
  # carries the same guard.
  icv_frm <- if (length(levels(dt$scanner)) > 1L) {
    icv_raw ~ s(age, k = 5, bs = "cs") + sex + s(scanner, bs = "re")
  } else {
    icv_raw ~ s(age, k = 5, bs = "cs") + sex
  }
  icv_fit <- mgcv::gam(
    icv_frm,
    data = dt[train_mask, ], method = "REML", na.action = na.exclude
  )
  icv_pred <- predict(icv_fit, newdata = dt, type = "response")
  dt$icv <- as.numeric(scale(dt$icv_raw - icv_pred, center = TRUE, scale = FALSE))
}

# Decide which modalities need ICV (default: SA + SubVol; CT typically does not).
icv_modalities <- cfg$normative$icv_modalities %||% c("SA", "SubVol")
needs_icv <- function(feat) {
  m <- fm$modality[match(feat, fm$feature)]
  !is.na(m) && (m %in% icv_modalities)
}

# -----------------------------
# Optionally save raw values
# -----------------------------
# A verbatim copy of the input feature table, and the largest file this stage
# can produce. Nothing in the pipeline reads it, so it is off by default; set
# outputs.write_raw_values: true if something downstream of yours consumes it.
# The skip is announced rather than silent, so a broken downstream step points
# straight back here.
write_raw_values <- isTRUE(cfg$outputs$write_raw_values %||% FALSE)
if (write_raw_values) {
  write_tsv(dt, out_raw)
  cat("Wrote: ", out_raw, "\n", sep = "")
} else {
  cat("[normative] skipping feats_rawVals.txt (outputs.write_raw_values is false)\n")
}

# -----------------------------
# Feature sizes (raw, pre-winsorization)
# -----------------------------
# Computed on the DISCOVERY subset when a split is active, so that `size` is a
# property of the sample the normative models were fitted on and not of the
# sample being scored.
feat_mean <- vapply(feature_cols,
                    function(f) mean(dt[[f]][train_mask], na.rm = TRUE),
                    numeric(1))
# `size` is the TISSUE VOLUME (mm^3) each feature describes:
#   CT / SA feature for region r -> mean thickness_r * mean area_r
#   SubVol structure s           -> mean volume_s
#
# Volume, not the mean of the feature itself. Because stage 03 pools globally,
# `size` is raised to beta and compared ACROSS modalities, so it must carry one
# unit for every feature. The feature means do not: a thickness is mm (~2.5),
# an area mm^2 (~2000), a volume mm^3 (~5000), and using them directly would
# suppress cortical thickness by an order of magnitude for reasons of units
# alone. At beta = 0 the size term is inert, but it must be coherent for beta
# to mean anything when turned up.
sizes_df <- fm %>% dplyr::mutate(raw_mean = as.numeric(feat_mean[feature]))

cort_vol <- sizes_df %>%
  dplyr::filter(modality %in% c("CT", "SA")) %>%
  dplyr::select(region, modality, raw_mean) %>%
  tidyr::pivot_wider(names_from = modality, values_from = raw_mean) %>%
  dplyr::mutate(cortical_volume = CT * SA) %>%
  dplyr::select(region, cortical_volume)

sizes_df <- sizes_df %>%
  dplyr::left_join(cort_vol, by = "region") %>%
  dplyr::mutate(size = ifelse(modality %in% c("CT", "SA"),
                              cortical_volume, raw_mean)) %>%
  dplyr::select(modality, region, size)

n_bad <- sum(!is.finite(sizes_df$size) | sizes_df$size <= 0)
if (n_bad)
  warning("feats_avgSize: ", n_bad, " feature(s) lack a usable size ",
          "(cortical region missing its CT or SA partner?); align_feature_set() ",
          "will drop them", call. = FALSE)

write_tsv(sizes_df, out_sizes)

# -----------------------------
# Winsorize feature values for modeling
# -----------------------------
dt_raw <- dt
dt[feature_cols] <- lapply(dt[feature_cols], winsorize, lower = 0.001, upper = 0.999)

# -----------------------------
# Fit normative models
# -----------------------------
# One payload per feature: the covariates plus that feature's winsorised and
# raw columns, and nothing else.
#
# The obvious formulation — future_lapply over feature names, with the function
# closing over dt and dt_raw — makes the whole subject x feature matrix a
# global. Every worker then receives a copy of both tables (at 95k subjects and
# 160 features, roughly 250 MB each way) in order to fit a model that reads one
# feature column and four covariates. It is also a scaling cliff, because
# future.globals.maxSize defaults to 500 MiB and a wider feature table trips it
# with an error that says nothing about the cause.
#
# Iterating over the payloads instead means future ships each worker only the
# chunk it is about to fit. It also makes newdata for predict() four columns
# wide rather than the full table, which mgcv would otherwise rebuild a model
# frame from on every one of the 160 calls.
covar_cols <- c("age", "sex", "scanner")
if ("icv" %in% names(dt)) covar_cols <- c(covar_cols, "icv")
covars <- dt[, covar_cols, drop = FALSE]
ids    <- dt$id

feat_payload <- lapply(feature_cols, function(f) list(
  feat     = f,
  y        = dt[[f]],
  y_raw    = dt_raw[[f]],
  with_icv = isTRUE(needs_icv(f))
))

# dt / dt_raw are not referenced past this point; dropping them keeps them out
# of the exported globals and halves peak memory in the parent process.
rm(dt, dt_raw); invisible(gc(FALSE))

cat(sprintf("Fitting %d features with %d workers...\n", length(feature_cols), n_cores))
res_list <- future_lapply(
  feat_payload,
  function(p) {
    d  <- covars
    d[[p$feat]] <- p$y
    dr <- stats::setNames(data.frame(p$y_raw, stringsAsFactors = FALSE), p$feat)
    fit_one_feature(
      feat           = p$feat,
      dt             = d,
      dt_raw         = dr,
      k_age          = k_age,
      use_site_RE    = use_site_RE,
      sex_interact   = sex_interact,
      needs_icv      = function(...) p$with_icv,
      train_mask     = train_mask,
      variance_model = variance_model,
      k_age_scale    = k_age_scale
    )
  },
  future.seed = TRUE
)

future::plan(future::sequential)
rm(feat_payload); invisible(gc(FALSE))

qc_df <- do.call(rbind, lapply(res_list, function(r) {
  if (!is.list(r)) return(NULL)
  data.frame(
    feature        = r$feat,
    ok             = isTRUE(r$ok),
    variance_model = r$variance_model %||% NA,
    sigma          = r$sigma %||% NA_real_,      # median sigma_i under gaulss
    sigma_p10      = r$sigma_p10 %||% NA_real_,
    sigma_p90      = r$sigma_p90 %||% NA_real_,
    sigma_ratio    = r$sigma_ratio %||% NA_real_, # p90/p10; 1.0 = homoscedastic
    edf_age        = r$edf_age %||% NA_real_,
    edf_age_scale  = r$edf_age_scale %||% NA_real_,
    r2             = r$r2 %||% NA_real_,
    n_train        = r$n_train %||% NA,
    msg            = r$msg %||% NA_character_,
    stringsAsFactors = FALSE
  )
}))

write_tsv(qc_df, out_qc)

# ---------------------------------------------------------------------------
# Loud accounting of the feature set.
#
# A feature that fails here never enters the z matrix, so it never reaches
# align_feature_set() in stage 03 either. Without this block its only trace
# would be one row in normative_qc.txt, and the number of features in the index
# would be an OUTCOME of model fitting rather than a specification — free to
# differ between runs, or between the hemisphere subsets used for the
# split-half. Stage 03 reconciles the full count against feats_avgSize.txt,
# which carries every feature in the map whether fitted or not; this block
# makes the same count visible at the point of failure.
# ---------------------------------------------------------------------------
n_fail <- sum(!qc_df$ok)
cat(sprintf("[normative] %d/%d features fitted; %d failed\n",
            sum(qc_df$ok), nrow(qc_df), n_fail))
if (n_fail) {
  cat("[normative] FEATURES DROPPED AT STAGE 02 (they will not reach the index):\n")
  print(qc_df[!qc_df$ok, c("feature", "edf_age", "r2", "msg")], row.names = FALSE)
}
if (variance_model) {
  fell_back <- sum(!is.na(qc_df$msg) & grepl("^gaulss failed", qc_df$msg))
  if (fell_back)
    warning("[normative] ", fell_back, " feature(s) fell back to the location-only ",
            "model because gaulss did not converge. The run is a MIXTURE of the two ",
            "models; see the msg column of normative_qc.txt.", call. = FALSE)
  sr <- qc_df$sigma_ratio[qc_df$ok & is.finite(qc_df$sigma_ratio)]
  if (length(sr))
    cat(sprintf("[normative] heteroscedasticity across subjects (sigma p90/p10): median %.2f, max %.2f\n",
                median(sr), max(sr)))
}

ok_res <- Filter(function(x) is.list(x) && isTRUE(x$ok), res_list)

if (!length(ok_res))
  stop("No feature fitted successfully at stage 02. See ", out_qc,
       " for the per-feature reason.", call. = FALSE)

# Assemble by building the column list first and converting once. Repeated
# `Z[[name]] <- value` on a data.table reallocates on each assignment, which at
# 160 features and 95k rows is 160 copies of the whole table.
Z  <- data.table::setDT(c(list(id = ids),
                          stats::setNames(lapply(ok_res, `[[`, "z"),
                                          vapply(ok_res, `[[`, character(1), "feat"))))
MU <- data.table::setDT(c(list(id = ids),
                          stats::setNames(lapply(ok_res, `[[`, "mu"),
                                          vapply(ok_res, `[[`, character(1), "feat"))))

write_tsv(Z,  out_z)
write_tsv(MU, out_mu)

# Long format using feature_map. melt() rather than pivot_longer(): this is the
# largest object the pipeline handles (subjects x features, ~15M rows at UK
# Biobank scale) and tidyr materialises considerably more of it at once.
zs_long <- data.table::melt(Z, id.vars = "id", variable.name = "feature",
                            value.name = "z", variable.factor = FALSE)
fm_dt   <- data.table::as.data.table(fm)[, .(feature, modality, region)]
zs_long <- merge(zs_long, fm_dt, by = "feature", all.x = TRUE, sort = FALSE)
zs_long <- zs_long[!is.na(modality) & !is.na(region)]
# Sort explicitly. Nothing downstream depends on row order, but a deterministic
# order makes two runs byte-comparable, which is the cheapest reproducibility
# check there is.
data.table::setorder(zs_long, id, feature)
zs_long <- zs_long[, .(id, modality, region, z)]

write_tsv(zs_long, out_long)
cat("Wrote: ", out_long, "\n")

write_run_info(out_stage, "02_run_normdevs", opt$config, params = list(
  k_age = k_age, k_age_scale = k_age_scale, variance_model = variance_model,
  use_site_RE = use_site_RE, sex_interact = sex_interact,
  icv_modalities = icv_modalities,
  n_features_in = length(feature_cols), n_features_fitted = length(ok_res),
  n_subjects = length(ids), train_ids = if (nzchar(train_ids_file)) train_ids_file else "none"
))
