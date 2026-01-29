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
# Optional: ICV residualization
#
# If icv_raw is present, compute residualized icv (centred residual) and use it
# as nuisance covariate for features that need ICV.
# -----------------------------
if ("icv_raw" %in% names(dt) && any(is.finite(dt$icv_raw))) {
  icv_fit <- mgcv::gam(
    icv_raw ~ s(age, k = 5, bs = "cs") + sex + s(scanner, bs = "re"),
    data = dt, method = "REML", na.action = na.exclude
  )
  # IMPORTANT: keep length == nrow(dt). With na.exclude, residuals align to input rows with NAs.
  dt$icv <- as.numeric(scale(stats::resid(icv_fit, type = "response"), center = TRUE, scale = FALSE))
}

# Decide which modalities need ICV (default: SA + SubVol; CT typically does not).
icv_modalities <- cfg$normative$icv_modalities %||% c("SA", "SubVol")
needs_icv <- function(feat) {
  m <- fm$modality[match(feat, fm$feature)]
  !is.na(m) && (m %in% icv_modalities)
}

# -----------------------------
# Save raw values (debugging / downstream)
# -----------------------------
write.table(dt, file = out_raw, sep = "\t", quote = FALSE, row.names = FALSE)

# -----------------------------
# Feature sizes (raw, pre-winsorization)
# -----------------------------
feat_mean <- vapply(feature_cols, function(f) mean(dt[[f]], na.rm = TRUE), numeric(1))
sizes_df <- fm %>%
  dplyr::mutate(size = as.numeric(feat_mean[feature])) %>%
  dplyr::select(modality, region, size)

write.table(sizes_df, file = out_sizes, sep = "\t", quote = FALSE, row.names = FALSE)

# -----------------------------
# Winsorize feature values for modeling
# -----------------------------
dt_raw <- dt
dt[feature_cols] <- lapply(dt[feature_cols], winsorize, lower = 0.001, upper = 0.999)

# -----------------------------
# Fit normative models
# -----------------------------
cat(sprintf("Fitting %d features with %d workers...\n", length(feature_cols), n_cores))
res_list <- future_lapply(
  feature_cols,
  function(f) fit_one_feature(
    feat = f,
    dt = dt,
    dt_raw = dt_raw,
    k_age = k_age,
    use_site_RE = use_site_RE,
    sex_interact = sex_interact,
    needs_icv = needs_icv
  ),
  future.seed = TRUE
)

future::plan(future::sequential)

qc_df <- do.call(rbind, lapply(res_list, function(r) {
  if (!is.list(r)) return(NULL)
  data.frame(
    feature = r$feat,
    ok      = isTRUE(r$ok),
    sigma   = r$sigma %||% NA_real_,
    edf_age = r$edf_age %||% NA_real_,
    r2      = r$r2 %||% NA_real_,
    n_train = r$n_train %||% NA,
    msg     = r$msg %||% NA_character_,
    stringsAsFactors = FALSE
  )
}))

write.table(qc_df, file = out_qc, sep = "\t", quote = FALSE, row.names = FALSE)

ok_res <- Filter(function(x) is.list(x) && isTRUE(x$ok), res_list)

Z <- data.table(id = dt$id)
MU <- data.table(id = dt$id)
for (r in ok_res) {
  Z[[r$feat]]  <- r$z
  MU[[r$feat]] <- r$mu
}

write.table(Z,  file = out_z,  sep = "\t", quote = FALSE, row.names = FALSE)
write.table(MU, file = out_mu, sep = "\t", quote = FALSE, row.names = FALSE)

# Long format using feature_map
Z_df <- as.data.frame(Z)
zs_long <- Z_df %>%
  tidyr::pivot_longer(
    cols = -id,
    names_to = "feature",
    values_to = "z"
  ) %>%
  dplyr::left_join(fm, by = "feature") %>%
  dplyr::filter(!is.na(modality) & !is.na(region)) %>%
  dplyr::select(id, modality, region, z)

write.table(zs_long, file = out_long, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Wrote: ", out_long, "\n")
