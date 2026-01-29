#!/usr/bin/env Rscript
rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})

source("helpers/run_helpers.R")
opt <- parse_common_args()
cfg <- read_config(opt$config)

# --- Inputs (FreeSurfer export) ---
fs_dir <- cfg$paths$freesurfer_stats_dir %||% NULL
demog_file <- cfg$inputs$demographics_file %||% NULL

if (is.null(fs_dir)) stop("Missing config: paths: freesurfer_stats_dir", call. = FALSE)
if (is.null(demog_file)) stop("Missing config: inputs: demographics_file", call. = FALSE)

fs_dir <- sub("/+$", "", fs_dir)
assert_file_exists(demog_file, "demographics_file")

# Expected file names can be overridden in config under inputs.freesurfer_files
fs_files <- cfg$inputs$freesurfer_files %||% list()
f_lh_thick <- fs_files$lh_thickness %||% "lh.thickness.txt"
f_rh_thick <- fs_files$rh_thickness %||% "rh.thickness.txt"
f_lh_area  <- fs_files$lh_area      %||% "lh.area.txt"
f_rh_area  <- fs_files$rh_area      %||% "rh.area.txt"
f_aseg     <- fs_files$aseg         %||% "subcorticalstats.txt"

assert_file_exists(file.path(fs_dir, f_lh_thick), "lh_thickness")
assert_file_exists(file.path(fs_dir, f_rh_thick), "rh_thickness")
assert_file_exists(file.path(fs_dir, f_lh_area),  "lh_area")
assert_file_exists(file.path(fs_dir, f_rh_area),  "rh_area")
assert_file_exists(file.path(fs_dir, f_aseg),     "aseg")

# --- Outputs ---
out_stage <- stage_dir(opt$outdir, "features")
out_wide  <- file.path(out_stage, "features_wide.tsv")
out_map   <- file.path(out_stage, "feature_map.tsv")

if ((file.exists(out_wide) || file.exists(out_map)) && !opt$overwrite) {
  stop("Outputs exist in ", out_stage, " (use --overwrite)", call. = FALSE)
}

# --- Load demographics ---
# Use config inputs.colmap to map columns; defaults assume already generic naming.
colmap <- cfg$inputs$colmap %||% list()
col_id      <- colmap$id      %||% "id"
col_age     <- colmap$age     %||% "age"
col_sex     <- colmap$sex     %||% "sex"
col_scanner <- colmap$scanner %||% "scanner"
col_icv     <- colmap$icv     %||% "icv"

# Allow icv to be explicitly disabled (null in YAML)
if (is.null(col_icv) || identical(col_icv, "null") || identical(col_icv, "NULL") || col_icv == "") {
  col_icv <- NA_character_
}

demog <- fread(demog_file, data.table = FALSE)

# If the expected colnames are not present, fall back to first 4 columns.
if (!all(c(col_id, col_age, col_sex, col_scanner) %in% names(demog))) {
  if (ncol(demog) < 4) stop("demographics_file must have >=4 columns", call. = FALSE)
  colnames(demog)[1:4] <- c(col_id, col_age, col_sex, col_scanner)
}

# Build generic covariate table
if (!is.na(col_icv) && (col_icv %in% names(demog))) {
  demog <- demog %>%
    dplyr::transmute(
      id      = as.character(.data[[col_id]]),
      age     = .data[[col_age]],
      sex     = .data[[col_sex]],
      scanner = .data[[col_scanner]],
      icv_raw = .data[[col_icv]]
    )
} else {
  demog <- demog %>%
    dplyr::transmute(
      id      = as.character(.data[[col_id]]),
      age     = .data[[col_age]],
      sex     = .data[[col_sex]],
      scanner = .data[[col_scanner]]
    )
}

# --- Load FreeSurfer stats ---
# We assume the first column is a subject identifier (MRID-like).
read_fs <- function(path) {
  x <- fread(path, header = TRUE, data.table = FALSE)
  colnames(x)[1] <- "MRID"
  x
}

lh_thick <- read_fs(file.path(fs_dir, f_lh_thick))
rh_thick <- read_fs(file.path(fs_dir, f_rh_thick))
thick <- dplyr::full_join(lh_thick, rh_thick, by = "MRID")

lh_area <- read_fs(file.path(fs_dir, f_lh_area))
rh_area <- read_fs(file.path(fs_dir, f_rh_area))
area <- dplyr::full_join(lh_area, rh_area, by = "MRID")

cortex <- dplyr::full_join(thick, area, by = "MRID")

aseg <- read_fs(file.path(fs_dir, f_aseg))

# Pragmatic aseg selection (matches your current PET/LUT conventions)
asegNonhemi <- c(
  "EstimatedTotalIntraCranialVol","3rd-Ventricle","4th-Ventricle","Brain-Stem",
  "CC_Posterior","CC_Mid_Posterior","CC_Central","CC_Mid_Anterior","CC_Anterior"
)
asegHemi <- c(
  "Lateral-Ventricle","Inf-Lat-Vent","Cerebellum-Cortex","Thalamus-Proper","Caudate",
  "Putamen","Pallidum","Hippocampus","Amygdala","Accumbens-area","VentralDC"
)

keep_aseg <- c("MRID", asegNonhemi, paste0("Left-",asegHemi), paste0("Right-",asegHemi))
keep_aseg <- keep_aseg[keep_aseg %in% names(aseg)]
aseg <- aseg[, keep_aseg, drop = FALSE]

# Suffix aseg features (except ICV) for modality mapping
if ("EstimatedTotalIntraCranialVol" %in% names(aseg)) {
  aseg_cols_to_suffix <- setdiff(names(aseg), c("MRID", "EstimatedTotalIntraCranialVol"))
  names(aseg)[match(aseg_cols_to_suffix, names(aseg))] <- paste0(aseg_cols_to_suffix, "_aseg")
}

all <- dplyr::full_join(aseg, cortex, by = "MRID")

# Normalize subject IDs: strip a common FS_ prefix if present, and sanitize column names.
all$MRID <- gsub("^FS_", "", all$MRID)
names(all) <- gsub("-", "_", names(all))
names(all) <- make.names(names(all), unique = TRUE)

# Join demographics: demog$id must match MRID after normalization.
dt <- dplyr::inner_join(demog, all, by = c("id" = "MRID"))

# If the FS aseg ICV exists but demog did not provide icv_raw, create it.
if (!("icv_raw" %in% names(dt)) && ("EstimatedTotalIntraCranialVol" %in% names(dt))) {
  dt <- dt %>% dplyr::rename(icv_raw = EstimatedTotalIntraCranialVol)
}

# Feature set: everything besides covariates
cov_cols <- c("id", "age", "sex", "scanner")
if ("icv_raw" %in% names(dt)) cov_cols <- c(cov_cols, "icv_raw")
feature_cols <- setdiff(names(dt), cov_cols)

# Feature map (suffix-based modality inference)
modality <- ifelse(grepl("_thickness$", feature_cols, ignore.case = TRUE), "CT",
            ifelse(grepl("_area$",     feature_cols, ignore.case = TRUE), "SA",
            ifelse(grepl("_aseg$",     feature_cols, ignore.case = TRUE), "SubVol", NA)))

region <- gsub("(_thickness|_area|_aseg)$", "", feature_cols, ignore.case = TRUE)

fm <- data.frame(
  feature = feature_cols,
  modality = modality,
  region = region,
  stringsAsFactors = FALSE
)
fm <- fm[!is.na(fm$modality), ]

# Restrict wide table to mapped features + covariates
# Note: keep covariates first for readability.
dt_wide <- dt[, c(cov_cols, fm$feature), drop = FALSE]

write.table(dt_wide, file = out_wide, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(fm,      file = out_map,  sep = "\t", quote = FALSE, row.names = FALSE)

cat("Wrote:\n- ", out_wide, "\n- ", out_map, "\n", sep = "")
