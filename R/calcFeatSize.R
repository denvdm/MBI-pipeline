# Computes feature-size or region-size metrics used for PET weighting or scaling.
# Supports volume- or area-based weighting strategies in downstream aggregation.

rm(list=ls())

###### load packages and set paths
pacman::p_load(dplyr,data.table,stringr,readr)

path <- "/cluster/projects/p33/users/dennisva/"
setwd(path)

# ===== Build cohort-derived region size table (per-feature average within cohort) =====

dt <- fread(paste0(path,"integrity/df/feats_rawVals.txt"),data.table=F)
feature_cols <- colnames(dt)

# exclude <- c("MRID","EstimatedTotalIntraCranialVol")
# feature_cols <- feature_cols[!feature_cols %in% exclude]

# Identify feature columns by suffix
ct_cols  <- grep("_thickness$", feature_cols, value = TRUE)
sa_cols  <- grep("_area$",      feature_cols, value = TRUE)
sv_cols  <- grep("_aseg$",      feature_cols, value = TRUE)

# Helper to compute cohort means safely
col_mean <- function(df, cols) vapply(cols, function(cn) mean(df[[cn]], na.rm = TRUE), numeric(1))

# --- Cortical parcel areas (used as size for CT and SA) ---
# Map thickness/area to a common "region stem" like 'lh-bankssts' / 'rh-superiorfrontal'
stem_from_col <- function(x) {
  # strip modality suffix then replace first "_" with "-" to mark hemisphere cleanly
  base <- sub("_(thickness|area)$", "", x)
  #sub("_", "-", base, fixed = TRUE)
}

# For SA we can read area directly; if some areas are missing, fall back to area derived from thickness stem
sa_means <- col_mean(dt, sa_cols)
sa_tbl <- tibble(
  col = names(sa_means),
  region = stem_from_col(names(sa_means)),
  size = as.numeric(sa_means)
) %>% select(region, size)

# Ensure CT has a matching area "size": use SA’s size by matching stems
ct_regions <- unique(stem_from_col(ct_cols))
ct_tbl <- tibble(region = ct_regions) %>%
  left_join(sa_tbl, by = "region")

# --- Subcortical structure volumes (aseg) ---
sv_means <- col_mean(dt, sv_cols)
sv_tbl <- tibble(
  region_raw = names(sv_means),
  region = sub("_aseg$", "", names(sv_means)),
  size = as.numeric(sv_means)
) %>% select(region, size)

# Optional: drop any undesired aseg (e.g., corpus callosum segments)
# drop_prefix <- c("CC_", "wm", "WM", "Unknown", "Optic", "hypointensities")
# sv_tbl <- sv_tbl %>%
#   filter(!str_starts(region, paste0("^(", paste(drop_prefix, collapse="|"), ")")))

# --- Assemble size_df (modality-specific) ---
size_df <- bind_rows(
  ct_tbl %>% mutate(modality = "CT") %>% select(region, modality, size),
  sa_tbl %>% mutate(modality = "SA") %>% select(region, modality, size),
  sv_tbl %>% mutate(modality = "SubVol") %>% select(region, modality, size)
)

# Save
write.table(size_df,file=paste0(path,"integrity/df/feats_avgSize.txt"),sep = "\t",quote=FALSE,row.names=FALSE)

