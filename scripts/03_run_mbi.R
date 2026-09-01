#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})

source("helpers/run_helpers.R")
opt <- parse_common_args()
cfg <- read_config(opt$config)

source_local("compute_mbi.R")

# in/out
normdev_long <- cfg$inputs$deviations_z_long %||% file.path(opt$outdir, "normdev", "deviations_z_long.txt")
sizes_file   <- cfg$inputs$feat_sizes        %||% file.path(opt$outdir, "normdev", "feats_avgSize.txt")
fdg_file     <- cfg$inputs$pet_weights       %||% file.path(opt$outdir, "extract_pet", "pet_fdg_featMedians.txt")

assert_file_exists(normdev_long, "deviations_z_long")
assert_file_exists(sizes_file,   "feat_sizes (feats_avgSize.txt)")
assert_file_exists(fdg_file,     "pet_weights (pet_fdg_featMedians.txt)")

out_stage    <- stage_dir(opt$outdir, "compute_mbi")
out_scores   <- file.path(out_stage, "mbi_flagship_scores.txt")
out_features <- file.path(out_stage, "mbi_feature_manifest.txt")
out_shares   <- file.path(out_stage, "mbi_influence_shares.txt")

if (file.exists(out_scores) && !opt$overwrite) {
  stop("Output exists: ", out_scores, " (use --overwrite)", call. = FALSE)
}

zs      <- fread(normdev_long, data.table = FALSE) %>%
  mutate(region = as.character(region), modality = as.character(modality))
size_df <- fread(sizes_file, data.table = FALSE) %>%
  mutate(region = as.character(region), modality = as.character(modality))
fdg     <- fread(fdg_file, data.table = FALSE) %>%
  mutate(region = as.character(region))

# Index parameters, overridable via config. See config_generic.yaml for what
# each one does and helpers/compute_mbi.R for the algebra behind them.
ratio        <- cfg$mbi$ratio        %||% "auto" # INERT for MBI; "auto" = read off the map
clamp_z      <- cfg$mbi$clamp_z      %||% 3      # safety rail on outlier parcels
emphasis     <- cfg$mbi$emphasis     %||% 1      # the only live weighting parameter
beta         <- cfg$mbi$beta         %||% 0      # 0 = no size weighting
coverage_thr <- cfg$mbi$coverage_thr %||% 0.80

# Subjects must contribute ALL of these modalities to receive a score. NULL
# scores anyone who clears coverage_thr on at least one modality and warns,
# because MBI is a contrast: a different feature set is a different linear
# functional, not a noisier estimate of the same one.
required_modalities <- cfg$mbi$required_modalities %||% NULL

# Keys the index does not read. A config carrying one of these describes a
# construction that is not the one being run, which is worth saying out loud.
unused <- intersect(c("pooling","map_method","alpha","alpha_theta","gamma","winsor_q"),
                    names(cfg$mbi))
if (length(unused))
  warning("config: the following mbi keys are not read by this pipeline and are ",
          "IGNORED: ", paste(unused, collapse = ", "),
          ". Remove them so the config documents what actually ran.",
          call. = FALSE)

# ---------------------------------------------------------------------------
# SubVol reliability gate.
#
# Very small aseg structures are segmented least reliably, so their deviations
# are dominated by measurement noise. Size stands in for reliability here; a
# per-feature test-retest ICC gate across ALL modalities is the better
# instrument where repeat scans exist, and the config hook for it is below.
#
# The threshold is ABSOLUTE (mm^3), deliberately not a quantile. A quantile cut
# would be recomputed over whatever SubVol set survived upstream filtering, so
# dropping a structure for an unrelated reason — no valid map weight, say —
# would silently shift which OTHER structures were gated out. Being absolute,
# this gate and the map-validity filter in align_feature_set() commute, and
# their order does not have to be defended.
#
# !! SET THIS FROM THE OBSERVED VOLUME DISTRIBUTION OF YOUR OWN SEGMENTATION.
# !! The default of 1000 mm^3 sits just below the amygdala in an aparc+aseg
# !! parcellation; it is not a universal constant.
# ---------------------------------------------------------------------------
subvol_min_mm3 <- cfg$mbi$subvol_min_mm3 %||% 1000

if (!is.null(cfg$mbi$subvol_gate_pct))
  warning("config: 'subvol_gate_pct' is not read by this pipeline and is ignored; ",
          "using absolute 'subvol_min_mm3' = ", subvol_min_mm3, call. = FALSE)

# Reserved for the ICC-based gate (all modalities). NULL = not in use.
icc_file    <- cfg$inputs$feat_icc  %||% NULL
icc_min     <- cfg$mbi$icc_min      %||% NULL

sub_sizes <- size_df %>% dplyr::filter(modality == "SubVol")
keep_sub  <- sub_sizes %>% dplyr::filter(size >= subvol_min_mm3) %>% dplyr::pull(region)
drop_sub  <- setdiff(sub_sizes$region, keep_sub)

message(sprintf("SubVol gate: keeping %d/%d structures at >= %.0f mm^3%s",
                length(keep_sub), nrow(sub_sizes), subvol_min_mm3,
                if (length(drop_sub))
                  paste0(" (dropped: ", paste(sort(drop_sub), collapse = ", "), ")") else ""))

# Feature provenance: record the full expected feature set BEFORE any gate, so
# the manifest can reconcile every feature in the map to its fate. feats_avgSize
# carries every feature the feature map defined, fitted or not, which makes it
# the only complete reference available at this stage.
prov <- size_df %>%
  dplyr::distinct(modality, region) %>%
  dplyr::mutate(status = NA_character_)

fitted_feats <- zs %>% dplyr::distinct(modality, region) %>% dplyr::mutate(.fit = TRUE)
prov <- prov %>%
  dplyr::left_join(fitted_feats, by = c("modality","region")) %>%
  dplyr::mutate(status = ifelse(is.na(.fit), "dropped: no normative fit (stage 02)", status)) %>%
  dplyr::select(-.fit)

# Two distinct reasons a SubVol feature can fail the gate, recorded separately:
# it can be genuinely too small, or it can have no usable size at all. Both
# leave it out of keep_sub, and lumping them together puts a wrong explanation
# in the one file whose purpose is to be trusted without re-derivation.
sub_nosize <- sub_sizes %>%
  dplyr::filter(!is.finite(as.numeric(size)) | as.numeric(size) <= 0) %>%
  dplyr::pull(region)
gated_small <- setdiff(setdiff(sub_sizes$region, keep_sub), sub_nosize)

prov <- prov %>%
  dplyr::mutate(status = dplyr::case_when(
    !is.na(status)                                        ~ status,
    modality == "SubVol" & region %in% sub_nosize         ~ "dropped: SubVol has no usable size",
    modality == "SubVol" & region %in% gated_small        ~ sprintf("dropped: SubVol < %.0f mm3", subvol_min_mm3),
    TRUE                                                  ~ status
  ))

zs <- zs %>% dplyr::filter(modality != "SubVol" | region %in% keep_sub)

if (!is.null(icc_file) && !is.null(icc_min)) {
  assert_file_exists(icc_file, "feat_icc")
  icc <- fread(icc_file, data.table = FALSE) %>%
    mutate(region = as.character(region), modality = as.character(modality))
  keep_icc <- icc %>% dplyr::filter(is.finite(icc), icc >= icc_min) %>%
    dplyr::select(region, modality)
  before_icc <- dplyr::distinct(zs, region, modality)
  zs <- zs %>% dplyr::semi_join(keep_icc, by = c("region","modality"))
  dropped_icc <- dplyr::anti_join(before_icc, dplyr::distinct(zs, region, modality),
                                  by = c("region","modality"))
  message(sprintf("ICC gate: %d -> %d features at ICC >= %.2f",
                  nrow(before_icc), nrow(dplyr::distinct(zs, region, modality)), icc_min))
  # The gate is per (modality, region), so it has to be recorded per (modality,
  # region) too. Without this the manifest labels an ICC-dropped feature by
  # whatever happened to its region, which is the wrong fate under the wrong
  # reason in the one file meant to be trusted without re-derivation.
  if (nrow(dropped_icc))
    prov <- prov %>%
      dplyr::left_join(dropped_icc %>% dplyr::mutate(.icc_out = TRUE),
                       by = c("modality","region")) %>%
      dplyr::mutate(status = ifelse(is.na(status) & !is.na(.icc_out),
                                    sprintf("dropped: ICC < %.2f", icc_min), status)) %>%
      dplyr::select(-.icc_out)
}

# ---------------------------------------------------------------------------
# IMPLAUSIBLE DEVIATIONS.
#
# Nothing upstream bounds z. Stage 02 fits on winsorised values but computes z
# from the RAW value, deliberately — a genuinely deviant subject should keep a
# large deviation — and `clamp_z` bounds the MAP, not the deviation. So a failed
# segmentation enters the contrast at full weight, and because MBI is a weighted
# sum over ~10^2 features with no single one dominating by design, one absurd
# value dominates by accident: on a 143-feature build at N_eff 52, a single
# feature at |z| = 56 moved that subject's index by 15 SD.
#
# SET TO MISSING, NOT CLIPPED. Clipping asserts a value — it says this region
# really is z_implausible SD deviant — and the index then propagates that
# assertion: on the same build, clipping at 5 still left a segmentation failure
# moving the subject by 1.7 SD, and at 10 by 3.0 SD. Setting the observation to
# NA asserts ignorance, which is the true state, and moves the subject by 0.36
# SD, the renormalisation from dropping one feature of 143. The existing rules
# then do the rest: coverage_thr decides whether the subject still contributes
# that modality, required_modalities whether they are scored at all.
#
# CHOOSE THE BOUND OUTSIDE THE BIOLOGICAL RANGE, not at its edge. Severe atrophy
# legitimately reaches |z| of 5 to 8 in a few regions, and those are the
# subjects the index exists to detect. On the reference cohort the 99.9th
# percentile of |z| is 4.17 and the 99.99th is 8.20, so a bound of 10 removes
# only the far tail. There is no universal value: read your own |z| distribution
# before setting one.
z_implausible <- cfg$mbi$z_implausible %||% NULL

if (!is.null(z_implausible)) {
  z_implausible <- as.numeric(z_implausible)
  if (!is.finite(z_implausible) || z_implausible <= 0)
    stop("config: mbi.z_implausible must be a positive number, or absent to ",
         "disable the rule.", call. = FALSE)
  bad <- which(is.finite(zs$z) & abs(zs$z) > z_implausible)
  if (length(bad)) {
    bad_tab <- zs[bad, ] %>% dplyr::count(modality, region, name = "n_obs") %>%
      dplyr::arrange(dplyr::desc(n_obs))
    message(sprintf(
      "z_implausible = %.2f: %d observation(s) set to NA, across %d subject(s) and %d feature(s) (%.5f%% of deviations)",
      z_implausible, length(bad), dplyr::n_distinct(zs$id[bad]),
      nrow(bad_tab), 100 * length(bad) / sum(is.finite(zs$z))))
    cat("\nFeatures contributing implausible deviations (top 10):\n")
    print(as.data.frame(utils::head(bad_tab, 10)), row.names = FALSE)
    write_tsv(zs[bad, ], file.path(out_stage, "mbi_implausible_z.txt"))
    cat("Wrote: ", file.path(out_stage, "mbi_implausible_z.txt"), "\n", sep = "")
    zs$z[bad] <- NA_real_
  } else {
    message(sprintf("z_implausible = %.2f: no deviation exceeds it.", z_implausible))
  }
} else {
  # No rule set. Say so if the data contains deviations no normative model can
  # plausibly have produced, rather than letting them into the index in silence.
  zmax <- suppressWarnings(max(abs(zs$z), na.rm = TRUE))
  if (is.finite(zmax) && zmax > 10)
    warning("mbi.z_implausible is not set and this cohort contains |z| up to ",
            sprintf("%.1f", zmax), ". Nothing bounds a deviation, so a failed ",
            "segmentation enters the contrast at full weight; see ",
            "mbi_subject_concentration.txt for whose index it carries.",
            call. = FALSE)
}

res <- compute_flagship_abg(
  zs = zs, fdg = fdg, sizes = size_df,
  ratio = ratio, clamp_z = clamp_z, emphasis = emphasis,
  beta = beta, coverage_thr = coverage_thr,
  required_modalities = required_modalities
)

mbi <- as.data.frame(res$scores)
write_tsv(mbi, out_scores)
cat("Wrote: ", out_scores, "\n")

# ---------------------------------------------------------------------------
# Feature manifest, FEATURE-LEVEL and complete.
#
# This file is the answer to "why does the index have fewer features than the
# feature map?", and it has to be derivable without hand-counting. It starts
# from the complete feature list and assigns every feature exactly one fate, so
# it sums to the input count by construction.
#
# The reason it starts from `prov` (built above, from feats_avgSize) rather
# than from what align_feature_set() saw: anything removed earlier — a feature
# whose normative fit failed, a SubVol structure below the size gate — is
# already absent from `zs`, so a manifest derived from `zs` would list it in
# neither the kept nor the dropped column. Keep the pre-gate snapshot.
# ---------------------------------------------------------------------------
# Reconcile against the FEATURES the index was actually built from, keyed on
# (modality, region). align_feature_set() decides by region, but every other
# gate in this stage decides by feature, so a region-level reconciliation
# mislabels any feature dropped for a modality-specific reason.
in_index <- res$mbi_details$weight_table %>%
  dplyr::distinct(modality, region) %>%
  dplyr::mutate(.in_index = TRUE)

prov <- prov %>%
  dplyr::left_join(in_index, by = c("modality","region")) %>%
  dplyr::mutate(status = dplyr::case_when(
    !is.na(status)      ~ status,
    !is.na(.in_index)   ~ "kept",
    TRUE                ~ "dropped: no valid map weight or size (align_feature_set)"
  )) %>%
  dplyr::select(-.in_index)

write_tsv(prov, out_features)
cat("Wrote: ", out_features, "\n")

cat("\nFeature provenance (every feature in the feature map, one fate each):\n")
prov_tab <- prov %>% dplyr::count(status, name = "n_features") %>%
  dplyr::arrange(dplyr::desc(n_features))
print(as.data.frame(prov_tab), row.names = FALSE)
cat(sprintf("TOTAL %d features in -> %d features in the index\n",
            nrow(prov), sum(prov$status == "kept")))

n_index <- sum(prov$status == "kept")
if (n_index != res$concentration$n_features)
  warning("manifest says ", n_index, " features kept but index_concentration() ",
          "counted ", res$concentration$n_features,
          ". These must agree; investigate before using this run.", call. = FALSE)

# Modality composition of the scored sample.
mc <- res$mbi_details$modality_composition
if (!is.null(mc) && nrow(mc)) {
  cat("\nModality composition of scored subjects:\n")
  print(as.data.frame(mc), row.names = FALSE)
  write_tsv(mc, file.path(out_stage, "mbi_modality_composition.txt"))
}

# Per-SUBJECT concentration. Answers "is anyone's index carried by one feature?",
# which nothing else in the pipeline reports and which nothing bounds: z is
# computed from the raw value, so a segmentation failure enters the contrast at
# full weight. Read top_share, not max_abs_z.
sc <- res$subject_concentration
write_tsv(sc, file.path(out_stage, "mbi_subject_concentration.txt"))
cat(sprintf("\nPer-subject concentration over %d subjects: top feature share median %.4f, p99 %.4f, max %.4f\n",
            nrow(sc), stats::median(sc$top_share),
            stats::quantile(sc$top_share, 0.99), max(sc$top_share)))
n_dom <- sum(sc$top_share > 0.25)
if (n_dom)
  cat(sprintf("  %d subject(s) (%.3f%%) have >25%% of their index carried by a single feature; see mbi_subject_concentration.txt\n",
              n_dom, 100 * n_dom / nrow(sc)))
cat(sprintf("  largest |z| entering any subject's index: %.2f\n", max(sc$max_abs_z)))

# Per-feature influence shares + N_eff, for QC and for the methods table.
write_tsv(res$concentration$shares, out_shares)
cat("Wrote: ", out_shares, "\n")
cat(sprintf("N_eff = %.2f over %d features, at the realised slope b = %.4f (%s)\n",
            res$concentration$n_eff, res$concentration$n_features,
            res$concentration$b, res$concentration$b_source))
cat(sprintf("        %.2f at the design-time cosine slope b = %.4f\n",
            res$concentration$n_eff_cosine, res$concentration$b_cosine))
cat(sprintf("ratio = %.3f [%s] (realised %.2f) | clamp_z = %.2f | emphasis = %.2f | beta = %.2f | R2(mbi_raw, gbi) = %.4f\n",
            res$diagnostics$ratio, res$diagnostics$ratio_source,
            res$diagnostics$realised_ratio,
            res$diagnostics$clamp_z, res$diagnostics$emphasis,
            res$diagnostics$beta, res$diagnostics$r2_raw_gbi))

cat("\nNOTE: `mbi` is the residual of mbi_raw on gbi WITHIN THIS SAMPLE. It is\n",
    "      orthogonal to gbi and standardised here and nowhere else. Analyses on\n",
    "      a subset of these subjects should re-residualise within that subset.\n", sep = "")

write_run_info(out_stage, "03_run_mbi", opt$config, params = list(
  emphasis = emphasis, clamp_z = clamp_z,
  ratio = round(res$diagnostics$ratio, 6), ratio_source = res$diagnostics$ratio_source,
  beta = beta,
  coverage_thr = coverage_thr, subvol_min_mm3 = subvol_min_mm3,
  required_modalities = required_modalities,
  n_features_in = nrow(prov), n_features_index = n_index,
  n_subjects = nrow(mbi),
  n_eff = round(res$concentration$n_eff, 4),
  n_eff_cosine = round(res$concentration$n_eff_cosine, 4),
  b_empirical = round(res$concentration$b, 6),
  b_cosine = round(res$concentration$b_cosine, 6),
  r2_mbiraw_gbi = round(res$diagnostics$r2_raw_gbi, 6),
  realised_ratio = round(res$diagnostics$realised_ratio, 4)
))
