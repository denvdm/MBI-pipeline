# MBI pipeline — rebuild changelog

Replaces 16 files. Unzip over the repo root; directory layout matches.
Files NOT in this package are unchanged from the original repo.

**Nothing here has been run.** No R in the authoring environment: syntax is
brace-balance-checked and the weight transform was verified by an independent
Python reimplementation, but the first real test is
`Rscript examples/toy/run_toy_example.R`.

---

## MUST DO BEFORE RUNNING

**Rerun stage 01 before stage 03.** The `weight` column in
`pet_fdg_featMedians.txt` changed meaning: it is now the RAW parcel median, not
a pre-transformed weight. Pointing the new code at an old file makes robust-z
transform an already-transformed vector. It will not error. It will be wrong.
(New files carry a `coverage` column; old ones do not.)

**Recalibrate.** `beta = 0.25` was selected against a cortex-only target under
a transform that no longer exists. `gamma` no longer exists at all. Run
`extras/calibrate_mbi_grid.R` before treating any result as final.

---

## Fix #1 — one feature set for both indices

Regions lacking a valid map weight were silently zeroed in MBI_raw but carried
at full weight in GBI, so their deviations loaded on MBI purely by
construction (CC_Posterior + CC_Anterior ≈ 4.8% of the index).

- `align_feature_set()` keeps only regions with a finite positive map weight
  AND deviation data AND a size entry; drops from BOTH; reports what went.
- `compute_flagship_abg()` refuses to residualise if the two weight tables
  disagree — the original bug was silent, and now cannot be.
- `index_concentration()` reports N_eff (effective number of contributing
  features) from weights alone, so candidate designs can be screened without
  subject data.

## Fix #2 — bounded, map-invariant weight transform

Old chain, split across two files: z-score → subtract min → normalise (stage
01), then winsorise → normalise → `^gamma` (compute_mbi). This pinned the whole
scale to the single coldest parcel (1e-6, → ~1e-15 after gamma: present in GBI,
absent from MBI_raw), depended on the map's distributional shape, and let one
outlier parcel move every weight.

New, entirely inside `transform_map_weights()`:

    u_r = clamp((x_r − median) / (clamp_z · MAD), ±1)
    w_r = ratio^(u_r / 2),  mean-normalised

`ratio` = MAX weight ratio, highest- vs lowest-demand region. `gamma` and
`winsor_q` are retired: `w^gamma == (ratio^gamma)^(u/2)`, so exponentiating IS
a change of ratio, and a clamped robust scale has no outliers to winsorise.
Both survive under `map_method = "legacy"`, which reproduces the old build.

Verified numerically (Python reimplementation):

| check | result |
|---|---|
| constant vector (GBI) → all ones | exact |
| gamma-absorption identity | max abs diff 2.2e-15 |
| one cold outlier → change in other weights | robust_z 0.004 vs legacy 10.7 |
| legacy realised weight ratio, clean data | ~6.7e16 |

Stage 01 now emits the raw parcel median plus `coverage`/`nvox_used`/
`nvox_total`. It no longer decides anything about weighting.

## SubVol gate: quantile → absolute

Was `quantile(size, 0.20)`, recomputed over whatever SubVol set survived
upstream filtering — so dropping a region for unrelated reasons silently
shifted which OTHER structures were gated out. Now `subvol_min_mm3: 1000`,
which is order-invariant, so the gate and map-validity filtering commute.

1000 mm³ sits just below amygdala (1590), the smallest structure surviving
map-validity filtering, and near the old effective cut (752). The principled
replacement is a per-feature test–retest ICC gate across all modalities; a
dormant config hook (`inputs$feat_icc`, `mbi$icm_min`) is in `03_run_mbi.R`.

## Parcel homogeneity (stage 01)

Stage 01 now carries `p05` / `p95` through from `pet_parcel_extract()` (both
were already computed and then discarded) and reports

    dispersion_r = (p95_r - p05_r) / median_r

per parcel, printed highest-first in the stage-01 log. High dispersion means
the parcel median is a poor summary of the parcel - a large heterogeneous
parcel can have a perfectly valid median that represents nothing in
particular. Outcome-blind and map-side, so it generalises to any prior.

Nothing is excluded on this basis automatically. The ranking is available
after stage 01 and before stage 03, so exclusions can be decided in-flight,
from the number rather than from region identity.

## Calibration grid

Sweeps **(ratio × beta)** instead of (gamma × beta) and reports per cell:
`mbi_map_r`, `n_eff`, `top_share`, `r2_raw_gbi`, `realised_ratio`. Prints the
unconstrained best (labelled do-not-adopt-blindly), the recommended cell
subject to `n_eff >= 30` AND `R2(raw,GBI) in [0.75, 0.95]`, and a concentration
frontier — best alignment achievable at each N_eff decile.

`mods_use` changed from `c("CT","SA")` to all three modalities. The r_tilt
calibration target was computed on CORTEX ONLY while the index it tunes is
overwhelmingly subcortical (SubVol r = 0.46 with MBI; CT 0.03, SA 0.01).

## De-duplication

`transform_map_weights()` is now the single definition; `compare_pet_funcs.R`
aliases it instead of carrying a copy whose comment said it "mirrors" the
index. `05_compare_pet_templates.R` gained `source_local("compute_mbi.R")`.

## General bug fixes

- `pet_parcel_extract.R` — winsorisation bounds were SWAPPED (`ql` took the
  99th percentile as the lower bound). Dormant (`winsor_q = NULL`) but the rest
  of the pipeline uses the 0.99 convention. Now accepts either.
- `compute_mbi.R` — `lm()` used default `na.omit`, so `residuals()` could
  return a short vector and misalign subjects on assignment. Now `na.exclude`.
- `compute_mbi.R` — supplied-`alpha` branch: YAML delivers a list, so
  `alpha[modality]` produced a list-column; and normalisation summed over
  modalities the subject lacks, shrinking scores for incomplete subjects.
  Both fixed, with named-vector validation.
- `01_extract_petvals.R` — a positional `colnames()<-` duplicated the
  `rename()` above it, relabelling BY POSITION. Deleted. `make.names(unique =
  TRUE)` silently renamed duplicate regions to `region.1`, matching nothing
  downstream; duplicates now error.
- `02_run_normdevs.R` — size weights computed on the full sample though
  Methods say discovery; now discovery. ICV GAM's `s(scanner, bs="re")` errored
  on single-site data; guarded, as `fit_one_feature()` already was.
- `run_helpers.R` — `start_logging()` was a no-op: its `on.exit` fired on
  return, closing both sinks immediately. Returns a closer instead.
- `make_split.R` — `cut()` errored on duplicated age-decile breaks.

## Known, deliberately NOT changed

1. `fit_one_feature()` drops features from the index via its `edf_age` / `r2`
   gate without recording it as a feature-set decision. Check
   `normative_qc.txt`; belongs in Methods either way.
2. Full-sample leakage: winsorisation quantiles and the GBI residualisation
   slope are estimated across discovery+replication. Negligible at n=83k, but
   nameable by a reviewer.
3. `collapse_zs()` in stage 04 averages z across CT/SA/SubVol per region before
   correlating with PET — different variances mixed 1:1. Fine as QC, not the index.
4. `prep_features.R` falls back to positional column naming if the demographics
   header does not match — same silent-relabelling pattern deleted in stage 01.
