# Selection rule for the MBI map weighting

This file is the parameter audit trail. The rule below was committed before the
parameter sweep was run. Amendments are recorded as dated decisions with the
grounds they were taken on, so that a reader can see which parameters were
pre-specified, which were changed afterwards, and why.

**On outcome blindness.** No outcome variable was consulted in decisions 1 to 5;
every parameter in the index was fixed on properties of the map and of the
weight vector alone. Decision 6 is the one exception to note explicitly: it was
taken on outcome-blind grounds, but an outcome panel had by then been run, and
the fact that it was materially unchanged is reported there as a check rather
than withheld. A reader should judge that decision knowing both things.

**Amend this file rather than changing a parameter silently.** A parameter that
moves without a dated entry here is indistinguishable, to a reader, from one
that was tuned against a result.

**On the figures quoted below.** Where a decision rests on a measurement, the
direction and the rough magnitude are given rather than the exact values of any
one cohort: the grounds are properties of the construction, not of a dataset.
Figures for a particular build come from that build's own stage outputs
(`run_info.txt`, `mbi_influence_shares.txt`, `extras/sweep_emphasis.R`).

---

## Rule as pre-specified

Among grid cells satisfying ALL of:

- `tail_bal` in [0.70, 1.40] — low- and high-CMRglc tails pull comparably
- `r2_raw_gbi` <= 0.90 — enough contrast left for the residual to be signal
- `n_features` unchanged against the linear / `emphasis = 1` reference — no region lost

select the cell with the highest `mbi_map_r` = cor(MBI, r_tilt). Ties, within
0.005, break toward the lower `ratio`, then the lower `emphasis`.

N_eff is REPORTED, not targeted.

---

## Decision 1 — 2026-08-13 — three selection criteria struck; parameters fixed directly

**Decided.** The three constraints above and the selection objective are void.
`ratio` and `emphasis` are set directly instead:

```
ratio    = 4    chosen for interpretability
emphasis = 2    a-priori quadratic emphasis on metabolic distance from the median
```

The emphasis grid is retained, and is run to REPORT sensitivity rather than to
select.

**Grounds.** `ratio` is algebraically inert for MBI under the additive map
transform: it scales a term that residualising on GBI removes, and it cancels
in the influence shares. Two of the three constraints therefore cannot
discriminate between cells.

- `r2_raw_gbi` is a property of MBI_raw, which varies with `ratio` while MBI is
  fixed. It cannot be a criterion for choosing MBI.
- The `tail_bal` band is governed by the skew of the FDG map itself, not by
  `ratio` or `emphasis`. No cell in the grid can satisfy it.
- `mbi_map_r` monotonically prefers `emphasis = 1`, because `r_tilt` is a
  full-gradient target and `emphasis` departs from a full gradient by design.
  Maximising it selects against the hypothesis the parameter encodes.

The replacement values were fixed before any further output was inspected.

**Effect on the index.** `emphasis = 2` is the operative choice; `ratio` is
inert, so the value carries labelling information only.

---

## Decision 2 — 2026-08-13 — clamp_z raised 2 → 3; log transform rejected; ratio relabelled

**Decided.** `clamp_z` = 3. `log(CMRglc)` is not used. `ratio` = 1.73.
`emphasis` remains 2 as pre-specified. Sensitivity to report: emphasis ∈ {1, 2, 3}.

**Grounds.**

- *clamp_z.* At `clamp_z = 2` the clamp engaged on 11 low-CMRglc regions and on
  0 high ones, winsorising one tail only — a one-sided distortion of the weight
  vector. At `clamp_z = 3` nothing on this map reaches the bound, and 3, 4 and 6
  give identical results, so the parameter functions as a safety rail rather
  than as a design choice. This ground is independent of any output metric.
- *log transform.* Considered as a symmetrising transform of the map. The map is
  left-skewed (skewness −0.69); logging worsens the skew to −0.93 and increases
  the low/high influence imbalance from 2.24 to 3.12. It does not do the thing
  it would have been adopted for.
- *ratio.* Set to the map's own observed max:min. Inert either way; this value
  makes the weight scale self-documenting.

**Effect on the index.** `clamp_z` changes which regions, if any, are
winsorised; on this map, none are. `ratio` is inert.

---

## Decision 3 — 2026-08-13 — configuration surface reduced. No parameter change

**Decided.** Every configuration option that was set to an inert value, and
every code branch retained only to reproduce an earlier build, is removed:
`pooling`, `map_method`, `alpha`, `alpha_theta`, `gamma`, `winsor_q`,
`per_modality_norm`, `use_size`, `use_global_norm`, `subvol_gate_pct`. Stage 03
warns if a config still carries any of them.

**Grounds.** Each removed option was already at the value the reported run used,
so no number moves. Earlier builds are reproducible from their git tags, which
is what tags are for; dead branches in current code are not a reproducibility
mechanism, and a reviewer reading five options where one is live has to be told
which four are decorative.

**Verified after the change**, on synthetic data with complete cases:

- N_eff identical at ratio 1.73 and ratio 64
- cor(MBI @1.73, MBI @64) = 1.000000000000 exactly
- cor(MBI, GBI) = 8e-17
- GBI bit-identical across emphasis
- emphasis 2 lowers N_eff relative to emphasis 1, as intended
- toy end-to-end smoke test passes

---

## Decision 4 — 2026-08-14 — three upstream and inclusion changes. CHANGES NUMBERS

Three changes outside the index formula. None touches `emphasis`, `ratio`,
`clamp_z`, `beta` or `subvol_min_mm3` — the pre-specified index parameters
stand — but all three require a rerun of the affected stages before any number
is quoted again.

### 4a. Normative model: location-only → location-scale (mgcv gaulss)

**Decided.** Stage 02 fits both the conditional mean and the conditional spread.
Config: `normative.variance_model` (default true), `normative.k_age_scale`
(default 5). Features that fail to converge fall back to the location-only fit
and record that in the `msg` column of `normative_qc.txt` rather than being
dropped.

**Grounds.** A single scalar sigma per feature assumes the spread around the
normative expectation is the same at every age. Regional morphometry becomes
more variable with age, so a single sigma inflates |z| at the ages where a
feature is most variable. On a simulated feature with a linearly age-increasing
SD, 3000 subjects:

|                     | tertile 1 | tertile 2 | tertile 3 |
|---|---|---|---|
| sd(z), single sigma | 0.675 | 0.981 | 1.261 |
| sd(z), gaulss       | 0.994 | 1.009 | 0.999 |
| rate \|z\| > 2, single sigma | 0.005 | 0.046 | 0.125 |
| rate \|z\| > 2, gaulss       | 0.045 | 0.053 | 0.040 |

This matters more for MBI than for a mean index. MBI is a spatial contrast, so
where the variance-by-age slope differs between regions the inflation is
spatially patterned — and spatially patterned inflation of z is arithmetically
indistinguishable from genuine deviation topography.

### 4b. Subject inclusion: required-modality rule added

**Decided.** `mbi.required_modalities`, set to `["CT","SA","SubVol"]`. `null`
restores unrestricted scoring with a warning. Requiring a modality absent from
the data is an error. `compute_flagship_abg()` additionally asserts that the MBI
and GBI passes scored the same subject set before residualising.

**Grounds.** `coverage_thr` decides whether a subject contributes a given
modality; nothing decided whether they contributed enough modalities for the
index to mean the same thing, so a subject with no subcortical data received a
score built from CT and SA alone. Because MBI is a contrast, that is a different
linear functional rather than a noisier estimate of the same one: on simulated
data, r = 0.92 against the full index when SubVol is absent, r = 0.77 when CT is.

### 4c. Feature provenance: the manifest accounts for every feature

**Decided.** Stage 03 starts from the complete feature list and assigns each
feature exactly one fate, and warns if the manifest count disagrees with
`index_concentration()`.

**Grounds.** The previous manifest was region-level and was derived from a `zs`
that had already been filtered, so a feature dropped by the normative fit or by
the SubVol gate appeared in neither the kept nor the dropped column. No file
reconciled the feature map to the index, and the feature count could only be
explained by re-deriving it by hand.

**Also corrected**, with no effect on the reported run: `r2` in
`normative_qc.txt` was NA for every feature, which disabled the low-edf arm of
the degenerate-fit gate, since `is.finite(NA)` is FALSE. It is now computed
explicitly rather than read from `summary()$r.sq`, which is NULL under gaulss in
any case.

**Verified after these changes**, on synthetic data:

- the guardrail drops exactly the incomplete subjects, and errors on an absent modality
- the MBI/GBI subject-set assertion fires correctly
- gaulss recovers a known age-varying sigma at cor = 0.9994
- toy end-to-end smoke test still passes 8/8

---

## Decision 5 — 2026-08-17 — spatial-null testing added. No parameter change

**Decided.** Whether the biological map is load-bearing is tested by permuting
which region receives which map value and rebuilding the contrast, under three
permutation schemes: unconstrained, homotopic (left–right homologous pairs
receive the same value), and variogram-matched surrogates in the manner of Burt
et al. This is an added analysis; no index parameter moves and no number in the
construction changes.

**Grounds.** The index is invariant across a wide range of builds, which is
reassuring but two-edged: every one of those builds is a mean-removed monotone
transform of the same map, so they all point the same way in feature space and
differ only in concentration. Invariance across them says nothing about whether
the map itself is doing work. A null over map *arrangements* does.

Three properties of the null govern how the result may be read:

- Label permutation destroys spatial autocorrelation, so a permutation null is
  ANTI-CONSERVATIVE. A non-significant result is decisive against the map; a
  significant one requires an autocorrelation-preserving null before any claim.
- Preserving left–right symmetry alone is not sufficient. What matters is
  smoothness: a smooth weight vector yields a more reliable index than a rough
  one irrespective of biological content, so both label-permutation schemes are
  beaten partly on smoothness rather than on metabolism.
- Variogram-matched surrogates are rank-reordered onto the original set of map
  values, so every surrogate carries the same value multiset as the real map and
  the weight transform behaves identically; only the spatial arrangement
  changes. Kernel bandwidth and blending weight are fitted by matching Moran's I,
  with variogram shape used to break ties. Report the surrogates' realised
  Moran's I against the real map's: under-matched autocorrelation leaves the null
  too easy to beat, and any p-value from such a run is a lower bound.

A surface spin test is not used as the primary null: it rotates a spherical
projection of the cortical surface, and several of the most influential features
are subcortical or brainstem, so a cortex-only spin would test a different index
from the one reported. `extras/region_centroids.R` supplies whole-brain centroid
distances instead, at the cost of treating each parcel as a point.

**Also decided.** A family-wise maximum statistic across an outcome panel must
be computed on standardised statistics, not raw ones. Null spreads differ
several-fold across outcomes, so a raw maximum is always whichever outcome has
the widest null, and the resulting p-value is uninformative.

**Note on scope.** The outcome-panel code implementing these nulls is
project-specific and is not part of this repository, which covers index
construction and its method-general diagnostics. This entry records the
decision, not the implementation.

---

## Decision 6 — 2026-08-18 — emphasis 2 → 1. CHANGES NUMBERS

**Decided.** `emphasis = 1` is the default and the reported setting. All other
parameters stand: `ratio` 1.73, `clamp_z` 3, `beta` 0, `coverage_thr` 0.80,
`subvol_min_mm3` 1000. Sensitivity to report remains emphasis ∈ {1, 2, 3}.

This supersedes the `emphasis = 2` fixed in Decision 1. Every index value moves,
so every downstream result must be recomputed before it is quoted again.

**Grounds.** Four, in order of weight.

1. *The index becomes parameter-free.* At `emphasis = 1` nothing in the
   construction is tunable: `ratio` is algebraically inert, `clamp_z` does not
   engage on this map, `beta` is 0. There is no quantity left that could have
   been chosen to produce a result. `emphasis = 2` leaves exactly one free
   parameter, and it is the one that concentrates influence.

2. *Concentration.* Influence is markedly more evenly spread. Measured on a
   reference cohort at `beta = 0`, the effective number of contributing features
   more than halves between `emphasis` 1 and 2 and falls by roughly three
   quarters between 1 and 3, while the largest single-feature share doubles and
   then quadruples and the low/high tail balance rises monotonically.
   `extras/sweep_emphasis.R` reports these for any build; the ordering is a
   property of the transform rather than of a cohort.

   The whole point of this construction is that no small set of regions should
   be able to drive the result. `emphasis = 2` more than halves the effective
   number of contributing features relative to the linear form, and doubles the
   share carried by the single largest one.

3. *Tail balance.* At `emphasis = 2` the low-CMRglc tail carries roughly four
   times the influence of the high tail, which makes "metabolic gradient" a
   strained description of what the index measures. At `emphasis = 1` it is
   closer to two, and the description is accurate.

4. *Construct validity.* Correlation with the assumption-free per-subject
   gradient target (`r_tilt`) is markedly higher at `emphasis = 1` than at 2.
   Decision 1 struck this quantity as a *selection* criterion, correctly, because
   it structurally favours the linear form. It remains admissible as evidence
   about which form better measures the stated construct.

**Invariance checks reported, not used to select.** The two settings produce
near-identical indices: subject-level cor(MBI at emphasis 1, MBI at emphasis 2)
is about 0.97. For scale, that exceeds MBI's own split-half reliability, so the
two settings agree with one another more closely than either
agrees with an independent half of the same brain — though both are built from
the same deviations and therefore share their noise, which inflates that
comparison. An outcome panel run under both settings differed negligibly. These
figures are reported to show the change is not a change of result; none of them
was the reason for making it.

**Effect.** `emphasis` is now at the value where the transform is inactive, so
the index is a plain linear contrast along the biological gradient and carries
no fitted quantity of any kind.

---

## Note on subject inclusion

`mbi.required_modalities` is set to `["CT", "SA", "SubVol"]` in the shipped
configuration. It is not optional dressing: a subject scored on a subset of the
modalities carries a different linear contrast, standardised into the same
column as everyone else. A configuration that leaves it unset will run, and
stage 03 will warn, but the resulting column mixes two different quantities.
Check the stage 03 log for that warning before quoting a sample size.

---

## Decision 7 — 2026-08-18 — `ratio` is derived, not set. MBI UNCHANGED

**Decided.** `mbi.ratio` defaults to `auto`: the nominal max:min of the weight
vector is read off the supplied map rather than written into the config. A
literal number is still accepted and is what `extras/verify_invariance.R` uses
to sweep the parameter.

**Grounds.** Decision 2 set `ratio = 1.73` because that was the observed max:min
of the map. It was therefore already a derived quantity, but derived by hand and
then frozen as a literal — which is the worst of both worlds. It reads as a
chosen number to anyone who has not read this file; it has to be recomputed by
hand whenever the map changes; and if nobody recomputes it, it goes stale in
silence. Deriving it removes the last number in the index block that a reader
could mistake for a tuning knob.

**Effect on the index: none.** `ratio` is algebraically inert for MBI, so this
changes only the scale of the discarded intermediate `MBI_raw`, and with it
`R²(MBI_raw, GBI)`. Verified against the implementation: on the reference map
`auto` resolves to 1.7300, weights are bit-identical to passing that number
explicitly, and the normalised contrast direction is unchanged between ratio
1.73, 4 and 64 (correlation 1.000000000000, max elementwise difference 4e-16).
Applied to a constant vector the transform still returns uniform weights, so GBI
is untouched.

**Consequence for reporting.** Quote the value stage 03 resolves and prints,
together with the REALISED contrast, which is lower whenever the clamp does not
engage on both tails.

---

## Decision 8 — 2026-08-31 — N_eff corrected to the index's own frame. NO INDEX VALUE MOVES

**Decided.** `index_concentration()` computes the contrast weights, `N_eff` and
the influence shares in the L1-normalised frame the index is built in, and the
empirical slope is converted onto that frame before it is used. The figures in
Decision 2 are superseded by the values stage 03 reports at pipeline v0.4.1.

**Grounds.** `mbi_raw` and `gbi` are weighted means with L1-normalised weights,
so the contrast is `e = W − b·U` with `W` and `U` on that scale and `b` the slope
between the two unstandardised means. Two things were off. The score columns are
standardised before the slope is fitted, so the fitted `b` was a correlation
rather than a slope on the weight scale; and `shares()` operated on L2-normalised
copies of `W` and `U`, rescaling the second term by `‖W‖/‖U‖`. Because `W` and
`U` are nearly parallel, `e` is a near-cancellation and a sub-percent error in
the effective slope moves the shares by several percent.

The design-time figure was never affected: under the cosine slope the two frames
are proportional, so `n_eff_cosine` and its share vector were already correct.
Only the figure quoted at the realised slope moves.

**Effect.** Diagnostic only. `mbi_raw`, `gbi` and `mbi` are bit-identical before
and after; verified on synthetic data, where the corrected `N_eff` and top share
now agree exactly with the contrast recomputed by hand from the index's own
weights, and sit within 0.15% of the cosine figure rather than several percent
above it. Stage 03 must be re-run once to regenerate `run_info.txt` and
`mbi_influence_shares.txt`; nothing downstream of the scores changes.

**Figures to quote:** those printed by stage 03 at v0.4.1 for the reported
build, and the corresponding sensitivity strip from `extras/sweep_emphasis.R`.
Decision 6's concentration table is to be reissued from that rerun; the ordering
it rests on — `N_eff` falling monotonically with `emphasis`, largest share
rising — is unchanged.

---

## Decision 9 — 2026-08-31 — the discovery/replication split is defined in the config. NO INDEX VALUE MOVES

**Decided.** The split is a deterministic function of the subject identifier and
a label, set in the `split:` block. `extras/make_split.R` writes the assignment
and its realised balance; stage 02 recomputes it rather than reading it.
`split.fit_on` defaults to `all`, so the normative models continue to be fitted
on every subject and no index value moves. A `TRAIN_IDS` path that does not
exist is an error rather than a silent fallback.

**Grounds.** The previous arrangement had two failure modes, both silent.

- *The split could simply not happen.* It was enabled by an environment
  variable, and an unset or mistyped `TRAIN_IDS` fitted the normative models on
  the whole cohort with no error and one line in a log.
- *A frozen split drifts.* Under an `id %in% train_ids` rule, every subject
  added to the cohort after the file was written falls to the held-out side
  without appearing in either list, so the realised discovery fraction declines
  with each cohort refresh and nothing reports it. On a synthetic cohort, a file
  covering a third of the subjects yielded a realised discovery fraction of
  0.332 against a nominal 0.6; that is now reported and warned about.

Regenerating the split each run would replace both with a worse problem:
`sample()` reassigns existing subjects whenever N changes, so a redrawn split
moves people between halves — measured at 23,559 of 95,000 subjects moving under
a cohort growth from 95,000 to 110,000 — and its random stream changed at R 3.6,
so a seed does not pin a split across environments. Hashing the identifier moves
nobody and is reproducible across versions and platforms; both are checked in
`extras/verify_split.R`.

What hashing gives up is stratification. `make_split.R` therefore reports the
realised discovery fraction by sex and age decile rather than assuming it: at
n = 95,000 the realised fraction sat within 0.004 of the requested one and every
decile within 0.003 of even. At small N, or where exact balance is required, an
external assignment supplied through `TRAIN_IDS` remains available.

**On `fit_on`.** The default is to fit the norms on the whole cohort. The
normative model contains no outcome — only age, sex, scanner and optionally ICV
— so fitting it on everyone is not leakage with respect to any later association
test, and it is the more precise estimate of the norms. The discovery/replication
distinction belongs where an outcome enters, which is downstream of this
pipeline; that is what the written assignment is for. `fit_on: discovery` is a
sensitivity arm, and it changes every index value.

**Effect on the index: none** at the shipped default. Verified on synthetic
data: `deviations_z_long.txt` is byte-identical with the `split:` block absent
and with a split defined at `fit_on: all`. Under `fit_on: discovery` the
deviations correlate 0.998 with the full-sample fit and the fitted feature set
differs — a reminder that the stage-02 feature set is an outcome of model
fitting, not a specification.

**On changing the label.** The label identifies the split. Changing it draws a new
one, which is its only purpose; a label changed after outcomes have been inspected
is a split chosen against a result. Any change is recorded here as a dated entry.

---

## Decision 10 — 2026-08-31 — low-r2 arm of the degenerate-fit gate removed; per-subject concentration reported. MAY CHANGE THE FEATURE SET

**Decided.** The degenerate-fit gate keeps its high-edf arm, now expressed as
`normative.edf_max_frac * k_age` per age smooth (0.8 * 10 = 8.0, unchanged at
the defaults). The low-r2 arm — drop a feature with `edf_age < 1.5` and
`r2 < 0.02` — is OFF unless `normative.min_r2` is set. Stage 03 additionally
writes `mbi_subject_concentration.txt`.

**Grounds, low-r2 arm.** A region whose morphometry barely depends on age still
yields a valid z-score; a low r2 is a fact about the covariates, not evidence
that the fit is degenerate. The arm therefore removed features for having a weak
age effect. That matters more here than it would for a mean index: the arm
selects on age-dependence, age-dependence is spatially patterned, and MBI is a
spatial contrast, so what it removed was a spatially patterned subset of the
feature set. It also made the feature set an OUTCOME of model fitting rather
than a specification — free to differ between the full-sample and
discovery-only fits of the same cohort, which is exactly the comparison the
split now invites.

**Grounds, per-subject concentration.** `index_concentration()` describes the
weights under a typical deviation; nothing described a given subject's own
score. Nothing bounds z either: stage 02 fits on winsorised values but computes
z from the raw value, deliberately, and `clamp_z` clamps the map rather than z.
On a synthetic build at N_eff 70, forcing a single feature to |z| = 8 moved that
subject's MBI by 1.6 SD and |z| = 20 moved it by 4.4 SD. The new file reports,
per subject, the largest single-feature share of their own index, the effective
number of features behind it, and the feature responsible, so a score carried by
one parcel is visible rather than inferred. No clamp is applied: the diagnostic
comes first, and any bound would be a further dated decision.

**Effect on the index.** The gate change can alter which features enter, and
therefore every score. On the reference cohort the low-r2 arm never fired —
every feature was fitted, with the largest age-smooth edf well inside the cap —
so removing it
leaves that run unchanged; it is recorded here because it will not necessarily
be dormant under the discovery-only fit. The concentration file is a report and
changes nothing.

---

## Decision 11 — 2026-08-31 — implausible deviations are set to missing, not clipped

**Decided.** `mbi.z_implausible` sets a magnitude above which a deviation is
treated as MISSING before aggregation. It is off unless set. `coverage_thr` and
`required_modalities` then decide whether the affected subject still contributes
that modality and whether they are scored at all. Stage 03 writes the removed
observations to `mbi_implausible_z.txt`, and warns when the rule is unset and the
cohort contains `|z|` above 10.

**Grounds.** Nothing bounded a deviation. Stage 02 computes z from the raw rather
than the winsorised value, deliberately — a genuinely deviant subject should keep
a large deviation — and `clamp_z` bounds the map, not z. On a reference cohort
the `|z|` distribution sits below 3 at the 99th percentile and runs past 50 at
the maximum. No morphometric quantity is fifty standard deviations from its
normative expectation; that is a segmentation failure. And because MBI is a
weighted sum over of order a hundred features with no single one dominating by
design, one absurd value dominates by accident: on a synthetic build of that
size at N_eff 52, a single feature at `|z|` = 56 moved that subject's index by
15 SD.

**Missing rather than clipped.** Clipping asserts a value the data does not
support — that the region really is z_implausible SD deviant — and the index then
propagates that assertion. On the same build, clipping at 5 still left the
artefact moving the subject by 1.7 SD and clipping at 10 by 3.0 SD, where setting
the observation to NA moved them by 0.36 SD, which is the renormalisation from
dropping one feature of the set. Ignorance is the true state of a failed
segmentation, and the existing inclusion rules already handle it.

**Bound used: `mbi.z_implausible = 10`.** Fixed on the sign structure of the
cohort's own tail, before any outcome was consulted.

For CT, SA and SubVol the orientation is +1, so real pathology — atrophy, tissue
loss — produces NEGATIVE deviations. The observed tail points the other way:

| bound | fraction of the tail that is negative |
|---|---|
| 5 | 0.42 |
| 6 | 0.33 |
| 8 | 0.12 |
| 10 | 0.02 |
| 15 | 0.00 |

A tail that is 98% one-sided above 10 and entirely one-sided above 15 is a
mechanism, not a distribution. Treating the real component as symmetric — the
assumption most favourable to a LOW bound, since real deviations are more likely
atrophy-dominant and therefore negative — the real part is about twice the
negative count, which gives the proportion of each cut that is artefact: about a
third at a bound of 6, three quarters at 8, and almost all of it at 10. A bound
of 6 would discard on the order of a thousand negative observations, plausibly
real severe atrophy in as many subjects, which is a directional bias against
exactly the subjects the index exists to detect. A bound of 10 removes well
under 0.01% of all deviations and costs a handful of subjects a modality under
`coverage_thr`.

**What the artefact is.** Not failed reconstructions. Above a bound of 10, no
(subject, region) pair is affected in both CT and SA, the median affected subject
has ONE affected feature, and every one of the eight worst features by rate is
surface area in a small medial parcel — posterior and isthmus cingulate,
parahippocampal, caudal anterior cingulate, paracentral — with the affected
parcels pairing bilaterally. That is a parcel-boundary problem in cortical
surface area near the medial
wall, not a subject-level QC failure. It is recorded rather than acted on beyond
the bound: excluding those parcels by identity is what Decision 5 rules out, and
within-parcel dispersion and the ICC gate remain the outcome-blind instruments if
a region-level exclusion is ever wanted.

There is no universal value and none is shipped in the generic config: read the
cohort's own distribution and its sign structure before setting one.

**Effect on the index.** None unless set. When set, it changes the scores of the
affected subjects only, and the count and identity of every removed observation
is written out rather than absorbed.

---

## Where the parameters now stand

| parameter | status |
|---|---|
| `emphasis` | 1 — transform inactive; the only parameter that changes MBI |
| `ratio` | derived from the map; algebraically inert |
| `clamp_z` | 3 — safety rail; check it does not engage on your map |
| `beta` | 0 — topography rather than burden |
| `coverage_thr`, `required_modalities`, `subvol_min_mm3` | subject and feature inclusion rules, not index parameters |

At these settings the index contains **no fitted or chosen quantity**. Every
number above is either inactive, derived from the data, or an inclusion rule
stated in advance.
