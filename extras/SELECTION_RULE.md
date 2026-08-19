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

2. *Concentration.* Influence is markedly more evenly spread. Measured on the
   reference cohort at `beta = 0`, over 143 features:

   | emphasis | N_eff | largest feature share | low/high tail balance | cor with the map-gradient target |
   |---|---|---|---|---|
   | **1** | **43.8** | **5.3%** | **2.16** | **0.926** |
   | 1.5 | 28.3 | 8.4% | 2.90 | 0.885 |
   | 2 | 20.6 | 11.8% | 3.97 | 0.843 |
   | 3 | 12.4 | 19.4% | 7.45 | 0.776 |

   The whole point of this construction is that no small set of regions should
   be able to drive the result. `emphasis = 2` more than halves the effective
   number of contributing features relative to the linear form, and doubles the
   share carried by the single largest one.

3. *Tail balance.* At `emphasis = 2` the low-CMRglc tail carries roughly four
   times the influence of the high tail, which makes "metabolic gradient" a
   strained description of what the index measures. At `emphasis = 1` the ratio
   is 2.16, and the description is accurate.

4. *Construct validity.* Correlation with the assumption-free per-subject
   gradient target (`r_tilt`) is 0.926 at `emphasis = 1` against 0.843 at 2.
   Decision 1 struck this quantity as a *selection* criterion, correctly, because
   it structurally favours the linear form. It remains admissible as evidence
   about which form better measures the stated construct.

**Invariance checks reported, not used to select.** The two settings produce
near-identical indices: subject-level cor(MBI at emphasis 1, MBI at emphasis 2)
= 0.9688 over 94,818 subjects. For scale, MBI's own split-half reliability is
0.805, so the two settings agree with one another more closely than either
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
