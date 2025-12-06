# MBI-pipeline

Tools and workflows for computing **Metabolic Brain Index (MBI)** metrics by combining MRI-based normative deviations with PET-derived metabolic priors.

The pipeline supports:
- Estimation of **GBI** (Global Brain Index): unweighted whole-brain deviation score
- Computation of **MBI_raw**: PET-weighted deviation index
- Computation of **MBI**: alignment metric derived as the residual of *MBI_raw ~ GBI*
- Scripts for PET preprocessing, PET extraction, normative modelling, aggregation, and outcome analyses
- A fully self-contained toy example

---

## Conceptual overview

### **1. GBI — Global Brain Index**
Unweighted whole-brain deviation score (mean deviation across regions).

### **2. MBI_raw — PET-weighted deviation index**
Weighted deviation score using region-specific PET metabolic weights.

### **3. MBI — alignment metric (formerly MDA)**
Residual of the model:

    MBI_raw ~ GBI

This isolates PET-weighted deviation **not explained by global deviation**, making MBI the primary interpretable alignment measure.

---

## Conceptual diagram (Mermaid)

```mermaid
flowchart LR
    A[Regional<br/>Deviations] --> B[GBI<br/>(unweighted mean)]
    A --> C[PET Weights]
    C --> D[MBI_raw<br/>(weighted sum)]
    B --> E[Residualization<br/>MBI_raw ~ GBI]
    D --> E
    E --> F[MBI<br/>(alignment metric)]
```

---

## Repository layout

- **R/** – core modules (normative modelling, PET weighting, MBI/GBI/MBI_raw computation, regression helpers)
- **scripts/** – pipeline runners (00_* PET prep; 01_* PET extraction; 02_* normative modelling etc.)
- **inst/templates/** – PET templates, lookup tables
- **inst/tools/** – auxiliary scripts (e.g., FSaverage → PET transform)
- **example/** – toy workflow: synthetic deviations, PET weights, MBI indices, plots

---

## Quick run example

From the repository root, run the toy example with:

```
Rscript example/run_toy_example.R
```

Outputs are written to:

- `example/data/`
- `example/plots/`

---

## Getting started

Clone the repository:

```
git clone https://github.com/denvdm/MBI-pipeline.git
cd MBI-pipeline
```

More documentation will follow in future releases (v0.1.0 planned).
