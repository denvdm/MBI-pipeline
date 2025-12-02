# MBI-pipeline

Pipeline and tools for computing Metabolic Brain Index (MBI) metrics by combining MRI-based normative deviations with PET-derived metabolic priors.

This repository will host:

- Core R functions for normative modelling, PET-weighted aggregation, and computation of MBI/GBI/MDA metrics  
- Runner scripts for full end-to-end pipelines  
- Example workflows with toy data  
- Utilities for relating MBI metrics to PET patterns and clinical or cognitive outcomes  

## Status

This is a work-in-progress reference implementation.  
Current priorities:

1. Cleaning and modularizing the existing internal scripts  
2. Providing a minimal end-to-end example with toy data  
3. Documenting interfaces for reuse across cohorts  

Expect breaking changes until the first tagged release (v0.1.0).

## Repository layout (planned)

- **R/** – core R modules (normative modelling, PET extraction, MBI/GBI/MDA computation, regression helpers)  
- **scripts/** – pipeline runner scripts (e.g., 00_* for PET prep, 01_* for extraction, 02_* for normative modelling, etc.)  
- **inst/templates/** – PET templates, lookup tables, region lists  
- **inst/tools/** – auxiliary scripts (e.g., FSaverage → PET transforms)  
- **example/** – toy dataset and end-to-end MBI workflow  
- **tests/** – minimal checks for core functionality  
- **.github/workflows/** – CI configuration (planned)  

## Getting started

Public API and detailed examples are under active development.

1. Clone the repository:

```bash
git clone https://github.com/denvdm/MBI-pipeline.git
cd MBI-pipeline
```

2. Inspect the R/ and scripts/ directories for the current implementation.

3. Once the example/ directory is populated, you will be able to run a small, self-contained MBI workflow on toy data.

More documentation will follow once the first stable interface is in place.