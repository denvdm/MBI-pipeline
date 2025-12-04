# Toy example for the MBI-pipeline
# --------------------------------
# This script simulates a small dataset of regional brain deviations and PET weights,
# then computes:
#   - GBI     : unweighted Global Brain Index (mean deviation across regions)
#   - MBI_raw: PET-weighted aggregate of deviations
#   - MBI    : residual of MBI_raw after regressing out GBI (alignment-style metric)
#
# The goal is to illustrate the MBI framework and the new naming convention:
#   MBI      = primary alignment metric (formerly MDA)
#   MBI_raw  = raw PET-weighted MBI aggregate (formerly "MBI")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
})

set.seed(42)

# -------------------------------------------------------------------------
# 1. Simulate toy deviations and PET weights
# -------------------------------------------------------------------------

n_subjects <- 200
n_regions  <- 30

subjects <- paste0("subj_", seq_len(n_subjects))
regions  <- paste0("region_", seq_len(n_regions))

# Simulate PET weights: positive, right-skewed, then normalized
toy_pet_weights <- tibble(
  region    = regions,
  pet_raw   = rexp(n_regions, rate = 1),          # skewed positive
  pet_weight = pet_raw / sum(pet_raw)             # normalize to sum 1
)

# Simulate subject-level "global vulnerability" (GBI driver)
subject_effect <- tibble(
  subject = subjects,
  g_global = rnorm(n_subjects, mean = 0, sd = 0.8)
)

# Simulate deviations: region-wise + subject-wise + noise
# Some regions are "hotter" (more strongly coupled to global)
region_sensitivity <- tibble(
  region = regions,
  sens   = runif(n_regions, min = 0.5, max = 1.5)
)

toy_deviations_long <- expand_grid(subject = subjects, region = regions) %>%
  left_join(subject_effect, by = "subject") %>%
  left_join(region_sensitivity, by = "region") %>%
  mutate(
    noise    = rnorm(n(), mean = 0, sd = 1),
    dev_z    = g_global * sens + noise            # this is our synthetic deviation z-score
  ) %>%
  select(subject, region, dev_z)

# -------------------------------------------------------------------------
# 2. Compute GBI, MBI_raw, and MBI
# -------------------------------------------------------------------------

# GBI: unweighted mean deviation per subject
gbi_df <- toy_deviations_long %>%
  group_by(subject) %>%
  summarize(
    GBI = mean(dev_z, na.rm = TRUE),
    .groups = "drop"
  )

# Join deviations with PET weights
dev_with_pet <- toy_deviations_long %>%
  left_join(toy_pet_weights %>% select(region, pet_weight),
            by = "region")

# MBI_raw: PET-weighted sum of deviations per subject
mbi_raw_df <- dev_with_pet %>%
  group_by(subject) %>%
  summarize(
    MBI_raw = sum(dev_z * pet_weight, na.rm = TRUE),
    .groups = "drop"
  )

# Combine GBI + MBI_raw
toy_idx <- gbi_df %>%
  left_join(mbi_raw_df, by = "subject")

# MBI: residual of MBI_raw after regressing out GBI
mbi_model <- lm(MBI_raw ~ GBI, data = toy_idx)
toy_idx$MBI <- resid(mbi_model)

# Optionally z-score indices for nicer scales
toy_idx <- toy_idx %>%
  mutate(
    GBI_z     = as.numeric(scale(GBI)),
    MBI_raw_z = as.numeric(scale(MBI_raw)),
    MBI_z     = as.numeric(scale(MBI))
  )

# -------------------------------------------------------------------------
# 3. Save outputs to example/data
# -------------------------------------------------------------------------

out_dir <- file.path("example", "data")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

readr::write_tsv(
  toy_idx %>% select(subject, GBI, MBI_raw, MBI, GBI_z, MBI_raw_z, MBI_z),
  file.path(out_dir, "toy_mbi_indices.tsv")
)

readr::write_tsv(
  toy_deviations_long,
  file.path(out_dir, "toy_deviations_long.tsv")
)

readr::write_tsv(
  toy_pet_weights,
  file.path(out_dir, "toy_pet_weights.tsv")
)

message("Toy example data written to example/data/")

# -------------------------------------------------------------------------
# 4. Simple visualization
# -------------------------------------------------------------------------

p1 <- ggplot(toy_idx, aes(x = GBI_z, y = MBI_raw_z)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "Relationship between GBI and MBI_raw (PET-weighted index)",
    x = "GBI (z-scored)",
    y = "MBI_raw (z-scored)"
  ) +
  theme_minimal()

p2 <- ggplot(toy_idx, aes(x = GBI_z, y = MBI_z)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "Relationship between GBI and MBI (residualized alignment metric)",
    x = "GBI (z-scored)",
    y = "MBI (z-scored; residual of MBI_raw ~ GBI)"
  ) +
  theme_minimal()

# Save plots
plots_dir <- file.path("example", "plots")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

ggsave(file.path(plots_dir, "toy_GBI_vs_MBI_raw.png"), p1, width = 6, height = 4, dpi = 300)
ggsave(file.path(plots_dir, "toy_GBI_vs_MBI.png"),     p2, width = 6, height = 4, dpi = 300)

message("Toy example plots written to example/plots/")

# -------------------------------------------------------------------------
# 5. Print a small preview
# -------------------------------------------------------------------------

print(head(toy_idx, 10))
message("Done. You can now inspect example/data/ and example/plots/.")
