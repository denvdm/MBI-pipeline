# Quantify similarity or alignment between MBI patterns and PET maps.
# Uses correlation, spatial matching, and region-weight comparisons.

###### load packages, set paths and load functions
library(dplyr,data.table,broom,tidyr)

path <- "/path/to/data"

source("R/compute_nmi.R")
source("R/relate_pet_funcs.R")

#=============================================================================
# calculate relationship (unweighted) integrity with metabolism
#=============================================================================
### read in z (from normative modeling), feat sizes (for scaling), and pet (median) values
zs <- fread(paste0(path,"df/deviations_z_long.txt"),data.table=F)
size_df <- fread(paste0(path,"df/feats_avgSize.txt"),data.table=F)
fdg <- fread(paste0(path,"df/pet_fdg_featMedians.txt"),data.table=F) 

mods_use <- c("CT","SA","SubVol")      # flagship comparison
subvol_gate_pct <- 0.20                # p20 gate (aligns with flagship)

fdg_use <- prep_fdg(fdg)
zs_gated <- gate_subvol(zs, size_df, pct = subvol_gate_pct)
Zc <- collapse_zs(zs_gated, fdg_use, mods_use = mods_use)

mods_tag <- paste(mods_use, collapse = "+")

# --- Run the three approaches and write to file
row_pc  <- pc_pet_alignment(Zc, fdg_use, mods_tag)
row_loo <- loo_pet_alignment(Zc, fdg_use)
align_tab <- bind_rows(row_pc, row_loo)
write.table(align_tab,file = paste0(path,"df/pet_alignment_summary.txt"),sep="\t",row.names=FALSE,quote=FALSE)

row_subj <- subject_pet_alignment(Zc, fdg_use)
write.table(row_subj,file = paste0(path,"df/subject_tilt_summary.txt"),sep="\t",row.names=FALSE,quote=FALSE)


### Grab input for brain maps

# 1) PC1 loadings (sign-aligned to PET)
pc1_tbl <- pc1_loadings_table(Zc, fdg_use, align_to_pet = TRUE)
write.table(pc1_tbl, paste0(path,"df/pc1_regional_z.txt"),sep = "\t", row.names = FALSE, quote = FALSE)

# 2) Mean z per region
meanZ <- mean_z_per_region(Zc)
write.table(meanZ, paste0(path,"df/mean_regional_z.txt"),sep = "\t", row.names = FALSE, quote = FALSE)

# 3) LOO beta map (optional)
looB <- loo_beta_map(Zc)
write.table(looB, paste0(path,"df/roi_beta_LOO.tsv"),sep = "\t", row.names = FALSE, quote = FALSE)

