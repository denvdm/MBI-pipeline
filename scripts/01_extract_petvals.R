# Extract regional PET values from templates or maps using a specified parcellation.
# Applies winsorization, coverage checks, and generates PET-derived weighting vectors.

###### load packages, set paths and load functions
library(RNifti)
path <- "/path/to/data"

### note: run 00_prep_petmap.sh first
source("R/pet_parcel_extract.R")

# inputs already in PET space (3mm); atlas NN-resampled, GM prob optional
pet_file   <- "source-castrillon2023_desc-cmrglc_space-MNI152_res-3mm_feature.nii.gz"
atlas_file <- "aparc+aseg_inPETspace_3mm.nii.gz"  # nearest-neighbour to PET grid
gm_file    <- "GMprob_inPETspace_3mm.nii.gz"          # liberal GM (e.g., p>=0.20)
lut <- read.csv("lut_aparc_aseg.csv")

res <- save_pet_parcel_csv(
  pet=pet_file, atlas=atlas_file, gm_prob=gm_file, output_csv="pet_fdg_parcellated.csv",
  gm_thr=0.20, use_median=TRUE, clamp_neg=TRUE, winsor_q=NULL,  # set winsor_q=0.005 if needed
  min_vox=50L, min_cov=0.30, lut=lut
)

### do some formatting to make suitable as input for NMI code
fdg_raw <- res$parcels
# align region naming
colnames(fdg_raw) <- c("index","fdg_median","fdg_mean","p05","p95","nvox_used","nvox_total","coverage","region")
fdg_raw$region <- gsub("^ctx-","",fdg_raw$region)
fdg_raw$region <- gsub("-","_",fdg_raw$region)
fdg_raw$region <- make.names(fdg_raw$region, unique = TRUE)

# z-scale FDG values, move to min 0, normalize to 1 and set column names
fdg_raw$fdg_z <- scale(fdg_raw$fdg_median)
fdg_raw$fdg_w <- fdg_raw$fdg_z - min(fdg_raw$fdg_z, na.rm = TRUE) + 1e-6
fdg_raw$fdg_w <- fdg_raw$fdg_w / mean(fdg_raw$fdg_w, na.rm = TRUE)
fdg <- fdg_raw %>% dplyr::transmute(region, weight = fdg_w)

write.table(fdg,file=paste0(path,"pet_fdg_featMedians.txt"),sep = "\t",qu