# Compute MBI, GBI, and MDA metrics from normative deviations and PET weighting matrices.
# Wrapper around the compute_mbi() functions in R/.

###### load packages, set paths and load functions
library(dplyr,data.table,stringr,tidyr)

path <- "/path/to/data"
source("R/compute_mbi.R")

#=============================================================================
# NMI FLAGSHIP: FDG × tissue-size weighting + sign-aware oriented z (see func)
#=============================================================================
### read in deviations (from normative modeling), feature sizes (for scaling), and processed pet (median) values
zs <- fread(paste0(path,"df/deviations_z_long.txt"),data.table=F)
size_df <- fread(paste0(path,"df/feats_avgSize.txt"),data.table=F)
fdg <- fread(paste0(path,"df/pet_fdg_featMedians.txt"),data.table=F)

## remove small subvol regions to lower noise
keep_sub <- size_df %>% filter(modality=="SubVol", size >= quantile(size, 0.20, na.rm=TRUE)) %>% pull(region)
zs <- zs %>% filter(modality!="SubVol" | region %in% keep_sub)

## Compute GBI, MBI and tilt (MDA)
res <- compute_flagship_abg(
  zs      = zs,
  fdg     = fdg,
  sizes   = size_df,
  beta    = 0.25,         # your β
  gamma   = 2.5,          # your γ
  alpha   = NULL,         # NULL = equal α across available modalities, alt e.g.: c(CT = 0.40, SA = 0.40, SubVol = 0.20)
  cap_frac     = 0.05,
  coverage_thr = 0.80
)
## extract metrics and write output
# nmi <- left_join(res_fdg$NMI[c("eid","NMI")],res_unw$NMI[c("eid","NMI")],by="eid")
# colnames(nmi) <- c("eid","nmi","gbi")
mbi <- as.data.frame(res$scores)
write.table(mbi,file=paste0(path,"df/mbi_flagship_scores.txt"),sep = "\t",quote=FALSE,row.names=FALSE)

## some checks
cor(mbi$gbi,mbi$mbi)
