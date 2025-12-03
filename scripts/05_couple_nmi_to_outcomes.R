# Fit linear/logistic models relating MBI metrics to clinical or cognitive outcomes.
# Wrapper around regression helpers in R/ to generate tables and summaries.
# NOTE: This section uses UKB-specific directories; generalization planned.

###### load packages, set paths and load functions
library(dplyr,data.table, purrr,pls,brglm2,pROC)
path <- "/path/to/data"

source("R/regression_funcs.R")

###### read in x,y, covs
nmi <- fread(paste0(path,"integrity/df/nmi_flagship_scores.txt"),data.table=F)

## Continuous traits
mdict <- fread(paste0(ukbdir,"metabo/df/metaboDict_240109.txt"),header=T,data.table = F,fill=T) 
metabo <- fread(paste0(ukbdir,"metabo/df/nmr_qced_long_240622.csv"),header=T,data.table = F,fill=T)
markers <- colnames(metabo)[colnames(metabo) %in% mdict$nmrName]
metabo <- metabo[metabo$visit_index==0,c("eid",markers)]
bmi <- fread(paste0(ukbdir,"metabo/df/bmi.csv"),header=T,data.table = F,fill=T)
bmi$bmi <- bmi$`21001-2.0`

## Covs
demog <- fread(paste0(ukbdir,"ukb/covariates/UKB54k_demogShort_220926.txt"),data.table=F)
colnames(demog) <- c("eid","age","sex","scanner")

dt <- reduce(list(nmi,metabo,bmi,demog),full_join,by = "eid")

## Dx
for(dx in c("SCZ","MDD","BIP","CAD","T2D")){
tmp <- fread(paste0(ukbdir,"metdx/df/dx/UKB500k_",dx,"_240906.txt"),header=T,data.table = F,fill=T)
dt$dx <- ifelse(dt$eid %in% tmp$eid,1,0)
colnames(dt)[colnames(dt)=="dx"] <- tolower(dx)
}

covars <- c("age","sex","scanner")
y_cont <- c("bmi",markers)
y_diag <- c("mdd","bip","scz","cad","t2d") 

#summary(lm(mdd~nmi_fdg+age+sex+scanner,data=dt))
dt <- make_tilt(dt, save_prov_path = "integrity/df/tilt_provenance.rds")

# Continuous outcomes
lin_tab <- as.data.frame(run_linear_panel(dt, y_vars=y_cont, covars=covars)); head(lin_tab)
write.table(lin_tab,file=paste0(path,"integrity/df/nmi_linreg_results.txt"),sep = "\t",quote=FALSE,row.names=FALSE)

# Diagnoses
log_tab <- as.data.frame(run_logistic_panel(dt, y_vars=y_diag, covars=covars, auc=TRUE)); head(log_tab)
write.table(log_tab,file=paste0(path,"integrity/df/nmi_logreg_results.txt"),sep = "\t",quote=FALSE,row.names=FALSE)

###### Confirm predictors impact low-FDG regions harder than high-FDG (given opposite signs)
zs <- fread(paste0(path,"integrity/df/deviations_z_long.txt"),data.table=F)
size_df <- fread(paste0(path,"integrity/df/feats_avgSize.txt"),data.table=F)
fdg <- fread(paste0(path,"integrity/df/pet_fdg_featMedians.txt"),data.table=F)

##### Regional maps for arbitrary outcomes 
y_vars <- c("bmi")  # or c("bmi","GlycA","Glucose","HbA1c", ...)
covars <- c("age","sex","scanner")

# Run the generic mapper directly on dt (the function selects needed cols)
reg_maps <- region_betas_for_outcomes(
  zs = zs,
  dt = dt,
  y_vars = y_vars,
  covars = covars,
  mods_use = c("CT","SA","SubVol"),
  size_df = size_df,
  subvol_gate_pct = 0.20,
  scale_y = TRUE  # effect is per 1-SD increase in outcome
)

# Write one combined table (long: outcome × region) that's ggseg-ready
outfile_maps <- if (length(y_vars) == 1L)paste0("roi_betas_by_", y_vars, ".txt") else "roi_betas_by_outcome.txt"
write.table(reg_maps, file = paste0(path,"integrity/df/",outfile_maps),sep = "\t", row.names = FALSE, quote = FALSE)

# PET-alignment summary per outcome (optional, useful in captions)
fdg_use <- if ("w_pet" %in% names(fdg)) dplyr::select(fdg, region, w_pet) else dplyr::transmute(fdg, region, w_pet = weight)

pet_align <- reg_maps |>
  dplyr::inner_join(fdg_use, by = "region") |>
  dplyr::group_by(outcome) |>
  dplyr::summarise(
    r_Pearson  = suppressWarnings(cor(beta, w_pet, use = "pairwise")),
    r_Spearman = suppressWarnings(cor(beta, w_pet, method = "spearman", use = "pairwise")),
    n_regions  = dplyr::n(), .groups = "drop"
  )
write.table(pet_align, file = paste0(path,"integrity/df/roi_betas_by_outcome__pet_alignment.txt"),sep = "\t", row.names = FALSE, quote = FALSE)

