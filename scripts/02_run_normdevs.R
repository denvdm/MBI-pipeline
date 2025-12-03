# Fit normative models for all structural MRI features and compute subject-level deviations.
# Uses GAM-based modelling with covariates, site effects, and standardized residual outputs.
# NOTE: This section uses UKB-specific directories; generalization planned.

###### load packages and set paths
library(dplyr,data.table,mgcv,future.apply,matrixStats)

path <- "/path/to/data"

# NOTE: The following section uses UK Biobank directory structure.
# Adjust to match your dataset layout.

## load functions ()
source("R/normative_modeling.R")

# Modeling options
k_age          <- 10         # basis dimension for s(age); mgcv will penalize
use_site_RE    <- TRUE       # include site as random effect
sex_interact   <- FALSE      # if TRUE adds s(age, by=sex) interaction smooth
family_gaussian <- gaussian()# keep gaussian; deviations are standardized residuals
n_cores <- max(1L, parallel::detectCores(logical = TRUE) - 1L)

###### Get brain features and covariates
### demog
demog <- fread(paste0(path,"UKB/covariates/UKB54k_demogShort_220926.txt"),data.table=F)
colnames(demog) <- c("eid","age","sex","scanner")
demog$eid <- as.character(demog$eid)

### cortical data (left, right, area and thickness)
lh_thick <- fread(paste0(ukbpath,"stats/batch_",batch,"/regularFSstats/lh.thickness.UKBB.txt"), header=T,data.table=F)
rh_thick <- fread(paste0(ukbpath,"stats/batch_",batch,"/regularFSstats/rh.thickness.UKBB.txt"), header=T,data.table=F)
colnames(lh_thick)[1] <- "MRID" -> colnames(rh_thick)[1]
thick <- full_join(lh_thick,rh_thick,by="MRID")
#thickFeats <- gsub("lh_","",colnames(lh_thick)[-1])

lh_area <- fread(paste0(ukbpath,"stats/batch_",batch,"/regularFSstats/lh.area.UKBB.txt"), header=T,data.table=F)
rh_area <- fread(paste0(ukbpath,"stats/batch_",batch,"/regularFSstats/rh.area.UKBB.txt"), header=T,data.table=F)
colnames(lh_area)[1] <- "MRID" -> colnames(rh_area)[1]
area <- full_join(lh_area,rh_area,by="MRID")
#areaFeats <- gsub("lh_","",colnames(lh_area)[-1])

cortex <- full_join(thick,area,by="MRID")

### aseg, do some feature selection
aseg <- fread(paste0(ukbpath,"stats/batch_",batch,"/regularFSstats/subcorticalstats.UKBB.txt"), header=T,data.table=F)
colnames(aseg)[1] <- "MRID"
#asegFeats <- substr(colnames(aseg)[grep("Left-",colnames(aseg))],6,30)
asegNonhemi <- c("EstimatedTotalIntraCranialVol","3rd-Ventricle","4th-Ventricle","Brain-Stem","CC_Posterior","CC_Mid_Posterior","CC_Central",
  "CC_Mid_Anterior","CC_Anterior")
asegHemi <- c("Lateral-Ventricle","Inf-Lat-Vent","Cerebellum-Cortex","Thalamus-Proper","Caudate",
  "Putamen","Pallidum","Hippocampus","Amygdala","Accumbens-area","VentralDC")
aseg <- aseg[,c("MRID",asegNonhemi,paste0("Left-",asegHemi),paste0("Right-",asegHemi))]
colnames(aseg)[-c(1,2)] <- paste0(colnames(aseg)[-c(1,2)],"_aseg")

### combine all features and make names syntactic
all <- full_join(aseg,cortex,by="MRID")
all$MRID <- gsub("^FS_","",all$MRID)
old_names <- names(all)
names(all) <- gsub("-","_",names(all))
names(all) <- make.names(names(all), unique = TRUE)

## Identify feature columns (also a place to exclude some, e.g. globals)
feature_cols <- setdiff(names(all), c("MRID","EstimatedTotalIntraCranialVol",
  "lh_MeanThickness_thickness","rh_MeanThickness_thickness","lh_WhiteSurfArea_area","rh_WhiteSurfArea_area"))

# write feature 'raw' values to file for input to calcFeatSize.R 
write.table(all[,feature_cols],file=paste0(path,"integrity/df/feats_rawVals.txt"),sep = "\t",quote=FALSE,row.names=FALSE)
dt <- inner_join(demog,all,by=c("eid"="MRID"))

# Tidy up common covariates, incl. residualizing ICV
dt <- dt[complete.cases(dt[, c("age", "sex", "scanner", "EstimatedTotalIntraCranialVol")]), ]
dt$sex <- as.factor(dt$sex)
dt$scanner <- as.factor(dt$scanner)

icv_fit <- mgcv::gam(EstimatedTotalIntraCranialVol ~ s(age, k = 5, bs = "cs") + sex + s(scanner, bs = "re"),data = dt, method = "REML")
dt$icv <- as.numeric(scale(resid(icv_fit, type = "response"), center = TRUE, scale = FALSE))

# Modality-specific covariates, by default we include ICV for all features. Could change this to _area and _aseg only. 
icv_feats <- grep("(_thickness|_area|_aseg)$", feature_cols,value=T)
needs_icv <- function(feat) feat %in% icv_feats

## Winsorizing
dt_raw <- dt
dt[feature_cols] <- lapply(dt[feature_cols], winsorize, lower = 0.001, upper = 0.999)

## Feature include/exclude patterns
# include_regex <- NULL  # e.g., "(_aseg|_area|_thickness)$"
# exclude_regex <- NULL  # e.g., "^DWI_FA_" to skip certain metrics
# if (!is.null(include_regex))  feature_cols <- feature_cols[grepl(include_regex, feature_cols)]
# if (!is.null(exclude_regex))  feature_cols <- setdiff(feature_cols, grep(exclude_regex, feature_cols, value = TRUE))

## Run normative modelling
cat(sprintf("Fitting %d features with %d workers...\n", length(feature_cols), n_cores))
res_list <- future_lapply(feature_cols, fit_one_feature, future.seed = TRUE)

## collect outputs
ok_res <- Filter(function(x) is.list(x) && isTRUE(x$ok), res_list)
fail_res <- Filter(function(x) is.list(x) && !isTRUE(x$ok), res_list)

# Deviations matrix (subjects × features)
Z <- data.table(eid = dt$eid)
MU <- data.table(eid = dt$eid)

for (r in ok_res) {
  Z[[r$feat]]  <- r$z
  MU[[r$feat]] <- r$mu
}

## Write outputs
write.table(Z,file=paste0(path,"integrity/df/deviations_z.txt"),sep = "\t",quote=FALSE,row.names=FALSE)
write.table(MU,file=paste0(path,"integrity/df/predicted_means.txt"),sep = "\t",quote=FALSE,row.names=FALSE)

## make z long form for NMI calculation (run_nmi.R)
#Z <- fread(paste0(path,"integrity/df/deviations_z.txt"),data.table=F)
Z <- as.data.frame(Z)

# Pick feature columns by suffix
fc   <- grep("(_thickness|_area|_aseg)$", names(Z), value = TRUE, ignore.case = TRUE)
# Flatten values + build keys
vals <- as.numeric(as.matrix(Z[fc]))
feat <- rep(fc, each = nrow(Z))
sid  <- rep(Z$eid, times = length(fc))
# Modality & region
mod  <- ifelse(grepl("_thickness$", feat, ignore.case=TRUE), "CT",
        ifelse(grepl("_area$",     feat, ignore.case=TRUE), "SA",
        ifelse(grepl("_aseg$",     feat, ignore.case=TRUE), "SubVol", NA)))
reg  <- sub("_thickness$|_area$|_aseg$", "", feat)

zs_long <- data.frame(eid = as.character(sid),modality = mod,region = reg,z = vals,stringsAsFactors = FALSE)
zs_long <- zs_long[!is.na(zs_long$modality), ]

write.table(zs_long,file=paste0(path,"integrity/df/deviations_z_long.txt"),sep = "\t",quote=FALSE,row.names=FALSE)
