# Calibration helpers for scaling or harmonizing feature distributions across sources.
# Intended for optional pre-processing before normative modelling or MBI computation.

###### load packages, set paths and load functions
pacman::p_load(dplyr,data.table,stringr,tidyr)

path <- "/cluster/projects/p33/users/dennisva/integrity/"
setwd(path)

source(paste0(path,"scripts/functions/compute_nmi.R"))

#=============================================================================
# NMI FLAGSHIP: FDG × tissue-size weighting + sign-aware oriented z (see func)
#=============================================================================
### read in deviations (from normative modeling), feature sizes (for scaling), and processed pet (median) values
zs <- fread(paste0(path,"df/deviations_z_long.txt"),data.table=F)
size_df <- fread(paste0(path,"df/feats_avgSize.txt"),data.table=F)
fdg <- fread(paste0(path,"df/pet_fdg_featMedians.txt"),data.table=F)

keep_sub <- size_df %>% filter(modality=="SubVol", size >= quantile(size, 0.20, na.rm=TRUE)) %>% pull(region)
zs <- zs %>% filter(modality!="SubVol" | region %in% keep_sub)


# pick z column once
z_col <- if ("z_oriented" %in% names(zs)) "z_oriented" else
         if ("z" %in% names(zs))          "z"          else stop("zs needs z_oriented or z")

# PET weights for alignment checks (accept 'w_pet' or 'weight')
fdg_use <- if ("w_pet" %in% names(fdg)) {
  fdg %>% select(region, w_pet) %>% filter(is.finite(w_pet))
} else if ("weight" %in% names(fdg)) {
  fdg %>% transmute(region, w_pet = weight) %>% filter(is.finite(w_pet))
} else stop("fdg must have column 'w_pet' or 'weight'")

# subject-level map×PET (computed once; constant across gamma/beta)
mods_use <- c("CT","SA")
r_tilt_df <- zs %>% filter(modality %in% mods_use) %>%
  select(eid, region, z = all_of(z_col)) %>%
  inner_join(fdg_use, by="region") %>%
  group_by(eid) %>%
  summarise(r_tilt = suppressWarnings(cor(z, w_pet, method="spearman", use="pairwise")),
            .groups="drop")

# PC1 loadings × PET (also constant across gamma/beta)
Z <- zs %>% filter(modality %in% mods_use) %>%
     select(eid, region, z = all_of(z_col)) %>%
     group_by(eid, region) %>% summarise(z = mean(z, na.rm=TRUE), .groups="drop") %>%
     semi_join(fdg_use, by="region")
M <- Z %>% pivot_wider(names_from=region, values_from=z) %>% arrange(eid) %>% select(-eid) %>% as.matrix()
X <- scale(M); X[is.na(X)] <- 0
p <- prcomp(X, center=FALSE, scale.=FALSE)
pc1_tbl <- tibble::tibble(region = colnames(X), loading = as.numeric(p$rotation[,1])) %>%
           inner_join(fdg_use, by="region")
pc1_pet_r_const <- suppressWarnings(cor(pc1_tbl$loading, pc1_tbl$w_pet, use="pairwise"))
sd_rtilt_const  <- sd(r_tilt_df$r_tilt, na.rm=TRUE)

nmi_grid_eval <- function(g, b) {
  # NMI (weighted vs unweighted)
  res_pet <- compute_nmi_flagship(zs, pet_w = fdg,  size_df = size_df, tune = list(gamma = g, beta = b))
  res_unw <- compute_nmi_flagship(zs, pet_w = NULL, size_df = size_df)

  nm <- merge(res_unw$NMI, res_pet$NMI, by="eid", suffixes=c("_unw","_pet"))
  corr_unw_pet <- suppressWarnings(cor(nm$NMI_unw, nm$NMI_pet, use="pairwise"))

  # PET_tilt (standardized residual of NMI_pet ~ NMI_unw)
  fit  <- lm(NMI_pet ~ NMI_unw, data = nm)
  tilt <- as.numeric(scale(residuals(fit)))
  nm$PET_tilt <- tilt

  # How much current tilt aligns with the subject-level map×PET correlation?
  tmp <- merge(nm[, c("eid","PET_tilt")], r_tilt_df, by="eid")
  tilt_map_r <- suppressWarnings(cor(tmp$PET_tilt, tmp$r_tilt, use="pairwise"))

  data.frame(gamma=g, beta=b,
             corr_unw_pet=corr_unw_pet,
             pc1_pet_r=pc1_pet_r_const,   # constant across rows (for reference)
             sd_rtilt=sd_rtilt_const,     # constant across rows (for reference)
             tilt_map_r=tilt_map_r)       # this SHOULD vary with gamma/beta
}


grid  <- expand.grid(gamma=c(2.2,2.3,2.4,2.5), beta=c(0.30,0.25,0.20), KEEP.OUT.ATTRS=FALSE)
evals <- do.call(rbind, mapply(function(ga, be) nmi_grid_eval(ga, be), grid$gamma, grid$beta, SIMPLIFY=FALSE))
tab <- evals[order(-evals$tilt_map_r), ]
write.table(tab,file=paste0(path,"df/nmi_calibration_runs.txt"),sep = "\t",quote=FALSE,row.names=FALSE)



check_concentration <- function(g, b, mods = c("CT","SA","SubVol")) {
  df <- dplyr::inner_join(fdg_use, size_df, by="region") |>
        dplyr::filter(modality %in% mods, is.finite(w_pet), is.finite(size))
  w  <- (df$w_pet^g) * (df$size^b)
  p  <- w / sum(w)
  neff <- 1 / sum(p^2)
  top10 <- sum(sort(p, decreasing=TRUE)[seq_len(ceiling(0.10*length(p)))])
  c(Neff = neff, Top10share = top10)
}

check_concentration(2.4,0.25)


