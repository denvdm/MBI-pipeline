#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# verify_invariance.R - numerical check that MBI is invariant to the shared
# variance of the two indices it is constructed from.
#
# R^2(MBI_raw, GBI) is a property of a discarded intermediate: it moves freely
# with `ratio` while MBI does not move at all. This script demonstrates that
# against the ACTUAL implementation in helpers/compute_mbi.R — not a
# re-implementation — on synthetic deviations, by sweeping ratio over three
# orders of magnitude and checking three propositions:
#
#   Proposition 1   MBI_raw = (GBI + kappa*T) / (1 + kappa*m)
#   Proposition 2   MBI     = standardised residual of T on GBI, for all kappa
#   Proposition 3   1 - R^2 = kappa^2 sT^2 (1-RT^2) / (sG^2 + 2 kappa C + kappa^2 sT^2)
#
# where T_i = sum_f U_f utilde_f y_if is the tilt score and contains no kappa.
#
# Run from the repo root:   Rscript extras/verify_invariance.R
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(dplyr))
source("helpers/compute_mbi.R")
set.seed(31)

n_ct <- 66; n_sub <- 14; n <- 8000
reg_c <- sprintf("ctx%02d", 1:n_ct); reg_s <- sprintf("sub%02d", 1:n_sub)
regions <- c(reg_c, reg_s)
xm  <- rbeta(length(regions), 5, 2); xm <- (xm - min(xm))/(max(xm) - min(xm))
map <- 1 + xm * 0.73                      # left-skewed, max:min = 1.73

featdf <- rbind(data.frame(region = reg_c, modality = "CT"),
                data.frame(region = reg_c, modality = "SA"),
                data.frame(region = reg_s, modality = "SubVol"))
Ff <- nrow(featdf)
Zm <- matrix(rnorm(n * Ff), n, Ff)
zs <- data.frame(id = rep(as.character(1:n), times = Ff),
                 region   = rep(featdf$region,   each = n),
                 modality = rep(featdf$modality, each = n),
                 z = as.vector(Zm), stringsAsFactors = FALSE)
sizes <- featdf %>% mutate(size = runif(n(), 1500, 9000))
fdg   <- data.frame(region = regions, weight = map, stringsAsFactors = FALSE)

beta <- 0                                  # size weighting off
# Deliberately run at emphasis = 2, which is NOT the default. The invariance to
# `ratio` is a property of the algebra and holds at every emphasis; testing it
# away from the default shows it is not an accident of the shipped setting.
a <- pmax(1e-12, sizes$size)^beta; A <- sum(a); U <- a / A

# utilde, exactly as transform_map_weights() constructs it
med <- median(map); s <- mad(map, center = med, constant = 1.4826)
uu  <- pmax(-1, pmin(1, (map - med)/(3 * s)))
ut  <- (sign(uu) * abs(uu)^2)[match(featdf$region, regions)]

m   <- sum(ut * a)/A
GBI <- as.vector(Zm %*% U)/sum(U)
Tt  <- as.vector(Zm %*% (U * ut))/sum(U)   # the tilt score; no kappa in it
resT <- as.numeric(scale(resid(lm(Tt ~ GBI))))

cat(sprintf("\n%7s %9s %16s %24s %15s %15s\n", "rho", "kappa",
            "max|Prop1 error|", "cor(MBI, std resid(T|GBI))",
            "1-R2 predicted", "1-R2 observed"))
for (rho in c(1.05, 1.2, 1.73, 4, 16, 64)) {
  k  <- (rho - 1)/(rho + 1)
  o  <- suppressMessages(compute_flagship_abg(zs, fdg, sizes, ratio = rho,
                          clamp_z = 3, emphasis = 2, beta = beta))
  s2 <- as.data.frame(o$scores); s2 <- s2[order(as.integer(s2$id)), ]

  d1 <- max(abs(s2$mbi_raw - as.numeric(scale((GBI + k*Tt)/(1 + k*m)))))   # Prop 1
  cc <- abs(cor(s2$mbi, resT))                                            # Prop 2
  sg <- var(GBI); st <- var(Tt); C <- cov(Tt, GBI)                        # Prop 3
  pred <- k^2 * st * (1 - C^2/(sg*st)) / (sg + 2*k*C + k^2*st)

  cat(sprintf("%7.2f %9.5f %16.1e %24.12f %15.3e %15.3e\n",
              rho, k, d1, cc, pred, 1 - o$diagnostics$r2_raw_gbi))
}
cat("\nExpected: column 3 at machine precision (~1e-15); column 4 exactly 1 at\n")
cat("every rho; columns 5 and 6 identical. The shared variance moves by three\n")
cat("orders of magnitude while the index does not change at all.\n")
