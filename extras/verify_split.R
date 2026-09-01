#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# verify_split.R - does the discovery/replication assignment do what it claims?
#
# The assignment is load-bearing for any replication claim, and it is the kind
# of component that fails quietly: a hash that does not avalanche still returns
# labels, still looks deterministic, and still writes a plausible file. It just
# does not give the requested fraction and does not spread subjects. That is why
# this runs against the ACTUAL implementation in helpers/split_assign.R rather
# than a description of it.
#
#   Check 1  the realised discovery fraction matches the requested one, for
#            sequential ids, integer ids and random string ids
#   Check 2  u is spread over [0, 1) rather than collapsing onto a lattice
#   Check 3  growing the cohort moves nobody between halves
#   Check 4  the same label reproduces the split; a different label draws an
#            independent one
#   Check 5  the assignment does not depend on the order of the input, and the
#            checksum identifies the SET rather than the sequence
#
# Run from the repo root:   Rscript extras/verify_split.R
# ---------------------------------------------------------------------------

source("helpers/run_helpers.R")
source_local("split_assign.R")

cat(sprintf("[verify_split] pipeline version %s\n\n", pipeline_version()))

ok  <- TRUE
chk <- function(cond, msg) {
  if (!isTRUE(cond)) { cat("  FAIL: ", msg, "\n", sep = ""); ok <<- FALSE }
  else               { cat("  ok  : ", msg, "\n", sep = "") }
}

set.seed(4)
N     <- 40000
FRAC  <- 0.6
LABEL <- "verify"

styles <- list(
  "sequential string ids" = sprintf("S%06d", seq_len(N)),
  "integer-like ids"      = as.character(seq.int(1000000, length.out = N)),
  "random string ids"     = replicate(N, paste(sample(c(letters, 0:9), 12, TRUE),
                                               collapse = ""))
)

# ---- Checks 1 and 2: fraction and spread ---------------------------------
cat("[1,2] realised fraction and spread\n")
for (nm in names(styles)) {
  a   <- assign_split(styles[[nm]], label = LABEL, discovery_frac = FRAC)
  fr  <- mean(a$split == "discovery")
  dec <- as.numeric(table(cut(a$u, seq(0, 1, 0.1))) / nrow(a))
  cat(sprintf("      %-22s frac %.5f | distinct u %d/%d | max decile dev %.5f\n",
              nm, fr, length(unique(a$u)), nrow(a), max(abs(dec - 0.1))))
  chk(abs(fr - FRAC) < 0.01,
      sprintf("%s: realised fraction within 0.01 of requested", nm))
  chk(length(unique(a$u)) > 0.99 * nrow(a),
      sprintf("%s: u does not collapse onto a lattice", nm))
  chk(max(abs(dec - 0.1)) < 0.01,
      sprintf("%s: u is spread evenly across deciles", nm))
}

# ---- Check 3: stability under cohort growth ------------------------------
ids   <- styles[[1]]
a     <- assign_split(ids, label = LABEL, discovery_frac = FRAC)
grown <- c(ids, sprintf("S%06d", seq.int(N + 1L, N + N %/% 4L)))
b     <- assign_split(grown, label = LABEL, discovery_frac = FRAC)
m     <- merge(a, b, by = "id")
moved <- sum(m$split.x != m$split.y)
cat(sprintf("\n[3] cohort %d -> %d: %d of %d existing subjects changed side\n",
            N, length(grown), moved, nrow(m)))
chk(moved == 0L, "growing the cohort moves nobody between halves")
chk(abs(mean(b$split == "discovery") - FRAC) < 0.01,
    sprintf("the fraction still holds after growth (%.5f)",
            mean(b$split == "discovery")))

# ---- Check 4: label semantics ---------------------------------------------
cat("\n[4] label\n")
chk(identical(assign_split(ids, label = LABEL, discovery_frac = FRAC)$split, a$split),
    "the same label reproduces the split exactly")
other <- assign_split(ids, label = paste0(LABEL, "-2"), discovery_frac = FRAC)
agree <- mean(a$split == other$split)
chance <- FRAC^2 + (1 - FRAC)^2
chk(abs(agree - chance) < 0.02,
    sprintf("a different label draws an independent split (agreement %.4f, chance %.4f)",
            agree, chance))

# ---- Check 5: order independence -----------------------------------------
cat("\n[5] order independence\n")
rev_a <- assign_split(rev(ids), label = LABEL, discovery_frac = FRAC)
chk(identical(rev_a$split[match(ids, rev_a$id)], a$split),
    "the assignment does not depend on the order of the input")
disc <- a$id[a$split == "discovery"]
chk(identical(split_checksum(disc), split_checksum(rev(disc))),
    "the checksum identifies the set, not the sequence")
chk(!identical(split_checksum(disc),
               split_checksum(b$id[b$split == "discovery"])),
    "the checksum distinguishes two different discovery sets")

# ---- Config guards --------------------------------------------------------
cat("\n[6] configuration is validated rather than approximated\n")
bad <- function(x) inherits(try(read_split_cfg(x), silent = TRUE), "try-error")
chk(bad(list(split = list(mode = "hash", label = ""))),
    "mode hash with an empty label is refused")
chk(bad(list(split = list(mode = "none", fit_on = "discovery"))),
    "fit_on discovery with no split defined is refused")
chk(bad(list(split = list(mode = "stratified"))), "an unknown mode is refused")
chk(bad(list(split = list(mode = "hash", label = "s", discovery_frac = 1.4))),
    "a fraction outside (0, 1) is refused")
chk(bad(list(split = list(mode = "hash", label = "s", seed = 42))),
    "an unknown key under split: is refused")
d <- read_split_cfg(list())
chk(identical(d$mode, "none") && identical(d$fit_on, "all"),
    "the default is no split, norms fitted on everyone")

cat("\n")
if (ok) {
  cat("[verify_split] ALL CHECKS PASSED\n")
} else {
  cat("[verify_split] CHECKS FAILED\n")
  quit(status = 1L)
}
