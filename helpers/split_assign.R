#!/usr/bin/env Rscript
#
# Deterministic discovery/replication assignment.
#
# WHAT THIS IS FOR
#
# A discovery/replication split has to satisfy two things that are easy to state
# and easy to get wrong:
#
#   1. It is fixed before any outcome is inspected.
#   2. It does not move. A subject assigned to replication stays in replication,
#      including after the cohort grows.
#
# A stratified draw satisfies neither reliably. sample() reassigns every subject
# whenever N changes, because adding one person to one stratum shifts the whole
# draw; and its random stream changed at R 3.6, so a seed does not pin a split
# across environments. Storing the result in a file fixes the movement but
# introduces staleness: subjects added after the file was written appear in
# neither half, and under a `id %in% train_ids` rule they are silently all
# assigned to the held-out side, so the realised discovery fraction drifts
# downward with every cohort refresh and nothing reports it.
#
# Making the assignment a property of the SUBJECT removes both failure modes.
# u(id) below depends on the subject identifier and a label, and on nothing else:
# not on who else is in the sample, not on how many, not on the order, not on
# the R version. Adding subjects leaves every existing assignment untouched and
# splits the new ones in the same proportion. The split is then fully specified
# by two config values, so there is no artefact to maintain and nothing to go
# stale.
#
# The cost is that it does not stratify. At cohort scale the strata balance
# themselves and split_balance() below reports by how much; at small N, or where
# exact balance is required, supply an external assignment instead.
#
# CHANGING THE LABEL DRAWS A NEW SPLIT. That is its only purpose, and it is not a
# tuning knob: a label changed after outcomes have been inspected is a split
# chosen against a result. Record any change in extras/SELECTION_RULE.md.

if (!exists("%||%")) `%||%` <- function(a, b) if (!is.null(a)) a else b

# ---- Exact 32-bit arithmetic on doubles ----------------------------------
# R has no unsigned 32-bit integer type and bitwXor() rejects anything at or
# above 2^31, so both operations are done on doubles in 16-bit halves. Every
# intermediate stays below 2^49, well inside the 2^53 range where double
# arithmetic is exact, which is what makes the result bit-identical on every
# platform and every R version with no package dependency.
.U32 <- 4294967296                     # 2^32

.mulmod32 <- function(a, b) {
  bh <- b %/% 65536; bl <- b %% 65536
  (((a * bh) %% .U32) * 65536 + a * bl) %% .U32
}

.xor32 <- function(a, b) {
  bitwXor(a %/% 65536, b %/% 65536) * 65536 + bitwXor(a %% 65536, b %% 65536)
}

# utf8ToInt() is not vectorised over a character vector; this is. Characters
# past the end of a shorter string read as 0, which is why the label is a PREFIX:
# it makes every hashed string differ from position 1 rather than only in a
# trailing region.
.chr_int <- function(x) {
  u <- unique(x)
  v <- vapply(u, function(s) if (nzchar(s)) as.numeric(utf8ToInt(s)) else 0,
              numeric(1), USE.NAMES = FALSE)
  v[match(x, u)]
}

# ---- Subject -> uniform [0, 1) -------------------------------------------
# FNV-1a over the bytes of "<label>|<id>", followed by a murmur3-style
# finaliser.
#
# The finaliser is not decoration. A plain polynomial or FNV hash does not
# avalanche in its low bits, so sequential identifiers — which is what cohort
# ids usually are — come out as a lattice rather than a spread: consecutive
# subjects land adjacent to one another, u takes few distinct values, and the
# realised discovery fraction misses the requested one badly. Measured on 95,000
# sequential ids, a small-multiplier polynomial hash gave a realised fraction of
# 0.684 against a requested 0.600, with u taking about 19 distinct values. The
# shift-multiply-shift finaliser is what turns the counter structure into a
# uniform draw; extras/verify_split.R checks that it did.
split_u <- function(id, label = "") {
  ids <- paste0(as.character(label), "|", as.character(id))
  if (!length(ids)) return(numeric(0))

  maxn <- max(nchar(ids))
  h    <- rep(2166136261, length(ids))              # FNV-1a offset basis
  for (k in seq_len(maxn)) {
    h <- .mulmod32(.xor32(h, .chr_int(substr(ids, k, k))), 16777619)
  }

  h <- .xor32(h, h %/% 65536)
  h <- .mulmod32(h, 2246822507)
  h <- .xor32(h, h %/% 8192)
  h <- .mulmod32(h, 3266489909)
  h <- .xor32(h, h %/% 65536)

  h / .U32
}

# ---- Assignment ----------------------------------------------------------
# Returns a data frame of id, split and u. `u` travels with the assignment so a
# reader can verify the labels rather than trust them.
assign_split <- function(id, label = "", discovery_frac = 0.6) {
  if (!is.finite(discovery_frac) || discovery_frac <= 0 || discovery_frac >= 1)
    stop("assign_split(): split.discovery_frac must be strictly between 0 and 1; got ",
         discovery_frac, call. = FALSE)
  id <- as.character(id)
  if (anyDuplicated(id))
    stop("assign_split(): duplicated subject ids; the assignment would be ",
         "ambiguous for: ",
         paste(utils::head(unique(id[duplicated(id)]), 5), collapse = ", "),
         call. = FALSE)
  u <- split_u(id, label)
  data.frame(id    = id,
             split = ifelse(u < discovery_frac, "discovery", "replication"),
             u     = u,
             stringsAsFactors = FALSE)
}

# ---- Balance report ------------------------------------------------------
# What stratification would have guaranteed, reported instead of assumed. Pass
# the covariates you would have stratified on; each is crossed with the split
# and the realised discovery fraction is reported per level.
split_balance <- function(assign_df, covars = NULL) {
  overall <- data.frame(
    variable = "(overall)", level = "(all)",
    n = nrow(assign_df),
    n_discovery = sum(assign_df$split == "discovery"),
    frac_discovery = mean(assign_df$split == "discovery"),
    stringsAsFactors = FALSE
  )
  if (is.null(covars) || !length(covars)) return(overall)

  rows <- lapply(names(covars), function(v) {
    lv <- as.character(covars[[v]])
    ok <- !is.na(lv)
    do.call(rbind, lapply(sort(unique(lv[ok])), function(l) {
      sel <- ok & lv == l
      data.frame(variable = v, level = l,
                 n = sum(sel),
                 n_discovery = sum(sel & assign_df$split == "discovery"),
                 frac_discovery = mean(assign_df$split[sel] == "discovery"),
                 stringsAsFactors = FALSE)
    }))
  })
  rbind(overall, do.call(rbind, rows))
}

# ---- Set checksum --------------------------------------------------------
# One short string identifying WHICH subjects are in the discovery half,
# independent of their order. Recorded in run_info.txt by every stage that
# touches the split, so "are these two runs using the same split?" is a string
# comparison rather than a set operation over two large files. A changed
# checksum with an unchanged label means the cohort changed, not the split rule.
split_checksum <- function(id) {
  if (!length(id)) return("empty")
  h   <- sort(floor(split_u(id, "checksum") * .U32))
  acc <- 2166136261
  for (k in seq_along(h)) acc <- .mulmod32(.xor32(acc, h[k]), 16777619)
  sprintf("%08x-%d", as.integer(acc %% 2147483647), length(id))
}

# ---- Config reader -------------------------------------------------------
# One place that interprets the `split:` block, so stage 02 and the split script
# cannot disagree about what the config asked for. Unknown keys and unknown
# values are errors: a split is not something to be silently approximated.
read_split_cfg <- function(cfg) {
  sp <- cfg$split
  if (is.null(sp)) sp <- list()

  known <- c("mode", "label", "discovery_frac", "fit_on")
  extra <- setdiff(names(sp), known)
  if (length(extra))
    stop("config: unknown key(s) under `split:`: ", paste(extra, collapse = ", "),
         ". Accepted: ", paste(known, collapse = ", "), ".", call. = FALSE)

  mode <- as.character(sp$mode %||% "none")
  if (!mode %in% c("none", "hash"))
    stop("config: split.mode must be \"none\" or \"hash\"; got \"", mode, "\".",
         call. = FALSE)

  fit_on <- as.character(sp$fit_on %||% "all")
  if (!fit_on %in% c("all", "discovery"))
    stop("config: split.fit_on must be \"all\" or \"discovery\"; got \"",
         fit_on, "\".", call. = FALSE)

  if (mode == "none" && fit_on == "discovery")
    stop("config: split.fit_on is \"discovery\" but split.mode is \"none\", ",
         "so there is no split to fit on. Set split.mode: hash, or supply an ",
         "external assignment through the TRAIN_IDS environment variable.",
         call. = FALSE)

  label <- as.character(sp$label %||% "")
  if (mode == "hash" && !nzchar(label))
    stop("config: split.mode is \"hash\" but split.label is empty. The label is ",
         "what identifies this split; set it and record it in ",
         "extras/SELECTION_RULE.md.", call. = FALSE)

  frac <- suppressWarnings(as.numeric(sp$discovery_frac %||% 0.6))
  if (mode == "hash" && (!is.finite(frac) || frac <= 0 || frac >= 1))
    stop("config: split.discovery_frac must be strictly between 0 and 1; got ",
         sp$discovery_frac, call. = FALSE)

  list(mode = mode, label = label, discovery_frac = frac, fit_on = fit_on)
}
