#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# region_centroids.R  -  region centroids in millimetres.
#
# Emits a region x (x, y, z) table for any analysis that needs inter-region
# distances: spatial-autocorrelation-preserving nulls, variogram matching,
# distance-weighted smoothing. It is derived from the same atlas stage 01 reads,
# so the geometry is exactly the geometry the map values were sampled on and
# region names come out in the same canonical form, ready to join.
#
# Centroids cover subcortical and brainstem regions as well as cortex, which a
# surface-based spin test cannot: a spin rotates a spherical projection of the
# cortical surface only. The cost is treating each parcel as a point.
#
# Run from the repo root:
#   Rscript extras/region_centroids.R --config config.yaml \
#     --out region_centroids.tsv
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(RNifti)
  library(optparse)
})
source("helpers/run_helpers.R")

opt <- parse_common_args(list(
  make_option("--out", type = "character", default = "region_centroids.tsv",
              help = "output TSV [default %default]")
))
cat(sprintf("[region_centroids] pipeline version %s\n", pipeline_version()))
cfg <- read_config(opt$config)

maps_dir   <- cfg$paths$maps_dir %||% NULL
atlas_file <- cfg$inputs$atlas_file %||% file.path(maps_dir, "aparc+aseg_inPETspace_3mm.nii.gz")
lut_file   <- cfg$inputs$lut_file   %||% file.path(maps_dir, "lut_aparc_aseg.csv")
assert_file_exists(atlas_file, "atlas_file")
assert_file_exists(lut_file,   "lut_file")

atl <- RNifti::readNifti(atlas_file)
xf  <- RNifti::xform(atl)
lut <- utils::read.csv(lut_file)

lab_col <- intersect(c("index","label","Label","ID","id"), names(lut))[1]
nam_col <- intersect(c("ROI_Name","name","Name","region","LabelName"), names(lut))[1]
if (is.na(lab_col) || is.na(nam_col))
  stop("could not find label and name columns in the LUT; it has: ",
       paste(names(lut), collapse = ", "), call. = FALSE)

# Same canonical form as scripts/01_extract_petvals.R, so names join downstream.
canon <- function(r) {
  r <- gsub("^ctx-",   "",    r)
  r <- gsub("-",       "_",   r)
  r <- gsub("^Left_",  "lh_", r)
  r <- gsub("^Right_", "rh_", r)
  make.names(r)
}

labs <- sort(unique(as.vector(atl)))
labs <- labs[labs > 0]
message(sprintf("atlas: %s\n  %d non-zero labels, dim %s",
                basename(atlas_file), length(labs),
                paste(dim(atl), collapse = " x ")))

vox <- which(atl > 0, arr.ind = TRUE)
val <- atl[atl > 0]

rows <- lapply(labs, function(L) {
  ix <- which(val == L)
  if (!length(ix)) return(NULL)
  # mean voxel index (1-based) -> 0-based -> mm via the affine
  ijk <- colMeans(vox[ix, , drop = FALSE]) - 1
  mm  <- xf %*% c(ijk, 1)
  nm  <- lut[[nam_col]][match(L, lut[[lab_col]])]
  if (is.na(nm)) return(NULL)
  data.frame(label = L, region = canon(as.character(nm)),
             x = mm[1], y = mm[2], z = mm[3], n_vox = length(ix),
             stringsAsFactors = FALSE)
})
cent <- do.call(rbind, rows)

if (any(duplicated(cent$region)))
  stop("duplicated region names after canonicalisation: ",
       paste(unique(cent$region[duplicated(cent$region)]), collapse = ", "),
       call. = FALSE)

dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
write_tsv(cent, opt$out)
message("Wrote: ", opt$out, "  (", nrow(cent), " regions)")

# Sanity checks a reader can act on.
lr <- ifelse(grepl("^lh_", cent$region), "L", ifelse(grepl("^rh_", cent$region), "R", "M"))
cat(sprintf("\nhemisphere counts: L %d, R %d, unpaired %d\n",
            sum(lr == "L"), sum(lr == "R"), sum(lr == "M")))
cat(sprintf("mean x by hemisphere: L %+.1f, R %+.1f mm  (should have opposite signs)\n",
            mean(cent$x[lr == "L"]), mean(cent$x[lr == "R"])))
D <- as.matrix(dist(cent[, c("x","y","z")]))
cat(sprintf("centroid distances: median %.1f mm, max %.1f mm\n",
            median(D[upper.tri(D)]), max(D)))
cat("\nSmallest parcels (check these are the ones you expect to be small):\n")
print(utils::head(cent[order(cent$n_vox), c("region","n_vox")], 6), row.names = FALSE)
