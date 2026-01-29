suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

# Assumption: run from repo root
repo_root <- function() normalizePath(getwd(), winslash = "/", mustWork = TRUE)

source_local <- function(filename) {
  root <- repo_root()
  candidates <- c(
    file.path(root, "R", filename),
    file.path(root, "funcs", filename),
    file.path(root, "helpers", filename)
  )
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit) || length(hit) == 0) {
    stop(
      "Cannot source ", filename, ". Looked in:\n- ",
      paste(candidates, collapse = "\n- "),
      call. = FALSE
    )
  }
  source(hit, local = FALSE)
  invisible(hit)
}

read_config <- function(path) {
  if (!file.exists(path)) stop("Config not found: ", path, call. = FALSE)
  cfg <- yaml::read_yaml(path)
  if (is.null(cfg$inputs)) stop("Config must contain top-level 'inputs:'", call. = FALSE)

  # Expand ${ENV_VAR} and ~ in any character fields that look like paths.
  # Keep this intentionally lightweight (no templating engine).
  expand_env <- function(x) {
  if (!is.character(x) || length(x) != 1L || is.na(x)) return(x)

  m <- gregexpr("\\$\\{[A-Za-z_][A-Za-z0-9_]*\\}", x, perl = TRUE)[[1]]
  if (length(m) == 1L && m[1] == -1L) return(x)

  hits <- regmatches(x, list(m))[[1]]
  if (length(hits) == 0) return(x)

  for (h in unique(hits)) {
    var <- sub("^\\$\\{", "", sub("\\}$", "", h))
    val <- Sys.getenv(var, unset = h)
    x <- gsub(h, val, x, fixed = TRUE)
  }
  x
}

  looks_like_path <- function(x) {
    is.character(x) && length(x) == 1L && !is.na(x) &&
      (grepl("/", x) || grepl("^~", x) || grepl("^\\.", x) || grepl("\\$\\{", x))
  }

  expand_paths_rec <- function(obj) {
    # Be robust to unexpected non-atomic objects that may appear in YAML-derived lists.
    if (is.function(obj)) return(obj)

    if (is.list(obj)) {
      return(lapply(obj, expand_paths_rec))
    }

    if (looks_like_path(obj)) {
      y <- expand_env(obj)
      y <- path.expand(y)

      # Resolve relative paths against repo root
      if (!grepl("^/", y) && !grepl("^[A-Za-z]:[/\\\\]", y)) {
        y <- file.path(repo_root(), y)
      }

      return(normalizePath(y, winslash = "/", mustWork = FALSE))
    }

    obj
  }

  cfg <- expand_paths_rec(cfg)
  cfg
}

assert_file_exists <- function(path, label = "file") {
  if (is.null(path) || !nzchar(path)) stop("Missing path for ", label, call. = FALSE)
  if (!file.exists(path)) stop("Not found (", label, "): ", path, call. = FALSE)
  invisible(TRUE)
}

stage_dir <- function(outdir, stage) {
  d <- file.path(outdir, stage)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

start_logging <- function(logfile) {
  con <- file(logfile, open = "wt")
  sink(con, type = "output")
  sink(con, type = "message")
  on.exit({
    sink(type = "message"); sink(type = "output")
    close(con)
  }, add = TRUE)
  invisible(TRUE)
}

parse_common_args <- function(extra_options = list(),
                              default_config = "config.yaml",
                              default_outdir = "out") {

  base <- list(
    make_option(c("-c", "--config"), type = "character", default = default_config,
                help = paste0("YAML config file [default: ", default_config, "]"), metavar = "FILE"),
    make_option(c("-o", "--outdir"), type = "character", default = default_outdir,
                help = paste0("Output directory [default: ", default_outdir, "]"), metavar = "DIR"),
    make_option(c("--prefix"), type = "character", default = "mbi",
                help = "Output prefix [default: %default]"),
    make_option(c("--overwrite"), action = "store_true", default = FALSE,
                help = "Overwrite outputs")
  )

  parser <- OptionParser(option_list = c(base, extra_options))
  opt <- parse_args(parser)

  # Resolve relative paths against repo root (so copy/paste is easy)
  if (!grepl("^/", opt$config)) opt$config <- file.path(repo_root(), opt$config)
  if (!grepl("^/", opt$outdir)) opt$outdir <- file.path(repo_root(), opt$outdir)

  opt
}
