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
  #
  # A variable that is not set is left in place VERBATIM rather than expanded to
  # an empty string or treated as fatal here. One config serves every stage, so
  # failing on any unresolved variable would stop stage 03 because a FreeSurfer
  # path it never reads is unset. Instead the unexpanded value travels through,
  # and assert_file_exists() raises the error at the point a stage actually asks
  # for that path — naming the variable, and only for paths that stage needs.
  unresolved <- character(0)

  expand_env <- function(x) {
  if (!is.character(x) || length(x) != 1L || is.na(x)) return(x)

  m <- gregexpr("\\$\\{[A-Za-z_][A-Za-z0-9_]*\\}", x, perl = TRUE)[[1]]
  if (length(m) == 1L && m[1] == -1L) return(x)

  hits <- regmatches(x, list(m))[[1]]
  if (length(hits) == 0) return(x)

  for (h in unique(hits)) {
    var <- sub("^\\$\\{", "", sub("\\}$", "", h))
    val <- Sys.getenv(var, unset = NA_character_)
    if (is.na(val) || !nzchar(val)) {
      unresolved <<- c(unresolved, var)
      next
    }
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

      # Still carrying ${VAR}? Then it is not a path yet. Return it untouched:
      # prefixing the repo root onto "${FEATURES_FILE}" would bury the variable
      # name inside a plausible-looking absolute path and make the eventual
      # error much harder to read.
      if (grepl("\\$\\{", y)) return(y)

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

  if (length(unresolved)) {
    u <- sort(unique(unresolved))
    message("config ", basename(path), ": ", length(u),
            " environment variable(s) not set and left unexpanded: ",
            paste0("$", u, collapse = ", "),
            ". This is only a problem if this stage needs one of them.")
  }
  cfg
}

# Version marker, so "which copy of the code am I actually running" is a
# one-line question. Recorded by write_run_info() into every stage directory.
pipeline_version <- function() {
  f <- file.path(repo_root(), "VERSION")
  if (file.exists(f)) trimws(readLines(f, warn = FALSE)[1]) else "unknown"
}

# ---- Provenance -----------------------------------------------------------
# Every stage writes run_info.txt beside its outputs: the code version, the
# config, the time, and the exact package versions in play. Without this a
# directory of results cannot be traced back to the code that produced it,
# which is the whole point of archiving a pipeline.
write_run_info <- function(stage_path, stage, config_path = NA_character_,
                           params = NULL) {
  si   <- utils::sessionInfo()
  pkgs <- character(0)
  if (length(si$otherPkgs))
    pkgs <- sort(vapply(si$otherPkgs,
                        function(p) paste0(p$Package, " ", p$Version),
                        character(1), USE.NAMES = FALSE))

  lines <- c(
    sprintf("stage             %s", stage),
    sprintf("pipeline_version  %s", pipeline_version()),
    sprintf("config            %s", as.character(config_path)),
    sprintf("timestamp         %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    sprintf("R                 %s", R.version.string),
    sprintf("platform          %s", R.version$platform)
  )
  if (length(params))
    lines <- c(lines, "", "parameters:",
               paste0("  ", names(params), " = ",
                      vapply(params, function(v)
                        if (is.null(v)) "NULL" else paste(v, collapse = ", "),
                        character(1))))
  if (length(pkgs))
    lines <- c(lines, "", "attached packages:", paste0("  ", pkgs))

  f <- file.path(stage_path, "run_info.txt")
  writeLines(lines, f)
  cat("Wrote: ", f, "\n", sep = "")
  invisible(f)
}

# ---- Table writer ---------------------------------------------------------
# data.table::fwrite where available, falling back to write.table. The
# arguments are pinned so the two produce byte-identical output: unquoted,
# no row names, literal "NA" for missing (fwrite defaults to an empty field,
# which fread would then read back as an empty string in a character column).
write_tsv <- function(x, file) {
  x <- as.data.frame(x)
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fwrite(x, file, sep = "\t", quote = FALSE, na = "NA",
                       row.names = FALSE)
  } else {
    utils::write.table(x, file, sep = "\t", quote = FALSE, row.names = FALSE)
  }
  invisible(file)
}

assert_file_exists <- function(path, label = "file") {
  if (is.null(path) || !nzchar(path)) stop("Missing path for ", label, call. = FALSE)

  # This is where an unset environment variable becomes an error: at the moment
  # a stage asks for a path it genuinely needs, naming the variable rather than
  # reporting a mangled path as "not found".
  if (grepl("\\$\\{", path)) {
    v <- regmatches(path, gregexpr("\\$\\{[A-Za-z_][A-Za-z0-9_]*\\}", path))[[1]]
    v <- unique(sub("^\\$\\{", "", sub("\\}$", "", v)))
    stop("Path for ", label, " still contains unset environment variable(s): ",
         paste0("$", v, collapse = ", "), "\n",
         "  value: ", path, "\n",
         "  Note that `VAR=value` alone is NOT visible to Rscript - it must be exported:\n",
         paste0("    export ", v, "=/path/to/...", collapse = "\n"),
         call. = FALSE)
  }

  if (!file.exists(path)) stop("Not found (", label, "): ", path, call. = FALSE)
  invisible(TRUE)
}

stage_dir <- function(outdir, stage) {
  d <- file.path(outdir, stage)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

start_logging <- function(logfile) {
  # Returns a closer rather than registering its own on.exit(): an on.exit()
  # here would fire when start_logging() returns and close both sinks
  # immediately, making the call a no-op. Caller does:
  #   stop_log <- start_logging(f); on.exit(stop_log())
  con <- file(logfile, open = "wt")
  sink(con, type = "output")
  sink(con, type = "message")
  invisible(function() {
    sink(type = "message"); sink(type = "output")
    close(con)
  })
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
