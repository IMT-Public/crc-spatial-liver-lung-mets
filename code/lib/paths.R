## ---------------------------------------------------------------------------
## lib/paths.R -- portable path resolution
##
## Works unchanged in three places:
##   * Code Ocean capsule   -> /data (read-only), /results (writable)
##   * local git clone      -> <repo>/data, <repo>/results
##   * custom locations     -> set CRC_DATA_DIR / CRC_RESULTS_DIR env vars
##
## All analysis scripts live in code/ and are run with code/ as the working
## directory, so the local fallback is "..".
## ---------------------------------------------------------------------------

DATA_DIR <- Sys.getenv(
  "CRC_DATA_DIR",
  unset = if (dir.exists("/data")) "/data" else file.path("..", "data")
)

RESULTS_DIR <- Sys.getenv(
  "CRC_RESULTS_DIR",
  unset = if (dir.exists("/results")) "/results" else file.path("..", "results")
)

data_path    <- function(...) file.path(DATA_DIR, ...)
results_path <- function(...) file.path(RESULTS_DIR, ...)

if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

## Results are grouped one directory per manuscript figure, so everything
## behind a given figure -- panels, score tables, p-values -- sits together:
##
##   results/figure1 .. figure4      panels and tables for main Figures 1-4
##   results/supplementary_figure1   Supplementary Figure 1
##
## fig_path("figure2", "Fig2a_mCAF.pdf") creates the directory on first use.
fig_path <- function(figure, ...) {
  d <- results_path(figure)
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  file.path(d, ...)
}

## Fail loudly and early if an expected input is missing, rather than deep
## inside an analysis chunk.
require_input <- function(...) {
  p <- data_path(...)
  if (!file.exists(p)) {
    stop("Required input file not found: ", p,
         "\nSee data/README.md for the expected data manifest.", call. = FALSE)
  }
  p
}

message("DATA_DIR    = ", normalizePath(DATA_DIR, mustWork = FALSE))
message("RESULTS_DIR = ", normalizePath(RESULTS_DIR, mustWork = FALSE))
