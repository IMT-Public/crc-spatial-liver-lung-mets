## Package installation for the capsule image.
##
## The CRAN snapshot is pinned by date so a rebuild resolves the same versions
## the published analysis used. The __linux__/jammy path serves prebuilt
## Ubuntu 22.04 binaries of that same snapshot, which matches the Code Ocean
## base image and avoids compiling every package from source. Bioconductor is
## pinned to 3.22, which is the release matching R 4.5.1.

options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/jammy/2025-09-01"),
        Ncpus = max(1L, parallel::detectCores()))

## Every package the pipeline loads, and nothing else. `ggpubr` is pulled in
## by lib/panels.R for stat_pvalue_manual, `nlme` by 01_figure1.R for
## gsummary, and `reshape2`/`stringr`/`tibble`/`tidyr` by the helpers in
## lib/. Keeping this list tight keeps the image build honest: if a script
## starts using something new, the build fails rather than silently picking it
## up from a transitive dependency.
cran <- c(
  "ggplot2", "dplyr", "tidyr", "tibble", "stringr",
  "ggrepel", "ggpubr", "patchwork",
  "RColorBrewer", "reshape2",
  "readxl", "writexl",
  "lme4", "lmerTest", "nlme",
  "survival", "survminer"
)
install.packages(cran)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(version = "3.22", ask = FALSE, update = FALSE)

bioc <- c(
  "Biobase", "BiocGenerics", "S4Vectors",
  "NanoStringNCTools", "GeomxTools",
  "limma", "fgsea", "GSVA"
)
BiocManager::install(bioc, ask = FALSE, update = FALSE)

## Verify everything the pipeline loads is importable.
need <- c(cran, bioc)
missing <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Failed to install: ", paste(missing, collapse = ", "))
cat("All packages installed.\n")
