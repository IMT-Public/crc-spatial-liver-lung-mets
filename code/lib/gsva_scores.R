## ---------------------------------------------------------------------------
## lib/gsva_scores.R -- per-ROI signature scores, averaged to one row per
##                      patient x site x area
##
## Most scores are read from GSVAscores_updated_03272025.csv. Two are scored
## here:
##   CD4.Tn.CCR7   50-gene CD4+CCR7+ Tcm list
##   Plasma.cells  9-gene plasma cell list
##
## gsva_scores()  -> per-ROI table
## average_rois() -> per patient/site/area means; every model uses this, so
##                   each patient contributes one value per group rather than
##                   three technical ROIs.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(GSVA)
  library(dplyr)
  library(stringr)
})

PLASMA_GENES <- c("MS4A1", "CD79A", "MZB1", "XBP1", "CD38",
                  "CD27", "IGHG1", "JCHAIN", "IGHA1")

gsva_scores <- function(d) {
  gs <- read.csv(require_input("data_boxplots", "GSVAscores_updated_03272025.csv"))
  load(require_input("data_boxplots", "CD4TnCCR7_top50genes.RDATA"))  # -> top.genes

  ## Scores are assigned positionally, so ROI order must agree.
  stopifnot(identical(sub("\\.dcc$", "", gs$ROI), colnames(d)))

  gs$CD4.Tn.CCR7 <- gsva(gsvaParam(d, list(CD4.Tn.CCR7 = top.genes)),
                         verbose = FALSE)[1, ]

  ## MS4A1 and CD38 are not on the GeoMx WTA panel; GSVA scores the
  ## remaining seven genes.
  gs$Plasma.cells <- gsva(gsvaParam(d, list(plasma = PLASMA_GENES)),
                          verbose = FALSE)[1, ]

  gs %>% mutate(Site = factor(clean_site(Site), levels = LEV_SITE))
}

## Collapse the (usually 3) ROIs per patient/site/area to their mean.
average_rois <- function(gs) {
  avg <- gs %>%
    mutate(combinedFactor = paste(patient, Site2, Area, sep = "_")) %>%
    select(-(1:5)) %>%
    group_by(combinedFactor) %>%
    summarize(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))

  out <- data.frame(
    matrix(unlist(strsplit(avg$combinedFactor, "_")), byrow = TRUE, ncol = 3),
    avg[, -1], check.names = FALSE)
  colnames(out)[1:3] <- c("patient", "Site", "Area")
  out
}

## Keep one copy of the scores per R session.
.score_cache <- new.env(parent = emptyenv())

roi_scores <- function() {
  if (is.null(.score_cache$roi))
    .score_cache$roi <- gsva_scores(geomx_expression()$d)
  .score_cache$roi
}

## Per-ROI scores for the Figure 2 panels: adjacent-normal and unused areas
## dropped, Site collapsed to Colon / Liver / Lung.
site_roi_scores <- function() {
  roi_scores() %>%
    filter(!Area %in% c("AN", "TBD", "TLS")) %>%
    mutate(Site = factor(Site2, levels = c("Colon", "Liver", "Lung")))
}

## One row per patient x site x area; every model uses this.
site_area_scores <- function() average_rois(roi_scores())
