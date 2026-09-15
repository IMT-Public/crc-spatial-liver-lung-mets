## ---------------------------------------------------------------------------
## lib/gsva_scores.R -- per-ROI signature scores, averaged to one row per
##                      patient x site x area
##
## Most scores were computed once from the Danaher/MCP-counter marker lists and
## ship as GSVAscores_updated_03272025.csv, a capsule input. Two are
## scored here because their gene lists arrived later:
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

  ## Scores are assigned positionally, so ROI order must agree. The CSV keeps
  ## the ".dcc" suffix that the expression matrix drops -- compare like with
  ## like, or the check silently passes over a real misalignment.
  stopifnot(identical(sub("\\.dcc$", "", gs$ROI), colnames(d)))

  gs$CD4.Tn.CCR7 <- gsva(gsvaParam(d, list(CD4.Tn.CCR7 = top.genes)),
                         verbose = FALSE)[1, ]

  ## 7 of the 9 plasma genes are on the GeoMx WTA panel (MS4A1 and CD38 are
  ## not); GSVA scores the intersection, as it did for the published figures.
  found <- intersect(PLASMA_GENES, rownames(d))
  message("plasma signature: ", length(found), "/", length(PLASMA_GENES),
          " genes found (missing: ",
          paste(setdiff(PLASMA_GENES, found), collapse = ", "), ")")
  stopifnot(length(found) > 0)

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

## Scoring takes a minute and Figure 2 needs both the per-ROI and the averaged
## table in the same session; keep one copy per R session.
.score_cache <- new.env(parent = emptyenv())

roi_scores <- function() {
  if (is.null(.score_cache$roi))
    .score_cache$roi <- gsva_scores(geomx_expression()$d)
  .score_cache$roi
}

## Per-ROI scores in the shape the panels want: adjacent-normal and unused
## areas dropped, Site collapsed to Colon / Liver / Lung.
##
## This is what the dots in the published Figure 2 are: one dot per ROI --
## 25 colon, 69 liver, 66 lung -- not one per patient, and not one per cell.
## The models still run on the averaged table below, so the panel draws the
## ROI-level spread over p-values fit on patient means. See lib/panels.R.
site_roi_scores <- function() {
  roi_scores() %>%
    filter(!Area %in% c("AN", "TBD", "TLS")) %>%
    mutate(Site = factor(Site2, levels = c("Colon", "Liver", "Lung")))
}

## The one call every model makes: one row per patient x site x area.
site_area_scores <- function() average_rois(roi_scores())
