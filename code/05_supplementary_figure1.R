## ---------------------------------------------------------------------------
## 05_supplementary_figure1.R -- Supplementary Figure 1
##
##   S1b-c  mCAF, tCAF, MRC1+CSF1R+ mac, SPP1+ mac,
##          core vs margin within each site      supplementary_figure1/SuppFig1bc_*.pdf
##   S1d    the same stromal signatures plus CXCL12,
##          metastasis vs its own adjacent normal  supplementary_figure1/SuppFig1d_*.pdf
##
## Supplementary Figure 1a has no code.
##
## S1b-c SHARES ITS ADJUSTMENT FAMILY WITH FIGURE 3b. `core_margin_family()`
## fits eight signatures across 15 pairwise contrasts under one
## Benjamini-Hochberg adjustment; 03_figure3.R takes cell types 5-8 and this
## script takes 1-4. The fit is repeated here so the figure stands alone.
##
## S1d IS A DIFFERENT FAMILY AND A DIFFERENT GROUPING. Five groups, not six:
##
##     Colon_prim   Liver_met   Liver_AN   Lung_met   Lung_AN
##
## i.e. each metastasis against the adjacent normal tissue of its own organ,
## with the primary as reference. CXCL12 is drawn from the Figure 4b-c gene
## family rather than from the signature family, which is why that panel's
## p-values come from a separate fit.
##
## As on Figure 3b, the boxes are per-ROI (208 here, including the 48
## adjacent-normal ROIs) while the models are fit on patient x group means.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(patchwork)
})

source("lib/paths.R")
source("lib/load_geomx.R")
source("lib/gsva_scores.R")
source("lib/models.R")
source("lib/panels.R")

gs.avg <- site_area_scores()
roi    <- roi_scores()
ex     <- geomx_expression()

## ===========================================================================
## S1b-c -- core vs margin
## ===========================================================================
contr.cm <- core_margin_family(gs.avg)$contrasts

SA6 <- paste0(rep(c("Colon", "Liver", "Lung"), each = 2), "_", c("In", "Ex"))
roi.cm <- roi %>%
  filter(Area %in% c("TIL", "T", "TB")) %>%
  mutate(Site_area = factor(
    paste0(Site2, "_", ifelse(Area %in% c("TIL", "T"), "In", "Ex")),
    levels = SA6))
stopifnot(nrow(roi.cm) == 160)

## Three of these four panels print every contrast reaching adj.P.Val < 0.05.
## SPP1+ mac prints a subset -- the In->In and Ex->Ex chains, 4 of its 8 --
## measured off the published figure.
SUPP1BC_BRACKETS <- list(
  Mph.SPP1 = c("Liver_In-Colon_In", "Lung_In-Liver_In",
               "Liver_Ex-Colon_In", "Lung_Ex-Liver_Ex"))

panels.bc <- list()
for (ct in CELL_TYPES[1:4]) {
  keep <- SUPP1BC_BRACKETS[[ct]]
  cn <- if (is.null(keep)) {
    contr.cm %>% filter(as.character(CellType) == ct, adj.P.Val < 0.05) %>%
      split_contrast(p_col = "adj.P.Val")
  } else {
    out <- contr.cm %>%
      filter(as.character(CellType) == ct, as.character(Contrast) %in% keep) %>%
      split_contrast(p_col = "adj.P.Val")
    stopifnot(nrow(out) == length(keep))
    out
  }

  panels.bc[[ct]] <- prism_panel(
    to_long(roi.cm, "Site_area", ct), cn,
    CELL_TYPE_LABELS[match(ct, CELL_TYPES)], YLAB_SCORE,
    cols = core_margin_cols, border = core_margin_border,
    ylim = if (ct %in% c("mCAF", "tCAF")) YLIM_SCORE else c(-1, 1),
    zero_line = TRUE, shape = "boxonly",
    bracket_order = "span", bracket_overflow = TRUE) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  save_panel(panels.bc[[ct]], "supplementary_figure1",
             sprintf("SuppFig1bc_%s.pdf", ct), w = 3.2)
  message(sprintf("  Supp 1b-c %-14s %d bracket(s)", ct, nrow(cn)))
}
save_panel(wrap_plots(panels.bc, nrow = 1), "supplementary_figure1",
           "SuppFig1b-c_all.pdf", w = 12.5, h = 4.6)

write.csv(roi.cm[c("ROI", "patient", "Site2", "Area", "Site_area", CELL_TYPES[1:8])],
          fig_path("supplementary_figure1", "SuppFig1bc_points_perROI.csv"),
          row.names = FALSE)
write.csv(contr.cm %>% mutate(Stars_adj.P.Val        = p_stars(adj.P.Val),
                              Stars_adj.P.Val.global = p_stars(adj.P.Val.global)),
          fig_path("supplementary_figure1", "SuppFig1bc_contrasts_with_stars.csv"),
          row.names = FALSE)

## ===========================================================================
## S1d -- metastasis vs its own adjacent normal
## ===========================================================================
## Bracket subsets measured off the published figure.
SUPP1D_BRACKETS <- list(
  tCAF     = c("Lung_met-Liver_met", "Lung_met-Liver_AN",
               "Liver_met-Liver_AN", "Lung_AN-Liver_met"),
  mCAF     = c("Lung_met-Liver_met", "Liver_met-Liver_AN", "Lung_AN-Liver_met"),
  Mph.PLTP = c("Liver_met-Colon_prim", "Lung_met-Liver_met", "Liver_met-Liver_AN"),
  aHSC     = c("Lung_met-Liver_met", "Liver_met-Liver_AN", "Liver_AN-Colon_prim"),
  CXCL12   = c("Liver-Colon", "Lung-Liver"))

supp1d_cells <- c(tCAF = "tCAF", mCAF = "mCAF", Mph.PLTP = "MRC1+CSF1R+ mac",
                  aHSC = "a-HSC", CXCL12 = "CXCL12")

roi.an <- roi %>%
  filter(!Area %in% c("TBD", "TLS")) %>%
  mutate(Site_area = factor(
    ifelse(Area == "AN", paste0(Site2, "_AN"),
           ifelse(Site2 == "Colon", "Colon_prim", paste0(Site2, "_met"))),
    levels = AN_LEVELS))
stopifnot(nrow(roi.an) == 208, !any(is.na(roi.an$Site_area)))

roi.an$CXCL12 <- as.numeric(ex$d["CXCL12",
  match(sub("[.]dcc$", "", roi.an$ROI), sub("[.]dcc$", "", ex$pd$Sample_ID))])
stopifnot(!any(is.na(roi.an$CXCL12)))

contr.an <- adjacent_normal_family(gs.avg)$by_site$contrasts
contr.cx <- site_expression_family(ex$pd, ex$d)$contrasts %>%
  filter(as.character(Gene) == "CXCL12") %>%
  dplyr::rename(CellType = Gene)

panels.d <- list()
for (ct in names(supp1d_cells)) {
  src  <- if (ct == "CXCL12") contr.cx else contr.an
  keep <- SUPP1D_BRACKETS[[ct]]
  cn <- src %>%
    filter(as.character(CellType) == ct, as.character(Contrast) %in% keep) %>%
    split_contrast(p_col = "adj.P.Val")
  stopifnot(nrow(cn) == length(keep))

  if (ct == "CXCL12") {
    ## The Figure 4b-c family names its groups Colon/Liver/Lung; this panel's
    ## axis uses the adjacent-normal naming. Without the remap ggpubr silently
    ## adds three phantom x categories.
    remap <- c(Colon = "Colon_prim", Liver = "Liver_met", Lung = "Lung_met")
    cn$group1 <- unname(remap[cn$group1])
    cn$group2 <- unname(remap[cn$group2])
    stopifnot(!any(is.na(c(cn$group1, cn$group2))))
  }

  panels.d[[ct]] <- prism_panel(
    to_long(roi.an, "Site_area", ct), cn, supp1d_cells[[ct]],
    if (ct == "CXCL12") YLAB_COUNT else YLAB_SCORE,
    cols = an_cols, border = an_border,
    ylim    = if (ct == "CXCL12") YLIM_CXCL12   else YLIM_SCORE,
    ybreaks = if (ct == "CXCL12") BREAKS_CXCL12 else waiver(),
    zero_line = (ct != "CXCL12"), shape = "boxonly",
    bracket_order = "span", bracket_overflow = TRUE) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  save_panel(panels.d[[ct]], "supplementary_figure1",
             sprintf("SuppFig1d_%s.pdf", ct), w = 3.4)
  message(sprintf("  Supp 1d   %-14s %d bracket(s)", ct, nrow(cn)))
}
save_panel(wrap_plots(panels.d, nrow = 1), "supplementary_figure1",
           "SuppFig1d_all.pdf", w = 16, h = 4.6)

write.csv(roi.an[c("ROI", "patient", "Site2", "Area", "Site_area",
                   names(supp1d_cells))],
          fig_path("supplementary_figure1", "SuppFig1d_points_perROI.csv"),
          row.names = FALSE)
write.csv(bind_rows(contr.an, contr.cx) %>%
            mutate(Stars_adj.P.Val        = p_stars(adj.P.Val),
                   Stars_adj.P.Val.global = p_stars(adj.P.Val.global)),
          fig_path("supplementary_figure1", "SuppFig1d_contrasts_with_stars.csv"),
          row.names = FALSE)

message("Supplementary Figure 1 done -> ",
        normalizePath(results_path("supplementary_figure1"), mustWork = FALSE))
