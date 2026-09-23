## ---------------------------------------------------------------------------
## 05_supplementary_figure1.R -- Supplementary Figure 1
##
##   S1b-c  mCAF, tCAF, MRC1+CSF1R+ mac, SPP1+ mac,
##          core vs margin within each site        supplementary_figure1/SuppFig1b-c_all.pdf
##   S1d    stromal signatures plus CXCL12,
##          metastasis vs its own adjacent normal  supplementary_figure1/SuppFig1d_all.pdf
##
## Supplementary Figure 1a has no code.
##
## S1b-c shares its limma family with Figure 3b (core_margin_family in
## lib/models.R). S1d compares five groups:
##
##     Colon_prim   Liver_met   Lung_met   Liver_AN   Lung_AN
##
## The CXCL12 p-values in S1d come from the Figure 4b-c gene family.
## The boxes are per-ROI; the models are fit on patient x group means.
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

## Brackets mark contrasts with adjusted p < 0.05, except for SPP1+ mac,
## where the contrasts shown are listed here.
SUPP1BC_BRACKETS <- list(
  Mph.SPP1 = c("Liver_In-Colon_In", "Lung_In-Liver_In",
               "Liver_Ex-Colon_In", "Lung_Ex-Liver_Ex"))

panels <- list()
for (ct in CELL_TYPES[1:4]) {
  keep <- SUPP1BC_BRACKETS[[ct]]
  cn <- contr.cm %>%
    filter(as.character(CellType) == ct,
           if (is.null(keep)) adj.P.Val < 0.05 else as.character(Contrast) %in% keep) %>%
    split_contrast()

  panels[[ct]] <- prism_panel(
    to_long(roi.cm, "Site_area", ct), cn,
    CELL_TYPE_LABELS[match(ct, CELL_TYPES)], YLAB_SCORE,
    ylim = if (ct %in% c("mCAF", "tCAF")) YLIM_SCORE else c(-1, 1),
    cols = core_margin_cols, border = core_margin_border,
    zero_line = TRUE, shape = "boxonly") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
save_panel(wrap_plots(panels, nrow = 1), "supplementary_figure1",
           "SuppFig1b-c_all.pdf", w = 12.5, h = 4.6)
message("  Supp 1b-c")

write.csv(roi.cm[c("ROI", "patient", "Site2", "Area", "Site_area", CELL_TYPES[1:4])],
          fig_path("supplementary_figure1", "SuppFig1bc_points_perROI.csv"),
          row.names = FALSE)
write.csv(contr.cm %>% filter(CellType %in% CELL_TYPES[1:4]) %>%
            mutate(Stars = p_stars(adj.P.Val)),
          fig_path("supplementary_figure1", "SuppFig1bc_contrasts.csv"),
          row.names = FALSE)

## ===========================================================================
## S1d -- metastasis vs its own adjacent normal
## ===========================================================================
## Contrasts shown on each panel.
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

roi.an$CXCL12 <- as.numeric(ex$d["CXCL12",
  match(sub("[.]dcc$", "", roi.an$ROI), sub("[.]dcc$", "", ex$pd$Sample_ID))])

contr.an <- adjacent_normal_family(gs.avg)$contrasts
contr.cx <- site_expression_family(ex$pd, ex$d)$contrasts %>%
  filter(as.character(Gene) == "CXCL12") %>%
  dplyr::rename(CellType = Gene)

panels <- list()
for (ct in names(supp1d_cells)) {
  src <- if (ct == "CXCL12") contr.cx else contr.an
  cn <- src %>%
    filter(as.character(CellType) == ct,
           as.character(Contrast) %in% SUPP1D_BRACKETS[[ct]]) %>%
    split_contrast()

  if (ct == "CXCL12") {
    ## Map the Colon/Liver/Lung group names onto this panel's axis.
    remap <- c(Colon = "Colon_prim", Liver = "Liver_met", Lung = "Lung_met")
    cn$group1 <- unname(remap[cn$group1])
    cn$group2 <- unname(remap[cn$group2])
  }

  panels[[ct]] <- prism_panel(
    to_long(roi.an, "Site_area", ct), cn, supp1d_cells[[ct]],
    if (ct == "CXCL12") YLAB_COUNT else YLAB_SCORE,
    ylim    = if (ct == "CXCL12") YLIM_CXCL12   else YLIM_SCORE,
    ybreaks = if (ct == "CXCL12") BREAKS_CXCL12 else waiver(),
    cols = an_cols, border = an_border,
    zero_line = (ct != "CXCL12"), shape = "boxonly") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
save_panel(wrap_plots(panels, nrow = 1), "supplementary_figure1",
           "SuppFig1d_all.pdf", w = 16, h = 4.6)
message("  Supp 1d")

write.csv(roi.an[c("ROI", "patient", "Site2", "Area", "Site_area",
                   names(supp1d_cells))],
          fig_path("supplementary_figure1", "SuppFig1d_points_perROI.csv"),
          row.names = FALSE)
write.csv(bind_rows(contr.an %>% filter(CellType %in% names(supp1d_cells)),
                    contr.cx) %>%
            mutate(Stars = p_stars(adj.P.Val)),
          fig_path("supplementary_figure1", "SuppFig1d_contrasts.csv"),
          row.names = FALSE)

message("Supplementary Figure 1 done -> ",
        normalizePath(results_path("supplementary_figure1"), mustWork = FALSE))
