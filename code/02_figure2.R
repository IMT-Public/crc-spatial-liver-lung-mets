## ---------------------------------------------------------------------------
## 02_figure2.R -- Figure 2a-d
##
##   2a  mCAF, tCAF                            figure2/Fig2_mCAF.pdf, _tCAF.pdf
##   2b  MRC1+CSF1R+ mac, SPP1+ mac            figure2/Fig2_Mph.*.pdf
##   2c  Proliferating B cells, plasma cells   figure2/Fig2_B.MKI67.pdf, ...
##   2d  CD4+ Tcm, CD8+ Tem                    figure2/Fig2_CD4*.pdf, _CD8*.pdf
##   all eight on one grid                     figure2/Fig2_all.pdf
##
## Figures 2e-g are multiplex IF images and have no code.
##
## Two things worth knowing about these panels:
##
## THE MODEL IS FIT ON PATIENT MEANS, THE BOXES SHOW ROIs. `site_gsva_family()`
## fits limma with duplicateCorrelation blocking on patient, over one row per
## patient x site x area. The published boxes are the 160 individual ROI
## scores. Both are correct and they are not the same distribution; the panel
## shows the data, the brackets show the model.
##
## THE FAMILY SPANS TWO FIGURES. Nine signatures x three site contrasts get a
## SINGLE Benjamini-Hochberg adjustment, and the ninth signature is aHSC,
## which is Figure 4a. 04_figure4.R refits this same family and selects the
## aHSC rows rather than adjusting aHSC on its own -- doing that would give a
## different p-value from the published one. The fit takes seconds; the
## guarantee is that either figure reproduces without the other.
##
## Brackets are every contrast reaching adj.P.Val < 0.05. That rule reproduces
## the published bracket set exactly on all eight panels; nothing is
## transcribed from the print here.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(patchwork)
})

source("lib/paths.R")
source("lib/load_geomx.R")
source("lib/gsva_scores.R")
source("lib/models.R")
source("lib/panels.R")

gs.avg <- site_area_scores()      # patient x site x area -- what the model sees
gs.roi <- site_roi_scores()       # 160 ROIs -- what the boxes show
fam    <- site_gsva_family(gs.avg)
contr  <- fam$contrasts

fig2_cells <- setNames(CELL_TYPE_LABELS[1:8], CELL_TYPES[1:8])

panels <- list()
for (ct in names(fig2_cells)) {
  cn <- contr %>%
    filter(as.character(CellType) == ct, adj.P.Val < 0.05) %>%
    split_contrast(p_col = "adj.P.Val")

  panels[[ct]] <- prism_panel(
    to_long(gs.roi, "Site", ct), cn, fig2_cells[[ct]], YLAB_SCORE,
    cols = site_cols, point_cols = site_point_cols,
    ## Proliferating B cells keep their own range; the rest share one.
    ylim = if (ct == "B.MKI67") c(-1, 1) else YLIM_SCORE,
    zero_line = TRUE, shape = "box",
    bracket_order = "span", bracket_overflow = TRUE)

  save_panel(panels[[ct]], "figure2", sprintf("Fig2_%s.pdf", ct))
  message(sprintf("  Fig 2   %-14s %d bracket(s)", ct, nrow(cn)))
}

save_panel(wrap_plots(panels, nrow = 2), "figure2", "Fig2_all.pdf", w = 11, h = 9)

write.csv(gs.roi[c("ROI", "patient", "Site", "Area", CELL_TYPES)],
          fig_path("figure2", "Fig2_points_perROI.csv"), row.names = FALSE)
write.csv(contr %>% mutate(Stars_adj.P.Val        = p_stars(adj.P.Val),
                           Stars_adj.P.Val.global = p_stars(adj.P.Val.global)),
          fig_path("figure2", "Fig2_contrasts_with_stars.csv"), row.names = FALSE)

message("Figure 2 done -> ", normalizePath(results_path("figure2"), mustWork = FALSE))
