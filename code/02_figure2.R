## ---------------------------------------------------------------------------
## 02_figure2.R -- Figure 2a-d
##
##   2a  mCAF, tCAF
##   2b  MRC1+CSF1R+ mac, SPP1+ mac
##   2c  Proliferating B cells, plasma cells
##   2d  CD4+ Tcm, CD8+ Tem
##   all eight panels in one file                figure2/Fig2_all.pdf
##
## Figures 2e-g are multiplex IF images and have no code.
##
## The boxes show the individual ROI scores; the p-values come from the limma
## model fit on patient x site x area means (lib/models.R). The nine
## signatures, including aHSC (Figure 4a), are adjusted together. Brackets
## mark contrasts with adjusted p < 0.05.
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
gs.roi <- site_roi_scores()       # ROIs -- what the boxes show
contr  <- site_gsva_family(gs.avg)$contrasts

fig2_cells <- setNames(CELL_TYPE_LABELS[1:8], CELL_TYPES[1:8])

panels <- list()
for (ct in names(fig2_cells)) {
  cn <- contr %>%
    filter(as.character(CellType) == ct, adj.P.Val < 0.05) %>%
    split_contrast()

  panels[[ct]] <- prism_panel(
    to_long(gs.roi, "Site", ct), cn, fig2_cells[[ct]], YLAB_SCORE,
    cols = site_cols, point_cols = site_point_cols,
    ylim = if (ct == "B.MKI67") c(-1, 1) else YLIM_SCORE,
    zero_line = TRUE, shape = "box")
}
save_panel(wrap_plots(panels, nrow = 2), "figure2", "Fig2_all.pdf", w = 11, h = 9)
message("  Fig 2a-d")

write.csv(gs.roi[c("ROI", "patient", "Site", "Area", CELL_TYPES[1:8])],
          fig_path("figure2", "Fig2_points_perROI.csv"), row.names = FALSE)
write.csv(contr %>% filter(CellType %in% CELL_TYPES[1:8]) %>%
            mutate(Stars = p_stars(adj.P.Val)),
          fig_path("figure2", "Fig2_contrasts.csv"), row.names = FALSE)

message("Figure 2 done -> ", normalizePath(results_path("figure2"), mustWork = FALSE))
