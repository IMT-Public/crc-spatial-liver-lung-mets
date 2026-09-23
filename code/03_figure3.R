## ---------------------------------------------------------------------------
## 03_figure3.R -- Figure 3
##
##   3b  B, plasma, CD4+ Tcm, CD8+ Tem, core vs margin   figure3/Fig3b_all.pdf
##   3c  Volcanoes, tumour core vs invasive margin       figure3/Fig3c_*.pdf
##
## Figure 3a is a schematic and 3d multiplex IF images; neither has code.
##
##     tumour core ("In")     = TIL + T
##     invasive margin ("Ex") = TB
##
##   3b  one limma fit over eight signature scores x 15 pairwise contrasts of
##       the six site x compartment groups, blocking on patient.
##   3c  one lmerTest fit per gene, within liver and within lung:
##           exprsn ~ Compartment + GeneDetectionRate + (1 | patient)
##       Satterthwaite df, BH within organ.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(patchwork); library(lmerTest)
})

source("lib/paths.R")
source("lib/load_geomx.R")
source("lib/gsva_scores.R")
source("lib/models.R")
source("lib/panels.R")
source("lib/de_volcano_gsea.R")

## ===========================================================================
## 3b -- signature scores, core vs margin within each site
## ===========================================================================
gs.avg <- site_area_scores()
contr  <- core_margin_family(gs.avg)$contrasts

## The boxes are per-ROI; the model above is fit on patient x group means.
SA6 <- paste0(rep(c("Colon", "Liver", "Lung"), each = 2), "_", c("In", "Ex"))
roi.cm <- roi_scores() %>%
  filter(Area %in% c("TIL", "T", "TB")) %>%
  mutate(Site_area = factor(
    paste0(Site2, "_", ifelse(Area %in% c("TIL", "T"), "In", "Ex")),
    levels = SA6))

## Contrasts shown on each panel.
FIG3B_BRACKETS <- list(
  B.MKI67      = character(0),
  Plasma.cells = "Lung_Ex-Liver_In",
  CD4.Tn.CCR7  = c("Liver_In-Colon_In", "Liver_Ex-Liver_In", "Liver_Ex-Colon_Ex",
                   "Lung_In-Colon_In", "Lung_Ex-Colon_Ex"),
  CD8.Tem.GZMK = c("Liver_Ex-Liver_In", "Lung_Ex-Lung_In", "Lung_In-Liver_Ex",
                   "Lung_Ex-Liver_In"))

FIG3B_YLIM <- list(B.MKI67     = c(-0.4, 0.6), Plasma.cells = c(-1.0, 1.5),
                   CD4.Tn.CCR7 = c(-1.0, 1.0), CD8.Tem.GZMK = c(-1.0, 1.0))

panels <- list()
for (ct in CELL_TYPES[5:8]) {
  cn <- contr %>%
    filter(as.character(CellType) == ct,
           as.character(Contrast) %in% FIG3B_BRACKETS[[ct]]) %>%
    split_contrast()

  panels[[ct]] <- prism_panel(
    to_long(roi.cm, "Site_area", ct), cn,
    CELL_TYPE_LABELS[match(ct, CELL_TYPES)], YLAB_SCORE,
    ylim = FIG3B_YLIM[[ct]], cols = core_margin_cols,
    border = core_margin_border, zero_line = TRUE, shape = "boxonly") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
save_panel(wrap_plots(panels, nrow = 1), "figure3", "Fig3b_all.pdf",
           w = 12.5, h = 4.6)
message("  Fig 3b")

write.csv(roi.cm[c("ROI", "patient", "Site2", "Area", "Site_area", CELL_TYPES[5:8])],
          fig_path("figure3", "Fig3b_points_perROI.csv"), row.names = FALSE)
write.csv(contr %>% filter(CellType %in% CELL_TYPES[5:8]) %>%
            mutate(Stars = p_stars(adj.P.Val)),
          fig_path("figure3", "Fig3b_contrasts.csv"), row.names = FALSE)

## ===========================================================================
## 3c -- per-gene core vs margin, within each organ
## ===========================================================================
ORGANS <- c("Liver", "Lung")
CORE   <- c("TIL", "T")
MARGIN <- "TB"

ex <- geomx_expression()
pd <- ex$pd
d  <- ex$d

fit_organ <- function(site) {
  keep <- pd$Site2 == site & pd$Area %in% c(CORE, MARGIN)
  sub  <- pd[keep, ]
  dm   <- d[, keep]

  ## Levels ordered margin-first; the estimate is reported as core - margin.
  sub$Compartment <- factor(ifelse(sub$Area %in% CORE, "TIL_T", "TB"),
                            levels = c("TB", "TIL_T"))

  out <- vector("list", nrow(dm))
  for (i in seq_len(nrow(dm))) {
    td <- data.frame(exprsn = dm[i, ],
                     sub[, c("GeneDetectionRate", "patient", "Compartment")])
    m <- suppressMessages(suppressWarnings(
      lmerTest::lmer(exprsn ~ Compartment + GeneDetectionRate + (1 | patient),
                     data = td)))
    lsm <- lmerTest::ls_means(m, which = "Compartment", pairwise = TRUE)
    stopifnot(identical(trimws(rownames(lsm)[1]),
                        "CompartmentTB - CompartmentTIL_T"))
    out[[i]] <- data.frame(Gene     = rownames(dm)[i],
                           Contrast = paste0("TIL_T - TB in ", site),
                           Estimate = -as.numeric(lsm[1, "Estimate"]),
                           p        = as.numeric(lsm[1, "Pr(>|t|)"]))
  }

  res <- do.call(rbind, out)
  names(res)[names(res) == "p"] <- "Pr(>|t|)"
  res$FDR <- p.adjust(res[["Pr(>|t|)"]], method = "BH")   # BH within organ
  res
}

for (org in ORGANS) {
  de <- fit_organ(org)
  ct <- paste0("TIL_T - TB in ", org)

  p <- volcano_plot(de, ct, x_lab = "Invasive margin <- log2(FC) -> Tumour core")
  ggsave(fig_path("figure3", sprintf("Fig3c_volcano_coreVsMargin_%s.pdf", org)),
         plot = p, width = 6, height = 5)
  write.csv(de[order(de$FDR), ],
            fig_path("figure3", sprintf("Fig3c_DE_%s.csv", org)),
            row.names = FALSE)
  message("  Fig 3c  ", org)
}

message("Figure 3 done -> ", normalizePath(results_path("figure3"), mustWork = FALSE))
