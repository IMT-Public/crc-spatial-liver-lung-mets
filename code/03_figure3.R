## ---------------------------------------------------------------------------
## 03_figure3.R -- Figure 3
##
##   3b  B, plasma, CD4+ Tcm, CD8+ Tem, core vs margin   figure3/Fig3b_*.pdf
##   3c  Volcanoes, tumour core vs invasive margin       figure3/Fig3c_*.pdf
##
## Figure 3a is a schematic and 3d multiplex IF images; neither has code.
##
## THE POOLING, which both panels share:
##
##     tumour core ("In")     = TIL + T
##     invasive margin ("Ex") = TB
##
## CORE VS MARGIN IS FIT TWICE HERE, at two different granularities, because
## the two panels ask different questions:
##
##   3b  one limma fit over eight signature scores x 15 pairwise contrasts of
##       the six site x compartment groups, blocking on patient. Seconds.
##   3c  one lmerTest fit PER GENE, within each organ separately:
##           exprsn ~ Compartment + GeneDetectionRate + (1 | patient)
##       6,188 genes x 3 organs, Satterthwaite df, BH within organ. ~5 min,
##       which is most of this script's runtime.
##
## The 3c fit runs all three organs -- the Source Data ships colon too, though
## only lung and liver are printed, and colon comes out essentially null,
## which is itself the point of the figure.
##
## Estimate is TIL_T - TB, so POSITIVE = higher in the tumour core. That is
## already the orientation the published panels are drawn in, so the volcanoes
## pass flip = FALSE.
##
## > The 3b fit also carries SUPPLEMENTARY FIGURE 1b-c: the CAF and macrophage
## > signatures are in the same Benjamini-Hochberg family as these four cell
## > types. 05_supplementary_figure1.R re-runs the same fit and selects its own
## > rows. Adjusting either set alone would give different p-values.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(patchwork)
  library(lmerTest); library(writexl)
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
fam    <- core_margin_family(gs.avg)
contr  <- fam$contrasts

## The boxes are per-ROI; the model above is fit on patient x group means.
SA6 <- paste0(rep(c("Colon", "Liver", "Lung"), each = 2), "_", c("In", "Ex"))
roi.cm <- roi_scores() %>%
  filter(Area %in% c("TIL", "T", "TB")) %>%
  mutate(Site_area = factor(
    paste0(Site2, "_", ifelse(Area %in% c("TIL", "T"), "In", "Ex")),
    levels = SA6))
stopifnot(nrow(roi.cm) == 160)

## Unlike Figure 2, Figure 3b does NOT draw every significant contrast: CD4+
## Tcm has 11 and the panel prints 5. The printed subsets were measured off
## the published figure and are listed here. Re-check them against the source
## figure files before any redeposit.
FIG3B_BRACKETS <- list(
  B.MKI67      = character(0),
  Plasma.cells = "Lung_Ex-Liver_In",
  CD4.Tn.CCR7  = c("Liver_In-Colon_In", "Liver_Ex-Liver_In", "Liver_Ex-Colon_Ex",
                   "Lung_In-Colon_In", "Lung_Ex-Colon_Ex"),
  CD8.Tem.GZMK = c("Liver_Ex-Liver_In", "Lung_Ex-Lung_In", "Lung_In-Liver_Ex",
                   "Lung_Ex-Liver_In"))

## Figure 3b does not share one y scale the way Figure 2 does -- each panel
## keeps its own published range.
FIG3B_YLIM <- list(B.MKI67     = c(-0.4, 0.6), Plasma.cells = c(-1.0, 1.5),
                   CD4.Tn.CCR7 = c(-1.0, 1.0), CD8.Tem.GZMK = c(-1.0, 1.0))

stopifnot(setequal(levels(roi.cm$Site_area), names(core_margin_cols)))

panels <- list()
for (ct in CELL_TYPES[5:8]) {
  keep <- FIG3B_BRACKETS[[ct]]
  cn <- contr %>%
    filter(as.character(CellType) == ct, as.character(Contrast) %in% keep) %>%
    split_contrast(p_col = "adj.P.Val")
  stopifnot(nrow(cn) == length(keep))

  panels[[ct]] <- prism_panel(
    to_long(roi.cm, "Site_area", ct), cn,
    CELL_TYPE_LABELS[match(ct, CELL_TYPES)], YLAB_SCORE,
    cols = core_margin_cols, border = core_margin_border,
    ylim = FIG3B_YLIM[[ct]], zero_line = TRUE, shape = "boxonly",
    bracket_order = "span", bracket_overflow = TRUE) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  save_panel(panels[[ct]], "figure3", sprintf("Fig3b_%s.pdf", ct), w = 3.2)
  message(sprintf("  Fig 3b  %-14s %d bracket(s)", ct, nrow(cn)))
}
save_panel(wrap_plots(panels, nrow = 1), "figure3", "Fig3b_all.pdf",
           w = 12.5, h = 4.6)

write.csv(roi.cm[c("ROI", "patient", "Site2", "Area", "Site_area", CELL_TYPES[1:8])],
          fig_path("figure3", "Fig3b_points_perROI.csv"), row.names = FALSE)
write.csv(contr %>% mutate(Stars_adj.P.Val        = p_stars(adj.P.Val),
                           Stars_adj.P.Val.global = p_stars(adj.P.Val.global)),
          fig_path("figure3", "Fig3b_contrasts_with_stars.csv"), row.names = FALSE)

## ===========================================================================
## 3c -- per-gene core vs margin, within each organ
## ===========================================================================
ORGANS <- c("Liver", "Lung", "Colon")
CORE   <- c("TIL", "T")
MARGIN <- "TB"

ex <- geomx_expression()
pd <- ex$pd
d  <- ex$d

fit_organ <- function(site) {
  keep <- pd$Site2 == site & pd$Area %in% c(CORE, MARGIN)
  sub  <- pd[keep, ]
  dm   <- d[, keep]

  ## Levels ordered margin-first, so the model coefficient is core - margin.
  sub$Compartment <- factor(ifelse(sub$Area %in% CORE, "TIL_T", "TB"),
                            levels = c("TB", "TIL_T"))
  message(sprintf("  Fig 3c  %-6s %3d ROIs, %d patients, %d core / %d margin",
                  site, nrow(sub), length(unique(sub$patient)),
                  sum(sub$Compartment == "TIL_T"), sum(sub$Compartment == "TB")))

  out <- vector("list", nrow(dm))
  for (i in seq_len(nrow(dm))) {
    td <- data.frame(exprsn = dm[i, ],
                     sub[, c("GeneDetectionRate", "patient", "Compartment")])
    ## Singular fits occur where between-patient variance is estimated at
    ## zero; the Satterthwaite test on the fixed effect is still what we want.
    m <- suppressMessages(suppressWarnings(
      lmerTest::lmer(exprsn ~ Compartment + GeneDetectionRate + (1 | patient),
                     data = td)))
    lsm <- lmerTest::ls_means(m, which = "Compartment", pairwise = TRUE)
    ## ls_means labels the row "CompartmentTB - CompartmentTIL_T"; negate for
    ## the TIL_T - TB convention. Assert the label, do not assume it.
    stopifnot(identical(trimws(rownames(lsm)[1]),
                        "CompartmentTB - CompartmentTIL_T"))
    out[[i]] <- data.frame(Gene     = rownames(dm)[i],
                           Subset   = site,
                           Contrast = paste0("TIL_T - TB in ", site),
                           Estimate = -as.numeric(lsm[1, "Estimate"]),
                           p        = as.numeric(lsm[1, "Pr(>|t|)"]))
    if (i %% 2000 == 0) message("            ", i, " / ", nrow(dm))
  }

  res <- do.call(rbind, out)
  names(res)[names(res) == "p"] <- "Pr(>|t|)"
  res$FDR <- p.adjust(res[["Pr(>|t|)"]], method = "BH")   # BH within organ
  de_bands(res)
}

fits   <- lapply(setNames(ORGANS, ORGANS), fit_organ)
all.de <- bind_rows(fits)

write.csv(all.de, fig_path("figure3", "Fig3c_DE_TIL_T-TB_all_organs.csv"),
          row.names = FALSE)
write_xlsx(setNames(fits, paste0("TIL_T - TB, ", ORGANS)),
           fig_path("figure3", "Fig3c_DE_coreMargin_lme.xlsx"))

for (org in ORGANS) {
  ct <- paste0("TIL_T - TB in ", org)

  ## Lung and liver have curated label sets; colon, which is not printed,
  ## falls back to its 20 most significant genes.
  labs <- VOLCANO_LABELS[[ct]]
  if (is.null(labs)) {
    hit  <- all.de[all.de$Contrast == ct & all.de$FDR < 0.05 &
                   abs(all.de$Estimate) > FC_CUTOFF, ]
    labs <- hit$Gene[order(hit$FDR)][seq_len(min(20, nrow(hit)))]
  }

  p <- volcano_plot(all.de, ct, label_genes = labs, flip = FALSE,
                    x_lab = "Invasive margin <- log2(FC) -> Tumour core")
  ggsave(fig_path("figure3", sprintf("Fig3c_volcano_coreVsMargin_%s.pdf", org)),
         plot = p, width = 6, height = 5)
  write.csv(all.de[all.de$Contrast == ct, ],
            fig_path("figure3", sprintf("Fig3c_points_%s.csv", org)),
            row.names = FALSE)
  message(sprintf("  Fig 3c  %-6s %d genes at FDR < 0.05 & |log2FC| > %.1f",
                  org, sum(all.de$Contrast == ct & all.de$FDR < 0.05 &
                           abs(all.de$Estimate) > FC_CUTOFF), FC_CUTOFF))
}

message("Figure 3 done -> ", normalizePath(results_path("figure3"), mustWork = FALSE))
