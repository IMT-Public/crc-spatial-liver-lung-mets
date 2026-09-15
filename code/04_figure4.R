## ---------------------------------------------------------------------------
## 04_figure4.R -- Figure 4
##
##   4a  aHSC signature score by site        figure4/Fig4a_aHSC.pdf
##   4b  CXCL12 expression by site           figure4/Fig4b_CXCL12.pdf
##   4c  CXCR4 expression by site            figure4/Fig4c_CXCR4.pdf
##   4d  PFS by aHSC score in TCGA-COAD/READ figure4/Fig4d_*.pdf
##
## Three separate adjustment families meet in this figure:
##
##   4a    the NINE-signature site family, shared with Figure 2a-d. aHSC is
##         the ninth. This script refits that whole family and selects the
##         aHSC rows -- adjusting aHSC alone would give a different p-value
##         from the published one.
##   4b-c  CXCL12 and CXCR4 are tested as two genes rather than signature
##         scores, so they form their own family.
##   4d    a Cox model on public bulk RNA-seq; no multiplicity at all.
##
## Panels 4a-c are published as violins, unlike the box plots of Figures 2-3,
## and are drawn over patient x site x area means (53 rows) rather than ROIs.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(patchwork)
  library(readxl); library(writexl)
  library(GSVA); library(survival); library(survminer)
})

source("lib/paths.R")
source("lib/load_geomx.R")
source("lib/gsva_scores.R")
source("lib/models.R")
source("lib/panels.R")

## ===========================================================================
## 4a-c -- aHSC, CXCL12, CXCR4 by site
## ===========================================================================
gs.avg   <- site_area_scores()
fam.site <- site_gsva_family(gs.avg)
contr4a  <- fam.site$contrasts %>% filter(CellType == "aHSC")

ex       <- geomx_expression()
fam.expr <- site_expression_family(ex$pd, ex$d)
contr4bc <- fam.expr$contrasts

gs.aHSC <- fam.site$scores[c("patient", "Site", "Area", "aHSC")]
d.avg   <- fam.expr$scores
gs.aHSC$Site <- factor(gs.aHSC$Site, levels = c("Colon", "Liver", "Lung"))
d.avg$Site   <- factor(d.avg$Site,   levels = c("Colon", "Liver", "Lung"))

## Brackets are every contrast reaching adj.P.Val < 0.05, which reproduces
## the published set on all three panels.
auto_brackets <- function(tbl, id_col, key) {
  tbl %>% filter(as.character(.data[[id_col]]) == key, adj.P.Val < 0.05) %>%
    split_contrast(p_col = "adj.P.Val")
}

p4a <- prism_panel(to_long(gs.aHSC, "Site", "aHSC"),
                   auto_brackets(contr4a, "CellType", "aHSC"),
                   "a-HSC", YLAB_SCORE, cols = violin_cols, ylim = YLIM_AHSC,
                   shape = "violin", bracket_order = "span",
                   bracket_overflow = TRUE)
save_panel(p4a, "figure4", "Fig4a_aHSC.pdf")
message("  Fig 4a  aHSC")

panels <- list(p4a)
for (g in c("CXCL12", "CXCR4")) {
  cn <- auto_brackets(contr4bc, "Gene", g)
  p <- prism_panel(to_long(d.avg, "Site", g), cn, g, YLAB_COUNT,
                   cols = violin_cols, shape = "violin",
                   ylim    = if (g == "CXCL12") YLIM_CXCL12   else YLIM_CXCR4,
                   ybreaks = if (g == "CXCL12") BREAKS_CXCL12 else BREAKS_CXCR4,
                   bracket_order = "span", bracket_overflow = TRUE)
  save_panel(p, "figure4",
             sprintf("Fig4%s_%s.pdf", c(CXCL12 = "b", CXCR4 = "c")[g], g))
  panels[[length(panels) + 1]] <- p
  message(sprintf("  Fig 4b-c %-8s %d bracket(s)", g, nrow(cn)))
}
save_panel(wrap_plots(panels, nrow = 1), "figure4", "Fig4a-c_all.pdf",
           w = 8.2, h = 4.4)

write_xlsx(list(GSVAscores_aHSC      = gs.aHSC,
                Fig4a.gsvaScore      = contr4a,
                exprn_data           = d.avg,
                Fig4bc.logExpression = contr4bc),
           fig_path("figure4", "contr-results-Fig4abc.xlsx"))
write.csv(gs.aHSC, fig_path("figure4", "Fig4a_points.csv"), row.names = FALSE)
write.csv(d.avg,   fig_path("figure4", "Fig4bc_points.csv"), row.names = FALSE)

## ===========================================================================
## 4d -- aHSC score and progression-free survival in TCGA-COAD/READ
## ===========================================================================
## Public bulk RNA-seq from cBioPortal, scored with the same aHSC gene list.
## some_cell_subtype_markers.xlsx holds several signatures; the sheet is
## selected BY NAME, not by position.
AHSC_SHEET <- "aHSC"
markers <- require_input("data_tcga", "some_cell_subtype_markers.xlsx")
stopifnot(AHSC_SHEET %in% excel_sheets(markers))
ahsc.genes <- read_xlsx(markers, AHSC_SHEET)[[1]]

clin <- read.delim(require_input("data_tcga", "data_clinical_patient.txt"), skip = 4)
d.cr <- read.delim(require_input(
  "data_tcga", "data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt"))

d.cr <- d.cr[!duplicated(d.cr[[1]]) & d.cr[[1]] != "", ]
colnames(d.cr) <- sub("\\.(01|02|03)$", "", colnames(d.cr))
rownames(d.cr) <- d.cr[[1]]
d.cr <- d.cr[, -(1:2)]
d.cr <- d.cr[apply(d.cr, 1, \(x) sum(is.na(x))) == 0, ]

clin[[1]] <- gsub("-", ".", clin[[1]])
clin <- clin[clin[[1]] %in% colnames(d.cr), ]
d.cr <- d.cr[, clin[[1]]]
stopifnot(all(clin[[1]] == colnames(d.cr)))

gsva.ahsc <- gsva(gsvaParam(as.matrix(d.cr), list(aHSC = ahsc.genes)),
                  verbose = FALSE)[1, ]
clin$PFS_STATUS <- as.numeric(sub("\\:.*", "", clin$PFS_STATUS))
clin$gsva <- gsva.ahsc

## Patients are stratified by the UPPER AND LOWER 40% of the score, which is
## what the published panel shows and where its n = 237 per group comes from:
##
##   patients with both expression and PFS data        592
##   Low  = score at or below the 40th percentile      0.40 x 592 = 237
##   High = score above the 60th percentile            0.40 x 592 = 237
##   the middle 20%, plotted in neither curve          118 dropped
##
## A true median split would put all 592 in the plot at 296 per group. The
## Methods say "upper and lower 35%", which gives 207 per group and does not
## match the panel; the split here reproduces the figure. The Cox model is fit
## on the CONTINUOUS score over all 592 patients, so it does not depend on the
## split at all.
q <- quantile(clin$gsva, c(.40, .60))
clin$GSVA.HL <- NA
clin$GSVA.HL[clin$gsva <= q[1]] <- "L"
clin$GSVA.HL[clin$gsva >  q[2]] <- "H"
n.hl <- table(clin$GSVA.HL)

surv <- with(clin, Surv(PFS_MONTHS, PFS_STATUS))
cox  <- coxph(surv ~ gsva, data = clin)
p.cox <- round(summary(cox)$coef[1, 5], 3)
hr    <- round(summary(cox)$coef[1, 2], 3)

## survfit orders strata by the sorted factor levels, so H comes first.
p4d <- ggsurvplot(survfit(surv ~ GSVA.HL, data = clin), clin, pval = TRUE,
                  legend.title = "",
                  legend.labs = sprintf("%s signature score (n=%d)",
                                        c("High", "Low"), n.hl[c("H", "L")]))$plot +
  labs(title = "PFS (TCGA-COAD+READ)",
       subtitle = paste0("CoxPH p = ", p.cox, ", HR = ", hr),
       x = "Time(months)") +
  theme(plot.title    = element_text(size = 14, face = "bold", hjust = .5),
        plot.subtitle = element_text(size = 11, hjust = .5))

## Wider than tall, as published -- the curves separate along the time axis
## and a square canvas squeezes them into the left third.
ggsave(fig_path("figure4", "Fig4d_aHSC_TCGA-COADREAD_PFS.pdf"), plot = p4d,
       width = 9, height = 5.5)
write.csv(clin[c(1, which(colnames(clin) %in%
                          c("PFS_MONTHS", "PFS_STATUS", "gsva", "GSVA.HL")))],
          fig_path("figure4", "Fig4d_tcga_patients.csv"), row.names = FALSE)
message(sprintf("  Fig 4d  N = %d; %d high / %d low; Cox p = %s, HR = %s",
                nrow(clin), n.hl["H"], n.hl["L"], p.cox, hr))

message("Figure 4 done -> ", normalizePath(results_path("figure4"), mustWork = FALSE))
