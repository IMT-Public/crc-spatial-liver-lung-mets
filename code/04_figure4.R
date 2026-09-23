## ---------------------------------------------------------------------------
## 04_figure4.R -- Figure 4
##
##   4a  aHSC signature score by site        figure4/Fig4a-c_all.pdf
##   4b  CXCL12 expression by site           figure4/Fig4a-c_all.pdf
##   4c  CXCR4 expression by site            figure4/Fig4a-c_all.pdf
##   4d  PFS by aHSC score in TCGA-COAD/READ figure4/Fig4d_*.pdf
##
##   4a    aHSC is adjusted together with the eight Figure 2 signatures
##         (site_gsva_family in lib/models.R).
##   4b-c  CXCL12 and CXCR4 form their own two-gene family.
##   4d    Cox model on TCGA-COAD/READ bulk RNA-seq.
##
## Panels 4a-c are drawn over patient x site x area means.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(patchwork)
  library(readxl); library(GSVA); library(survival); library(survminer)
})

source("lib/paths.R")
source("lib/load_geomx.R")
source("lib/gsva_scores.R")
source("lib/models.R")
source("lib/panels.R")

## ===========================================================================
## 4a-c -- aHSC, CXCL12, CXCR4 by site
## ===========================================================================
fam.site <- site_gsva_family(site_area_scores())
contr4a  <- fam.site$contrasts %>% filter(CellType == "aHSC")

ex       <- geomx_expression()
fam.expr <- site_expression_family(ex$pd, ex$d)
contr4bc <- fam.expr$contrasts

gs.aHSC <- fam.site$scores[c("patient", "Site", "Area", "aHSC")]
d.avg   <- fam.expr$scores

## Brackets mark contrasts with adjusted p < 0.05.
sig_brackets <- function(tbl, id_col, key) {
  tbl %>% filter(as.character(.data[[id_col]]) == key, adj.P.Val < 0.05) %>%
    split_contrast()
}

p4a <- prism_panel(to_long(gs.aHSC, "Site", "aHSC"),
                   sig_brackets(contr4a, "CellType", "aHSC"),
                   "a-HSC", YLAB_SCORE, ylim = YLIM_AHSC,
                   cols = violin_cols, shape = "violin")

panels <- list(p4a)
for (g in c("CXCL12", "CXCR4")) {
  panels[[g]] <- prism_panel(
    to_long(d.avg, "Site", g), sig_brackets(contr4bc, "Gene", g), g, YLAB_COUNT,
    ylim    = if (g == "CXCL12") YLIM_CXCL12   else YLIM_CXCR4,
    ybreaks = if (g == "CXCL12") BREAKS_CXCL12 else BREAKS_CXCR4,
    cols = violin_cols, shape = "violin")
}
save_panel(wrap_plots(panels, nrow = 1), "figure4", "Fig4a-c_all.pdf",
           w = 8.2, h = 4.4)
message("  Fig 4a-c")

write.csv(gs.aHSC, fig_path("figure4", "Fig4a_points.csv"), row.names = FALSE)
write.csv(contr4a %>% mutate(Stars = p_stars(adj.P.Val)),
          fig_path("figure4", "Fig4a_contrasts.csv"), row.names = FALSE)
write.csv(d.avg, fig_path("figure4", "Fig4bc_points.csv"), row.names = FALSE)
write.csv(contr4bc %>% mutate(Stars = p_stars(adj.P.Val)),
          fig_path("figure4", "Fig4bc_contrasts.csv"), row.names = FALSE)

## ===========================================================================
## 4d -- aHSC score and progression-free survival in TCGA-COAD/READ
## ===========================================================================
## Public bulk RNA-seq from cBioPortal (coadread_tcga_pan_can_atlas_2018),
## scored with the aHSC gene list used on the GeoMx data.
##
## The TCGA expression matrix is not redistributed. The capsule includes only
## the processed per-patient table this panel needs:
##
##   data_tcga/TCGA_COADREAD_aHSC_GSVA_PFS.csv
##     PATIENT_ID   TCGA barcode
##     PFS_MONTHS   progression-free survival time, months
##     PFS_STATUS   1 = progressed, 0 = censored
##     gsva         aHSC GSVA score
##
## To regenerate it from the TCGA files:
##   1. cd code && ./fetch_tcga.sh     downloads the two cBioPortal files
##   2. uncomment the block below and run 04_figure4.R once.

## --- scoring from TCGA (commented out for the reproducible run) ------------
## AHSC_SHEET <- "aHSC"
## markers <- require_input("data_tcga", "some_cell_subtype_markers.xlsx")
## stopifnot(AHSC_SHEET %in% excel_sheets(markers))
## ahsc.genes <- read_xlsx(markers, AHSC_SHEET)[[1]]
##
## clin <- read.delim(require_input("data_tcga", "data_clinical_patient.txt"), skip = 4)
## d.cr <- read.delim(require_input(
##   "data_tcga", "data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt"))
##
## d.cr <- d.cr[!duplicated(d.cr[[1]]) & d.cr[[1]] != "", ]
## colnames(d.cr) <- sub("\\.(01|02|03)$", "", colnames(d.cr))
## rownames(d.cr) <- d.cr[[1]]
## d.cr <- d.cr[, -(1:2)]
## d.cr <- d.cr[apply(d.cr, 1, \(x) sum(is.na(x))) == 0, ]
##
## clin[[1]] <- gsub("-", ".", clin[[1]])
## clin <- clin[clin[[1]] %in% colnames(d.cr), ]
## d.cr <- d.cr[, clin[[1]]]
## stopifnot(all(clin[[1]] == colnames(d.cr)))
##
## gsva.ahsc <- gsva(gsvaParam(as.matrix(d.cr), list(aHSC = ahsc.genes)),
##                   verbose = FALSE)[1, ]
## clin$PFS_STATUS <- as.numeric(sub("\\:.*", "", clin$PFS_STATUS))
## clin$gsva <- gsva.ahsc
##
## write.csv(clin[c("PATIENT_ID", "PFS_MONTHS", "PFS_STATUS", "gsva")],
##           data_path("data_tcga", "TCGA_COADREAD_aHSC_GSVA_PFS.csv"),
##           row.names = FALSE)
## ---------------------------------------------------------------------------

clin <- read.csv(require_input("data_tcga", "TCGA_COADREAD_aHSC_GSVA_PFS.csv"))

## High = score above the 60th percentile, Low = at or below the 40th
## percentile. The Cox model is fit on the continuous score over all patients.
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

ggsave(fig_path("figure4", "Fig4d_aHSC_TCGA-COADREAD_PFS.pdf"), plot = p4d,
       width = 9, height = 5.5)
write.csv(clin[c("PATIENT_ID", "PFS_MONTHS", "PFS_STATUS", "gsva", "GSVA.HL")],
          fig_path("figure4", "Fig4d_tcga_patients.csv"), row.names = FALSE)
message("  Fig 4d")

message("Figure 4 done -> ", normalizePath(results_path("figure4"), mustWork = FALSE))
