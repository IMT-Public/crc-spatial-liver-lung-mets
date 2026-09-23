## ---------------------------------------------------------------------------
## 01_figure1.R -- Figure 1
##
##   1b  CIBERSORTx immune composition by site        figure1/Fig1b_*.pdf
##   1c  Volcano plots, Liver-Colon and Lung-Colon    figure1/Fig1c_*.pdf
##   1d  Hallmark GSEA on the same two contrasts      figure1/Fig1d_*.pdf
##
## Figure 1a is a study schematic and has no code.
##
## Figures 1c and 1d use the per-gene linear mixed-effects results (patient
## random intercept, GeneDetectionRate covariate) in
## DE_res_combinedROIs_lme.xlsx.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(nlme); library(RColorBrewer)
})

source("lib/paths.R")
source("lib/load_geomx.R")
source("lib/de_volcano_gsea.R")

## ---------------------------------------------------------------------------
## 1b -- CIBERSORTx immune composition
## ---------------------------------------------------------------------------
## The CIBERSORTx table is row-aligned to the raw object, so this panel uses
## load_geomx() rather than the filtered geomx_expression().
targets <- load_geomx()
ct.frac <- read.csv(require_input("data_cibersortx",
                                  "CIBERSORTx_Job14_Results_02152024.csv"),
                    as.is = TRUE)
ct.frac <- ct.frac[, !colnames(ct.frac) %in% c("P.value", "Correlation", "RMSE")]
stopifnot(all(rownames(pData(targets)) == ct.frac$Mixture))

df <- data.frame(pData(targets)[c("patient", "Site2", "Area")], ct.frac[, -1])
colnames(df)[1:2] <- c("Patient", "Site")

## Collapse the 18 CIBERSORTx classes to the 12 the panel shows.
df$B.cells       <- df$B.naive + df$B.memory
df$CD4.T.cells   <- df$"T.CD4.naive" + df$"T.CD4.memory"
df$CD8.T.cells   <- df$"T.CD8.naive" + df$"T.CD8.memory"
df$Other.myeloid <- df$pDC + df$mDC + df$"monocytes.C" + df$"monocytes.NC.I"
df <- df[, !grepl("mDC|pDC|naive|memory|monocyte", colnames(df))]
colnames(df)[colnames(df) == "mast"]   <- "mast.cells"
colnames(df)[colnames(df) == "plasma"] <- "plasma.cells"
colnames(df)[colnames(df) == "NK"]     <- "NK.cells"
colnames(df)[colnames(df) == "Treg"]   <- "Tregs"

cells.ordered <- c("neutrophils", "macrophages", "Other.myeloid", "CD4.T.cells",
                   "CD8.T.cells", "B.cells", "plasma.cells", "mast.cells",
                   "NK.cells", "Tregs", "endothelial.cells", "fibroblasts")
df <- data.frame(df[, 1:3], df[, cells.ordered])

## Immune fractions only, renormalised to sum to 1 within each ROI.
df.imm <- df[, !colnames(df) %in% c("fibroblasts", "endothelial.cells")]
df.imm[, -(1:3)] <- sweep(df.imm[, -(1:3)], 1, rowSums(df.imm[, -(1:3)]), "/")

## Average within site x area first, then sum over areas, so an over-sampled
## area does not dominate its site.
prop.site <- df.imm[df.imm$Area != "TLS", colnames(df.imm) != "Patient"] %>%
  gsummary(form = ~ 1 | Site / Area) %>%
  gsummary(form = ~ 1 | Site, FUN = sum) %>%
  reshape2::melt(id.vars = c("Site", "Area"), variable.name = "celltype",
                 value.name = "proportion") %>%
  mutate(Site = factor(Site, levels = c("Colon", "Liver", "Lung"))) %>%
  select(-Area)

p1b <- ggplot(prop.site, aes(x = Site, y = proportion, fill = celltype)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = c(brewer.pal(8, "Dark2")[1:2], brewer.pal(8, "Set3")[1:8])) +
  labs(x = "", y = "Proportion", fill = "Cell type") +
  theme_bw() +
  theme(text         = element_text(colour = "black"),
        axis.text    = element_text(colour = "black"),
        axis.text.x  = element_text(face = "bold", size = 8, colour = "black"),
        axis.title.y = element_text(face = "bold", size = 8, colour = "black"),
        legend.text  = element_text(colour = "black"),
        legend.title = element_text(face = "bold", size = 8, colour = "black"))

ggsave(fig_path("figure1", "Fig1b_cibersortx_stacked_barplot.pdf"), p1b,
       width = 4, height = 5)
write.csv(prop.site, fig_path("figure1", "Fig1b_proportions.csv"),
          row.names = FALSE)
message("  Fig 1b")

## ---------------------------------------------------------------------------
## 1c -- volcanoes
## ---------------------------------------------------------------------------
## Genes of interest labelled on each panel.
FIG1C_LABELS <- list(
  `Liver - Colon` = c("APOC1", "ARID1B", "CCL24", "CD79A", "CTNNB1", "CXCL14",
                      "FN1", "FZD1", "GATA2", "IFI30", "ITGB1", "MZT2B",
                      "SERPINA1", "TIMP1", "WNT3", "ZFP57"),
  `Lung - Colon`  = c("ACTA2", "AGBL5", "ARID1B", "CCL24", "CD163", "CD274",
                      "COL1A1", "CXCL14", "FZD1", "HLA-E", "IGHA1", "IGHM",
                      "IL17RA", "MZT2B", "PFN1", "RNASE1", "SFTPB", "TGFB1",
                      "VEGFA", "ZFP57", "WNT3", "CTNNB1", "GATA2"))

de <- load_de_results()

for (ct in names(FIG1C_LABELS)) {
  nm <- sub(" - ", "_vs_", ct)
  ggsave(fig_path("figure1", sprintf("Fig1c_volcano_%s.pdf", nm)),
         volcano_plot(de, ct, label_genes = FIG1C_LABELS[[ct]]),
         width = 6, height = 5)

  d <- de[de$Contrast == ct, c("Gene", "Contrast", "Estimate", "Pr(>|t|)", "FDR")]
  write.csv(d[order(d$FDR), ],
            fig_path("figure1", sprintf("Fig1c_DE_%s.csv", nm)),
            row.names = FALSE)
  message("  Fig 1c  ", ct)
}

## ---------------------------------------------------------------------------
## 1d -- Hallmark GSEA
## ---------------------------------------------------------------------------
## Every Hallmark set reaching nominal p < 0.05 is shown; sets that also reach
## FDR < 0.05 are starred.
for (ct in names(FIG1C_LABELS)) {
  g  <- hallmark_gsea(de, ct)
  nm <- sub(" - ", "_vs_", ct)
  ggsave(fig_path("figure1", sprintf("Fig1d_gsea_hallmark_%s.pdf", nm)),
         gsea_dotplot(g, ct), width = 6, height = 5)
  write.csv(gsea_table(g), fig_path("figure1", sprintf("Fig1d_gsea_%s.csv", nm)),
            row.names = FALSE)
  message("  Fig 1d  ", ct)
}

message("Figure 1 done -> ", normalizePath(results_path("figure1"), mustWork = FALSE))
