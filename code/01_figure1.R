## ---------------------------------------------------------------------------
## 01_figure1.R -- Figure 1
##
##   1b  CIBERSORTx immune composition by site        figure1/Fig1b_*.pdf
##   1c  Volcano plots, Liver-Colon and Lung-Colon    figure1/Fig1c_*.pdf
##   1d  Hallmark GSEA on the same two contrasts      figure1/Fig1d_*.pdf
##
## Figure 1a is a study schematic and has no code.
##
## FIGURE 1c IS A MIXED MODEL
## --------------------------
## Figure 1c is a per-gene linear mixed-effects fit -- patient random
## intercept, GeneDetectionRate covariate -- read from
## `DE_res_combinedROIs_lme.xlsx`, the same table Figure 1d is ranked from. It
## takes hours over ~18k genes and is not re-run here; it ships as an input.
##
## The panels as originally SUBMITTED were something else: a per-gene Welch
## two-sample t-test on the `log_q` assay, adjacent-normal ROIs dropped, the 11
## TLS ROIs KEPT (171 ROIs), Benjamini-Hochberg within each 6,188-gene
## contrast, which reproduces the journal's Source Data exactly, to 5.6e-17.
## That fit treats 171 ROIs from 8 patients as independent observations, which
## they are not, and it is kept only so the submitted panel can be regenerated
## and checked:
##
##   CRC_FIG1C_METHOD=lme     (default)  the mixed model -- the published fit
##   CRC_FIG1C_METHOD=welch              the originally submitted panel
##
## The two agree on direction (1.000), on effect size (r = 0.986 / 0.987) and
## on 34 of the 36 labelled genes; IGHM and IGHA1 clear FDR 0.05 under the
## t-test (0.014, 0.050) and not under the mixed model (0.069, 0.223).
##
## FIGURE 1d IS ALWAYS THE MIXED MODEL, whatever CRC_FIG1C_METHOD says.
## Ranking fgsea on load_de_results() reproduces the printed dot plots set for
## set, in order; ranking on the t-test gives seven sets and three instead of
## six and five. The switch above must not reach 1d.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(nlme); library(RColorBrewer)
})

source("lib/paths.R")
source("lib/load_geomx.R")
source("lib/de_volcano_gsea.R")

FIG1C_METHOD <- match.arg(Sys.getenv("CRC_FIG1C_METHOD", unset = "lme"),
                          c("lme", "welch"))

## ---------------------------------------------------------------------------
## 1b -- CIBERSORTx immune composition
## ---------------------------------------------------------------------------
## The CIBERSORTx table is row-aligned to the raw object, so this one panel
## uses load_geomx() rather than the filtered geomx_expression().
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
message("  Fig 1b  immune fractions")

## ---------------------------------------------------------------------------
## 1c -- volcanoes
## ---------------------------------------------------------------------------
## The published panels label a hand-picked set of genes of interest, not an
## algorithmic top-N. The lists were lifted as text from the Source Data vector
## PDFs. Note the lung panel's label reads IL17RA, not IL7R.
FIG1C_PUBLISHED_LABELS <- list(
  `Liver - Colon` = c("APOC1", "ARID1B", "CCL24", "CD79A", "CTNNB1", "CXCL14",
                      "FN1", "FZD1", "GATA2", "IFI30", "ITGB1", "MZT2B",
                      "SERPINA1", "TIMP1", "WNT3", "ZFP57"),
  `Lung - Colon`  = c("ACTA2", "AGBL5", "ARID1B", "CCL24", "CD163", "CD274",
                      "COL1A1", "CXCL14", "FZD1", "HLA-E", "IGHA1", "IGHM",
                      "IL17RA", "MZT2B", "PFN1", "RNASE1", "SFTPB", "TGFB1",
                      "VEGFA", "ZFP57"))

## The Wnt/beta-catenin and GATA2 genes are labelled on the published
## Liver-Colon panel but not on Lung-Colon, where all four are also
## significant. Adding them makes the Lung-Colon panel a superset of the print;
## the two sets are kept apart so the published one stays auditable.
FIG1C_EXTRA_LABELS <- c("ARID1B", "WNT3", "CTNNB1", "GATA2")
FIG1C_LABELS <- lapply(FIG1C_PUBLISHED_LABELS,
                       \(x) union(x, FIG1C_EXTRA_LABELS))

de <- switch(FIG1C_METHOD,
             welch = welch_de_results(names(FIG1C_LABELS)),
             lme   = load_de_results())
message("  Fig 1c  differential expression: ", FIG1C_METHOD)

for (ct in names(FIG1C_LABELS)) {
  labs <- FIG1C_LABELS[[ct]]
  miss <- setdiff(labs, de$Gene[de$Contrast == ct])
  if (length(miss))
    message("    !! label absent from the DE table: ", paste(miss, collapse = ", "))

  nm <- sub(" - ", "_vs_", ct)
  ggsave(fig_path("figure1", sprintf("Fig1c_volcano_%s.pdf", nm)),
         volcano_plot(de, ct, label_genes = labs), width = 6, height = 5)

  ## The differential expression behind the panel: every gene, with the
  ## estimate, the nominal p, the adjusted p and the band it falls in.
  d <- de[de$Contrast == ct, c("Gene", "Contrast", "Estimate", "Pr(>|t|)",
                               "FDR", "Color")]
  write.csv(d[order(d$FDR), ],
            fig_path("figure1", sprintf("Fig1c_DE_%s.csv", nm)),
            row.names = FALSE)

  ## Every label with the numbers behind it, so the panel can be checked
  ## without re-reading the figure.
  st <- de[de$Contrast == ct & de$Gene %in% labs,
           c("Gene", "Estimate", "Pr(>|t|)", "FDR", "Color")]
  st$published  <- st$Gene %in% FIG1C_PUBLISHED_LABELS[[ct]]
  st$passes_FDR <- st$FDR < 0.05
  st$passes_FC  <- abs(st$Estimate) > FC_CUTOFF
  write.csv(st[order(st$FDR), ],
            fig_path("figure1", sprintf("Fig1c_label_status_%s.csv", nm)),
            row.names = FALSE)
  message(sprintf("    %-14s %d labels, %d below FDR 0.05, %d below |log2FC| %.1f",
                  nm, nrow(st), sum(!st$passes_FDR), sum(!st$passes_FC), FC_CUTOFF))
}

## ---------------------------------------------------------------------------
## 1d -- Hallmark GSEA
## ---------------------------------------------------------------------------
## Always the mixed model, never FIG1C_METHOD: see the header. The panel shows
## every Hallmark set reaching the nominal p < 0.05 -- six on Liver - Colon,
## five on Lung - Colon -- with an asterisk on those that also survive BH.
de_1d <- if (FIG1C_METHOD == "lme") de else load_de_results()

for (ct in names(FIG1C_LABELS)) {
  g  <- hallmark_gsea(de_1d, ct)
  nm <- sub(" - ", "_vs_", ct)
  ggsave(fig_path("figure1", sprintf("Fig1d_gsea_hallmark_%s.pdf", nm)),
         gsea_dotplot(g, ct), width = 6, height = 5)

  tb <- gsea_table(g)
  write.csv(tb, fig_path("figure1", sprintf("Fig1d_gsea_%s.csv", nm)),
            row.names = FALSE)
  message(sprintf("  Fig 1d  %-14s %d sets plotted (p < 0.05), %d starred (FDR < 0.05)",
                  nm, sum(tb$plotted), sum(tb$starred)))
}

message("Figure 1 done -> ", normalizePath(results_path("figure1"), mustWork = FALSE))
