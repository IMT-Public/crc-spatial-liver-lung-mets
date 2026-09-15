## ---------------------------------------------------------------------------
## lib/de_volcano_gsea.R -- volcano plots and Hallmark GSEA dot plots
##
## Two differential-expression paths feed the Figure 1c volcanoes, both
## returning the same schema:
##
##   welch_de_results()  per-gene Welch t-test, computed here. This is the fit
##                       the published Figure 1c was drawn from.
##   load_de_results()   the precomputed per-gene linear mixed-effects results
##                       in DE_res_combinedROIs_lme.xlsx. That fit takes hours
##                       over ~18k genes and is not re-run; it ships as an
##                       input.
##
## 01_figure1.R picks between them on CRC_FIG1C_METHOD; see the README section
## "Figure 1c uses a t-test".
##
## The spreadsheet has one sheet per contrast. NOTE: the sheet names and the
## values in the `Contrast` column differ for the last one -- sheet
## "TIL.TB - T__inColon" carries Contrast "TIL_TB - T". Everything here keys
## off the Contrast column, whose verified values are in DE_CONTRASTS.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(fgsea)
})

## Verified against the file itself; see DE_CONTRASTS use in each figure script.
DE_CONTRASTS <- c(liver_colon = "Liver - Colon",
                  liver_lung  = "Liver - Lung",
                  lung_colon  = "Lung - Colon",
                  til_tb_t    = "TIL_TB - T")

## The published volcanoes (Fig 1c, Fig 3c) are a single blue, drawn as open
## circles. The four-way significance banding below is Wenbin's original colour
## scheme; it is still computed because it is useful in the exported table, but
## it is no longer mapped to the plot.
VOLCANO_BLUE <- "#4A7EBB"

## The fold-change cutoff is 0.5, not the 1 in Wenbin's script. Measured off the
## published Fig 1c: the dashed verticals sit 57.25 px either side of x = 0 on
## an axis whose left edge is labelled -5 at 587.5 px from x = 0, i.e. +/-0.49.
## The lower panel, which has a fully ticked axis (-4, -2, 0, 2 at 132 px per
## unit), gives +/-0.506 independently. A cutoff of 1 would need the axis to
## reach -9 to place those genes, which contradicts its own -5 label.
FC_CUTOFF <- 0.5

NS_BAND <- paste0("NS or FC < ", FC_CUTOFF)
VOLCANO_COLORS <- setNames(
  c("dodgerblue", "lightblue", "orange2", "gray"),
  c("FDR < 0.001", "FDR < 0.05", "P < 0.05", NS_BAND))

## All four contrasts in one table, with the significance banding the volcano
## plots colour by. A gene needs |log2 FC| > FC_CUTOFF to be called at all.
load_de_results <- function() {
  f <- require_input("data_volcanoplot_pathway", "DE_res_combinedROIs_lme.xlsx")
  res <- do.call(rbind, lapply(seq_along(DE_CONTRASTS), \(i) read_xlsx(f, i)))

  de_bands(res)
}

## The significance banding and the signed-significance column used to rank
## candidate labels. Shared by both differential-expression paths.
de_bands <- function(res) {
  res$Color <- NS_BAND
  res$Color[res$`Pr(>|t|)` < 0.05] <- "P < 0.05"
  res$Color[res$FDR < 0.05]        <- "FDR < 0.05"
  res$Color[res$FDR < 0.001]       <- "FDR < 0.001"
  res$Color[abs(res$Estimate) < FC_CUTOFF] <- NS_BAND
  res$Color <- factor(res$Color, levels = names(VOLCANO_COLORS)[c(4, 3, 2, 1)])
  res$invert_P <- (-log10(res$`Pr(>|t|)`)) * sign(res$Estimate)
  res
}

## ---------------------------------------------------------------------------
## The published Figure 1c is NOT the mixed model above.
##
## `Fig 1C.xlsx` (the journal's Source Data) is reproduced to 6e-17 -- exactly,
## up to floating point -- by a per-gene Welch two-sample t-test on the `log_q`
## assay: adjacent-normal ROIs excluded, the 11 TLS ROIs RETAINED (171 ROIs),
## Benjamini-Hochberg applied separately within each 6,188-gene contrast.
## The audit that establishes this lives in the working repository; the short
## version is in the README under "Figure 1c uses a t-test".
##
## Two differences from the LME in the shipped spreadsheet: no random effect
## for patient (ROIs from one patient are treated as independent) and no
## GeneDetectionRate covariate. Which of the two is the right analysis is a
## scientific decision, not a reproducibility one. This function reproduces
## what was printed.
##
## Returns the same schema as load_de_results(), so volcano_plot() and
## hallmark_gsea() take it unchanged.
welch_de_results <- function(contrasts = c("Liver - Colon", "Lung - Colon")) {
  targets <- load_geomx()
  pd <- pData(targets)
  m  <- assayDataElement(targets, "log_q")
  stopifnot(!is.null(pd$Site2), !is.null(pd$Area))

  one <- function(contrast) {
    target <- trimws(sub("-.*", "", contrast))
    keep <- pd$Area != "AN" & pd$Site2 %in% c("Colon", target)
    grp  <- factor(as.character(pd$Site2[keep]), levels = c("Colon", target))
    stopifnot(nlevels(droplevels(grp)) == 2)
    z <- t(apply(m[, keep, drop = FALSE], 1, function(y) {
      tt <- t.test(y ~ grp)                     # Welch is R's default
      c(Estimate = unname(diff(tt$estimate)), P = tt$p.value)
    }))
    data.frame(Gene       = rownames(m),
               Subset     = "combinedROIs",
               Contrast   = contrast,
               Estimate   = z[, "Estimate"],
               `Pr(>|t|)` = z[, "P"],
               FDR        = p.adjust(z[, "P"], method = "BH"),
               check.names = FALSE, stringsAsFactors = FALSE)
  }

  res <- do.call(rbind, lapply(contrasts, one))
  rownames(res) <- NULL
  de_bands(res)
}

## The y axis carries -log10(adjusted p) but is labelled in powers of ten
## (10^0, 10^-20, ...), so the reader sees the p-value itself on a log scale --
## the way the published panels are drawn. Breaks are chosen to give ~4 labels
## whatever the contrast's dynamic range: 0-70 for liver-vs-colon, 0-35 for
## liver-vs-lung.
volcano_breaks <- function(y) {
  m <- max(y[is.finite(y)], na.rm = TRUE)
  steps <- c(1, 2, 5, 10, 20, 25, 50, 100)
  step  <- steps[which(steps >= m / 4)][1]
  seq(0, ceiling(m / step) * step, by = step)
}

volcano_labels <- function(breaks) {
  parse(text = ifelse(breaks == 0, "10^0", paste0("10^-", breaks)))
}

## The published panels label a hand-picked set of genes of interest, NOT the
## algorithmic top 15 per side that Wenbin's script produces -- only a handful
## of the printed labels come out of that rule. No such list exists anywhere in
## the original code; the single manual addition in it (CD79A) is commented
## out, so the lists below were read off the published figure.
##
## Both lists match the print label for label -- 16 on the top panel, 20 on
## the bottom -- confirmed from the Source Data vector PDFs and again from a
## 600 dpi render of page 26 of the accepted manuscript. Two earlier readings were wrong; this list carries the fixes:
##   * the bottom-panel label is IL17RA, not "IL7RA"/IL7R. IL17RA is on the
##     panel and is strongly colon-enriched (-0.67, FDR 1.2e-10); IL7R is a
##     different gene and is not significant in that contrast (FDR 0.64). The
##     Results text names IL17RA too.
##   * MZT2B and ZFP57 are labelled on the bottom panel as well, not only the
##     top one.
## Figure 1c is drawn from the MIXED MODEL (load_de_results), the analysis the
## paper stands behind; welch_de_results() is kept because it is what the
## originally submitted panels were drawn from. Three of the 36 labels sit
## outside the coloured band under the mixed model and stay labelled anyway,
## as genes of interest:
##   * CTNNB1  Liver panel, FDR 4.1e-4 but |log2FC| 0.46, inside the 0.5 line
##             (it is outside the band under the t-test too, at 0.43).
##   * IGHM    Lung panel, FDR 0.069.
##   * IGHA1   Lung panel, FDR 0.223.
## IGHM and IGHA1 are the only two labels the choice of fit moves: under the
## t-test they reach FDR 0.014 and 0.050. All 36 labels are tabulated under
## both fits in the working repository's REPRODUCIBILITY notes.
## A contrast with no entry here falls back to the original top-15 rule.
VOLCANO_LABELS <- list(
  `Liver - Colon` = c(
    ## higher in colon
    "MZT2B", "ZFP57", "ARID1B", "WNT3", "FZD1", "CXCL14", "CTNNB1", "CD79A", "GATA2",
    ## higher in liver -- the immunosuppressive / fibrotic program in the Results
    "CCL24", "ITGB1", "TIMP1", "SERPINA1", "FN1", "APOC1", "IFI30"),
  `Lung - Colon` = c(
    ## higher in colon
    "MZT2B", "ZFP57", "ARID1B", "FZD1", "CD274", "IL17RA", "ACTA2", "CXCL14",
    "VEGFA", "COL1A1", "CD163", "TGFB1",
    ## higher in lung
    "AGBL5", "CCL24", "SFTPB", "RNASE1", "PFN1", "HLA-E", "IGHM", "IGHA1"),

  ## The two panels of Figure 3c: tumour core (TIL + T) against invasive
  ## margin (TB), within lung (top) and within liver (bottom). Both contrasts
  ## are refit by 03_figure3.R. Every one of these 31 labels comes out of that
  ## refit with the sign the published panel shows, and all 31 reach
  ## FDR < 0.05. The published panels ring CXCL13/IGHA1/IGHM/JCHAIN (lung) and
  ## CXCR4/CD3E/CD28 (liver) in red; the rings are decoration and are not
  ## reproduced, but every ringed gene is in the lists below.
  `TIL_T - TB in Lung`  = c("CD74", "HLA-DRA", "CXCL13", "CCL19", "IGHA1",
                            "IGHM", "JCHAIN", "CDX2", "S100A6", "KRT8"),
  `TIL_T - TB in Liver` = c("ALB", "APOA1", "CCL19", "IL7R", "CXCR4", "CXCL9",
                            "MMP9", "CCR7", "CD3E", "PTPRC", "HLA-DRA", "CD28",
                            "CD74", "ALDOA", "EPCAM", "KRT8", "MIF", "CEACAM1",
                            "S100A10", "IFI27", "TGFBI")
)

## The Fig 3c label readings that were settled from the data rather than from
## the print, kept here because the published panel is hard to read:
##   * the core-side label on the lung panel is CDX2, not COX2 -- PTGS2/COX2 is
##     not on the GeoMx WTA panel at all, while CDX2 is, and is strongly
##     core-enriched (FDR 1.7e-6 in the colon core/margin fit).
##   * the core-side label is TGFB-I, not TGFB-1. Both are on the panel, but
##     TGFBI is the core-enriched one (-0.65) while TGFB1 goes the other way
##     (+0.90), and the label sits on the core side.
##   * "IF123" is IFI27 (-0.68, core).
## One label on the liver panel, immediately right of PTPRC and partly hidden
## under a red ring, is UNREAD: a 5-character symbol, with MS4A1/CD20 ruled out
## because it is not on the panel. Read it off the source figure and add it to
## `TIL_T - TB in Liver` above.

## One volcano. Genes are plotted against their adjusted p-value; the horizontal
## dash is the FDR < 0.05 threshold and the vertical dashes the log2
## fold-change cutoff, FC_CUTOFF.
##
## label_genes  a character vector of genes to label; defaults to the published
##              list for this contrast. Pass NULL for Wenbin's original rule --
##              the 15 most significant genes on each side that reach
##              FDR < 0.001.
## title  off by default. The contrast name ("Liver - Colon") reads backwards
##        against the axis, where the LEFT side is higher in colon, so printing
##        it over the panel invites exactly that misreading. The x-axis label
##        already names both directions. Pass title = contrast to bring it back.
## flip   negate Estimate before plotting. The core/margin sheets store
##        TIL_TB - T (positive = higher at the margin), while the published
##        panels put the margin on the left. 03_figure3.R fits TIL_T - TB
##        directly, so it does not need this.
## x_lab  override the automatic "Higher in A <- log2(FC) -> Higher in B".
volcano_plot <- function(res, contrast, title = NULL,
                         label_genes = VOLCANO_LABELS[[contrast]],
                         flip = FALSE, x_lab = NULL) {
  d <- res[res$Contrast == contrast, ]
  stopifnot(nrow(d) > 0)

  if (flip) d$Estimate <- -d$Estimate
  d$neglog10_FDR <- -log10(d$FDR)

  if (is.null(label_genes)) {
    top_g <- c(d$Gene[order(d$invert_P, decreasing = TRUE)][1:15],
               d$Gene[order(d$invert_P, decreasing = FALSE)][1:15])
    labelled <- filter(d, Gene %in% top_g & FDR < 0.001)
  } else {
    absent <- setdiff(label_genes, d$Gene)
    if (length(absent))
      warning("not in the '", contrast, "' results, so not labelled: ",
              paste(absent, collapse = ", "), call. = FALSE)
    labelled <- filter(d, Gene %in% label_genes)
  }

  if (is.null(x_lab)) {
    sides <- c(sub(".*\\s", "", contrast), sub("\\s-.*", "", contrast))
    if (flip) sides <- rev(sides)
    x_lab <- paste("Higher in", sides[1], "<- log2(FC) ->", "Higher in", sides[2])
  }

  brk <- volcano_breaks(d$neglog10_FDR)

  p <- ggplot(d, aes(x = Estimate, y = neglog10_FDR, label = Gene)) +
    geom_vline(xintercept = c(FC_CUTOFF, -FC_CUTOFF), lty = "dashed",
               colour = "grey40") +
    geom_hline(yintercept = -log10(0.05), lty = "dashed", colour = "grey40") +
    geom_point(shape = 21, colour = VOLCANO_BLUE, fill = NA,
               stroke = 0.35, size = 1.5) +
    ## Anchor at 0 so the 10^0 tick is drawn, as on the published panels.
    scale_y_continuous(breaks = brk, labels = volcano_labels(brk),
                       limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    labs(x = x_lab, y = expression(-log[10](italic(p) * "-adjusted"))) +
    geom_text_repel(data = labelled, size = 3, point.padding = 0.15,
                    color = "black", min.segment.length = .1,
                    box.padding = .2, max.overlaps = 50) +
    theme_bw(base_size = 11) +
    ## Everything that carries text is forced to black; theme_bw draws axis
    ## text in grey30, which prints faint.
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          text = element_text(colour = "black"),
          axis.text = element_text(colour = "black"),
          axis.title = element_text(colour = "black"),
          plot.title = element_text(face = "bold", hjust = .5,
                                    colour = "black"))

  if (!is.null(title)) p <- p + labs(title = title)
  p
}

## fgsea against MSigDB Hallmark, ranking genes by signed -log(p).
##
## BOTH HALVES OF FIGURE 1 ARE THE MIXED MODEL. Figure 1d always was; Figure
## 1c is now too, so the same load_de_results() table feeds both and the panel
## pair is one analysis shown two ways.
##
## The published Fig 1d dot plots were checked set by set against page 26 of
## the accepted manuscript PDF: Liver-Colon prints six sets (Allograft
## rejection, Myc targets v1, EMT, Unfolded protein response, Complement,
## KRAS signaling dn) and Lung-Colon five (Myc targets v1, EMT, Myogenesis,
## Bile acid metabolism, KRAS signaling dn). Ranking on load_de_results()
## reproduces both lists exactly, in the printed order; ranking on
## welch_de_results() does not -- it gives seven and three. Never pass the
## t-test results here, whatever CRC_FIG1C_METHOD is set to.
##
## nPermSimple is raised from the fgsea default of 1000 because two Hallmark
## sets on the Lung - Colon ranks have unbalanced gene-level statistics and
## come back with pval/NES = NA below 50000. At 50000 no pathway is dropped on
## either contrast, with either differential-expression path.
hallmark_gsea <- function(res, contrast, seed = 54321, nPermSimple = 50000) {
  load(require_input("data_volcanoplot_pathway", "pathwayInfo.RDATA"))  # -> pathway_info
  d <- res[res$Contrast == contrast, ]
  stopifnot(nrow(d) > 0)

  ranks <- -log(d$`Pr(>|t|)`) * sign(d$Estimate)
  names(ranks) <- d$Gene
  ranks <- ranks[order(names(ranks), decreasing = TRUE)]

  set.seed(seed)
  g <- fgsea(pathways = pathway_info[["HALLMARK"]], stats = ranks,
             eps = 0.0, minSize = 15, maxSize = 500, nPermSimple = nPermSimple)
  if (any(is.na(g$pval)))
    warning(sum(is.na(g$pval)), " Hallmark set(s) returned NA on '", contrast,
            "'; raise nPermSimple.", call. = FALSE)
  g
}

## The Figure 1d dot plot.
##
## SELECTION RULE -- every Hallmark set that reaches the NOMINAL p < 0.05, in
## either direction, ordered by NES. That is six sets on Liver - Colon and
## five on Lung - Colon, which is the published panel exactly. Stating it as a
## threshold rather than as a top-N is deliberate: the panel shows everything
## that clears nominal significance and nothing else, so there is no arbitrary
## "why six" to defend. `max_per_side` is a safety cap for other contrasts and
## does not bind on either Figure 1 contrast.
##
## The sets that also survive Benjamini-Hochberg across the 50 Hallmark sets
## are starred at the upper right of their dot: Allograft rejection and Myc
## targets v1 on Liver - Colon, Myc targets v1 on Lung - Colon. The dot colour
## stays on the nominal p, so the asterisk is the only place the adjusted p
## reaches the panel; the caption states both thresholds.
##
## title  off by default. The contrast name reads backwards against the x
##        axis, where the LEFT side is higher in the second organ, so the
##        subtitle names both directions instead. Pass title = contrast to
##        bring it back.
gsea_dotplot <- function(gsea, contrast, title = NULL,
                         p_cutoff = 0.05, padj_cutoff = 0.05,
                         max_per_side = 10) {
  parts <- strsplit(contrast, "-")[[1]]
  subtitle <- paste("higher in", trimws(parts[2]),
                    "<-> higher in", trimws(parts[1]))

  d <- rbind(gsea %>% filter(NES > 0, pval < p_cutoff) %>% arrange(NES) %>%
               tail(max_per_side),
             gsea %>% filter(NES < 0, pval < p_cutoff) %>% arrange(NES) %>%
               head(max_per_side)) %>%
    arrange(NES)
  if (nrow(d) == 0)
    stop("no Hallmark set reaches p < ", p_cutoff, " on '", contrast, "'",
         call. = FALSE)
  d$index <- as.factor(seq_len(nrow(d)))
  starred <- d[d$padj < padj_cutoff, ]

  ggplot(d, aes(x = NES, y = index, color = -log10(pval))) +
    geom_point(alpha = .8, size = 5) +
    geom_text(data = starred, aes(x = NES, y = index), label = "*",
              nudge_x = .28, nudge_y = .13, size = 7, color = "black",
              inherit.aes = FALSE, show.legend = FALSE) +
    scale_y_discrete(labels = d$pathway) +
    scale_color_gradient2(low = "lightgray", mid = "pink", high = "red2") +
    xlim(c(-3.5, 3.5)) +
    geom_vline(aes(xintercept = 0), lty = "dashed") +
    labs(title = title, subtitle = subtitle, x = "NES", y = "",
         caption = sprintf("all sets shown reach nominal p < %s;  * FDR < %s",
                           format(p_cutoff), format(padj_cutoff))) +
    theme_bw() +
    ## Every piece of text is black -- pathway names, axis, legend, caption.
    ## theme_bw's default grey30 on the axis text prints faint next to the
    ## pale end of the colour scale.
    theme(text = element_text(colour = "black"),
          plot.title = element_text(hjust = .5, size = 12, face = "bold",
                                    colour = "black"),
          plot.subtitle = element_text(hjust = .5, size = 11, colour = "black"),
          plot.caption = element_text(hjust = 0, size = 9, colour = "black"),
          axis.text = element_text(size = 11, colour = "black"),
          axis.title = element_text(size = 12, colour = "black"),
          legend.text = element_text(size = 11, colour = "black"),
          legend.title = element_text(size = 12, face = "bold",
                                      colour = "black"))
}

## Which Hallmark sets the panel above draws, as a table: the full fgsea
## result with the two flags the panel encodes, so the selection can be
## checked without re-reading the figure.
gsea_table <- function(gsea, p_cutoff = 0.05, padj_cutoff = 0.05) {
  gsea %>%
    as.data.frame() %>%
    transmute(pathway, size, NES, ES,
              pval, padj,
              plotted   = pval < p_cutoff,
              starred   = padj < padj_cutoff) %>%
    arrange(pval)
}
