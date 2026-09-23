## ---------------------------------------------------------------------------
## lib/de_volcano_gsea.R -- volcano plots and Hallmark GSEA dot plots
##
## load_de_results() reads the per-gene linear mixed-effects results in
## DE_res_combinedROIs_lme.xlsx (one sheet per contrast). Everything keys off
## the `Contrast` column, whose values are listed in DE_CONTRASTS.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(fgsea)
})

DE_CONTRASTS <- c(liver_colon = "Liver - Colon",
                  liver_lung  = "Liver - Lung",
                  lung_colon  = "Lung - Colon",
                  til_tb_t    = "TIL_TB - T")

VOLCANO_BLUE <- "#4A7EBB"

## log2 fold-change cutoff drawn on the volcano plots.
FC_CUTOFF <- 0.5

load_de_results <- function() {
  f <- require_input("data_volcanoplot_pathway", "DE_res_combinedROIs_lme.xlsx")
  do.call(rbind, lapply(seq_along(DE_CONTRASTS), \(i) read_xlsx(f, i)))
}

## The y axis carries -log10(adjusted p), labelled in powers of ten, with
## ~4 breaks whatever the contrast's dynamic range.
volcano_breaks <- function(y) {
  m <- max(y[is.finite(y)], na.rm = TRUE)
  steps <- c(1, 2, 5, 10, 20, 25, 50, 100)
  step  <- steps[which(steps >= m / 4)][1]
  seq(0, ceiling(m / step) * step, by = step)
}

volcano_labels <- function(breaks) {
  parse(text = ifelse(breaks == 0, "10^0", paste0("10^-", breaks)))
}

## Genes of interest labelled on the Figure 3c panels: tumour core (TIL + T)
## against invasive margin (TB), within lung and within liver.
VOLCANO_LABELS <- list(
  `TIL_T - TB in Lung`  = c("CD74", "HLA-DRA", "CXCL13", "CCL19", "IGHA1",
                            "IGHM", "JCHAIN", "CDX2", "S100A6", "KRT8"),
  `TIL_T - TB in Liver` = c("ALB", "APOA1", "CCL19", "IL7R", "CXCR4", "CXCL9",
                            "MMP9", "CCR7", "CD3E", "PTPRC", "HLA-DRA", "CD28",
                            "CD74", "ALDOA", "EPCAM", "KRT8", "MIF", "CEACAM1",
                            "S100A10", "IFI27", "TGFBI")
)

## One volcano. Genes are plotted against their adjusted p-value; the horizontal
## dash is the FDR < 0.05 threshold and the vertical dashes the log2
## fold-change cutoff, FC_CUTOFF.
##
## label_genes  genes to label; defaults to VOLCANO_LABELS for this contrast.
## x_lab        overrides the automatic "Higher in A <- log2(FC) -> Higher in B".
volcano_plot <- function(res, contrast,
                         label_genes = VOLCANO_LABELS[[contrast]],
                         x_lab = NULL) {
  d <- res[res$Contrast == contrast, ]
  stopifnot(nrow(d) > 0)
  d$neglog10_FDR <- -log10(d$FDR)

  absent <- setdiff(label_genes, d$Gene)
  if (length(absent))
    warning("not in the '", contrast, "' results, so not labelled: ",
            paste(absent, collapse = ", "), call. = FALSE)
  labelled <- filter(d, Gene %in% label_genes)

  if (is.null(x_lab)) {
    sides <- c(sub(".*\\s", "", contrast), sub("\\s-.*", "", contrast))
    x_lab <- paste("Higher in", sides[1], "<- log2(FC) ->", "Higher in", sides[2])
  }

  brk <- volcano_breaks(d$neglog10_FDR)

  ggplot(d, aes(x = Estimate, y = neglog10_FDR, label = Gene)) +
    geom_vline(xintercept = c(FC_CUTOFF, -FC_CUTOFF), lty = "dashed",
               colour = "grey40") +
    geom_hline(yintercept = -log10(0.05), lty = "dashed", colour = "grey40") +
    geom_point(shape = 21, colour = VOLCANO_BLUE, fill = NA,
               stroke = 0.35, size = 1.5) +
    scale_y_continuous(breaks = brk, labels = volcano_labels(brk),
                       limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    labs(x = x_lab, y = expression(-log[10](italic(p) * "-adjusted"))) +
    geom_text_repel(data = labelled, size = 3, point.padding = 0.15,
                    color = "black", min.segment.length = .1,
                    box.padding = .2, max.overlaps = 50, seed = 1) +
    theme_bw(base_size = 11) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          text = element_text(colour = "black"),
          axis.text = element_text(colour = "black"),
          axis.title = element_text(colour = "black"),
          plot.title = element_text(face = "bold", hjust = .5,
                                    colour = "black"))
}

## fgsea against MSigDB Hallmark, ranking genes by signed -log(p).
## nPermSimple is raised from the fgsea default so that no Hallmark set
## returns pval/NES = NA.
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

## The Figure 1d dot plot: every Hallmark set reaching nominal p < p_cutoff,
## ordered by NES; sets with FDR < padj_cutoff are starred.
gsea_dotplot <- function(gsea, contrast, p_cutoff = 0.05, padj_cutoff = 0.05) {
  parts <- strsplit(contrast, "-")[[1]]
  subtitle <- paste("higher in", trimws(parts[2]),
                    "<-> higher in", trimws(parts[1]))

  d <- gsea %>% filter(pval < p_cutoff) %>% arrange(NES)
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
    labs(subtitle = subtitle, x = "NES", y = "",
         caption = sprintf("all sets shown reach nominal p < %s;  * FDR < %s",
                           format(p_cutoff), format(padj_cutoff))) +
    theme_bw() +
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

## The fgsea result behind the dot plot.
gsea_table <- function(gsea) {
  gsea %>%
    as.data.frame() %>%
    transmute(pathway, size, NES, ES, pval, padj) %>%
    arrange(pval)
}
