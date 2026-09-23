## ---------------------------------------------------------------------------
## lib/models.R -- limma models behind the box and violin panels
##
## Each family is one limma fit with duplicateCorrelation blocking on patient,
## all pairwise contrasts between the group levels, and Benjamini-Hochberg
## adjustment within each contrast across the signatures (or genes) of the
## family (`adj.P.Val`). Two families span more than one figure:
##
##   site_gsva_family()    Fig 2a-d and Fig 4a   (aHSC is signature 9)
##   core_margin_family()  Fig 3b and Supp Fig 1b-c
##
## Each family is fit whole and the figure scripts select their rows.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(limma)
  library(dplyr)
  library(stringr)
})

## The 9 signatures, in manuscript order.
CELL_TYPES <- c("mCAF", "tCAF", "Mph.PLTP", "Mph.SPP1", "B.MKI67",
                "Plasma.cells", "CD4.Tn.CCR7", "CD8.Tem.GZMK", "aHSC")

CELL_TYPE_LABELS <- c("mCAF", "tCAF", "MRC1+CSF1R+ mac", "SPP1+ mac",
                      "Proliferating B cells", "Plasma cells", "CD4+ Tcm",
                      "CD8+ Tem", "a-HSC")

## Supplementary Figure 1d group order: the primary and the two metastases,
## then the adjacent-normal tissue of each metastasis.
AN_LEVELS <- c("Colon_prim", "Liver_met", "Lung_met", "Liver_AN", "Lung_AN")

## All pairwise contrasts between levels, later minus earlier: with levels
## c("Colon","Liver","Lung") this is Liver-Colon, Lung-Colon, Lung-Liver.
make_contr_names <- function(lev) {
  out <- NULL
  for (i in seq_along(lev))
    for (j in 2:length(lev))
      if (i < j) out <- c(out, paste0(lev[j], "-", lev[i]))
  out
}

## The fit shared by all families.
##   values    data frame, one column per signature (or gene), one row per sample
##   group     factor giving the comparison groups, one per row of `values`
##   block     patient id, one per row of `values`
##   label     name for the row-identifier column ("CellType" or "Gene")
##   moderate  empirical-Bayes variance moderation across the columns of
##             `values` (signature families). FALSE gives the ordinary t
##             (Figure 4b-c, two genes).
fit_pairwise_contrasts <- function(values, group, block, label = "CellType",
                                   moderate = TRUE) {
  group  <- factor(group)
  values <- as.matrix(values)

  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)

  ## limma wants samples in columns.
  corfit <- limma::duplicateCorrelation(t(values), design, block = block)
  fit    <- lmFit(t(values), design, block = block,
                  correlation = corfit$consensus.correlation)

  contr.names <- make_contr_names(levels(group))
  fit2 <- contrasts.fit(fit, makeContrasts(contrasts = contr.names,
                                           levels = design))

  if (moderate) {
    fit2 <- eBayes(fit2)
    per.contrast <- lapply(seq_along(contr.names), function(i) {
      topTable(fit2, coef = i, number = ncol(values), sort.by = "none") %>%
        tibble::rownames_to_column(label) %>%
        mutate(Contrast = contr.names[i]) %>%
        select(all_of(label), Contrast, P.Value, adj.P.Val)
    })
  } else {
    ## Ordinary t on the same generalised-least-squares fit; BH within each
    ## contrast, as topTable() does it.
    tstat <- fit2$coefficients / (fit2$stdev.unscaled * fit2$sigma)
    pval  <- 2 * pt(-abs(tstat), df = fit2$df.residual)
    per.contrast <- lapply(seq_along(contr.names), function(i) {
      data.frame(rownames(pval), contr.names[i], pval[, i],
                 p.adjust(pval[, i], "BH"), row.names = NULL) %>%
        setNames(c(label, "Contrast", "P.Value", "adj.P.Val"))
    })
  }

  do.call(rbind, per.contrast) %>%
    mutate(!!label := factor(.data[[label]], levels = colnames(values)),
           Contrast = factor(Contrast, levels = contr.names)) %>%
    arrange(.data[[label]], Contrast)
}

## --- Family A : signature score by site --------------------------------------
## Colon / Liver / Lung, areas pooled, adjacent-normal excluded.
## Feeds Fig 2a-d (8 cell types) and Fig 4a (aHSC).
site_gsva_family <- function(gs.avg) {
  scores <- gs.avg %>%
    mutate(Site = factor(sub("_(met|prim)$", "", Site),
                         levels = c("Colon", "Liver", "Lung"))) %>%
    filter(Area != "AN") %>%
    select(patient, Site, Area, all_of(CELL_TYPES))

  contrasts <- fit_pairwise_contrasts(
    scores[, CELL_TYPES], scores$Site, scores$patient) %>%
    mutate(Renamed = CELL_TYPE_LABELS[match(CellType, CELL_TYPES)], .after = CellType)

  list(scores = scores, contrasts = contrasts)
}

## --- Family B : signature score, tumour core vs invasive margin --------------
## TIL + T pooled as "In", TB as "Ex", within each site.
## Feeds Fig 3b (cell types 5-8) and Supp Fig 1b-c (cell types 1-4).
core_margin_family <- function(gs.avg) {
  cells <- CELL_TYPES[1:8]   # aHSC is not part of this figure family
  scores <- gs.avg %>%
    mutate(Area2 = case_when(Area %in% c("T", "TIL") ~ "In",
                             Area == "TB"            ~ "Ex",
                             TRUE                    ~ NA_character_),
           Site = factor(sub("_(met|prim)$", "", Site),
                         levels = c("Colon", "Liver", "Lung")),
           Site_area = factor(paste(Site, Area2, sep = "_"),
                              levels = paste0(rep(levels(Site), each = 2), "_",
                                              rep(c("In", "Ex"), 3)))) %>%
    filter(!is.na(Area2)) %>%
    select(patient, Site_area, all_of(cells))

  contrasts <- fit_pairwise_contrasts(
    scores[, cells], scores$Site_area, scores$patient) %>%
    mutate(Renamed = CELL_TYPE_LABELS[match(CellType, CELL_TYPES)], .after = CellType)

  list(scores = scores, contrasts = contrasts)
}

## --- Family C : CXCL12 / CXCR4 expression by site ----------------------------
## Feeds Fig 4b-c and the CXCL12 panel of Supp Fig 1d.
site_expression_family <- function(pd, d, genes = c("CXCL12", "CXCR4")) {
  expr <- data.frame(pd[c("Site2", "patient", "Area")],
                     t(d[genes, ]), check.names = FALSE) %>%
    mutate(Site = Site2) %>%
    filter(Area != "AN") %>%
    mutate(combinedFactor = paste(patient, Site, Area, sep = "_")) %>%
    group_by(combinedFactor) %>%
    summarize(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))

  expr <- data.frame(
    matrix(unlist(strsplit(expr$combinedFactor, "_")), byrow = TRUE, ncol = 3),
    expr[, -1], check.names = FALSE)
  colnames(expr)[1:3] <- c("patient", "Site", "Area")
  expr$Site <- factor(expr$Site, levels = c("Colon", "Liver", "Lung"))

  ## Two genes: ordinary t, no variance moderation.
  list(scores = expr,
       contrasts = fit_pairwise_contrasts(expr[, genes], expr$Site,
                                          expr$patient, label = "Gene",
                                          moderate = FALSE))
}

## --- Family D : metastasis vs its own adjacent-normal tissue -----------------
## Colon_prim / Liver_met / Liver_AN / Lung_met / Lung_AN, areas pooled.
## Feeds Supp Fig 1d.
adjacent_normal_family <- function(gs.avg) {
  scores <- gs.avg %>%
    select(patient, Site, Area, all_of(CELL_TYPES)) %>%
    mutate(Site_area = paste(Site, Area, sep = "_"),
           Site_area = str_replace(Site_area, "Colon", "Colon_prim"),
           Site_area = str_replace(Site_area, "Liver", "Liver_met"),
           Site_area = str_replace(Site_area, "Lung", "Lung_met"),
           Site_area = str_replace(Site_area, "met_AN", "AN"),
           Site_area = str_replace(Site_area, "Lung_met_(T|TB|TIL)$", "Lung_met"),
           Site_area = str_replace(Site_area, "Liver_met_(T|TB|TIL)$", "Liver_met"),
           Site_area = str_replace(Site_area, "Colon_prim_(T|TB|TIL)$", "Colon_prim"),
           .after = Area)

  list(scores    = scores,
       contrasts = fit_pairwise_contrasts(scores[, CELL_TYPES],
                                          scores$Site_area, scores$patient))
}
