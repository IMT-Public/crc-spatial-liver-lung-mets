## ---------------------------------------------------------------------------
## lib/models.R -- the mixed-model families behind every box plot
##
## READ THIS BEFORE SPLITTING ANYTHING FURTHER.
##
## Each function below fits ONE family: a single limma fit with
## duplicateCorrelation blocking on patient, all pairwise contrasts between the
## group levels, and ONE Benjamini-Hochberg adjustment (`adj.P.Val.global`)
## across the whole family. Two of these families cross figure boundaries:
##
##   site_gsva_family()    Fig 2a-d   AND  Fig 4a   (aHSC is cell type 9)
##   core_margin_family()  Fig 3b     AND  Supp Fig 1b-c
##
## The p-values printed on Fig 4a and Supp Fig 1b-c are therefore adjusted
## together with Fig 2 and Fig 3b respectively. Re-fitting a figure on its own
## and adjusting only its own rows gives DIFFERENT numbers from the published
## ones. That is why each family is fit whole here and the figure scripts only
## select rows from the result -- never re-adjust.
##
## `adj.P.Val` (per contrast, across signatures) also comes back from limma and
## is kept in the tables, but the panels use `adj.P.Val.global`.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(limma)
  library(dplyr)
  library(stringr)
})

## The 9 signatures, in the order the manuscript uses them. aHSC is scored
## alongside the Figure 2 cell types because it belongs to that BH family.
CELL_TYPES <- c("mCAF", "tCAF", "Mph.PLTP", "Mph.SPP1", "B.MKI67",
                "Plasma.cells", "CD4.Tn.CCR7", "CD8.Tem.GZMK", "aHSC")

CELL_TYPE_LABELS <- c("mCAF", "tCAF", "MRC1+CSF1R+ mac", "SPP1+ mac",
                      "Proliferating B cells", "Plasma cells", "CD4+ Tcm",
                      "CD8+ Tem", "a-HSC")

## Supplementary Figure 1d group order, as published: the primary and the two
## metastases first, then the adjacent-normal tissue of each metastasis.
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

## The fit shared by all four families.
##   values    data frame, one column per signature (or gene), one row per sample
##   group     factor giving the comparison groups, one per row of `values`
##   block     patient id, one per row of `values`
##   label     name for the row-identifier column ("CellType" or "Gene")
##   moderate  empirical-Bayes variance moderation across the columns of
##             `values`. TRUE for the signature families, which have eight or
##             nine columns to borrow across. FALSE for Figure 4b-c, which has
##             TWO genes: moderating across two columns borrows almost nothing
##             and is what separates our numbers from the published Source
##             Data, which used the ordinary t. See the note on
##             site_expression_family().
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
    ## Ordinary t on the same generalised-least-squares fit: the residual
    ## standard deviation of each column on its own, no prior. BH is applied
    ## within each contrast, exactly as topTable() does it.
    tstat <- fit2$coefficients / (fit2$stdev.unscaled * fit2$sigma)
    pval  <- 2 * pt(-abs(tstat), df = fit2$df.residual)
    per.contrast <- lapply(seq_along(contr.names), function(i) {
      data.frame(rownames(pval), contr.names[i], pval[, i],
                 p.adjust(pval[, i], "BH"), row.names = NULL) %>%
        setNames(c(label, "Contrast", "P.Value", "adj.P.Val"))
    })
  }

  ## One BH adjustment over the whole family -- see the header.
  do.call(rbind, per.contrast) %>%
    mutate(!!label := factor(.data[[label]], levels = colnames(values)),
           Contrast = factor(Contrast, levels = contr.names),
           adj.P.Val.global = p.adjust(P.Value, "BH")) %>%
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
## Feeds Fig 4b-c only.
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

  ## moderate = FALSE: two genes is not a family to borrow variance across,
  ## and the published Source Data (Fig4abc.xlsx, sheet Fig4bc.logExpression)
  ## is the ordinary t -- reproduced here to the last digit. Fig 4a, with nine
  ## cell types, is moderated and matches its own Source Data as well.
  list(scores = expr,
       contrasts = fit_pairwise_contrasts(expr[, genes], expr$Site,
                                          expr$patient, label = "Gene",
                                          moderate = FALSE))
}

## Supplementary Figure 1d shows CXCL12 alongside the four signature panels, so
## the gene expression has to be summarised over the same five groups. Only the
## display is new: the brackets on the published panel are the Colon/Liver/Lung
## p-values from site_expression_family() above (the same ones on Fig 4b), so
## nothing is refit over the AN groups here.
site_area_expression <- function(pd, d, genes = "CXCL12") {
  stopifnot(all(genes %in% rownames(d)))

  expr <- data.frame(pd[c("Site2", "patient", "Area")],
                     t(d[genes, , drop = FALSE]), check.names = FALSE) %>%
    mutate(Site_area = case_when(Area == "AN"      ~ paste0(Site2, "_AN"),
                                 Site2 == "Colon"  ~ "Colon_prim",
                                 TRUE              ~ paste0(Site2, "_met")),
           combinedFactor = paste(patient, Site_area, sep = "|")) %>%
    group_by(combinedFactor) %>%
    summarize(across(all_of(genes), \(x) mean(x, na.rm = TRUE)), .groups = "drop")

  parts <- do.call(rbind, strsplit(expr$combinedFactor, "|", fixed = TRUE))
  data.frame(patient   = parts[, 1],
             Site_area = factor(parts[, 2], levels = AN_LEVELS),
             expr[, genes, drop = FALSE], check.names = FALSE)
}

## --- Family D : metastasis vs its own adjacent-normal tissue -----------------
## Adjacent normal (AN) exists for liver and lung only. Requested by a reviewer;
## feeds Supp Fig 1d. Two groupings, each its own BH family:
##   by_site   Colon_prim / Liver_AN / Liver_met / Lung_AN / Lung_met
##   by_area   the same, split further into core ("In") and margin ("Ex")
adjacent_normal_family <- function(gs.avg) {
  by_site_scores <- gs.avg %>%
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

  by_area_cells  <- CELL_TYPES[c(1:4, 9)]
  by_area_scores <- by_site_scores %>%
    mutate(Area2 = case_when(Area %in% c("T", "TIL") ~ "In",
                             Area == "TB"            ~ "Ex",
                             Area == "AN"            ~ "AN"),
           Site_area = paste(Site_area, Area2, sep = "_"),
           Site_area = factor(sub("_AN_AN$", "_AN", Site_area),
                              levels = c("Colon_prim_In", "Colon_prim_Ex",
                                         "Liver_AN", "Liver_met_In", "Liver_met_Ex",
                                         "Lung_AN", "Lung_met_In", "Lung_met_Ex"))) %>%
    select(patient, Site_area, all_of(by_area_cells))

  list(
    by_site = list(
      scores    = by_site_scores,
      contrasts = fit_pairwise_contrasts(by_site_scores[, CELL_TYPES],
                                         by_site_scores$Site_area,
                                         by_site_scores$patient)),
    by_area = list(
      scores    = by_area_scores,
      contrasts = fit_pairwise_contrasts(by_area_scores[, by_area_cells],
                                         by_area_scores$Site_area,
                                         by_area_scores$patient))
  )
}
