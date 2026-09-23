## ---------------------------------------------------------------------------
## lib/panels.R -- box and violin panels
##
## Figures 2a-d, 3b, 4a-c and Supplementary Figures 1b-d, drawn from the score
## tables and adjusted p-values in lib/models.R. prism_panel() builds one
## panel; save_panel() writes it to its figure folder.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggpubr)
})

## --- Prism-like theme -------------------------------------------------------
## Prism defaults: white panel, black axis lines only, ticks pointing out,
## no grid, sans-serif labels, no box around the plot.
theme_prism_like <- function(base_size = 11) {
  theme_classic(base_size = base_size, base_family = "") +
    theme(
      axis.line        = element_line(colour = "black", linewidth = 0.5),
      axis.ticks       = element_line(colour = "black", linewidth = 0.5),
      axis.ticks.length = unit(3, "pt"),
      text             = element_text(colour = "black"),
      axis.text        = element_text(colour = "black", size = base_size),
      axis.title       = element_text(colour = "black", size = base_size, face = "bold"),
      plot.title       = element_text(hjust = 0.5, face = "bold", colour = "black",
                                      size = base_size + 1),
      legend.position  = "none",
      panel.grid       = element_blank(),
      plot.margin      = margin(6, 8, 6, 8)
    )
}

## --- Palettes ----------------------------------------------------------------
## Each figure has its own palette; some groups are drawn as open boxes
## (coloured outline, white fill):
##
##   Fig 2, Supp Fig 1d   Colon black, Liver pink, Lung green; AN groups open
##   Fig 3b, Supp Fig 1bc Colon salmon, Liver blue, Lung olive; "Ex" open
##   Fig 4a-c             Colon grey, Liver orange, Lung olive (violins)
site_cols <- c(Colon = "#000000", Liver = "#FFBEBC", Lung = "#7FC77F")

## Fig 2 jitters the ROIs in the site colour, not in black.
site_point_cols <- c(Colon = "#FF9300", Liver = "#FF2600", Lung = "#018F00")

## Fig 3b and Supp Fig 1b-c: tumour core ("In") filled, invasive margin ("Ex")
## drawn as an open box in the same colour.
core_margin_base <- c(Colon = "#FF7E79", Liver = "#7A81FF", Lung = "#478A00")
core_margin_cols <- setNames(
  c(rbind(core_margin_base, "white")),
  paste0(rep(names(core_margin_base), each = 2), c("_In", "_Ex")))
core_margin_border <- setNames(
  c(rbind(core_margin_base, core_margin_base)),
  names(core_margin_cols))

## Supp Fig 1d: the metastases filled, their adjacent normal open.
an_cols   <- c(Colon_prim = "#000000", Liver_met = "#FFBEBC", Lung_met = "#7FC77F",
               Liver_AN = "white", Lung_AN = "white")
an_border <- c(Colon_prim = "#000000", Liver_met = "#FFBEBC", Lung_met = "#7FC77F",
               Liver_AN = "#FFBEBC", Lung_AN = "#7FC77F")

## Fig 4a-c violins.
violin_cols <- c(Colon = "#A9A9A9", Liver = "#FF9300", Lung = "#4F8F00")

## y scales. The cell-type-score panels share a common range; the two
## expression panels keep their own.
YLIM_SCORE   <- c(-1, 1.5)
YLIM_AHSC    <- c(-2, 2)
YLIM_CXCL12  <- c(0, 8)
YLIM_CXCR4   <- c(0, 10)
BREAKS_CXCL12 <- seq(0, 8, 2)
BREAKS_CXCR4  <- seq(0, 10, 2)
YLAB_SCORE   <- "cell type score"
YLAB_COUNT   <- "Log2 Count"

## Significance stars.
p_stars <- function(p) {
  ifelse(is.na(p), "",
  ifelse(p < 0.0001, "****",
  ifelse(p < 0.001,  "***",
  ifelse(p < 0.01,   "**",
  ifelse(p < 0.05,   "*", "ns")))))
}

## --- Core panel builder -----------------------------------------------------
## df          long data: one row per plotted observation, columns `group`, `value`
## contrasts   data.frame with group1, group2, p (already multiplicity-adjusted);
##             only significant brackets are drawn
## shape       "box"     box + whiskers with the observations jittered over it
##                       (Fig 2a-d)
##             "boxonly" box + whiskers, no points (Fig 3b, Supp Fig 1b-d)
##             "violin"  violin, no box and no points (Fig 4a-c)
## cols        fill colour per group
## border      outline colour per group; defaults to black. Pass a named vector
##             to draw open boxes.
## point_cols  jitter colour per group (Fig 2); NULL keeps the white/black dot.
## ylim        y-axis range. Significance brackets are stacked above it,
##             narrowest at the bottom.
## zero_line   dotted line at a cell-type score of 0.
prism_panel <- function(df, contrasts, title, ylab, ylim,
                        cols = site_cols, border = NULL, point_cols = NULL,
                        ybreaks = waiver(), zero_line = FALSE,
                        shape = c("box", "boxonly", "violin")) {
  shape <- match.arg(shape)
  if (is.null(border))
    border <- setNames(rep("black", length(cols)), names(cols))

  df <- df %>% filter(!is.na(value))

  glev <- levels(factor(df$group))
  sig <- contrasts %>%
    mutate(label = p_stars(p)) %>%
    filter(label != "ns") %>%
    arrange(abs(match(group2, glev) - match(group1, glev)),
            pmin(match(group1, glev), match(group2, glev)))

  if (nrow(sig) > 0) {
    sig$y.position <- ylim[2] + diff(ylim) * 0.062 * seq_len(nrow(sig))
    y.scale <- scale_y_continuous(breaks = ybreaks, expand = expansion(mult = 0))
  } else {
    y.scale <- scale_y_continuous(limits = ylim, breaks = ybreaks,
                                  expand = expansion(mult = 0))
  }

  p <- ggplot(df, aes(x = group, y = value))

  if (zero_line)
    p <- p + geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.4,
                        colour = "black")

  if (shape == "violin") {
    ## Full density (untrimmed), equal maximum width per group, bandwidth
    ## 4.2% of the y-axis range.
    p <- p + geom_violin(aes(fill = group, colour = group), width = 0.85,
                         trim = FALSE, linewidth = 0.4, scale = "width",
                         bw = diff(ylim) * 0.042) +
      scale_colour_manual(values = cols, na.value = "grey70", guide = "none")
  } else {
    p <- p + geom_boxplot(aes(fill = group, colour = group), width = 0.6,
                          outlier.shape = NA, linewidth = 0.6) +
      scale_colour_manual(values = border, na.value = "black", guide = "none")
    if (shape == "box") {
      dot_col <- if (is.null(point_cols)) "black"
                 else unname(point_cols[as.character(df$group)])
      dot_fill <- if (is.null(point_cols)) "white" else dot_col
      p <- p + geom_point(position = position_jitter(width = 0.3, height = 0,
                                                     seed = 1),
                          size = 0.9, shape = 21, stroke = 0.25,
                          colour = dot_col, fill = dot_fill)
    }
  }

  p <- p +
    scale_fill_manual(values = cols, na.value = "grey70") +
    y.scale +
    labs(title = title, x = NULL, y = ylab) +
    theme_prism_like()

  if (nrow(sig) > 0) {
    ## Brackets sit above the axis range; clip is turned off and the title
    ## is lifted by the height of the bracket stack.
    lift <- 0.075 * nrow(sig) + 0.02
    p <- p + ggpubr::stat_pvalue_manual(
      sig, label = "label", tip.length = 0.012,
      bracket.size = 0.4, size = 3.4, inherit.aes = FALSE) +
      coord_cartesian(ylim = ylim, clip = "off") +
      theme(plot.margin = unit(c(0.02, 0.02, 0.02, 0.02), "npc"),
            plot.title = element_text(margin = margin(b = lift, unit = "npc")))
  }
  p
}

## Reshape a wide score table (one column per cell type) into the long form
## prism_panel() expects.
to_long <- function(wide, group_col, value_col) {
  wide %>%
    dplyr::select(group = all_of(group_col), value = all_of(value_col)) %>%
    dplyr::mutate(group = factor(group, levels = levels(factor(wide[[group_col]]))))
}

## Split "Lung-Liver" style contrast names into the two group labels.
split_contrast <- function(tbl, contrast_col = "Contrast", p_col = "adj.P.Val") {
  parts <- strsplit(as.character(tbl[[contrast_col]]), "-", fixed = TRUE)
  data.frame(
    group1 = vapply(parts, `[`, character(1), 2),
    group2 = vapply(parts, `[`, character(1), 1),
    p      = as.numeric(tbl[[p_col]]),
    stringsAsFactors = FALSE
  )
}

## Panels are written to the folder of the figure they belong to.
save_panel <- function(p, figure, file, w = 2.6, h = 3.2) {
  ggsave(fig_path(figure, file), p, width = w, height = h, useDingbats = FALSE)
  invisible(p)
}
