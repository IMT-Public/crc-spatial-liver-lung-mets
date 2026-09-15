## ---------------------------------------------------------------------------
## lib/panels.R -- the Prism-like box-plot panels
##
## Figures 2a-d, 3b, 4a-c and Supplementary Figures 1b-d were published from
## GraphPad Prism 10.6.1, drawn from the score tables and adjusted p-values in
## lib/models.R. These builders redraw the same panels in R from the same
## numbers, so every published panel is reproducible from code. The published
## panels are the Prism versions; these are visually equivalent, not
## pixel-identical.
##
## prism_panel(df, contrasts, title, ylab) is the single panel builder; every
## figure script calls it and then save_panel() to write into its own folder.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggpubr)
  library(patchwork)
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

## --- Published palettes ------------------------------------------------------
## Sampled straight out of the manuscript PDF. The paper does NOT use one
## palette throughout -- each figure has its own, and Fig 3b / Supp Fig 1d draw
## some groups as open boxes (coloured outline, white fill):
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

## Published y scales. The cell-type-score panels are pinned to a common range
## so panels can be compared by eye; the two expression panels are not scores
## and keep their own.
YLIM_SCORE   <- c(-1, 1.5)
YLIM_AHSC    <- c(-2, 2)
YLIM_CXCL12  <- c(0, 8)
YLIM_CXCR4   <- c(0, 10)
## The published count panels tick every 2, not on R's default pretty scale.
BREAKS_CXCL12 <- seq(0, 8, 2)
BREAKS_CXCR4  <- seq(0, 10, 2)
YLAB_SCORE   <- "cell type score"
YLAB_COUNT   <- "Log2 Count"

## Prism significance convention.
p_stars <- function(p) {
  ifelse(is.na(p), "",
  ifelse(p < 0.0001, "****",
  ifelse(p < 0.001,  "***",
  ifelse(p < 0.01,   "**",
  ifelse(p < 0.05,   "*", "ns")))))
}

## --- Core panel builder -----------------------------------------------------
## df        long data: one row per plotted observation, columns `group`, `value`
## contrasts data.frame with group1, group2, p  (already multiplicity-adjusted)
## Only significant brackets are drawn, as in the published panels; set
## show_ns = TRUE to draw "ns" as well.
##
## `shape` picks which published style to draw. The paper is not uniform here:
##   "box"     box + whiskers with the individual observations jittered over it
##             (Fig 2a-d)
##   "boxonly" box + whiskers, no points (Fig 3b, Supp Fig 1b-d)
##   "violin"  violin, no box and no points (Fig 4a-c)
##
## NOTE the distinction between what is plotted and what is tested. `df` is the
## distribution the panel draws; `contrasts` carries p-values from the model in
## lib/models.R, which is always fit on the patient x site x area means. Figure
## 2 draws per-ROI scores over p-values fit on the means -- that is what the
## published Prism panels do, not an oversight.
## border      outline colour per group; defaults to black, which is what a
##             filled Prism box uses. Pass a named vector to draw open boxes.
## point_cols  jitter colour per group (Fig 2); NULL keeps the white/black dot.
## ylim        fixed y range. Significance brackets are then laid out inside the
##             headroom between the data and ylim[2] rather than pushing the
##             axis upwards, so every panel keeps the published scale.
## zero_line   the dashed line at a cell-type score of 0 that Figs 2, 3b and
##             Supp Fig 1d all carry (and their legends call out).
prism_panel <- function(df, contrasts, title, ylab,
                        cols = site_cols, border = NULL, point_cols = NULL,
                        show_ns = FALSE, ylim = NULL, ybreaks = waiver(),
                        zero_line = FALSE,
                        shape = c("box", "boxonly", "violin"),
                        bracket_order = c("p", "span"),
                        bracket_overflow = FALSE) {
  shape <- match.arg(shape)
  bracket_order <- match.arg(bracket_order)
  if (is.null(border))
    border <- setNames(rep("black", length(cols)), names(cols))

  df <- df %>% filter(!is.na(value))
  rng  <- range(df$value, na.rm = TRUE)
  span <- diff(rng)
  if (!is.finite(span) || span == 0) span <- 1

  ## Bracket stacking order. "p" (the default, and what results/ has always
  ## used) puts the smallest p at the bottom. "span" puts the NARROWEST
  ## bracket at the bottom, which is what Prism does and what the published
  ## Figure 3b shows.
  glev <- levels(factor(df$group))
  sig <- contrasts %>%
    mutate(label = p_stars(p)) %>%
    filter(if (show_ns) TRUE else label != "ns")
  sig <- if (bracket_order == "span") {
    sig %>% arrange(abs(match(group2, glev) - match(group1, glev)),
                    pmin(match(group1, glev), match(group2, glev)))
  } else {
    sig %>% arrange(p)
  }

  if (is.null(ylim)) {
    if (nrow(sig) > 0) {
      sig$y.position <- rng[2] + span * (0.08 * seq_len(nrow(sig)))
      top <- max(sig$y.position) + span * 0.06
    } else {
      top <- rng[2] + span * 0.05
    }
    y.scale <- scale_y_continuous(expand = expansion(mult = c(0.06, 0.02)),
                                  limits = c(rng[1] - span * 0.06, top))
  } else {
    ## Fixed scale: fit the brackets into whatever headroom is left above the
    ## data, so they never push past the published axis.
    if (nrow(sig) > 0) {
      if (bracket_overflow) {
        ## Published Prism panels keep the axis at its stated range and stack
        ## the brackets ABOVE it, outside the plotting box. Squeezing five
        ## brackets into the headroom below ylim[2] instead (the default) runs
        ## them into the data. Needs clip = "off" and top margin, added below.
        sig$y.position <- ylim[2] + diff(ylim) * 0.062 * seq_len(nrow(sig))
      } else {
        base <- min(rng[2], ylim[2] - diff(ylim) * 0.04)
        step <- (ylim[2] - base) / (nrow(sig) + 1)
        sig$y.position <- base + step * seq_len(nrow(sig))
      }
    }
    ## Drop the scale limits ONLY when brackets are actually being placed
    ## above them -- otherwise a panel with no significant contrast would
    ## free-scale and lose its published axis.
    y.scale <- if (bracket_overflow && nrow(sig) > 0)
      scale_y_continuous(breaks = ybreaks, expand = expansion(mult = 0))
    else
      scale_y_continuous(limits = ylim, breaks = ybreaks,
                         expand = expansion(mult = 0))
  }

  p <- ggplot(df, aes(x = group, y = value))

  if (zero_line)
    p <- p + geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.4,
                        colour = "black")

  if (shape == "violin") {
    ## No transformation of the scores -- a violin is a kernel density of the
    ## same numbers the box plots show.
    ##
    ## Three rendering choices, all measured off the printed Fig 4a-c panels
    ## at 600 dpi (nine violin outlines traced pixel by pixel and fitted):
    ##
    ##   scale = "width"  All nine printed violins are exactly 244 px wide at
    ##     their widest, so Prism scales every group to the same maximum width
    ##     rather than to equal area.
    ##   trim = FALSE     The printed violins run well past the smallest and
    ##     largest observed value, so Prism draws the whole smoothed density,
    ##     tails and all. trim = TRUE cuts it flat at the data range -- the
    ##     honest outline, but not the printed one.
    ##   bw = 4.2% of the axis  Prism's bandwidth is NOT Silverman's rule. It
    ##     tracks the plotted axis range, not the spread of the data: across
    ##     the nine violins bw/axis-range has a coefficient of variation of
    ##     0.21, against 0.44 for bw/sd, 0.44 for bw/data-range and 0.66 for
    ##     bw/bw.nrd0. That is why no single `adjust` multiplier can match all
    ##     three panels -- the adjust the print implies runs from 0.6 on Fig 4a
    ##     Liver to 3.5 on Fig 4b Liver. Setting the bandwidth in axis units
    ##     instead lands within 0.01 RMSE of the per-violin optimum on seven of
    ##     the nine. The kernel shape is irrelevant: gaussian, Epanechnikov,
    ##     triangular and biweight all fit the printed outlines equally well.
    ##
    ## A free-scaled panel (ylim = NULL) has no axis to key off, so it keeps
    ## the old multiplier.
    vio <- c(list(aes(fill = group, colour = group), width = 0.85,
                  trim = FALSE, linewidth = 0.4, scale = "width"),
             if (is.null(ylim)) list(adjust = 0.6)
             else list(bw = diff(ylim) * 0.042))
    p <- p + do.call(geom_violin, vio) +
      scale_colour_manual(values = cols, na.value = "grey70", guide = "none")
  } else {
    p <- p + geom_boxplot(aes(fill = group, colour = group), width = 0.6,
                          outlier.shape = NA, linewidth = 0.6) +
      scale_colour_manual(values = border, na.value = "black", guide = "none")
    if (shape == "box") {
      ## Wide enough that 66-69 ROIs in one column stay individually visible;
      ## at the old 0.14 they piled into a line.
      ## The dot colour is passed as a plain vector rather than mapped, because
      ## `fill` is already spoken for by the box.
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
    p <- p + ggpubr::stat_pvalue_manual(
      sig, label = "label", tip.length = 0.012,
      bracket.size = 0.4, size = 3.4, inherit.aes = FALSE)
  }

  if (!is.null(ylim) && bracket_overflow && nrow(sig) > 0) {
    ## The brackets are drawn above the panel, where the title normally sits,
    ## so the title is pushed up by the height of the stack and the plot
    ## margin grows to make room for both.
    ## The bracket stack is offset in DATA units (a fraction of the y range),
    ## so the space to reserve for it scales with the panel's height -- which
    ## differs between a single 3.2in panel and a cell of a combined grid.
    ## Expressing the margin in "npc" (a fraction of the plot) tracks that;
    ## a fixed pt lift collides with the title in the taller combined figures.
    lift <- 0.075 * nrow(sig) + 0.02
    p <- p + coord_cartesian(ylim = ylim, clip = "off") +
      theme(plot.margin = unit(c(0.02, 0.02, 0.02, 0.02), "npc"),
            ## ggplot stacks margin -> title -> panel, so the brackets (drawn
            ## above the panel with clip off) land in the title's space; the
            ## title's bottom margin is what actually clears them.
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
split_contrast <- function(tbl, contrast_col = "Contrast", p_col = "adj.P.Val.global") {
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
