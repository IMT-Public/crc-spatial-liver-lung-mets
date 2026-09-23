## ---------------------------------------------------------------------------
## lib/load_geomx.R -- the GeoMx object every figure starts from
##
## targets.new is a NanoStringGeoMxSet: 219 ROIs from 8 patients, counts
## normalised housekeeping-then-Q3 (assay slots `log_q` and `q_norm`).
##
## Two consumers with different needs:
##   load_geomx()        raw object + phenotype table, ROI order untouched
##                       (Fig 1b: the CIBERSORTx table is row-aligned to it)
##   geomx_expression()  log_q matrix + phenotype table with Site relabelled
##                       and TBD/TLS dropped (every model downstream)
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(GeomxTools)
  library(dplyr)
})

GEOMX_RDATA <- c("data_volcanoplot_pathway",
                 "targets.new-219-ROIs-newHKthenQ3.RDATA")

## Keep one copy of the object per R session.
.geomx_cache <- new.env(parent = emptyenv())

load_geomx <- function() {
  if (is.null(.geomx_cache$targets)) {
    load(do.call(require_input, as.list(GEOMX_RDATA)))  # -> targets.new
    .geomx_cache$targets <- targets.new
  }
  .geomx_cache$targets
}

## Normalise the free-text site labels ("Liver mets (met)", "Colon tumor").
clean_site <- function(x) {
  x |>
    stringr::str_replace("\\s", "_") |>
    stringr::str_replace_all("\\(|\\)", "") |>
    stringr::str_replace("mets", "met") |>
    stringr::str_replace("tumor", "prim")
}

LEV_AREA <- c("TIL", "TB", "T", "AN")
LEV_SITE <- c("Liver_AN", "Liver_met", "Lung_AN", "Lung_met", "Colon_prim")

## Returns list(pd, d): phenotype rows and the log_q matrix, same order,
## with the TBD and TLS areas dropped.
geomx_expression <- function() {
  targets.new <- load_geomx()

  pd <- pData(targets.new) %>%
    mutate(Site = clean_site(Site)) %>%
    filter(!Area %in% c("TBD", "TLS")) %>%
    mutate(Area = factor(Area, levels = LEV_AREA),
           Site = factor(Site, levels = LEV_SITE))

  d <- assayDataElement(targets.new, "log_q")
  colnames(d) <- sub("\\.dcc", "", colnames(d))
  d <- d[, pd$Sample_ID]
  stopifnot(all(colnames(d) == pd$Sample_ID))

  list(pd = pd, d = d)
}
