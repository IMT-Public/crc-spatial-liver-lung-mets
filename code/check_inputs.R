## ---------------------------------------------------------------------------
## check_inputs.R -- preflight for the data directory.
##
## Verifies that every file the pipeline reads is present AND has the structure
## the code assumes (object names inside .RDATA, sheet names and columns in the
## spreadsheets, row-order agreement between the CIBERSORTx table and the GeoMx
## object). Run this before code/run.
##
##   cd code && Rscript check_inputs.R
##   CRC_DATA_DIR=/scratch/.../data Rscript check_inputs.R
## ---------------------------------------------------------------------------

source("lib/paths.R")

PASS <- 0L; FAIL <- 0L; WARN <- 0L

ok   <- function(msg) { PASS <<- PASS + 1L; cat("  [ ok ] ", msg, "\n", sep = "") }
bad  <- function(msg) { FAIL <<- FAIL + 1L; cat("  [FAIL] ", msg, "\n", sep = "") }
warn <- function(msg) { WARN <<- WARN + 1L; cat("  [warn] ", msg, "\n", sep = "") }

human <- function(b) {
  u <- c("B", "KB", "MB", "GB"); i <- 1
  while (b >= 1024 && i < length(u)) { b <- b / 1024; i <- i + 1 }
  sprintf("%.1f %s", b, u[i])
}

## Returns the path if present, else NULL (and records the failure).
have <- function(..., required = TRUE) {
  p <- data_path(...)
  rel <- file.path(...)
  if (file.exists(p)) {
    ok(sprintf("%-52s %s", rel, human(file.info(p)$size)))
    invisible(p)
  } else {
    if (required) bad(sprintf("%-52s MISSING", rel)) else warn(sprintf("%-52s absent (optional)", rel))
    invisible(NULL)
  }
}

## Load a .RDATA into a private env and report which objects it carries.
peek <- function(path) {
  e <- new.env()
  tryCatch({ load(path, envir = e); e }, error = function(err) {
    bad(paste("could not load:", conditionMessage(err))); NULL
  })
}

cat("\n=========================================================\n")
cat(" Input check for DATA_DIR = ", normalizePath(DATA_DIR, mustWork = FALSE), "\n", sep = "")
cat("=========================================================\n")

## --- 1. GeoMx object -------------------------------------------------------
cat("\n[1/4] data_volcanoplot_pathway/  (every figure)\n")

p <- have("data_volcanoplot_pathway", "targets.new-219-ROIs-newHKthenQ3.RDATA")
targets.new <- NULL
if (!is.null(p)) {
  e <- peek(p)
  if (!is.null(e)) {
    if ("targets.new" %in% ls(e)) {
      targets.new <- get("targets.new", envir = e)
      nf <- length(Biobase::featureNames(targets.new))
      ns <- length(Biobase::sampleNames(targets.new))
      ok(sprintf("  object 'targets.new' present: %d genes x %d ROIs", nf, ns))
      assays <- tryCatch(Biobase::assayDataElementNames(targets.new), error = function(z) character())
      for (a in c("log_q", "q_norm")) {
        if (a %in% assays) ok(sprintf("  assay '%s' present", a))
        else bad(sprintf("  assay '%s' MISSING (needed by every figure)", a))
      }
      pdcols <- tryCatch(colnames(Biobase::pData(targets.new)), error = function(z) character())
      need <- c("patient", "Area", "Site2", "GeneDetectionRate", "slide", "Sample_ID")
      miss <- setdiff(need, pdcols)
      if (!length(miss)) ok("  pData has patient/Area/Site2/GeneDetectionRate/slide/Sample_ID")
      else bad(paste("  pData missing columns:", paste(miss, collapse = ", ")))
      ## PHI sweep -- nothing here should look like an identifier or a date.
      susp <- grep("mrn|name|dob|birth|ssn|date", pdcols, ignore.case = TRUE, value = TRUE)
      if (length(susp)) warn(paste("  review before publishing, pData columns:",
                                   paste(susp, collapse = ", ")))
      else ok("  no obviously identifying pData column names")
      ## Slides must carry arbitrary labels, never pathology accession numbers.
      acc <- grepl("S[0-9]{2}-[0-9]{4,}", as.character(Biobase::pData(targets.new)$slide))
      if (any(acc)) bad(sprintf("  pData$slide holds %d accession-like values; use the recoded object",
                                sum(acc)))
      else ok("  pData$slide carries no pathology accession numbers")
    } else {
      bad(paste("  expected object 'targets.new'; file contains:",
                paste(ls(e), collapse = ", ")))
    }
  }
}

p <- have("data_volcanoplot_pathway", "pathwayInfo.RDATA")
if (!is.null(p)) {
  e <- peek(p)
  if (!is.null(e)) {
    if ("pathway_info" %in% ls(e)) {
      pi <- get("pathway_info", envir = e)
      if ("HALLMARK" %in% names(pi))
        ok(sprintf("  pathway_info$HALLMARK present (%d gene sets)", length(pi$HALLMARK)))
      else bad("  pathway_info has no 'HALLMARK' element (Fig 1d GSEA needs it)")
    } else bad(paste("  expected 'pathway_info'; file contains:", paste(ls(e), collapse = ", ")))
  }
}

## 01_figure1.R reads this for Figure 1c (the default, CRC_FIG1C_METHOD=lme)
## and always for Figure 1d.
p <- have("data_volcanoplot_pathway", "DE_res_combinedROIs_lme.xlsx")
if (!is.null(p) && requireNamespace("readxl", quietly = TRUE)) {
  sh <- readxl::excel_sheets(p)
  if (length(sh) >= 4) ok(sprintf("  %d sheets: %s", length(sh), paste(sh, collapse = ", ")))
  else bad(sprintf("  Figures 1c-d expect 4 contrast sheets but found %d", length(sh)))
  d1 <- suppressMessages(readxl::read_xlsx(p, 1, n_max = 5))
  need <- c("Gene", "Contrast", "Estimate", "Pr(>|t|)", "FDR")
  miss <- setdiff(need, names(d1))
  if (!length(miss)) ok("  columns Gene/Contrast/Estimate/Pr(>|t|)/FDR present")
  else bad(paste("  missing columns:", paste(miss, collapse = ", ")))
}

## --- 2. GSVA score inputs --------------------------------------------------
cat("\n[2/4] data_boxplots/  (02, 03, 04, 05)\n")

p <- have("data_boxplots", "GSVAscores_updated_03272025.csv")
if (!is.null(p)) {
  g <- utils::read.csv(p, nrows = 5)
  need <- c("ROI", "patient", "Site", "Site2", "Area")
  miss <- setdiff(need, names(g))
  if (!length(miss)) ok("  annotation columns ROI/patient/Site/Site2/Area present")
  else bad(paste("  missing annotation columns:", paste(miss, collapse = ", ")))
  cells <- c("mCAF", "tCAF", "Mph.PLTP", "Mph.SPP1", "B.MKI67", "CD8.Tem.GZMK", "aHSC")
  miss <- setdiff(cells, names(g))
  if (!length(miss)) ok("  all 7 signature columns present")
  else bad(paste("  missing signature columns:", paste(miss, collapse = ", ")))
}

p <- have("data_boxplots", "CD4TnCCR7_top50genes.RDATA")
if (!is.null(p)) {
  e <- peek(p)
  if (!is.null(e)) {
    if ("top.genes" %in% ls(e))
      ok(sprintf("  'top.genes' present (%d genes)", length(get("top.genes", envir = e))))
    else bad(paste("  expected 'top.genes'; file contains:", paste(ls(e), collapse = ", ")))
  }
}

## --- 3. CIBERSORTx ---------------------------------------------------------
cat("\n[3/4] data_cibersortx/  (01_figure1.R -> Fig 1b)\n")

p <- have("data_cibersortx", "CIBERSORTx_Job14_Results_02152024.csv")
if (!is.null(p)) {
  ct <- utils::read.csv(p, as.is = TRUE)
  if ("Mixture" %in% names(ct)) ok(sprintf("  'Mixture' column present, %d rows", nrow(ct)))
  else bad("  no 'Mixture' column")
  need <- c("B.naive", "B.memory", "T.CD4.naive", "T.CD8.naive", "macrophages",
            "neutrophils", "fibroblasts", "endothelial.cells", "Treg", "plasma", "mast", "NK")
  miss <- setdiff(need, names(ct))
  if (!length(miss)) ok("  all cell-type columns 01_figure1.R merges are present")
  else bad(paste("  missing cell-type columns:", paste(miss, collapse = ", ")))
  ## Script 03 asserts row-order agreement; catch a mismatch here, not there.
  if (!is.null(targets.new) && "Mixture" %in% names(ct)) {
    rn <- rownames(Biobase::pData(targets.new))
    if (length(rn) == nrow(ct) && all(rn == ct$Mixture))
      ok("  row order matches pData(targets.new) exactly")
    else
      bad(sprintf("  ROW ORDER MISMATCH vs targets.new (%d vs %d rows) - 01_figure1.R asserts this and will stop",
                  length(rn), nrow(ct)))
  }
}

## --- 4. TCGA ---------------------------------------------------------------
cat("\n[4/4] data_tcga/  (04_figure4.R -> Fig 4d)\n")

p <- have("data_tcga", "some_cell_subtype_markers.xlsx")
if (!is.null(p) && requireNamespace("readxl", quietly = TRUE)) {
  sh <- readxl::excel_sheets(p)
  ok(sprintf("  %d sheets: %s", length(sh), paste(sh, collapse = ", ")))
  ## 04_figure4.R selects the aHSC sheet by name, so sheet order does not
  ## matter -- but the name has to exist.
  if ("aHSC" %in% sh) ok("  sheet 'aHSC' present -> 04_figure4.R can build Fig 4d")
  else bad("  no sheet named 'aHSC'; 04_figure4.R selects it by name and will stop")
}

p <- have("data_tcga", "data_clinical_patient.txt")
if (!is.null(p)) {
  cl <- utils::read.delim(p, skip = 4, nrows = 5)
  need <- c("PFS_MONTHS", "PFS_STATUS")
  miss <- setdiff(need, names(cl))
  if (!length(miss)) ok("  PFS_MONTHS / PFS_STATUS present after skip=4")
  else bad(paste("  missing after skip=4 (header offset wrong?):", paste(miss, collapse = ", ")))
}

p <- have("data_tcga", "data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt")
if (!is.null(p)) {
  hdr <- utils::read.delim(p, nrows = 2)
  ok(sprintf("  %d columns (2 annotation + ~%d samples)", ncol(hdr), ncol(hdr) - 2))
  ahsc <- c("COL3A1", "COL1A1", "COL1A2", "ACTA2", "TAGLN")
  genes <- utils::read.delim(p, colClasses = c("character", rep("NULL", ncol(hdr) - 1)))[[1]]
  miss <- setdiff(ahsc, genes)
  if (!length(miss)) ok("  all 5 aHSC signature genes present in the expression matrix")
  else bad(paste("  aHSC genes absent from the matrix (Fig 4d will be wrong):",
                 paste(miss, collapse = ", ")))
}

## --- summary ---------------------------------------------------------------
cat("\n=========================================================\n")
cat(sprintf(" %d passed, %d failed, %d warnings\n", PASS, FAIL, WARN))
if (FAIL > 0) {
  cat(" NOT ready to run. Resolve the [FAIL] items above.\n")
  cat("=========================================================\n")
  quit(status = 1)
}
cat(" All required inputs present and structurally valid.\n")
cat(" Next: ./run\n")
cat("=========================================================\n")
