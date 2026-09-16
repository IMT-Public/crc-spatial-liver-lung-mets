# Spatial profiling of primary CRC and synchronous lung and liver metastases

## Abstract

**Purpose:** Liver metastases from microsatellite-stable colorectal cancer (MSS
CRC) demonstrate profound resistance to immune checkpoint therapy (ICT) compared
with extrahepatic disease, but the tissue-specific mechanisms remain unclear. We
sought to determine how organ-specific microenvironments shape immune exclusion
in MSS CRC metastases.

**Experimental Design:** We performed spatial transcriptomic profiling and
multiplex immunofluorescence on synchronous liver and lung metastases, along
with matched primary CRC tumors, from patients with MSS CRC. Spatially resolved
regions of tumor core and invasive margin were analyzed for immune and stromal
composition. Clinical relevance was assessed using survival analyses in the
TCGA-COAD/READ cohort.

**Results:** Liver metastases exhibited a dense, fibrotic microenvironment
enriched with activated hepatic stellate cells (aHSCs), MRC1+/SPP1+
immunosuppressive macrophages, and elevated CXCL12–CXCR4 signaling. Spatial
analyses revealed stromal-associated effector T-cell exclusion from the tumor
core in liver metastases, with T cells restricted to the invasive margin.
Conversely, synchronous lung metastases showed sparse stroma and uniform
infiltration of proliferating B and T cells. Coordinated enrichment of aHSCs and
macrophages defined a liver-specific stromal–myeloid niche associated with
immune exclusion. High aHSC gene signature expression correlated with
significantly worse progression-free survival in the TCGA cohort.

**Conclusions:** These findings identify liver-specific stromal programming,
associated with aHSCs and CXCL12-CXCR4 signaling, as a candidate contributor of
immune exclusion in MSS CRC. Targeting liver-specific stromal–immune
interactions may represent a therapeutic strategy to enhance ICT efficacy in
liver-dominant metastatic CRC.

## Code and Data info

This archive reproduces the results of "Spatial profiling of primary CRC and
synchronous lung and liver metastases identifies a liver-specific stromal niche
associated with immune exclusion". To run everything, press **Run**. This
executes `code/run`, a shell script that runs one R script per figure. Click
the Environment icon to see the installed packages (R 4.5.1, Bioconductor
3.22). Complete legends from the current manuscript are provided below. The
code-generated panels are mapped separately under **Included Materials**.

### Complete figure captions

**Figure 1. Comparative transcriptomic and cellular landscapes of primary colon
tumors and metastatic lesions in the liver and lung**

**(a)** Workflow of the research design. Liver and lung metastasis samples from
patients who underwent synchronous resections were chosen. ROIs of tumor
(PanCK, green) and CD45 infiltration (red) were chosen for sequential RNA
expression analysis using GeoMx DSP. Archival samples of primary colorectal
tumors (*n* = 3) were also included in the analysis. Scale bar = 1 mm. Shaded
areas within the images indicate the distinct ROIs selected for downstream
analysis, highlighting representative regions across the tissue samples.
**(b)** Cell type composition in primary colon tumors (*n* = 3), liver
metastases (*n* = 8), and lung metastases (*n* = 8) using CIBERSORT. **(c)**
Volcano plot showing differentially expressed genes (DEGs) between **(top)**
colon and liver tumors and **(bottom)** colon and lung tumors. Vertical dashed
lines indicate the log2 fold-change cutoff, and the horizontal dashed line
indicates the adjusted *p*-value threshold. **(d)** Enrichment of biological
pathways associated with DEGs across metastatic sites compared to primary
tumors, illustrating site-specific molecular programs. Statistical significance
was assessed using linear mixed-effects models with patient as a random effect,
followed by Benjamini–Hochberg adjustment for multiple comparisons. Pathways
with nominal *p* < 0.05 are shown; asterisks indicate pathways that remained
significant after Benjamini–Hochberg correction (adjusted *p* < 0.05).

**Figure 2. Differential enrichment of stromal, myeloid, B cell, and T cell
subpopulations across primary colon tumors and matched liver and lung
metastases.**

**(a)** Fibroblast subtypes mCAFs and tumor-associated tCAFs. **(b)**
MRC1⁺CSF1R⁺ macrophages and SPP1⁺ macrophages exhibit site-specific enrichment,
with higher scores in liver metastases highlighting a myeloid immunosuppressive
phenotype. **(c)** Proliferating B cells and plasma cells show greater abundance
in lung metastases, suggesting differential humoral immune activation. **(d)**
Memory CD4⁺ Tcm and CD8⁺ T Tem cells are variably enriched across sites, with
increased scores in lung metastases indicating site-dependent adaptive immune
engagement. The horizontal dashed line indicates a cell-type score of 0.
**(e)** Primary CRC and its metastases to the liver and lung are shown using
multiplex immunofluorescence. Colon cancer shows the presence of fibroblasts
(cyan, white), CD8 cells (yellow), and plasma cells (orange). **(f)** Liver
metastases show the presence of fibroblasts (cyan, white), CD4 T cells (red),
and macrophages (green, blue). The inset shows higher magnification of
fibroblasts (cyan) and macrophages (green). **(g)** Lung metastases show the
presence of CD8 cells (yellow) and plasma cells (orange). The inset shows higher
magnification of proliferating B cells (pink, white arrow). Scale bar = 100 μm.
Statistical significance was assessed using linear mixed-effects models with
patient as a random effect, followed by Benjamini–Hochberg adjustment for
multiple comparisons. Adjusted *p* values are shown: \*\*\*\* *p* < 0.0001,
\*\*\* *p* < 0.001, \*\* *p* < 0.01, and \* *p* < 0.05.

**Figure 3. Spatial immune profiling reveals site- and region-specific
enrichment of immune populations in primary and metastatic colorectal cancer**

**(a)** Morphology marker-stained representative section of a liver metastasis
for digital spatial analysis. Green, PanCK; red, CD45; blue, nuclear stain.
Representative ROI polygons (tumor core and invasive margin) were chosen for
subsequent gene-expression analysis. **(b)** Cell type score comparisons across
tumor core (In) and invasive margin (Ex) regions in primary colon tumors, liver
metastases, and lung metastases. Select populations shown include proliferating
B cells, plasma cells, CD4⁺ central memory T cells (Tcm), and CD8⁺ effector
memory T cells. The horizontal dashed line indicates a cell-type score of 0.
**(c)** Volcano plot showing differentially expressed genes **(top)** in the
tumor core versus invasive margin of lung metastases and **(bottom)** in the
tumor core versus invasive margin of liver metastases. Vertical dashed lines
indicate the log2 fold-change cutoff, and the horizontal dashed line indicates
the adjusted *p*-value threshold. **(d)** **Left:** Primary CRC shows tumor cells
(purple) and immune CD4 and CD8 T cells (red, yellow), plasma cells (orange),
and macrophages (green) at the tumor margin. **Middle:** Liver metastases show
tumor cells (purple), immune T cells (red, yellow), and B cells (pink) separated
by a rim of fibroblasts (white, cyan). The inset shows higher magnification of T
cells and fibroblasts at the tumor border. **Right:** Lung metastases show immune
T cells (red, yellow) and B cells (pink) infiltrating the metastatic tumor-cell
(purple) lesion. Fewer fibroblasts (white, cyan) are seen in the surrounding
invasive margin. The inset shows B cells (pink) and plasma cells (orange). Scale
bar = 100 μm. Statistical significance was assessed using linear mixed-effects
models with patient as a random effect, followed by Benjamini–Hochberg
adjustment for multiple comparisons. Adjusted *p* values are shown:
\*\*\*\* *p* < 0.0001, \*\*\* *p* < 0.001, \*\* *p* < 0.01, and \* *p* < 0.05.

**Figure 4. Hepatic stellate cell activation and CXCL12-CXCR4 signaling are
associated with an immunosuppressive niche in liver metastases.**

**(a)** Cell type scores for aHSCs show selective enrichment in liver metastases
compared to primary colon tumors and lung metastases. **(b)** *CXCL12*
expression is enriched in liver metastases relative to primary CRC and lung
metastases, supporting a *CXCL12*-rich liver metastatic niche. **(c)**
Expression of *CXCR4*, the cognate receptor for *CXCL12*, is also elevated in
liver lesions, supporting a *CXCL12–CXCR4* signaling milieu that may contribute
to immune exclusion. **(d)** In the TCGA-COAD+READ cohort, which predominantly
comprises primary colorectal tumors, the aHSC/fibrotic stromal gene signature
showed a significant association with PFS (Cox proportional hazards *p* =
0.022, HR = 1.426; log-rank *p* = 0.024), supporting the prognostic relevance
of this stromal program in CRC rather than direct validation of the
liver-specific metastatic niche. Statistical significance in (a–c) was
assessed using linear mixed-effects models with patient as a random effect,
followed by Benjamini–Hochberg adjustment for multiple comparisons; adjusted
*p* values are shown (\*\*\* *p* < 0.001, \*\* *p* < 0.01, \* *p* < 0.05).
Survival analyses in (d) were evaluated by Cox proportional hazards and
log-rank tests and were not subject to multiple-comparison adjustment.

**Supplementary Figure 1a. Graphical abstract**

Highlights:

- CRC liver metastases exhibit a dense immunosuppressive microenvironment
  characterized by activated hepatic stellate cells (aHSCs) and fibrotic
  stroma.
- CRC liver metastases show increased immunosuppressive macrophage infiltration
  in association with an activated hepatic stellate cell/myofibroblast program
  and CXCL12/CXCR4 signaling.
- T cells are spatially excluded and primarily restricted to the invasive tumor
  margin in liver metastases, highlighting a strong immune-exclusion mechanism.
- In contrast, synchronous lung metastases demonstrate sparse fibroblasts and
  uniformly infiltrating proliferating B cells and T cells, indicating an
  immunologically permissive microenvironment.
- Primary CRC tumors exhibit a fibroblast-rich, mesenchymal microenvironment
  with moderate infiltration of T cells at the invasive margin and
  immunosuppressive macrophages, representing an intermediate microenvironment.

**Supplementary Figure 1b-c. Spatial immune differences across metastatic
sites**

**(b)** Cell type scores for mCAFs and tCAFs and **(c)** MRC1⁺CSF1R⁺
macrophages and SPP1⁺ macrophages across tumor core (In) and invasive margin
(Ex) regions in primary colon tumors, liver metastases, and lung metastases.
Statistical significance was assessed using linear mixed-effects models with
patient as a random effect, followed by Benjamini–Hochberg adjustment for
multiple comparisons. Adjusted *p* values are shown: \*\*\*\* *p* < 0.0001,
\*\* *p* < 0.01, and \* *p* < 0.05.

**Supplementary Figure 1d. Adjacent non-malignant liver and lung ROIs compared
to metastatic tumors**

Cell type scores for tCAFs and mCAFs, MRC1⁺CSF1R⁺ macrophages, aHSCs, and
*CXCL12* expression across primary colon tumors, liver metastases, lung
metastases, and adjacent non-malignant liver and lung tissue. Statistical
significance was assessed using linear mixed-effects models with patient as a
random effect, followed by Benjamini–Hochberg adjustment for multiple
comparisons. Adjusted *p* values are shown: \*\*\* *p* < 0.001,
\*\* *p* < 0.01, and \* *p* < 0.05.

**Supplementary Figure 2. Multiplex IF quantification for biomarkers across
primary and metastatic sites**

**(a)** Multiplex IF quantification of proliferating B (CD20⁺Ki-67⁺) cells,
plasma cells (CD79a⁺), and CD4⁺ and CD8⁺ T-cell densities across tumor core and
tumor border in primary colon tumors, liver metastases, and lung metastases.
**(b)** Multiplex IF images show individual channels for CD4, SMA, CDX2, CD79a,
CD20, CD8, CD68, CD206, FAP, and Ki-67 in primary CRC, liver metastases, and
lung metastases. Statistical significance was assessed using paired *t* tests.
*p* values are shown: \*\*\* *p* < 0.001 and \* *p* < 0.05.

The captions above are complete. The public analysis code generates Figure
1b–d, Figure 2a–d, Figure 3b–c, Figure 4a–d, and Supplementary Figure 1b–d, as
mapped below. The remaining schematic and multiplex immunofluorescence panels
are included for caption completeness but are not generated by this archive.


## Included Materials

```
code/          one R script per figure, plus shared code in code/lib; code/run is the entry point
data/          read-only inputs (not in git; see data/README.md)
environment/   Dockerfile + install.R defining the compute environment
results/       figures and tables generated by code/run (not in git)
```

| Script | Figure panels | Output |
|---|---|---|
| `01_figure1.R` | 1b, 1c, 1d | `results/figure1/` |
| `02_figure2.R` | 2a–d | `results/figure2/` |
| `03_figure3.R` | 3b, 3c | `results/figure3/` |
| `04_figure4.R` | 4a–d | `results/figure4/` |
| `05_supplementary_figure1.R` | Supp 1b–c, Supp 1d | `results/supplementary_figure1/` |

Each script writes its panels next to a CSV of the values plotted, so results
can be checked without re-running. Each script can also be run on its own.

Shared code in `code/lib/`:

```
paths.R            /data and /results resolution
load_geomx.R       the GeoMx object and its normalised expression matrix
gsva_scores.R      per-ROI signature scores
models.R           mixed-model fits and Benjamini–Hochberg adjustment
panels.R           box and violin panels
de_volcano_gsea.R  differential expression, volcano plots, Hallmark GSEA
```

The box and violin panels in the published figures were drawn in GraphPad
Prism from the values these scripts produce; the R panels show the same values
but are not pixel-identical.

## Reproducing the analyses

**On Code Ocean.** Press **Run**. Results are written to `/results`.

**Locally.**

1. Install R 4.5.1 and the packages in `environment/install.R` (Bioconductor
   3.22), or build the image from `environment/Dockerfile`.
2. Place the input files in `data/` as listed in `data/README.md`. The two TCGA
   files are public and can be downloaded with `cd code && ./fetch_tcga.sh`.
3. Check the inputs and run:

```bash
cd code
Rscript check_inputs.R     # exits 1 if anything is missing or malformed
./run                      # every figure (about 7 minutes)
./run 02                   # Figure 2 only
```

Paths resolve to `/data` and `/results` when those exist (Code Ocean),
otherwise to `../data` and `../results`. Override them with `CRC_DATA_DIR` and
`CRC_RESULTS_DIR`.

## Data availability

Raw GeoMx DSP data are deposited in the European Genome-phenome Archive under
accession **EGA00002691488** (controlled access). The processed, de-identified
ROI-level data needed to reproduce the figures are distributed with the Code
Ocean capsule; see `data/README.md`. No patient-level data are stored in this
repository.

TCGA-COAD/READ expression and clinical data are public, from cBioPortal study
`coadread_tcga_pan_can_atlas_2018`.

## Citation

> Basu S, Anandappa G, *et al.* Spatial profiling of primary CRC and
> synchronous lung and liver metastases identifies a liver-specific stromal
> niche associated with immune exclusion. *Clin Cancer Res* 2026.
> DOI assigned at publication.

See also `CITATION.cff`.

Analysis code: Wenbin Liu, Alejandro Jiménez-Sánchez, Yulong Chen,
Jiun-Sheng Chen.

## License

Code is released under the MIT License (`LICENSE`).
