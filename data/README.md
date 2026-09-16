# Data manifest

The study files here are processed, de-identified, ROI-level matrices — not
raw sequencing data — and are tracked in git. The two public TCGA files are
not tracked; download them with `code/fetch_tcga.sh` (see below). All nine
files are bundled in the Code Ocean capsule.

Raw GeoMx DSP data are under controlled access at the European
Genome-phenome Archive, accession **EGA00002691488**.

Verify the directory at any time:

```bash
cd code && Rscript check_inputs.R
```

## What has to be here

| File | Size | Needed by |
|---|---|---|
| `data_volcanoplot_pathway/targets.new-219-ROIs-newHKthenQ3.RDATA` | 8.9 MB | all five scripts |
| `data_volcanoplot_pathway/pathwayInfo.RDATA` | 2.6 MB | `01_figure1.R` (Fig 1d) |
| `data_volcanoplot_pathway/DE_res_combinedROIs_lme.xlsx` | 1.4 MB | `01_figure1.R` (Fig 1c–d) |
| `data_boxplots/GSVAscores_updated_03272025.csv` | 69 KB | 02–05 |
| `data_boxplots/CD4TnCCR7_top50genes.RDATA` | 292 B | 02–05 |
| `data_cibersortx/CIBERSORTx_Job14_Results_02152024.csv` | 59 KB | `01_figure1.R` (Fig 1b) |
| `data_tcga/some_cell_subtype_markers.xlsx` | 13 KB | `04_figure4.R` (Fig 4d) |
| `data_tcga/data_clinical_patient.txt` | 171 KB | `04_figure4.R` (Fig 4d) |
| `data_tcga/data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt` | 83 MB | `04_figure4.R` (Fig 4d) |

## The two TCGA files are public — download them

Both come from cBioPortal study **`coadread_tcga_pan_can_atlas_2018`**
(TCGA-COAD and TCGA-READ, PanCancer Atlas processing). Fetch them with:

```bash
cd code && ./fetch_tcga.sh
```

which downloads both files from the cBioPortal datahub repository
(`public/coadread_tcga_pan_can_atlas_2018/`). To do it by hand instead, go to
<https://www.cbioportal.org/study/summary?id=coadread_tcga_pan_can_atlas_2018>,
use the download link, and take these two files from the archive:

```
data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt
data_clinical_patient.txt
```

A Code Ocean Reproducible Run has **no network access**, so both are bundled
into `/data` for the capsule rather than downloaded at runtime.

## What each study-specific file is

### `data_volcanoplot_pathway/`

| File | Description |
|---|---|
| `targets.new-219-ROIs-newHKthenQ3.RDATA` | `NanoStringGeoMxSet` object: 219 ROIs from 8 patients, housekeeping- then Q3-normalised counts (`log_q` and `q_norm` assay slots) plus the ROI phenotype table. 208 ROIs survive the `Area != "TBD" & Area != "TLS"` filter; 160 remain after dropping adjacent normal. |
| `pathwayInfo.RDATA` | `pathway_info`, a list of MSigDB gene sets. Only the `HALLMARK` element is read. |
| `DE_res_combinedROIs_lme.xlsx` | Per-gene linear mixed-effects differential expression, 4 sheets (Liver–Colon, Liver–Lung, Lung–Colon, TIL_TB–T in colon). Columns `Gene, Subset, Contrast, Estimate, Pr(>\|t\|), FDR`. The fit takes hours over ~18k genes, so this ships as an input rather than being re-run. |

### `data_boxplots/`

| File | Description |
|---|---|
| `GSVAscores_updated_03272025.csv` | Per-ROI GSVA scores for the single-cell-derived signatures (mCAF, tCAF, `Mph.PLTP`, `Mph.SPP1`, `B.MKI67`, `aHSC`, …) with `ROI, patient, Site, Site2, Area` annotation columns. |
| `CD4TnCCR7_top50genes.RDATA` | `top.genes`: the 50-gene CD4+CCR7+ Tcm signature. |

### `data_cibersortx/`

| File | Description |
|---|---|
| `CIBERSORTx_Job14_Results_02152024.csv` | CIBERSORTx output (safeTME signature matrix × Q3-normalised mixture matrix): 18 cell-type fractions per ROI. **Row order matches `pData(targets.new)`**, and `01_figure1.R` asserts it. |

### `data_tcga/`

| File | Description | Source |
|---|---|---|
| `data_clinical_patient.txt` | TCGA-COAD/READ patient clinical table (`PFS_MONTHS`, `PFS_STATUS`, …). | cBioPortal |
| `data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt` | RSEM z-scores against all samples. | cBioPortal |
| `some_cell_subtype_markers.xlsx` | One sheet per signature. The `aHSC` sheet holds `COL3A1, COL1A1, COL1A2, ACTA2, TAGLN` (Correia et al., *Nature* 2021) and is selected **by name**. | this study |

## De-identification

Patients are labelled `1`–`8`, and slides `slide_01`–`slide_27`. No figure
script reads the slide label. ROIs are labelled `ROI_001`–`ROI_219` (with the
`.dcc` suffix where the GeoMx object uses it) in place of GeoMx DSP scan IDs,
consistently across the GeoMx object, the GSVA score table and the CIBERSORTx
table. Lab-internal run metadata (plate, well, run date, internal sample code)
has been removed from the object's `protocolData`. The key linking these codes
to the original identifiers is held privately by the study team.
`check_inputs.R` fails if a slide label looks like a pathology accession
number or an ROI still carries a DSP scan ID.
