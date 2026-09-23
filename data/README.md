# Data

Processed, de-identified ROI-level GeoMx DSP data (219 ROIs from 8 patients)
used to reproduce the figures. Raw data are under controlled access at the
European Genome-phenome Archive, accession **EGA00002691488**.

| File | Contents | Used by |
|---|---|---|
| `data_volcanoplot_pathway/targets.new-219-ROIs-newHKthenQ3.RDATA` | GeoMx object: normalised expression and ROI annotation | all figures |
| `data_volcanoplot_pathway/pathwayInfo.RDATA` | MSigDB Hallmark gene sets | Fig 1d |
| `data_volcanoplot_pathway/DE_res_combinedROIs_lme.xlsx` | Mixed-model differential expression results | Fig 1c–d |
| `data_boxplots/GSVAscores_updated_03272025.csv` | Per-ROI signature scores | Figs 2–4, Supp Fig 1 |
| `data_boxplots/CD4TnCCR7_top50genes.RDATA` | CD4+ CCR7+ T-cell gene signature | Figs 2–4, Supp Fig 1 |
| `data_cibersortx/CIBERSORTx_Job14_Results_02152024.csv` | CIBERSORTx cell-type fractions per ROI | Fig 1b |
| `data_tcga/TCGA_COADREAD_aHSC_GSVA_PFS.csv` | Processed TCGA table: aHSC GSVA score and PFS per patient (592 patients) | Fig 4d |
| `data_tcga/some_cell_subtype_markers.xlsx` | Cell-type marker gene lists | regenerating Fig 4d scores |

**TCGA data.** To avoid redistributing third-party data, the TCGA-COAD/READ
expression matrix and clinical file are not in this repository or the Code
Ocean capsule. Only the processed table needed to reproduce Fig 4d is
included. Its columns are `PATIENT_ID` (TCGA barcode), `PFS_MONTHS`,
`PFS_STATUS` (1 = progressed, 0 = censored) and `gsva` (the aHSC GSVA score).

To regenerate the table from the source data:

1. Download the two public files from cBioPortal (study
   `coadread_tcga_pan_can_atlas_2018`,
   <https://www.cbioportal.org/study/summary?id=coadread_tcga_pan_can_atlas_2018>),
   subject to the TCGA and cBioPortal terms of use:
   `data_clinical_patient.txt` and
   `data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt`. `cd code && ./fetch_tcga.sh`
   downloads both into `data/data_tcga/`.
2. In `code/04_figure4.R`, uncomment the block marked
   "scoring from raw TCGA" and run the script once. It rewrites
   `TCGA_COADREAD_aHSC_GSVA_PFS.csv`.
