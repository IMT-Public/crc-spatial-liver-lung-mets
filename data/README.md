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
| `data_tcga/some_cell_subtype_markers.xlsx` | Cell-type marker gene lists | Fig 4d |
| `data_tcga/data_clinical_patient.txt` | TCGA clinical data (not in git) | Fig 4d |
| `data_tcga/data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt` | TCGA expression z-scores (not in git) | Fig 4d |

**TCGA files.** The last two files are public (cBioPortal study
`coadread_tcga_pan_can_atlas_2018`) and are not stored in git. Download them
with `cd code && ./fetch_tcga.sh`. The Code Ocean capsule already includes
them.
