# Data availability

The study reuses public datasets from the Gene Expression Omnibus (GEO) and PRIDE. Accession identifiers, analytical populations, and source publications are listed in metadata/datasets.csv and metadata/references.csv.

The repository includes processed program scores, cohort membership tables, cellular and spatial summaries, proteomic enrichment summaries, figure inputs, and Supplementary Tables S1–S9. Large source expression, sequencing, and proteomics data are available from the repositories hosting the original studies. GSE282617 sample-level scores can be reconstructed from its GEO processed FPKM matrix. For GSE193725, the analyzed CT-solidity effects are included; lesion-level scores needed to refit that contrast are not distributed here. The PXD046637 per-gene ranking table is not included; see docs/DATA_PROVENANCE.md for its source and reconstruction route.

GSE68465 is retained for program-coverage screening and is shown in Figure S3A. V2–V4 coverage was below 80%; the cohort was therefore excluded from outcome-score analyses.
