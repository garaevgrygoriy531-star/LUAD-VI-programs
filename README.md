# Vascular invasion-associated program changes concentrate in invasive lung adenocarcinoma and span multiple tissue compartments

## Study overview

This repository contains analyses of four transcriptional programs defined in an independent study of microscopic vascular invasion. We examined how these programs vary across histologic states and invasive LUAD in pathology cohorts, then assessed their associations with vascular invasion, CT solidity, paired tumor regions, proteomic enrichment, and overall survival in stage I disease. Single-nucleus analyses localized program-score differences across epithelial, fibroblast/stromal, endothelial, and myeloid lineages. Independent single-cell and spatial datasets provide cellular and regional concordance.

The analyses characterize cohort-level associations and the distribution of program-score differences across broad tissue compartments.

## Repository contents

- code/ contains analysis, figure-rendering, and table-export scripts.
- data_processed/ contains processed cohort scores, analysis membership, figure inputs, static figure assets, and table sources.
- figures/ contains the five main and five supplementary publication figures.
- metadata/ describes source datasets, program definitions, panels, variables, and file provenance.
- results/ contains model estimates and analysis summaries.
- docs/ describes data access, licensing, provenance, and reproduction.
- environment/ records software versions and package requirements.

## Reproducing the analyses

Run the commands in docs/REPRODUCIBILITY.md from the repository root. The included processed data support the pathology, microscopic vascular-invasion, stage I survival, stage-architecture, and cellular analyses. Analyses that require upstream source matrices are identified separately.

## Reproducing figures and supplementary tables

Prepare the derived Figure 3C and Supplementary Figure S4D inputs, render the figures to a new directory, and export Supplementary Tables S1–S9 using the commands in the reproducibility guide. Figure 1A uses the supplied native SVG artwork; the UMAP point layer in Figure 5A is supplied as a PNG asset. The publication PDFs in figures/ are the corresponding final display files.

## Data sources

The study reuses public datasets from GEO and PRIDE. Accessions, analysis populations, and cohort roles are listed in metadata/datasets.csv. GSE68465 was used for program-coverage screening; V2–V4 coverage did not meet the prespecified 80% threshold, so this cohort was not included in outcome-score analyses. Upstream expression matrices not included here can be obtained from their source repositories.

## Software environment

The analysis and rendering scripts use R and Python. Tested versions and package requirements are listed under environment/.

## Citation

Associated manuscript: Vascular invasion-associated program changes concentrate in invasive lung adenocarcinoma and span multiple tissue compartments.

Formal bibliographic and archival DOI information will be added when available.

## License

LICENSE-CODE applies to original code only. Source datasets and derived data remain subject to their source repositories, publications, and data-use conditions. See docs/DATA_LICENSE_AND_SOURCES.md.

