# Reproducibility

Run commands from the repository root.

## Figures and supplementary tables

Prepare the derived panel inputs from the included result tables:

    python code/09_figures/prepare_f3c_data.py
    python code/09_figures/prepare_s4d_data.py

Render the figures to a separate output directory. The script uses supplied estimates and confidence limits; it does not refit upstream models.

    Rscript --vanilla code/09_figures/render_figures.R . results/rendered

Figure 1A is rendered from the supplied native SVG asset in data_processed/figure_assets/. Figure 5A uses the supplied UMAP point layer in data_processed/plot_ready/.

Export Supplementary Tables S1–S9 from the included workbook:

    python code/10_tables/export_supplementary_tables.py

Run the repository checks:

    python tests/run_repository_checks.py

## Analyses supported by included processed data

    Rscript --vanilla code/02_stage_architecture/fit_stage_architecture.R .
    Rscript --vanilla code/02_stage_architecture/fit_stage_ia_contrasts.R .
    Rscript --vanilla code/03_pathology_replication/fit_v1_adjusted_stage_models.R .
    Rscript --vanilla code/04_proliferation_and_VI/fit_vi_associations.R .
    Rscript --vanilla code/04_cross_context/run_cross_context_pathology.R .
    Rscript --vanilla code/06_stageI_outcome/fit_stage_i_survival.R .
    Rscript --vanilla code/07_cellular_discovery/analyze_cellular_processed_scores.R .
    Rscript --vanilla code/07_cellular_discovery/analyze_fine_celltype_scores.R .

Stage-architecture models use a patient random intercept and maximum likelihood; model support is compared using AICc. The separate GSE319666 Stage IA contrast uses REML with t-based 95% confidence intervals. Figure 3A contrasts use a patient-cluster bootstrap with 1,000 resamples and seed 20261016. Each resampled patient contributes all of that patient's specimens; duplicated patient clusters are treated as independent groups.

## Analyses requiring source matrices

GSE282617 sample scores are included. To reconstruct them from the GEO processed FPKM matrix, download the matrix and run:

    python code/04_cross_context/reconstruct_gse282617_scores.py --root . --matrix <path-to-matrix.csv.gz>

Public metadata do not provide a patient crosswalk for GSE282617, so its contrast uses samples as the inferential units. For GSE193725, the repository includes the analyzed CT-solidity effects and cohort counts; the lesion-level scores required to refit that model are not included.

Fine-cell-type models can be refit from the included processed scores. Rebuilding program scores from raw expression or count matrices requires source data and processing details from the original studies. Accession identifiers and available processed objects are listed in metadata/datasets.csv and metadata/provenance_manifest.csv.

Proteomic enrichment and sensitivity summaries are study-derived outputs. The per-gene rank table is not included; see docs/DATA_PROVENANCE.md for its source and reconstruction route.

## Software

R and Python versions and package requirements are listed under environment/.
