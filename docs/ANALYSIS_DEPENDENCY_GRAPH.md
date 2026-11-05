# Analysis dependencies

metadata/analysis_dependency_graph.csv lists source cohorts, processed inputs, analysis scripts, outputs, and inferential units for the principal analyses. metadata/figure_table_map.csv maps each figure panel and table to its inputs, output, and script. metadata/panel_manifest.csv describes each panel's scientific role and statistical representation.

Stage-architecture models account for repeated specimens within patients. Stage I survival models use patients as the analysis unit. Cross-context, regional, cellular, and spatial analyses use the inferential unit recorded for each cohort in metadata/datasets.csv.

Figure 3C is calculated from the unadjusted and V1-adjusted contrasts supplied for Figure 3A. Supplementary Figure S4D is assembled from the V1-adjusted decomposition result table in results/cellular/; the same result is represented in Section F of Supplementary Table S8.
