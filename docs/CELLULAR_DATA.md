# Cellular and spatial data

The GSE308103 discovery atlas contains 798,100 nuclei from 75 specimens and 23 patients. The repository includes nucleus and specimen annotations, fine-cell-type program scores, cell-eligibility summaries, and lineage-level results. Quantitative analyses use 25 eligible fine cell types across epithelial, fibroblast/stromal, endothelial, and myeloid lineages; other populations remain visible as atlas context.

The primary cellular decomposition partitions invasive-LUAD-minus-earlier-state contrasts into within-cell-type program-score differences, fine-cell-type composition, covariance, and total effect within each broad lineage. The V1-adjusted sensitivity analysis recalculates these components after adjustment for the proliferation-related V1 program. Its 5,000 patient-cluster bootstrap refits the V1 adjustment in each replicate. A patient-only random-intercept model is included as a model-specification sensitivity.

External concordance is assessed in GSE189357, which includes 115,253 cells from nine treatment-naive patients. GSE189487 provides spatial support from six specimens in the same study population and is not an independent patient cohort. The included scripts refit analyses from processed scores; rebuilding scores from original count matrices requires the source studies' data and upstream annotations.
