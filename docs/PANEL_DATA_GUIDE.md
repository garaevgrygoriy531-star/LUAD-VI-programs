# Panel and table inputs

metadata/figure_table_map.csv maps figure panels and tables to their processed inputs, outputs, and rendering or export scripts. The figure renderer reads the tables in data_processed/plot_ready/. Figure 1A uses the native SVG asset in data_processed/figure_assets/F1a.svg; Figure 5A uses the supplied UMAP point-layer image.

The supplementary workbook is data_processed/supplementary_tables/supplementary_tables_s1_s9.xlsx. Its worksheets are exported to data_processed/supplementary_tables/csv/ by code/10_tables/export_supplementary_tables.py. Figure 3C inputs are derived from the integrated pathology contrasts, and Supplementary Figure S4D inputs are derived from the included V1-adjusted decomposition result table.
