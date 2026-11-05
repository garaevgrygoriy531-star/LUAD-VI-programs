# Environment

The included analyses were run with R 4.5.1 and Python 3.12.14. Package versions for the R scripts are listed in R_packages.csv; Python requirements are listed in requirements.txt. R package versions in this table identify the running environment, not the build version of each Windows binary.

Install the R packages in R_packages.csv from a CRAN mirror. Install Python packages with:

    python -m pip install -r environment/requirements.txt

Run the commands in docs/REPRODUCIBILITY.md from the repository root. The figure renderer also requires a graphics device with Cairo support and an Arial-compatible font.
