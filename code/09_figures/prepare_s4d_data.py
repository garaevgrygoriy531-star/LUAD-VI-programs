#!/usr/bin/env python3
"""Build Figure S4D inputs from the supplied V1-adjusted decomposition results."""
from pathlib import Path
import csv
ROOT = Path(__file__).resolve().parents[2]
source = ROOT / "results/cellular/v1_adjusted_decomposition_sensitivity.csv"
target = ROOT / "data_processed/plot_ready/S4d.csv"
with source.open(encoding="utf-8-sig", newline="") as f:
    rows = list(csv.DictReader(f))
out = []
for r in rows:
    for component, stem in (("Within-cell-type", "V1_adjusted_within"),
                            ("Fine-cell-type composition", "V1_adjusted_composition")):
        out.append({"lineage": r["lineage"], "program": r["program"], "component": component,
                    "estimate": r[stem], "ci_low": r[stem + "_CI_low"],
                    "ci_high": r[stem + "_CI_high"], "bootstrap_replicates": r["n_bootstrap_valid"],
                    "analysis": "V1-adjusted primary decomposition"})
with target.open("w", encoding="utf-8-sig", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(out[0]))
    w.writeheader()
    w.writerows(out)
