#!/usr/bin/env python3
"""Build Figure 3C attenuation values from the supplied pathology coefficients."""
from pathlib import Path
import csv
ROOT = Path(__file__).resolve().parents[2]
source = ROOT / "data_processed/plot_ready/F3a_integrated.csv"
target = ROOT / "data_processed/plot_ready/F3c.csv"
with source.open(encoding="utf-8-sig", newline="") as f:
    rows = list(csv.DictReader(f))
out = []
for dataset in dict.fromkeys(r["dataset"] for r in rows):
    for program in ("V2", "V3", "V4"):
        group = [r for r in rows if r["dataset"] == dataset and r["program"] == program]
        unadjusted = next(r for r in group if r["series"] == "Unadjusted")
        adjusted = next(r for r in group if r["series"] == "V1-adjusted")
        pct = 100 * (1 - abs(float(adjusted["estimate"])) / abs(float(unadjusted["estimate"])))
        out.append({"dataset": dataset, "program": program, "attenuation_percent": f"{pct:.1f}",
                    "unadjusted_coefficient": unadjusted["estimate"], "adjusted_coefficient": adjusted["estimate"]})
with target.open("w", encoding="utf-8-sig", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(out[0]))
    w.writeheader()
    w.writerows(out)
