from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

PROGRAMS = ("V1", "V2", "V3", "V4")


def read_members(root: Path) -> dict[str, list[str]]:
    members = {}
    for program in PROGRAMS:
        path = root / "metadata" / "program_members" / f"{program}.txt"
        members[program] = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    return members


def main() -> None:
    parser = argparse.ArgumentParser(description="Reconstruct sample-level GSE282617 expression scores.")
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2], help="Repository root")
    parser.add_argument("--matrix", type=Path, required=True, help="GEO processed FPKM matrix (.csv or .csv.gz)")
    args = parser.parse_args()
    root = args.root.resolve()
    manifest_path = root / "metadata" / "gse282617_sample_manifest.csv"
    manifest = pd.read_csv(manifest_path, dtype={"sample_id": str, "GSM": str})
    included = manifest.loc[manifest["in_primary_contrast"].astype(bool)].copy()
    if len(included) != 60 or included["pathology_group"].value_counts().to_dict() != {"IAC": 23, "AIS": 20, "MIA": 17}:
        raise ValueError("The eligible GSE282617 sample set does not match the documented 20 AIS, 17 MIA, and 23 IAC samples.")

    matrix = pd.read_csv(args.matrix, compression="infer")
    if "gene_name" not in matrix.columns:
        raise ValueError("Expected a gene_name column in the GEO processed matrix.")
    column_by_sample = {column[5:]: column for column in matrix.columns if column.startswith("FPKM.")}
    missing_samples = sorted(set(included["sample_id"]) - set(column_by_sample))
    if missing_samples:
        raise ValueError(f"Missing eligible samples in expression matrix: {missing_samples}")
    sample_columns = [column_by_sample[sid] for sid in included["sample_id"]]
    log_expression = np.log2(matrix[sample_columns].astype(float).to_numpy() + 1.0)
    symbols = matrix["gene_name"].fillna("").astype(str).str.strip().to_numpy()

    # Resolve duplicate gene symbols by their mean log2(FPKM + 1) across the 60 eligible samples.
    chosen: dict[str, tuple[float, int]] = {}
    for i, symbol in enumerate(symbols):
        if not symbol:
            continue
        mean_value = float(np.nanmean(log_expression[i]))
        if symbol not in chosen or mean_value > chosen[symbol][0]:
            chosen[symbol] = (mean_value, i)
    expression = {symbol: log_expression[i] for symbol, (_, i) in chosen.items()}
    memberships = read_members(root)

    score_rows = []
    symbol_rows = []
    score_data: dict[str, np.ndarray] = {}
    for program, genes in memberships.items():
        present = [gene for gene in genes if gene in expression]
        if len(present) / len(genes) < 0.80:
            raise ValueError(f"{program} has less than 80% exact-symbol coverage: {len(present)}/{len(genes)}")
        values = np.vstack([expression[gene] for gene in present])
        means = values.mean(axis=1, keepdims=True)
        sds = values.std(axis=1, ddof=1, keepdims=True)
        if np.any(sds == 0):
            raise ValueError(f"Zero-variance genes found in {program}.")
        raw_score = ((values - means) / sds).mean(axis=0)
        standardized_score = (raw_score - raw_score.mean()) / raw_score.std(ddof=1)
        score_data[program] = standardized_score
        for gene in genes:
            alias = gene.replace(".", "-") if "." in gene else ""
            symbol_rows.append({
                "program": program,
                "signature_symbol": gene,
                "exact_match": gene in expression,
                "hyphenated_candidate": alias,
                "hyphenated_candidate_present": bool(alias and alias in expression),
                "alias_applied": False,
            })

    out = included[["sample_id", "GSM", "pathology_group", "analysis_role", "in_primary_contrast"]].copy()
    for program in PROGRAMS:
        out[program] = score_data[program]
    out.to_csv(root / "data_processed" / "cross_context" / "gse282617_sample_scores.csv", index=False, float_format="%.15g")
    pd.DataFrame(symbol_rows).to_csv(root / "metadata" / "gse282617_symbol_matching.csv", index=False)

    coverage = []
    for program, genes in memberships.items():
        matched = sum(gene in expression for gene in genes)
        coverage.append({"program": program, "signature_genes": len(genes), "exactly_measured_genes": matched, "coverage_fraction": matched / len(genes)})
    coverage_path = root / "results" / "cross_context" / "gse282617_program_coverage.csv"
    coverage_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(coverage).to_csv(coverage_path, index=False)
    print(pd.DataFrame(coverage).to_string(index=False))
    print(f"Wrote {len(out)} eligible samples; patient IDs are not available in the GEO sample metadata.")


if __name__ == "__main__":
    main()
