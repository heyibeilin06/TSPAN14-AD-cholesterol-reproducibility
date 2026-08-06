#!/usr/bin/env python3
"""Quantify donor reuse across GTEx neural-tissue junction-count contexts."""

from __future__ import annotations

import argparse
from itertools import combinations
from pathlib import Path

import pandas as pd


FOCAL_TISSUES = [
    "Brain - Anterior cingulate cortex (BA24)",
    "Brain - Hippocampus",
    "Brain - Putamen (basal ganglia)",
    "Brain - Spinal cord (cervical c-1)",
]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t")
    focal = df[df["tissue"].isin(FOCAL_TISSUES)].copy()
    donors = {tissue: set(group["donor_id"].dropna()) for tissue, group in focal.groupby("tissue")}

    tissue_summary = (
        focal.groupby("tissue", as_index=False)
        .agg(samples=("sample_id", "nunique"), donors=("donor_id", "nunique"))
    )
    overlap_rows = []
    for first, second in combinations(FOCAL_TISSUES, 2):
        a, b = donors.get(first, set()), donors.get(second, set())
        union = a | b
        overlap_rows.append(
            {
                "tissue_1": first,
                "tissue_2": second,
                "donors_tissue_1": len(a),
                "donors_tissue_2": len(b),
                "shared_donors": len(a & b),
                "jaccard_index": len(a & b) / len(union) if union else 0.0,
                "fraction_tissue_1_shared": len(a & b) / len(a) if a else 0.0,
                "fraction_tissue_2_shared": len(a & b) / len(b) if b else 0.0,
            }
        )

    membership = (
        focal[["donor_id", "tissue"]]
        .drop_duplicates()
        .assign(present=1)
        .pivot(index="donor_id", columns="tissue", values="present")
        .fillna(0)
        .astype(int)
        .reset_index()
    )
    membership["n_focal_tissues"] = membership[FOCAL_TISSUES].sum(axis=1)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    tissue_summary.to_csv(args.output_dir / "gtex_focal_tissue_donor_counts.tsv", sep="\t", index=False)
    pd.DataFrame(overlap_rows).to_csv(args.output_dir / "gtex_focal_tissue_donor_overlap.tsv", sep="\t", index=False)
    membership.to_csv(args.output_dir / "gtex_focal_tissue_donor_membership.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
