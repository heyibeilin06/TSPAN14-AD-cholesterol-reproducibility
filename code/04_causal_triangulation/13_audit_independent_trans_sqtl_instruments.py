#!/usr/bin/env python3
"""Audit INTERVAL trans-sQTL evidence for the exact TSPAN14 exon5-6 event.

The public INTERVAL file is significance-filtered at P < 1e-5. Therefore this
script treats absent lipid instruments as censored, not as null associations.
"""

from __future__ import annotations

import argparse
import gzip
from collections import Counter
from pathlib import Path

import pandas as pd


TSPAN14 = "ENSG00000108219"
EXACT_START = 80509471
EXACT_END = 80512144


def normalize_gene(value: str) -> str:
    return value.split(".", 1)[0]


def parse_feature(feature: str) -> tuple[int | None, int | None]:
    parts = feature.split(":")
    try:
        return int(parts[1]), int(parts[2])
    except (IndexError, ValueError):
        return None, None


def is_exact_feature(feature: str) -> bool:
    start, end = parse_feature(feature)
    return start == EXACT_START and abs(end - EXACT_END) <= 1


def complement(allele: str) -> str:
    return allele.translate(str.maketrans("ACGT", "TGCA"))


def align_interval_effect(row: pd.Series, ea: str, oa: str) -> tuple[float | None, str]:
    iea, ioa = str(row.effect_allele).upper(), str(row.other_allele).upper()
    ea, oa = ea.upper(), oa.upper()
    if (iea, ioa) == (ea, oa):
        return float(row.b), "direct"
    if (iea, ioa) == (oa, ea):
        return -float(row.b), "swapped"
    if (complement(iea), complement(ioa)) == (ea, oa):
        return float(row.b), "strand"
    if (complement(iea), complement(ioa)) == (oa, ea):
        return -float(row.b), "strand_swapped"
    return None, "incompatible"


def read_tspan14_hits(path: Path) -> pd.DataFrame:
    rows: list[list[str]] = []
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        gene_idx = header.index("gene_id")
        for line in handle:
            values = line.rstrip("\n").split("\t")
            if normalize_gene(values[gene_idx]) == TSPAN14:
                rows.append(values)
    frame = pd.DataFrame(rows, columns=header)
    for column in ("pos_b38", "pos_b37", "af", "b", "b_se", "pval"):
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame["exact_exon5_6"] = frame.phenotype_id.map(is_exact_feature)
    return frame


def extract_lipid_instruments(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t")
    frame = frame[(frame.exposure.isin(["TC", "LDL", "nonHDL"])) & (frame.outcome == "AD")]
    instruments = frame[["SNP", "effect_allele", "other_allele", "beta_x", "se_x", "p_x", "exposure"]].copy()
    return instruments.drop_duplicates(["SNP", "exposure"])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--interval", type=Path, required=True)
    parser.add_argument("--mr-input", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    hits = read_tspan14_hits(args.interval)
    hits.to_csv(args.outdir / "09_interval_tspan14_trans_sqtl_hits.tsv", sep="\t", index=False)

    feature_summary = (
        hits.groupby("phenotype_id", as_index=False)
        .agg(n_hits=("variant_id", "size"), min_p=("pval", "min"), n_chromosomes=("chr", "nunique"))
    )
    feature_summary["exact_exon5_6"] = feature_summary.phenotype_id.map(is_exact_feature)
    feature_summary = feature_summary.sort_values(["exact_exon5_6", "min_p"], ascending=[False, True])
    feature_summary.to_csv(args.outdir / "10_interval_tspan14_feature_summary.tsv", sep="\t", index=False)

    exact = hits[hits.exact_exon5_6].copy()
    exact.to_csv(args.outdir / "11_interval_exact_exon5_6_trans_sqtl.tsv", sep="\t", index=False)

    instruments = extract_lipid_instruments(args.mr_input)
    all_overlap = instruments.merge(hits, left_on="SNP", right_on="variant_id", how="inner")
    all_overlap.to_csv(args.outdir / "16_lipid_iv_any_tspan14_trans_sqtl_overlap.tsv", sep="\t", index=False)
    if all_overlap.empty:
        all_overlap_summary = pd.DataFrame(
            columns=["exposure", "n_rows", "n_lipid_instruments", "n_splicing_features", "minimum_p"]
        )
    else:
        all_overlap_summary = (
            all_overlap.groupby("exposure", as_index=False)
            .agg(
                n_rows=("SNP", "size"),
                n_lipid_instruments=("SNP", "nunique"),
                n_splicing_features=("phenotype_id", "nunique"),
                minimum_p=("pval", "min"),
            )
        )
    all_overlap_summary.to_csv(
        args.outdir / "17_lipid_iv_any_tspan14_trans_sqtl_summary.tsv", sep="\t", index=False
    )

    overlap = instruments.merge(exact, left_on="SNP", right_on="variant_id", how="left", indicator=True)
    aligned = overlap.apply(
        lambda row: align_interval_effect(row, row.effect_allele_x, row.other_allele_x)
        if row._merge == "both"
        else (None, "not_reported_p_lt_1e-5"),
        axis=1,
        result_type="expand",
    )
    overlap[["sqtl_beta_aligned_to_lipid_effect_allele", "alignment"]] = aligned
    overlap.to_csv(args.outdir / "12_lipid_iv_exact_trans_sqtl_overlap.tsv", sep="\t", index=False)

    counts = Counter(
        overlap.assign(reported=overlap["_merge"].eq("both").astype(int))
        .groupby("exposure")["reported"]
        .sum()
        .astype(int)
        .to_dict()
    )
    audit = pd.DataFrame(
        [
            {"metric": "TSPAN14 trans-sQTL rows", "value": len(hits)},
            {"metric": "TSPAN14 splicing features", "value": hits.phenotype_id.nunique()},
            {"metric": "Exact exon5-6 trans-sQTL rows", "value": len(exact)},
            {"metric": "TC instruments reported for exact event", "value": counts["TC"]},
            {"metric": "LDL instruments reported for exact event", "value": counts["LDL"]},
            {"metric": "nonHDL instruments reported for exact event", "value": counts["nonHDL"]},
            {"metric": "Public-file selection threshold", "value": "P < 1e-5"},
            {"metric": "Missing association interpretation", "value": "censored; cannot assume beta=0"},
        ]
    )
    audit.to_csv(args.outdir / "13_interval_trans_sqtl_audit.tsv", sep="\t", index=False)

    print(audit.to_string(index=False))


if __name__ == "__main__":
    main()
