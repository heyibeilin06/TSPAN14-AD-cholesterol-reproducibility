"""Assess public-brain co-usage of two TSPAN14 junctions.

This replaces inaccessible controlled individual-level RNA-seq with the open
GTEx v8 exon-exon junction-count matrix. It tests co-usage only; it does not
recompute an sQTL and must not be interpreted as transcript-level PSI.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import os
from pathlib import Path

import pandas as pd
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[2]
DATA = Path(os.environ.get("SLM_DATA_ROOT", ROOT / "data" / "raw"))
OUT = ROOT / "outputs" / "p1_public_junction_cousage"
DEFAULT_JUNCTIONS = {
    # STAR IDs encode intron boundaries; see the explicit mapping in the output manifest.
    "ex5_6": "chr10_80509472_80512143",
    "ex6_7": "chr10_80512270_80514018",
}
SOURCE_URLS = {
    "junction_matrix": "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_STARv2.5.3a_junctions.gct.gz",
    "sample_attributes": "https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
}


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if path.suffix == ".gz" else path.open(encoding="utf-8", errors="replace")


def extract_brain_junction_counts(matrix_path: Path, attributes_path: Path, junctions: dict[str, str]) -> pd.DataFrame:
    """Return selected open-GTEx junction counts for samples annotated as brain."""
    attributes = pd.read_csv(attributes_path, sep="\t", dtype=str, usecols=["SAMPID", "SMTSD"])
    brain = attributes.loc[attributes["SMTSD"].str.startswith("Brain", na=False)].copy()
    tissue_by_sample = dict(zip(brain["SAMPID"], brain["SMTSD"], strict=True))
    requested = {junction: label for label, junction in junctions.items()}
    rows: dict[str, list[str]] = {}

    with open_text(matrix_path) as handle:
        first = handle.readline().rstrip("\n")
        if first.startswith("#"):
            handle.readline()  # GCT dimension line
        header = handle.readline().rstrip("\n").split("\t")
        samples = header[2:]
        selected_indices = [(index + 2, sample) for index, sample in enumerate(samples) if sample in tissue_by_sample]
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            label = requested.get(fields[0])
            if label:
                rows[label] = [fields[index] for index, _ in selected_indices]
                if len(rows) == len(junctions):
                    break

    missing = sorted(set(junctions) - set(rows))
    if missing:
        raise ValueError(f"Requested junction(s) absent from matrix: {', '.join(missing)}")
    result = pd.DataFrame({
        "sample_id": [sample for _, sample in selected_indices],
        "donor_id": ["-".join(sample.split("-")[:2]) for _, sample in selected_indices],
        "tissue": [tissue_by_sample[sample] for _, sample in selected_indices],
        **rows,
    })
    for label in junctions:
        result[label] = pd.to_numeric(result[label], errors="raise")
    return result


def summarize_cousage(frame: pd.DataFrame, first: str, second: str) -> pd.DataFrame:
    """Summarize zero-aware junction-count co-usage within each tissue and overall."""
    rows: list[dict[str, object]] = []
    for tissue, subset in [("All brain tissues", frame), *frame.groupby("tissue", sort=True)]:
        x, y = subset[first].map(math.log1p), subset[second].map(math.log1p)
        if len(subset) >= 3:
            test = spearmanr(x, y)
            rho, p_value = test.statistic, test.pvalue
        else:
            rho, p_value = float("nan"), float("nan")
        rows.append({
            "tissue": tissue,
            "n_samples": len(subset),
            "n_donors": subset["donor_id"].nunique() if "donor_id" in subset else len(subset),
            "n_both_nonzero": int(((subset[first] > 0) & (subset[second] > 0)).sum()),
            "spearman_rho": rho,
            "spearman_p": p_value,
            "metric_definition": "Spearman correlation of log1p exon-exon junction read counts; not transcript-level PSI.",
            "interpretation_boundary": "Open normal-brain RNA-seq co-usage only; no genotype, AD-state, or causal inference.",
        })
    return pd.DataFrame(rows)


def write_tsv(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--junction-matrix", type=Path, required=True)
    parser.add_argument("--sample-attributes", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, default=OUT)
    args = parser.parse_args()
    counts = extract_brain_junction_counts(args.junction_matrix, args.sample_attributes, DEFAULT_JUNCTIONS)
    summary = summarize_cousage(counts, "ex5_6", "ex6_7")
    write_tsv(args.out_dir / "gtex_v8_tspan14_brain_junction_counts.tsv", counts)
    write_tsv(args.out_dir / "gtex_v8_tspan14_brain_cousage_summary.tsv", summary)
    manifest = pd.DataFrame([{
        "resource": "GTEx v8 open STAR exon-exon junction read counts",
        "junction_matrix": str(args.junction_matrix),
        "sample_attributes": str(args.sample_attributes),
        "junction_matrix_url": SOURCE_URLS["junction_matrix"],
        "sample_attributes_url": SOURCE_URLS["sample_attributes"],
        "junction_matrix_bytes": args.junction_matrix.stat().st_size,
        "sample_attributes_bytes": args.sample_attributes.stat().st_size,
        "junction_ex5_6": DEFAULT_JUNCTIONS["ex5_6"],
        "junction_ex6_7": DEFAULT_JUNCTIONS["ex6_7"],
        "coordinate_mapping": "MiGA project chr10:80509471-80512144 -> GTEx STAR chr10_80509472_80512143; external chr10:80512269-80514018 -> GTEx STAR chr10_80512270_80514018. All four IDs have the TSPAN14 Ensembl gene ID in their respective resources.",
        "scope": "Public normal-brain junction co-usage replacement for inaccessible controlled individual-level RNA-seq; not PSI, sQTL, or AD-state evidence.",
    }])
    write_tsv(args.out_dir / "gtex_v8_public_resource_manifest.tsv", manifest)


if __name__ == "__main__":
    main()
