#!/usr/bin/env python3
"""Prepare donor- and state-level microglial pseudobulk data from GSE243292."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import date
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


REACTIVE_STATES = {"ARM", "Dystrophic", "Ex_microglia"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def aggregate(matrix, groups: pd.Series) -> tuple[pd.DataFrame, pd.DataFrame]:
    categories = sorted(groups.unique().tolist())
    rows = []
    metadata = []
    for category in categories:
        mask = groups.to_numpy() == category
        summed = np.asarray(matrix[mask].sum(axis=0)).ravel().astype(np.int64)
        rows.append(summed)
        metadata.append({"sample": category, "n_cells": int(mask.sum())})
    counts = pd.DataFrame(np.vstack(rows).T, columns=categories)
    return counts, pd.DataFrame(metadata)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--minimum-state-cells", type=int, default=10)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    dataset = ad.read_h5ad(args.h5ad)
    matrix = dataset.X.tocsr() if sparse.issparse(dataset.X) else sparse.csr_matrix(dataset.X)
    obs = dataset.obs.copy()
    obs["sampleID"] = obs["sampleID"].astype(str)
    obs["pathology_stage"] = obs["atscore"].map({"A-T-": 0, "A+T-": 1, "A+T+": 2})
    obs["ad_binary"] = (obs["atscore"] == "A+T+").astype(int)
    obs["apoe4_dosage"] = obs["apoe"].map({"E3/E3": 0, "E3/E4": 1, "E4/E4": 2})
    obs["trem2_r47h"] = (obs["trem2"] == "R47H").astype(int)
    obs["state_group"] = obs["microglia_subpopulations"].map(
        lambda state: "Reactive_associated" if state in REACTIVE_STATES else "Other_microglia"
    )
    total_umis = np.asarray(matrix.sum(axis=1)).ravel()
    obs["total_umis"] = total_umis

    donor_counts, donor_metadata = aggregate(matrix, obs["sampleID"])
    donor_counts.index = dataset.var_names.astype(str)
    donor_counts.index.name = "gene"
    donor_rows = []
    for sample_id, group in obs.groupby("sampleID", observed=True):
        first = group.iloc[0]
        donor_rows.append(
            {
                "sample": sample_id,
                "atscore": first["atscore"],
                "pathology_stage": int(first["pathology_stage"]),
                "ad_binary": int(first["ad_binary"]),
                "apoe": first["apoe"],
                "apoe4_dosage": int(first["apoe4_dosage"]),
                "trem2": first["trem2"],
                "trem2_r47h": int(first["trem2_r47h"]),
                "n_cells": int(len(group)),
                "median_umis": float(group["total_umis"].median()),
                "reactive_associated_fraction": float(
                    (group["state_group"] == "Reactive_associated").mean()
                ),
            }
        )
    donor_metadata = pd.DataFrame(donor_rows).sort_values(
        "sample", key=lambda series: series.astype(int)
    )
    donor_counts = donor_counts[donor_metadata["sample"].tolist()]

    state_labels = obs["sampleID"] + "__" + obs["state_group"]
    state_counts, state_metadata = aggregate(matrix, state_labels)
    state_counts.index = dataset.var_names.astype(str)
    state_counts.index.name = "gene"
    state_metadata[["donor", "state_group"]] = state_metadata["sample"].str.split(
        "__", n=1, expand=True
    )
    state_metadata = state_metadata[
        state_metadata["n_cells"] >= args.minimum_state_cells
    ].copy()
    state_metadata = state_metadata.sort_values(
        ["donor", "state_group"], key=lambda series: series.astype(int) if series.name == "donor" else series
    )
    state_counts = state_counts[state_metadata["sample"].tolist()]

    tspan_index = np.where(dataset.var_names.astype(str) == "TSPAN14")[0]
    if len(tspan_index) != 1:
        raise ValueError(f"Expected one TSPAN14 feature, observed {len(tspan_index)}")
    tspan = np.asarray(matrix[:, tspan_index[0]].toarray()).ravel()
    obs_summary = obs[
        ["sampleID", "atscore", "microglia_subpopulations", "state_group"]
    ].copy()
    obs_summary["TSPAN14_count"] = tspan
    state_summary = (
        obs_summary.groupby(["microglia_subpopulations", "atscore"], observed=True)
        .agg(
            n_nuclei=("TSPAN14_count", "size"),
            mean_count=("TSPAN14_count", "mean"),
            detection_fraction=("TSPAN14_count", lambda values: (values > 0).mean()),
        )
        .reset_index()
    )

    donor_counts.to_csv(
        args.output_dir / "gse243292_microglia_donor_pseudobulk_counts.tsv", sep="\t"
    )
    donor_metadata.to_csv(
        args.output_dir / "gse243292_microglia_donor_metadata.tsv",
        sep="\t",
        index=False,
    )
    state_counts.to_csv(
        args.output_dir / "gse243292_microglia_state_pseudobulk_counts.tsv", sep="\t"
    )
    state_metadata.to_csv(
        args.output_dir / "gse243292_microglia_state_metadata.tsv",
        sep="\t",
        index=False,
    )
    state_summary.to_csv(
        args.output_dir / "gse243292_tspan14_state_expression.tsv",
        sep="\t",
        index=False,
    )
    pd.DataFrame(
        [
            {
                "accession": "GSE243292",
                "resource": "Microglia processed H5AD",
                "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243292",
                "local_source": str(args.h5ad),
                "bytes": args.h5ad.stat().st_size,
                "sha256": sha256(args.h5ad),
                "access_date": date.today().isoformat(),
            }
        ]
    ).to_csv(
        args.output_dir / "gse243292_microglia_source_manifest.tsv",
        sep="\t",
        index=False,
    )
    audit = {
        "n_donors": int(obs["sampleID"].nunique()),
        "n_nuclei": int(dataset.n_obs),
        "n_genes": int(dataset.n_vars),
        "pathology_counts": {
            str(key): int(value)
            for key, value in donor_metadata["atscore"].value_counts().items()
        },
        "state_definition": {
            "Reactive_associated": sorted(REACTIVE_STATES),
            "Other_microglia": sorted(
                set(obs["microglia_subpopulations"].astype(str)) - REACTIVE_STATES
            ),
        },
        "minimum_state_cells": args.minimum_state_cells,
    }
    (args.output_dir / "gse243292_microglia_pseudobulk_audit.json").write_text(
        json.dumps(audit, indent=2), encoding="utf-8"
    )
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
