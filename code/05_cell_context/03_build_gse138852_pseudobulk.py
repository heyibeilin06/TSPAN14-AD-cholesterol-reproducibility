#!/usr/bin/env python3
"""Build pool-level pseudobulk counts for the public GSE138852 snRNA-seq cohort.

The released object pools two donors per library and therefore cannot support
donor-level inference. We aggregate only prespecified broad cell types and keep
the paired library block for a conservative sensitivity analysis.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd


CELL_TYPE_MAP = {
    "mg": "Microglia",
    "neuron": "Neurons",
    "oligo": "Oligodendrocytes",
    "astro": "Astrocytes",
    "OPC": "Oligodendrocyte_precursor_cells",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def pool_from_cell(cell: str) -> str:
    match = re.search(r"_((?:AD|Ct)\d+_(?:AD|Ct)\d+)$", cell)
    if not match:
        raise ValueError(f"Cannot recover pooled library from cell ID: {cell}")
    return match.group(1)


def pair_from_pool(pool: str) -> str:
    match = re.search(r"(\d+)_", pool)
    if not match:
        raise ValueError(f"Cannot recover paired block from pool: {pool}")
    return {"1": "pair_1", "3": "pair_2", "5": "pair_3"}[match.group(1)]


def load_metadata(path: Path, minimum_cells: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    metadata = pd.read_csv(path).rename(
        columns={
            "Unnamed: 0": "cell",
            "oupSample.batchCond": "condition_raw",
            "oupSample.cellType": "cell_type_raw",
        }
    )
    metadata["pool"] = metadata["cell"].map(pool_from_cell)
    metadata["pair"] = metadata["pool"].map(pair_from_pool)
    metadata["condition"] = metadata["condition_raw"].map(
        {"AD": "AD", "ct": "Control"}
    )
    metadata["cell_type"] = metadata["cell_type_raw"].map(CELL_TYPE_MAP)
    metadata = metadata[metadata["cell_type"].notna()].copy()
    metadata["sample"] = metadata["cell_type"] + "__" + metadata["pool"]

    sample_metadata = (
        metadata.groupby(
            ["sample", "cell_type", "pool", "pair", "condition"], as_index=False
        )
        .size()
        .rename(columns={"size": "n_cells"})
    )
    valid = sample_metadata.groupby("cell_type").filter(
        lambda group: (
            group["condition"].value_counts().min() >= 3
            and group["n_cells"].min() >= minimum_cells
        )
    )
    keep_samples = set(valid["sample"])
    metadata = metadata[metadata["sample"].isin(keep_samples)].copy()
    sample_metadata = valid.sort_values(["cell_type", "pair", "condition"])
    return metadata, sample_metadata


def aggregate_counts(
    count_path: Path, metadata: pd.DataFrame, sample_metadata: pd.DataFrame
) -> pd.DataFrame:
    samples = sample_metadata["sample"].tolist()
    sample_index = {sample: index for index, sample in enumerate(samples)}

    with gzip.open(count_path, "rt", encoding="utf-8", newline="") as handle:
        header = handle.readline().rstrip("\r\n")
        cells = [value.strip('"') for value in header.split(",")[1:]]
        cell_to_position = {cell: position for position, cell in enumerate(cells)}
        missing = sorted(set(metadata["cell"]) - set(cell_to_position))
        if missing:
            raise ValueError(f"{len(missing)} metadata cells are absent from count matrix")

        selected_positions = np.array(
            [cell_to_position[cell] for cell in metadata["cell"]], dtype=np.int64
        )
        selected_groups = np.array(
            [sample_index[sample] for sample in metadata["sample"]], dtype=np.int64
        )

        genes: list[str] = []
        rows: list[np.ndarray] = []
        for line_number, line in enumerate(handle, start=2):
            comma = line.find(",")
            if comma < 0:
                raise ValueError(f"Malformed count row at line {line_number}")
            gene = line[:comma].strip('"')
            values = np.fromstring(line[comma + 1 :], sep=",", dtype=np.int64)
            if len(values) != len(cells):
                raise ValueError(
                    f"Count width mismatch for {gene}: {len(values)} != {len(cells)}"
                )
            aggregate = np.bincount(
                selected_groups,
                weights=values[selected_positions],
                minlength=len(samples),
            ).astype(np.int64)
            genes.append(gene)
            rows.append(aggregate)

    frame = pd.DataFrame(np.vstack(rows), index=genes, columns=samples)
    frame.index.name = "gene"
    return frame


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--counts", type=Path, required=True)
    parser.add_argument("--covariates", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--minimum-cells", type=int, default=20)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    metadata, sample_metadata = load_metadata(args.covariates, args.minimum_cells)
    pseudobulk = aggregate_counts(args.counts, metadata, sample_metadata)

    pseudobulk.to_csv(
        args.output_dir / "gse138852_pool_pseudobulk_counts.tsv", sep="\t"
    )
    sample_metadata.to_csv(
        args.output_dir / "gse138852_pool_pseudobulk_metadata.tsv",
        sep="\t",
        index=False,
    )
    manifest = pd.DataFrame(
        [
            {
                "accession": "GSE138852",
                "file_role": "raw gene-by-nucleus count matrix",
                "local_source": str(args.counts),
                "bytes": args.counts.stat().st_size,
                "sha256": sha256(args.counts),
                "access_date": date.today().isoformat(),
            },
            {
                "accession": "GSE138852",
                "file_role": "released cell annotations and pooled condition labels",
                "local_source": str(args.covariates),
                "bytes": args.covariates.stat().st_size,
                "sha256": sha256(args.covariates),
                "access_date": date.today().isoformat(),
            },
        ]
    )
    manifest.to_csv(
        args.output_dir / "gse138852_source_manifest.tsv", sep="\t", index=False
    )
    audit = {
        "n_cells_retained": int(len(metadata)),
        "n_genes": int(pseudobulk.shape[0]),
        "n_pseudobulk_samples": int(pseudobulk.shape[1]),
        "cell_types": sorted(sample_metadata["cell_type"].unique().tolist()),
        "minimum_cells_per_pool": args.minimum_cells,
        "inference_level": "pooled library; two donors per released pool",
    }
    (args.output_dir / "gse138852_pool_pseudobulk_audit.json").write_text(
        json.dumps(audit, indent=2), encoding="utf-8"
    )
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
