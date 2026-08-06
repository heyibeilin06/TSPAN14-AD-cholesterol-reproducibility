"""Prepare allele-aligned regional inputs for PC-GMM cis-MR and cis-MVMR."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def is_palindromic(a1: str, a2: str) -> bool:
    return {str(a1).upper(), str(a2).upper()} in ({"A", "T"}, {"C", "G"})


def allele_factor(source_a1: str, source_a2: str, target_a1: str, target_a2: str) -> int | None:
    source = (str(source_a1).upper(), str(source_a2).upper())
    target = (str(target_a1).upper(), str(target_a2).upper())
    if is_palindromic(*source) or is_palindromic(*target):
        return None
    if source == target:
        return 1
    if source == target[::-1]:
        return -1
    complement = (COMPLEMENT.get(source[0], ""), COMPLEMENT.get(source[1], ""))
    if complement == target:
        return 1
    if complement == target[::-1]:
        return -1
    return None


def read_selected(path: Path, snps: set[str]) -> pd.DataFrame:
    compression = "gzip" if path.suffix == ".gz" else "infer"
    chunks = []
    for chunk in pd.read_csv(path, sep="\t", compression=compression, chunksize=400_000):
        selected = chunk.loc[chunk["SNP"].astype(str).isin(snps)]
        if not selected.empty:
            chunks.append(selected)
    if not chunks:
        return pd.DataFrame(columns=["SNP", "A1", "A2", "BETA", "SE", "N"])
    return pd.concat(chunks, ignore_index=True).drop_duplicates("SNP")


def harmonize_to_ld(
    qtl: pd.DataFrame,
    outcomes: dict[str, pd.DataFrame],
    ld_meta: pd.DataFrame,
) -> pd.DataFrame:
    frame = qtl[["SNP", "A1", "A2", "b", "SE", "p"]].drop_duplicates("SNP").copy()
    frame = frame.merge(ld_meta[["SNP", "A1_LD", "A2_LD"]], on="SNP", how="inner")
    qtl_factor = [allele_factor(a1, a2, la1, la2) for a1, a2, la1, la2 in frame[["A1", "A2", "A1_LD", "A2_LD"]].itertuples(index=False, name=None)]
    frame["qtl_factor"] = qtl_factor
    frame = frame.loc[frame["qtl_factor"].notna()].copy()
    frame["beta_splice"] = frame["b"].astype(float) * frame["qtl_factor"].astype(float)
    frame["se_splice"] = frame["SE"].astype(float)

    for label, outcome in outcomes.items():
        columns = outcome[["SNP", "A1", "A2", "BETA", "SE", "N"]].drop_duplicates("SNP").copy()
        columns = columns.rename(columns={c: f"{c}_{label}" for c in ["A1", "A2", "BETA", "SE", "N"]})
        frame = frame.merge(columns, on="SNP", how="inner")
        factors = [
            allele_factor(a1, a2, la1, la2)
            for a1, a2, la1, la2 in frame[[f"A1_{label}", f"A2_{label}", "A1_LD", "A2_LD"]].itertuples(index=False, name=None)
        ]
        frame[f"factor_{label}"] = factors
        frame = frame.loc[frame[f"factor_{label}"].notna()].copy()
        frame[f"beta_{label}"] = frame[f"BETA_{label}"].astype(float) * frame[f"factor_{label}"].astype(float)
        frame[f"se_{label}"] = frame[f"SE_{label}"].astype(float)
        frame[f"n_{label}"] = pd.to_numeric(frame[f"N_{label}"], errors="coerce")

    keep = ["SNP", "A1_LD", "A2_LD", "p", "beta_splice", "se_splice"]
    for label in outcomes:
        keep.extend([f"beta_{label}", f"se_{label}", f"n_{label}"])
    return frame[keep].reset_index(drop=True)


def greedy_ld_prune(frame: pd.DataFrame, ld: np.ndarray, *, r2_threshold: float) -> list[int]:
    order = frame.sort_values("p").index.tolist()
    selected: list[int] = []
    for index in order:
        if all(float(ld[index, prior]) ** 2 < r2_threshold for prior in selected):
            selected.append(index)
    return sorted(selected)


def nearest_correlation(matrix: np.ndarray, floor: float = 1e-8) -> np.ndarray:
    symmetric = (matrix + matrix.T) / 2
    values, vectors = np.linalg.eigh(symmetric)
    repaired = (vectors * np.maximum(values, floor)) @ vectors.T
    scale = np.sqrt(np.diag(repaired))
    repaired = repaired / np.outer(scale, scale)
    return (repaired + repaired.T) / 2


def effect_geometry(splice: np.ndarray, lipid: np.ndarray) -> dict[str, float | bool]:
    effects = np.column_stack([np.asarray(splice, dtype=float), np.asarray(lipid, dtype=float)])
    correlation = float(np.corrcoef(effects.T)[0, 1])
    scaled = (effects - effects.mean(axis=0)) / effects.std(axis=0, ddof=1)
    singular_values = np.linalg.svd(scaled, compute_uv=False)
    ratio = float(singular_values[1] / singular_values[0]) if singular_values[0] else 0.0
    return {
        "effect_correlation": correlation,
        "absolute_effect_correlation": abs(correlation),
        "second_to_first_singular_value_ratio": ratio,
        "separable": bool(abs(correlation) < 0.95 and ratio > 0.1),
    }


def parse_label_path(value: str) -> tuple[str, Path]:
    label, separator, path = value.partition("=")
    if not separator:
        raise argparse.ArgumentTypeError("Use LABEL=PATH for --outcome.")
    return label, Path(path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--qtl", required=True, type=Path)
    parser.add_argument("--ld-store", required=True, type=Path)
    parser.add_argument("--outcome", action="append", required=True, type=parse_label_path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--max-r2", type=float, default=0.95)
    args = parser.parse_args()

    from magenpy import LDMatrix  # Imported only in the dedicated Python 3.12 runtime.

    args.output_dir.mkdir(parents=True, exist_ok=True)
    qtl = pd.read_csv(args.qtl, sep="\t")
    store = LDMatrix.from_path(str(args.ld_store))
    metadata = store.to_snp_table().copy()
    metadata = metadata.rename(columns={"A1": "A1_LD", "A2": "A2_LD"})
    candidate_snps = set(qtl["SNP"].astype(str)) & set(metadata["SNP"].astype(str))
    outcomes = {label: read_selected(path, candidate_snps) for label, path in args.outcome}
    harmonized = harmonize_to_ld(qtl, outcomes, metadata)

    lookup = {str(snp): index for index, snp in enumerate(np.asarray(store.snps).astype(str))}
    global_indices = np.asarray([lookup[snp] for snp in harmonized["SNP"]], dtype=int)
    start, end = int(global_indices.min()), int(global_indices.max()) + 1
    regional = store.load_data(
        start_row=start,
        end_row=end,
        dtype="float64",
        return_square=True,
        return_symmetric=True,
        return_as_csr=True,
    ).toarray()
    ld = regional[np.ix_(global_indices - start, global_indices - start)]
    retained = greedy_ld_prune(harmonized, ld, r2_threshold=args.max_r2)
    harmonized = harmonized.iloc[retained].reset_index(drop=True)
    ld = ld[np.ix_(retained, retained)]
    minimum_eigenvalue_before = float(np.linalg.eigvalsh((ld + ld.T) / 2).min())
    repaired = nearest_correlation(ld)
    minimum_eigenvalue_after = float(np.linalg.eigvalsh(repaired).min())
    relative_repair = float(np.linalg.norm(repaired - ld, ord="fro") / np.linalg.norm(ld, ord="fro"))

    harmonized.to_csv(args.output_dir / "01_pc_mvmr_harmonized.tsv", sep="\t", index=False)
    pd.DataFrame(repaired, index=harmonized["SNP"], columns=harmonized["SNP"]).to_csv(
        args.output_dir / "02_pc_mvmr_ld.tsv", sep="\t"
    )
    geometry = []
    for label in outcomes:
        if label == "AD":
            continue
        geometry.append({"lipid": label, **effect_geometry(harmonized["beta_splice"], harmonized[f"beta_{label}"])})
    pd.DataFrame(geometry).to_csv(args.output_dir / "03_exposure_effect_geometry.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "qtl_variants": len(qtl),
                "ld_overlap": len(candidate_snps),
                "complete_allele_harmonized": len(retained),
                "maximum_pairwise_r2": args.max_r2,
                "minimum_ld_eigenvalue_before_repair": minimum_eigenvalue_before,
                "minimum_ld_eigenvalue_after_repair": minimum_eigenvalue_after,
                "relative_frobenius_ld_repair": relative_repair,
            }
        ]
    ).to_csv(args.output_dir / "00_pc_mvmr_input_qc.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
