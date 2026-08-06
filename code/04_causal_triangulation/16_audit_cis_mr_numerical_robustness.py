#!/usr/bin/env python3
"""Audit LD-aware cis-MR against LD regularization and instrument redundancy."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from statistics import NormalDist

import numpy as np
import pandas as pd


NORMAL = NormalDist()


def generalized_ivw(bx, sx, by, sy, r, ridge=0.0, eigen_floor=None):
    estimate = float(np.median(by / bx))
    for _ in range(100):
        covariance = r * (np.outer(sy, sy) + estimate**2 * np.outer(sx, sx))
        covariance = (covariance + covariance.T) / 2
        if ridge:
            covariance += np.eye(len(bx)) * ridge * float(np.mean(np.diag(covariance)))
        values, vectors = np.linalg.eigh(covariance)
        if eigen_floor is not None:
            floor = eigen_floor * float(values.max())
            values = np.maximum(values, floor)
        inverse = vectors @ np.diag(1 / values) @ vectors.T
        denominator = float(bx @ inverse @ bx)
        updated = float((bx @ inverse @ by) / denominator)
        if abs(updated - estimate) < 1e-12:
            estimate = updated
            break
        estimate = updated
    se = math.sqrt(1 / denominator)
    return estimate, se, 2 * NORMAL.cdf(-abs(estimate / se))


def greedy_clump(order, r, threshold):
    selected = []
    for i in order:
        if all(r[i, j] ** 2 < threshold for j in selected):
            selected.append(i)
    return selected


def nearest_correlation(matrix, floor=1e-6):
    """Project a symmetric perturbation to a positive-definite correlation matrix."""
    matrix = (matrix + matrix.T) / 2
    values, vectors = np.linalg.eigh(matrix)
    matrix = vectors @ np.diag(np.maximum(values, floor)) @ vectors.T
    scale = np.sqrt(np.diag(matrix))
    matrix = matrix / np.outer(scale, scale)
    np.fill_diagonal(matrix, 1.0)
    return matrix


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--atlas", type=Path, required=True)
    parser.add_argument("--ld", type=Path, required=True)
    parser.add_argument("--bim", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    atlas = pd.read_csv(args.atlas, sep="\t")
    bim = pd.read_csv(
        args.bim,
        sep=r"\s+",
        header=None,
        names=["chromosome", "plink_id", "cm", "BP", "plink_a1", "plink_a2"],
    )
    r = np.loadtxt(args.ld)
    qtl = atlas[atlas["layer"] == "Exact exon5-6 sQTL"].drop_duplicates("SNP")
    mapping = bim.merge(
        qtl[["SNP", "BP", "sqtl_effect_allele", "sqtl_other_allele", "orientation", "p_value"]],
        on="BP",
        how="left",
    )
    if mapping["SNP"].isna().any():
        raise ValueError("At least one PLINK variant could not be mapped by GRCh37 position.")
    mapping["risk_aligned_allele"] = np.where(
        mapping["orientation"] == 1,
        mapping["sqtl_effect_allele"],
        mapping["sqtl_other_allele"],
    )
    mapping["plink_to_risk_sign"] = np.where(
        mapping["plink_a1"] == mapping["risk_aligned_allele"], 1, -1
    )
    if not np.all(
        (mapping["plink_a1"] == mapping["risk_aligned_allele"])
        | (mapping["plink_a2"] == mapping["risk_aligned_allele"])
    ):
        raise ValueError("Risk-aligned allele is absent from a PLINK allele pair.")
    signed_r = r * np.outer(mapping["plink_to_risk_sign"], mapping["plink_to_risk_sign"])

    eigenvalues = np.linalg.eigvalsh(signed_r)
    proportions = eigenvalues / eigenvalues.sum()
    diagnostics = pd.DataFrame(
        [
            {"metric": "minimum_eigenvalue", "value": eigenvalues.min()},
            {"metric": "maximum_eigenvalue", "value": eigenvalues.max()},
            {"metric": "condition_number", "value": np.linalg.cond(signed_r)},
            {"metric": "effective_rank_entropy", "value": np.exp(-np.sum(proportions * np.log(proportions)))},
            {"metric": "effective_rank_eigenvalue_gt_0.1", "value": np.sum(eigenvalues > 0.1)},
        ]
    )

    results = []
    perturbation_results = []
    rng = np.random.default_rng(20260807)
    order = np.argsort(mapping["p_value"].to_numpy())
    for outcome in ["AD", "TC", "LDL", "nonHDL"]:
        exposure = atlas[atlas["layer"] == "Exact exon5-6 sQTL"].set_index("SNP")
        target = atlas[atlas["outcome"] == outcome].set_index("SNP")
        names = mapping["SNP"].tolist()
        bx = exposure.loc[names, "risk_aligned_beta"].to_numpy(float)
        sx = exposure.loc[names, "se"].to_numpy(float)
        by = target.loc[names, "risk_aligned_beta"].to_numpy(float)
        sy = target.loc[names, "se"].to_numpy(float)

        for ridge in [0.0, 1e-8, 1e-6, 1e-4, 1e-2, 1e-1]:
            est, se, p = generalized_ivw(bx, sx, by, sy, signed_r, ridge=ridge)
            results.append({"outcome": outcome, "analysis": "ridge", "parameter": ridge, "n_instruments": len(names), "estimate": est, "se": se, "p_value": p})
        for floor in [1e-8, 1e-6, 1e-4, 1e-2, 5e-2, 1e-1]:
            est, se, p = generalized_ivw(bx, sx, by, sy, signed_r, eigen_floor=floor)
            results.append({"outcome": outcome, "analysis": "relative_eigenvalue_floor", "parameter": floor, "n_instruments": len(names), "estimate": est, "se": se, "p_value": p})
        for shrinkage in [0.01, 0.05, 0.1, 0.2, 0.5]:
            shrunk = (1 - shrinkage) * signed_r + shrinkage * np.eye(len(names))
            est, se, p = generalized_ivw(bx, sx, by, sy, shrunk)
            results.append({"outcome": outcome, "analysis": "LD_shrinkage_to_identity", "parameter": shrinkage, "n_instruments": len(names), "estimate": est, "se": se, "p_value": p})
        for threshold in [0.2, 0.4, 0.6, 0.8]:
            selected = greedy_clump(order, signed_r, threshold)
            est, se, p = generalized_ivw(bx[selected], sx[selected], by[selected], sy[selected], signed_r[np.ix_(selected, selected)])
            results.append({"outcome": outcome, "analysis": "pairwise_r2_clump", "parameter": threshold, "n_instruments": len(selected), "estimate": est, "se": se, "p_value": p})
        for omitted in range(len(names)):
            selected = [i for i in range(len(names)) if i != omitted]
            est, se, p = generalized_ivw(bx[selected], sx[selected], by[selected], sy[selected], signed_r[np.ix_(selected, selected)])
            results.append({"outcome": outcome, "analysis": "leave_one_variant_out", "parameter": names[omitted], "n_instruments": len(selected), "estimate": est, "se": se, "p_value": p})
        lead = int(np.argmin(exposure.loc[names, "p_value"].to_numpy(float)))
        ratio = by[lead] / bx[lead]
        ratio_se = math.sqrt(sy[lead] ** 2 / bx[lead] ** 2 + (by[lead] ** 2 * sx[lead] ** 2) / bx[lead] ** 4)
        results.append({"outcome": outcome, "analysis": "lead_variant_Wald", "parameter": names[lead], "n_instruments": 1, "estimate": ratio, "se": ratio_se, "p_value": 2 * NORMAL.cdf(-abs(ratio / ratio_se))})

        # Quantify sensitivity to plausible uncertainty in the external LD panel.
        # Noise is added to off-diagonal correlations and each draw is projected
        # back to the space of positive-definite correlation matrices.
        for noise_sd in [0.02, 0.05]:
            draw_estimates = []
            draw_ses = []
            for _ in range(1000):
                noise = rng.normal(0, noise_sd, size=signed_r.shape)
                noise = np.triu(noise, 1)
                noise = noise + noise.T
                perturbed = nearest_correlation(signed_r + noise)
                est, se, _ = generalized_ivw(bx, sx, by, sy, perturbed)
                draw_estimates.append(est)
                draw_ses.append(se)
            estimates = np.asarray(draw_estimates)
            ses = np.asarray(draw_ses)
            perturbation_results.append(
                {
                    "outcome": outcome,
                    "off_diagonal_noise_sd": noise_sd,
                    "n_draws": len(estimates),
                    "estimate_median": np.median(estimates),
                    "estimate_p025": np.quantile(estimates, 0.025),
                    "estimate_p975": np.quantile(estimates, 0.975),
                    "se_median": np.median(ses),
                    "proportion_positive": np.mean(estimates > 0),
                }
            )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    mapping.to_csv(args.output_dir / "cis_mr_ld_allele_mapping.tsv", sep="\t", index=False)
    pd.DataFrame(signed_r, index=mapping["SNP"], columns=mapping["SNP"]).to_csv(args.output_dir / "cis_mr_signed_ld_matrix.tsv", sep="\t")
    diagnostics.to_csv(args.output_dir / "cis_mr_ld_diagnostics.tsv", sep="\t", index=False)
    pd.DataFrame(results).to_csv(args.output_dir / "cis_mr_numerical_sensitivity.tsv", sep="\t", index=False)
    pd.DataFrame(perturbation_results).to_csv(
        args.output_dir / "cis_mr_ld_perturbation_sensitivity.tsv", sep="\t", index=False
    )


if __name__ == "__main__":
    main()
