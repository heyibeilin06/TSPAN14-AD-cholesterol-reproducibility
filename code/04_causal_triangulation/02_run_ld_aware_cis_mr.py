"""Run LD-aware generalized IVW and MR-Egger cis-MR for an exact sQTL feature."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from statistics import NormalDist

import numpy as np
import pandas as pd
from scipy.stats import chi2


NORMAL = NormalDist()


def ld_clump(variants: list[str], ld: np.ndarray, *, r2_threshold: float) -> list[str]:
    selected: list[int] = []
    for index, _ in enumerate(variants):
        if all(float(ld[index, chosen]) ** 2 < r2_threshold for chosen in selected):
            selected.append(index)
    return [variants[index] for index in selected]


def orient_ld(ld: np.ndarray, signs: np.ndarray) -> np.ndarray:
    return ld * np.outer(signs, signs)


def generalized_ivw(beta_x: np.ndarray, se_x: np.ndarray, beta_y: np.ndarray, se_y: np.ndarray, ld: np.ndarray) -> dict[str, float]:
    estimate = float(np.median(beta_y / beta_x))
    for _ in range(20):
        covariance = ld * (np.outer(se_y, se_y) + estimate**2 * np.outer(se_x, se_x))
        covariance = (covariance + covariance.T) / 2 + np.eye(len(beta_x)) * 1e-12
        inverse = np.linalg.pinv(covariance)
        denominator = float(beta_x @ inverse @ beta_x)
        updated = float((beta_x @ inverse @ beta_y) / denominator)
        if abs(updated - estimate) < 1e-12:
            estimate = updated
            break
        estimate = updated
    covariance = ld * (np.outer(se_y, se_y) + estimate**2 * np.outer(se_x, se_x))
    covariance = (covariance + covariance.T) / 2 + np.eye(len(beta_x)) * 1e-12
    inverse = np.linalg.pinv(covariance)
    denominator = float(beta_x @ inverse @ beta_x)
    standard_error = math.sqrt(1 / denominator)
    z_score = estimate / standard_error
    residual = beta_y - estimate * beta_x
    q_statistic = float(residual @ inverse @ residual)
    degrees_freedom = max(1, len(beta_x) - 1)
    return {
        "estimate": estimate,
        "se": standard_error,
        "pvalue": 2 * NORMAL.cdf(-abs(z_score)),
        "q_statistic": q_statistic,
        "q_df": degrees_freedom,
        "q_pvalue": float(chi2.sf(q_statistic, degrees_freedom)),
    }


def generalized_egger(beta_x: np.ndarray, se_x: np.ndarray, beta_y: np.ndarray, se_y: np.ndarray, ld: np.ndarray) -> dict[str, float]:
    slope = float(np.median(beta_y / beta_x))
    design = np.column_stack([np.ones(len(beta_x)), beta_x])
    for _ in range(20):
        covariance = ld * (np.outer(se_y, se_y) + slope**2 * np.outer(se_x, se_x))
        covariance = (covariance + covariance.T) / 2 + np.eye(len(beta_x)) * 1e-12
        inverse = np.linalg.pinv(covariance)
        covariance_coefficients = np.linalg.pinv(design.T @ inverse @ design)
        coefficients = covariance_coefficients @ design.T @ inverse @ beta_y
        if abs(float(coefficients[1]) - slope) < 1e-12:
            break
        slope = float(coefficients[1])
    intercept, slope = (float(coefficients[0]), float(coefficients[1]))
    intercept_se, slope_se = (math.sqrt(float(covariance_coefficients[0, 0])), math.sqrt(float(covariance_coefficients[1, 1])))
    return {
        "egger_intercept": intercept,
        "egger_intercept_se": intercept_se,
        "egger_intercept_pvalue": 2 * NORMAL.cdf(-abs(intercept / intercept_se)),
        "egger_slope": slope,
        "egger_slope_se": slope_se,
        "egger_slope_pvalue": 2 * NORMAL.cdf(-abs(slope / slope_se)),
    }


def alignment_factor(qtl_a1: str, qtl_a2: str, outcome_a1: str, outcome_a2: str) -> int | None:
    if qtl_a1 == outcome_a1 and qtl_a2 == outcome_a2:
        return 1
    if qtl_a1 == outcome_a2 and qtl_a2 == outcome_a1:
        return -1
    return None


def parse_outcome(value: str) -> tuple[str, Path]:
    label, separator, path = value.partition("=")
    if not separator:
        raise argparse.ArgumentTypeError("Each --outcome must use LABEL=PATH.")
    return label, Path(path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--qtl", required=True)
    parser.add_argument("--ld", required=True)
    parser.add_argument("--allele-orientation", required=True, help="ALT-dosage to sQTL-effect-allele signs for selected instruments.")
    parser.add_argument("--outcome", action="append", type=parse_outcome, required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--r2-threshold", type=float, default=0.6)
    parser.add_argument("--instrument-p-threshold", type=float, default=5e-8)
    args = parser.parse_args()

    qtl = pd.read_csv(args.qtl, sep="\t")
    ld = pd.read_csv(args.ld, sep="\t", index_col=0)
    orientation = pd.read_csv(args.allele_orientation, sep="\t")
    strong = qtl[(qtl["p"] <= args.instrument_p_threshold) & qtl["SNP"].isin(ld.index)].sort_values("p").drop_duplicates("SNP")
    initial_ld = ld.loc[strong["SNP"], strong["SNP"]].to_numpy(dtype=float)
    selected_names = ld_clump(strong["SNP"].tolist(), initial_ld, r2_threshold=args.r2_threshold)
    selected_qtl = strong.set_index("SNP").loc[selected_names].reset_index()
    selected_qtl = selected_qtl.merge(orientation[["SNP", "alt_to_sqtl_effect_sign"]], on="SNP", how="left")
    if selected_qtl["alt_to_sqtl_effect_sign"].isna().any():
        raise ValueError("Every LD-aware instrument requires a verified ALT-to-sQTL-effect sign.")
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    selected_qtl.to_csv(output_dir / "p1_exact_ba24_ld_aware_cis_mr_instruments.tsv", sep="\t", index=False)

    results: list[dict[str, object]] = []
    for trait, path in args.outcome:
        outcome = pd.read_csv(path, sep="\t")
        frame = selected_qtl.merge(outcome, left_on="SNP", right_on="snp", how="inner")
        frame["factor"] = [alignment_factor(*row) for row in frame[["A1", "A2", "ad_a1", "ad_a2"]].itertuples(index=False, name=None)]
        frame = frame[frame["factor"].notna()].copy()
        frame = frame[~((frame["A1"].isin(["A", "T"])) & (frame["A2"].isin(["A", "T"])))].copy()
        frame = frame[~((frame["A1"].isin(["C", "G"])) & (frame["A2"].isin(["C", "G"])))].copy()
        if len(frame) < 3:
            raise ValueError(f"Fewer than three non-palindromic LD-aware instruments remain for {trait}.")
        if trait == "AD":
            beta_y = frame["ad_beta"].to_numpy(dtype=float) * frame["factor"].to_numpy(dtype=float)
            se_y = np.sqrt(frame["ad_varbeta"].to_numpy(dtype=float))
        else:
            beta_y = frame["trait_beta_aligned_to_ad_a1"].to_numpy(dtype=float) * frame["factor"].to_numpy(dtype=float)
            se_y = np.sqrt(frame["trait_varbeta"].to_numpy(dtype=float))
        frame = frame.set_index("SNP").loc[[name for name in selected_names if name in frame["SNP"].tolist()]].reset_index()
        beta_x = frame["b"].to_numpy(dtype=float)
        selected_ld = orient_ld(
            ld.loc[frame["SNP"], frame["SNP"]].to_numpy(dtype=float),
            frame["alt_to_sqtl_effect_sign"].to_numpy(dtype=float),
        )
        if trait == "AD":
            beta_y = frame["ad_beta"].to_numpy(dtype=float) * frame["factor"].to_numpy(dtype=float)
            se_y = np.sqrt(frame["ad_varbeta"].to_numpy(dtype=float))
        else:
            beta_y = frame["trait_beta_aligned_to_ad_a1"].to_numpy(dtype=float) * frame["factor"].to_numpy(dtype=float)
            se_y = np.sqrt(frame["trait_varbeta"].to_numpy(dtype=float))
        se_x = frame["SE"].to_numpy(dtype=float)
        ivw = generalized_ivw(beta_x, se_x, beta_y, se_y, selected_ld)
        egger = generalized_egger(beta_x, se_x, beta_y, se_y, selected_ld)
        results.append(
            {
                "outcome": trait,
                "n_ld_aware_instruments": len(frame),
                "instrument_p_threshold": args.instrument_p_threshold,
                "pairwise_ld_r2_threshold": args.r2_threshold,
                "minimum_instrument_f": float(np.min((beta_x / frame["SE"].to_numpy(dtype=float)) ** 2)),
                "method": "LD-aware generalized IVW cis-MR",
                **ivw,
                **egger,
                "interpretation": "LD-aware cis-MR estimate; correlated instruments were retained only with covariance-matrix correction.",
                "causal_boundary": "Does not identify a lipid-to-AD mediation effect because sQTL-to-lipid and sQTL-to-AD estimates share one cis configuration.",
            }
        )
    pd.DataFrame(results).to_csv(output_dir / "p1_exact_ba24_ld_aware_cis_mr.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
