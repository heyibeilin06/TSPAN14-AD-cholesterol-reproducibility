"""Run bidirectional, multivariable, and mediation-oriented MR audits.

The genome-wide component uses one or more LD-pruned instruments from European
GWAS summary statistics. The TSPAN14 component is explicitly cis-local and
tests whether splice and lipid effects can be separated within the locus.
"""

from __future__ import annotations

import argparse
import gzip
import math
from pathlib import Path
from statistics import NormalDist

import numpy as np
import pandas as pd
from scipy.stats import chi2


NORMAL = NormalDist()
PALINDROMIC = {frozenset(("A", "T")), frozenset(("C", "G"))}
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def p_value(z_value: float) -> float:
    return 2 * NORMAL.cdf(-abs(float(z_value)))


def harmonize(exposure: pd.DataFrame, outcome: pd.DataFrame) -> pd.DataFrame:
    merged = exposure.merge(outcome, on="SNP", suffixes=("_x", "_y"), how="inner")
    rows = []
    for row in merged.itertuples(index=False):
        a1x, a2x = str(row.A1_x).upper(), str(row.A2_x).upper()
        a1y, a2y = str(row.A1_y).upper(), str(row.A2_y).upper()
        if frozenset((a1x, a2x)) in PALINDROMIC:
            continue
        factor = None
        if (a1x, a2x) == (a1y, a2y):
            factor = 1
        elif (a1x, a2x) == (a2y, a1y):
            factor = -1
        elif (a1x, a2x) == (COMPLEMENT.get(a1y), COMPLEMENT.get(a2y)):
            factor = 1
        elif (a1x, a2x) == (COMPLEMENT.get(a2y), COMPLEMENT.get(a1y)):
            factor = -1
        if factor is None:
            continue
        rows.append(
            {
                "SNP": row.SNP,
                "CHR": int(row.CHR_x),
                "BP": int(row.BP_x),
                "effect_allele": a1x,
                "other_allele": a2x,
                "beta_x": float(row.BETA_x),
                "se_x": float(row.SE_x),
                "p_x": float(row.P_x),
                "eaf_x": float(row.EAF_x),
                "n_x": float(row.N_x),
                "beta_y": float(row.BETA_y) * factor,
                "se_y": float(row.SE_y),
                "p_y": float(row.P_y),
                "eaf_y": float(row.EAF_y) if pd.notna(row.EAF_y) else np.nan,
                "n_y": float(row.N_y),
            }
        )
    return pd.DataFrame(rows)


def ivw_mr(beta_x: np.ndarray, beta_y: np.ndarray, se_y: np.ndarray) -> dict[str, float]:
    weights = 1 / np.square(se_y)
    denominator = float(np.sum(weights * np.square(beta_x)))
    estimate = float(np.sum(weights * beta_x * beta_y) / denominator)
    fixed_se = math.sqrt(1 / denominator)
    residual = beta_y - estimate * beta_x
    q_statistic = float(np.sum(weights * np.square(residual)))
    q_df = max(1, len(beta_x) - 1)
    phi = max(1.0, q_statistic / q_df)
    se = fixed_se * math.sqrt(phi)
    return {
        "estimate": estimate,
        "se": se,
        "pvalue": p_value(estimate / se),
        "q_statistic": q_statistic,
        "q_df": q_df,
        "q_pvalue": float(chi2.sf(q_statistic, q_df)),
        "random_effect_scale": phi,
    }


def weighted_median_mr(beta_x: np.ndarray, beta_y: np.ndarray, se_y: np.ndarray) -> float:
    ratio = beta_y / beta_x
    weights = np.square(beta_x) / np.square(se_y)
    order = np.argsort(ratio)
    cumulative = np.cumsum(weights[order]) / weights.sum()
    return float(ratio[order][np.searchsorted(cumulative, 0.5)])


def egger_mr(beta_x: np.ndarray, beta_y: np.ndarray, se_y: np.ndarray) -> dict[str, float]:
    design = np.column_stack((np.ones(len(beta_x)), beta_x))
    weights = np.diag(1 / np.square(se_y))
    covariance = np.linalg.pinv(design.T @ weights @ design)
    coefficients = covariance @ design.T @ weights @ beta_y
    residual = beta_y - design @ coefficients
    scale = max(1.0, float(residual.T @ weights @ residual) / max(1, len(beta_x) - 2))
    covariance *= scale
    intercept, slope = map(float, coefficients)
    intercept_se, slope_se = np.sqrt(np.diag(covariance))
    return {
        "egger_intercept": intercept,
        "egger_intercept_se": float(intercept_se),
        "egger_intercept_pvalue": p_value(intercept / intercept_se),
        "egger_estimate": slope,
        "egger_se": float(slope_se),
        "egger_pvalue": p_value(slope / slope_se),
    }


def multivariable_ivw(
    exposure_effects: np.ndarray,
    outcome_effects: np.ndarray,
    outcome_se: np.ndarray,
    covariance: np.ndarray | None = None,
) -> dict[str, object]:
    if covariance is None:
        covariance = np.diag(np.square(outcome_se))
    inverse = np.linalg.pinv(covariance)
    information = exposure_effects.T @ inverse @ exposure_effects
    coefficient_covariance = np.linalg.pinv(information)
    coefficients = coefficient_covariance @ exposure_effects.T @ inverse @ outcome_effects
    residual = outcome_effects - exposure_effects @ coefficients
    q_statistic = float(residual.T @ inverse @ residual)
    q_df = max(1, len(outcome_effects) - exposure_effects.shape[1])
    correlation = float(np.corrcoef(exposure_effects.T)[0, 1])
    standardized = exposure_effects / np.maximum(
        np.linalg.norm(exposure_effects, axis=0, keepdims=True), 1e-15
    )
    condition_number = float(np.linalg.cond(standardized))
    identifiable = bool(
        np.linalg.matrix_rank(exposure_effects) == exposure_effects.shape[1]
        and abs(correlation) < 0.9
        and condition_number < 10
    )
    standard_errors = np.sqrt(np.diag(coefficient_covariance))
    return {
        "coefficients": coefficients,
        "standard_errors": standard_errors,
        "pvalues": np.array(
            [p_value(beta / se) for beta, se in zip(coefficients, standard_errors)]
        ),
        "q_statistic": q_statistic,
        "q_df": q_df,
        "q_pvalue": float(chi2.sf(q_statistic, q_df)),
        "exposure_effect_correlation": correlation,
        "condition_number": condition_number,
        "matrix_rank": int(np.linalg.matrix_rank(exposure_effects)),
        "identifiable": identifiable,
    }


def read_significant(path: Path, threshold: float = 5e-8) -> pd.DataFrame:
    parts = []
    for chunk in pd.read_csv(path, sep="\t", compression="gzip", chunksize=400_000):
        selected = chunk.loc[
            (pd.to_numeric(chunk["P"], errors="coerce") <= threshold)
            & chunk["SNP"].astype(str).str.startswith("rs")
        ]
        if not selected.empty:
            parts.append(selected)
    if not parts:
        return pd.DataFrame(columns=["SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "EAF", "N"])
    return pd.concat(parts, ignore_index=True).drop_duplicates("SNP").sort_values("P")


def read_selected(path: Path, snps: set[str]) -> pd.DataFrame:
    parts = []
    for chunk in pd.read_csv(path, sep="\t", compression="gzip", chunksize=400_000):
        selected = chunk.loc[chunk["SNP"].isin(snps)]
        if not selected.empty:
            parts.append(selected)
    if not parts:
        return pd.DataFrame(columns=["SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "EAF", "N"])
    return pd.concat(parts, ignore_index=True).drop_duplicates("SNP")


def ld_clump_by_blocks(frame: pd.DataFrame, ld_root: Path, r2_threshold: float) -> pd.DataFrame:
    from magenpy import LDMatrix  # imported only in the Python 3.12 analysis runtime

    retained = []
    for chromosome, candidates in frame.groupby("CHR"):
        path = ld_root / f"chr_{int(chromosome)}"
        if not path.exists():
            continue
        store = LDMatrix.from_path(str(path))
        snps = np.asarray(store.snps).astype(str)
        lookup = {snp: i for i, snp in enumerate(snps)}
        candidates = candidates.loc[candidates["SNP"].isin(lookup)].sort_values("P")
        if candidates.empty:
            continue
        blocks = store.zarr_group.attrs["Estimator properties"]["LD blocks"]
        bp = np.asarray(store.bp_position)
        block_members: dict[int, list[tuple[pd.Series, int]]] = {}
        for _, row in candidates.iterrows():
            index = lookup[row["SNP"]]
            block_id = next(
                i for i, (start, end) in enumerate(blocks) if start <= bp[index] < end
            )
            block_members.setdefault(block_id, []).append((row, index))
        for block_id, members in block_members.items():
            start_bp, end_bp = blocks[block_id]
            global_indices = np.where((bp >= start_bp) & (bp < end_bp))[0]
            start, end = int(global_indices[0]), int(global_indices[-1]) + 1
            matrix = store.load_data(
                start_row=start,
                end_row=end,
                dtype="float64",
                return_square=True,
                return_symmetric=True,
                return_as_csr=True,
            ).toarray()
            selected_indices = []
            for row, global_index in sorted(members, key=lambda item: float(item[0]["P"])):
                local_index = global_index - start
                if all(matrix[local_index, prior] ** 2 < r2_threshold for prior in selected_indices):
                    selected_indices.append(local_index)
                    retained.append(row)
    if not retained:
        return frame.iloc[0:0].copy()
    return pd.DataFrame(retained).drop_duplicates("SNP").sort_values(["CHR", "BP"])


def summarize_univariable(exposure: str, outcome: str, harmonized: pd.DataFrame, analysis: str) -> dict:
    bx = harmonized["beta_x"].to_numpy(float)
    by = harmonized["beta_y"].to_numpy(float)
    sey = harmonized["se_y"].to_numpy(float)
    if len(harmonized) == 1:
        estimate = float(by[0] / bx[0])
        se = float(abs(sey[0] / bx[0]))
        return {
            "analysis": analysis,
            "exposure": exposure,
            "outcome": outcome,
            "n_instruments": 1,
            "method": "Wald ratio",
            "estimate": estimate,
            "se": se,
            "pvalue": p_value(estimate / se),
            "q_statistic": np.nan,
            "q_pvalue": np.nan,
            "weighted_median": np.nan,
            "egger_estimate": np.nan,
            "egger_pvalue": np.nan,
            "egger_intercept_pvalue": np.nan,
        }
    ivw = ivw_mr(bx, by, sey)
    egger = egger_mr(bx, by, sey) if len(harmonized) >= 3 else {}
    return {
        "analysis": analysis,
        "exposure": exposure,
        "outcome": outcome,
        "n_instruments": len(harmonized),
        "method": "multiplicative random-effects IVW",
        **ivw,
        "weighted_median": weighted_median_mr(bx, by, sey),
        "egger_estimate": egger.get("egger_estimate", np.nan),
        "egger_pvalue": egger.get("egger_pvalue", np.nan),
        "egger_intercept_pvalue": egger.get("egger_intercept_pvalue", np.nan),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--ld-root", type=Path, required=True)
    parser.add_argument("--qtl", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--r2-threshold", type=float, default=0.001)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    paths = {
        "AD": args.data_root / "ad_bellenguez_2022_harmonized.tsv.gz",
        "TC": args.data_root / "glgc_tc_eur_harmonized.tsv.gz",
        "LDL": args.data_root / "glgc_ldl_eur_harmonized.tsv.gz",
        "nonHDL": args.data_root / "glgc_nonhdl_eur_harmonized.tsv.gz",
        "HDL": args.data_root / "glgc_hdl_eur_harmonized.tsv.gz",
        "TG": args.data_root / "glgc_tg_eur_harmonized.tsv.gz",
    }
    instrument_sets = {}
    selection_rows = []
    for trait, path in paths.items():
        significant = read_significant(path)
        clumped = ld_clump_by_blocks(significant, args.ld_root, args.r2_threshold)
        instrument_sets[trait] = clumped
        selection_rows.append(
            {
                "trait": trait,
                "genome_wide_significant_before_ld": len(significant),
                "ukb_ld_available_and_clumped": len(clumped),
                "p_threshold": 5e-8,
                "r2_threshold": args.r2_threshold,
                "ld_reference": "UK Biobank European LD matrix",
            }
        )
    pd.DataFrame(selection_rows).to_csv(
        args.output_dir / "01_genomewide_instrument_selection.tsv", sep="\t", index=False
    )

    union = set().union(*(set(frame["SNP"]) for frame in instrument_sets.values()))
    lookup = {trait: read_selected(path, union) for trait, path in paths.items()}
    results = []
    harmonized_parts = []
    for exposure in paths:
        for outcome in paths:
            if exposure == outcome or (exposure != "AD" and outcome != "AD"):
                continue
            h = harmonize(instrument_sets[exposure], lookup[outcome])
            h["exposure"] = exposure
            h["outcome"] = outcome
            harmonized_parts.append(h)
            if not h.empty:
                results.append(summarize_univariable(exposure, outcome, h, "genomewide_bidirectional"))
    pd.concat(harmonized_parts, ignore_index=True).to_csv(
        args.output_dir / "02_genomewide_bidirectional_harmonized.tsv", sep="\t", index=False
    )
    pd.DataFrame(results).to_csv(
        args.output_dir / "03_genomewide_bidirectional_mr.tsv", sep="\t", index=False
    )

    # APOE-region exclusion is a prespecified sensitivity analysis for all genome-wide models.
    no_apoe_results = []
    for exposure in paths:
        for outcome in paths:
            if exposure == outcome or (exposure != "AD" and outcome != "AD"):
                continue
            instruments = instrument_sets[exposure].loc[
                ~(
                    instrument_sets[exposure]["CHR"].eq(19)
                    & instrument_sets[exposure]["BP"].between(40_000_000, 50_000_000)
                )
            ]
            h = harmonize(instruments, lookup[outcome])
            if not h.empty:
                no_apoe_results.append(summarize_univariable(exposure, outcome, h, "genomewide_APOE_40_50Mb_excluded"))
    pd.DataFrame(no_apoe_results).to_csv(
        args.output_dir / "04_genomewide_bidirectional_apoe_sensitivity.tsv", sep="\t", index=False
    )

    # Global lipid MVMR uses the union of instruments for the conventional LDL/HDL/TG model.
    lipid_model = ("LDL", "HDL", "TG")
    lipid_union = set().union(*(set(instrument_sets[x]["SNP"]) for x in lipid_model))
    base = lookup["AD"].loc[lookup["AD"]["SNP"].isin(lipid_union)].copy()
    exposure_frames = []
    for trait in lipid_model:
        frame = lookup[trait].loc[lookup[trait]["SNP"].isin(lipid_union)].copy()
        frame = frame.rename(columns={column: f"{column}_{trait}" for column in frame.columns if column != "SNP"})
        exposure_frames.append(frame)
    wide = base.rename(columns={column: f"{column}_AD" for column in base.columns if column != "SNP"})
    for frame in exposure_frames:
        wide = wide.merge(frame, on="SNP", how="inner")
    aligned_rows = []
    for row in wide.itertuples(index=False):
        record = {"SNP": row.SNP}
        ad_a1, ad_a2 = row.A1_AD, row.A2_AD
        if frozenset((ad_a1, ad_a2)) in PALINDROMIC:
            continue
        valid = True
        for trait in lipid_model:
            a1, a2 = getattr(row, f"A1_{trait}"), getattr(row, f"A2_{trait}")
            if (a1, a2) == (ad_a1, ad_a2):
                factor = 1
            elif (a1, a2) == (ad_a2, ad_a1):
                factor = -1
            else:
                valid = False
                break
            record[f"beta_{trait}"] = getattr(row, f"BETA_{trait}") * factor
        if valid:
            record.update(beta_AD=row.BETA_AD, se_AD=row.SE_AD)
            aligned_rows.append(record)
    global_mvmr_data = pd.DataFrame(aligned_rows)
    global_mvmr_data.to_csv(args.output_dir / "05_global_lipid_mvmr_harmonized.tsv", sep="\t", index=False)
    if len(global_mvmr_data) >= 4:
        fit = multivariable_ivw(
            global_mvmr_data[[f"beta_{x}" for x in lipid_model]].to_numpy(float),
            global_mvmr_data["beta_AD"].to_numpy(float),
            global_mvmr_data["se_AD"].to_numpy(float),
        )
        rows = []
        for i, trait in enumerate(lipid_model):
            rows.append(
                {
                    "model": "LDL+HDL+TG -> AD",
                    "exposure": trait,
                    "n_instruments": len(global_mvmr_data),
                    "direct_estimate": fit["coefficients"][i],
                    "se": fit["standard_errors"][i],
                    "pvalue": fit["pvalues"][i],
                    "q_statistic": fit["q_statistic"],
                    "q_pvalue": fit["q_pvalue"],
                    "condition_number": fit["condition_number"],
                    "matrix_rank": fit["matrix_rank"],
                }
            )
        pd.DataFrame(rows).to_csv(args.output_dir / "06_global_lipid_mvmr.tsv", sep="\t", index=False)

    # Local exact-splice reciprocal and MVMR analyses.
    qtl = pd.read_csv(args.qtl, sep="\t").drop_duplicates("SNP")
    local_snps = set(qtl["SNP"])
    local = {trait: read_selected(path, local_snps) for trait, path in paths.items()}
    qtl_std = qtl.rename(columns={"b": "BETA", "p": "P", "Freq": "EAF"}).copy()
    qtl_std["EAF"] = pd.to_numeric(qtl_std["EAF"], errors="coerce").fillna(0.25)
    qtl_std["N"] = 147
    qtl_std = qtl_std[["SNP", "Chr", "BP", "A1", "A2", "BETA", "SE", "P", "EAF", "N"]].rename(columns={"Chr": "CHR"})
    local_results = []
    local_harmonized = []
    strict_qtl = ld_clump_by_blocks(qtl_std.loc[qtl_std["P"] <= 5e-8], args.ld_root, 0.1)
    for outcome in paths:
        h = harmonize(strict_qtl, local[outcome])
        if not h.empty:
            h["exposure"], h["outcome"] = "exact_exon5_6_sQTL", outcome
            local_harmonized.append(h)
            local_results.append(summarize_univariable("exact_exon5_6_sQTL", outcome, h, "strict_local_forward"))
    for exposure in paths:
        local_exposure = ld_clump_by_blocks(
            local[exposure].loc[local[exposure]["P"] <= 5e-8], args.ld_root, 0.1
        )
        h = harmonize(local_exposure, qtl_std)
        if not h.empty:
            h["exposure"], h["outcome"] = exposure, "exact_exon5_6_sQTL"
            local_harmonized.append(h)
            local_results.append(summarize_univariable(exposure, "exact_exon5_6_sQTL", h, "locus_restricted_reverse"))
    pd.concat(local_harmonized, ignore_index=True).to_csv(
        args.output_dir / "07_local_reciprocal_harmonized.tsv", sep="\t", index=False
    )
    pd.DataFrame(local_results).to_csv(
        args.output_dir / "08_local_reciprocal_mr.tsv", sep="\t", index=False
    )

    relaxed = ld_clump_by_blocks(qtl_std.loc[qtl_std["P"] <= 1e-3], args.ld_root, 0.1)
    local_mvmr_rows = []
    local_mvmr_data = []
    for lipid in ("TC", "LDL", "nonHDL"):
        hs = harmonize(relaxed, local[lipid])
        ha = harmonize(relaxed, local["AD"])
        merged = hs[["SNP", "beta_x", "se_x", "beta_y", "se_y"]].rename(
            columns={"beta_x": "splice_beta", "se_x": "splice_se", "beta_y": "lipid_beta", "se_y": "lipid_se"}
        ).merge(
            ha[["SNP", "beta_y", "se_y"]].rename(columns={"beta_y": "ad_beta", "se_y": "ad_se"}),
            on="SNP",
        )
        merged["lipid"] = lipid
        local_mvmr_data.append(merged)
        if len(merged) < 3:
            continue
        fit = multivariable_ivw(
            merged[["splice_beta", "lipid_beta"]].to_numpy(float),
            merged["ad_beta"].to_numpy(float),
            merged["ad_se"].to_numpy(float),
        )
        for i, exposure in enumerate(("exact_exon5_6_sQTL", lipid)):
            local_mvmr_rows.append(
                {
                    "model": f"exact_exon5_6_sQTL + {lipid} -> AD",
                    "exposure": exposure,
                    "n_instruments": len(merged),
                    "instrument_selection": "sQTL P<=1e-3; UKB EUR LD r2<0.1",
                    "direct_estimate": fit["coefficients"][i],
                    "se": fit["standard_errors"][i],
                    "pvalue": fit["pvalues"][i],
                    "exposure_effect_correlation": fit["exposure_effect_correlation"],
                    "condition_number": fit["condition_number"],
                    "matrix_rank": fit["matrix_rank"],
                    "q_statistic": fit["q_statistic"],
                    "q_pvalue": fit["q_pvalue"],
                    "identifiable": fit["identifiable"],
                }
            )
    pd.concat(local_mvmr_data, ignore_index=True).to_csv(
        args.output_dir / "09_local_mvmr_harmonized.tsv", sep="\t", index=False
    )
    pd.DataFrame(local_mvmr_rows).to_csv(
        args.output_dir / "10_local_splice_lipid_mvmr.tsv", sep="\t", index=False
    )

    coverage = pd.DataFrame(
        [
            {
                "analysis": "genome-wide lipid <-> AD bidirectional MR",
                "status": "completed",
                "scope": "Five lipid traits; IVW, weighted median, MR-Egger, heterogeneity, APOE exclusion",
            },
            {
                "analysis": "exact sQTL -> AD/lipids",
                "status": "completed_strict_independent_and_existing_LD_aware",
                "scope": "Strict r2<0.1 local instruments plus prior correlated-instrument cis-MR",
            },
            {
                "analysis": "AD/lipids -> exact sQTL",
                "status": "completed_locus_restricted_only",
                "scope": "Cannot represent genome-wide liability because released sQTL associations are cis-only",
            },
            {
                "analysis": "exact sQTL + lipid -> AD MVMR",
                "status": "fitted_with_identifiability_diagnostics",
                "scope": "Relaxed sQTL instruments used only to test separability; inference conditional on diagnostics",
            },
            {
                "analysis": "two-step lipid -> sQTL -> AD mediation",
                "status": "not_point_identified",
                "scope": "Genome-wide lipid instruments lack released trans associations with the BA24 junction; local shared instruments violate pathway separation",
            },
        ]
    )
    coverage.to_csv(args.output_dir / "11_mr_mediation_coverage_decision.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
