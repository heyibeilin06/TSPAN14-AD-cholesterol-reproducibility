#!/usr/bin/env python3
"""Targeted, hypothesis-led audit of TSPAN14 expression context.

This analysis is deliberately limited to pre-specified genes and modules linked
to the TSPAN14-ADAM10 axis or membrane-lipid biology. It is an exploratory
context analysis, not evidence for splice-event mediation.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.io import mmread


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs" / "manuscript_refinement_v2"
OUT.mkdir(parents=True, exist_ok=True)

GENES = ["ADAM10", "TREM2", "APP", "NOTCH1", "APOE", "LPL", "ABCA1", "TYROBP", "PLCG2"]
MODULES = {
    "ADAM10_substrate_context": ["ADAM10", "TREM2", "APP", "NOTCH1"],
    "membrane_lipid_context": ["APOE", "LPL", "ABCA1"],
    "microglial_signalling_context": ["TYROBP", "PLCG2"],
}


def bh_adjust(pvalues):
    p = np.asarray(pvalues, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.minimum.accumulate((ranked * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    out = np.empty_like(adjusted)
    out[order] = np.minimum(adjusted, 1.0)
    return out


def residualize(y, covariates):
    x = pd.get_dummies(covariates, drop_first=True, dtype=float)
    x = x.replace([np.inf, -np.inf], np.nan)
    keep = x.notna().all(axis=1) & np.isfinite(y)
    x = x.loc[keep]
    y = np.asarray(y)[keep.to_numpy()]
    x = np.column_stack([np.ones(len(x)), x.to_numpy(dtype=float)])
    beta = np.linalg.lstsq(x, y, rcond=None)[0]
    residuals = y - x @ beta
    return residuals, keep


def log_cpm(counts):
    lib = counts.sum(axis=1)
    return np.log2((counts / np.maximum(lib[:, None], 1)) * 1e6 + 0.5)


def analyse_dataset(name, expression, metadata, covariate_columns, cell_context):
    present = [g for g in ["TSPAN14", *GENES] if g in expression.columns]
    expression = expression[present].copy()
    covariates = metadata[covariate_columns].copy()
    results = []

    t_resid, keep = residualize(expression["TSPAN14"].to_numpy(), covariates)
    kept_expression = expression.loc[keep].copy()
    kept_covariates = covariates.loc[keep]

    for gene in GENES:
        if gene not in kept_expression:
            continue
        y_resid, gene_keep = residualize(kept_expression[gene].to_numpy(), kept_covariates)
        if not gene_keep.all():
            continue
        rho, p = stats.spearmanr(t_resid, y_resid)
        results.append({
            "dataset": name,
            "cell_context": cell_context,
            "feature_type": "gene",
            "feature": gene,
            "n": len(t_resid),
            "partial_spearman_rho": rho,
            "p_value": p,
        })

    z = kept_expression.apply(lambda x: (x - x.mean()) / x.std(ddof=1), axis=0)
    for module, genes in MODULES.items():
        available = [g for g in genes if g in z and z[g].notna().all()]
        if len(available) < 2:
            continue
        score = z[available].mean(axis=1).to_numpy()
        score_resid, score_keep = residualize(score, kept_covariates)
        if not score_keep.all():
            continue
        rho, p = stats.spearmanr(t_resid, score_resid)
        results.append({
            "dataset": name,
            "cell_context": cell_context,
            "feature_type": "module",
            "feature": module,
            "n": len(t_resid),
            "partial_spearman_rho": rho,
            "p_value": p,
        })

    result = pd.DataFrame(results)
    result["fdr_within_dataset"] = bh_adjust(result["p_value"])
    return result


def load_seaad():
    base = ROOT / "outputs" / "p1_cell_context"
    matrix = mmread(base / "p1_seaad_AnG_pseudobulk_umi_counts.mtx")
    if sparse.issparse(matrix):
        matrix = matrix.tocsr()
    genes = pd.read_csv(base / "p1_seaad_AnG_pseudobulk_genes.tsv", sep="\t")["gene"].tolist()
    metadata = pd.read_csv(base / "p1_seaad_AnG_pseudobulk_metadata.tsv", sep="\t")
    index = {gene: i for i, gene in enumerate(genes)}
    selected = [g for g in ["TSPAN14", *GENES] if g in index]
    cols = [index[g] for g in selected]
    values = matrix[:, cols].toarray() if sparse.issparse(matrix) else np.asarray(matrix)[:, cols]
    expression = pd.DataFrame(log_cpm(values), columns=selected)
    metadata = metadata.reset_index(drop=True)
    metadata["high_ad"] = (metadata["adnc"] == "High").astype(int)
    metadata["sex_male"] = (metadata["sex"] == "Male").astype(int)
    metadata["log_median_umis"] = np.log1p(metadata["median_umis"])
    metadata["log_n_cells"] = np.log1p(metadata["n_cells"])
    return expression, metadata


def load_gse243292():
    base = ROOT / "outputs" / "p1_cell_context_v16"
    counts = pd.read_csv(base / "gse243292_microglia_donor_pseudobulk_counts.tsv", sep="\t").set_index("gene")
    metadata = pd.read_csv(base / "gse243292_microglia_donor_metadata.tsv", sep="\t")
    selected = [g for g in ["TSPAN14", *GENES] if g in counts.index]
    values = counts.loc[selected, metadata["sample"].astype(str)].T.to_numpy(dtype=float)
    expression = pd.DataFrame(log_cpm(values), columns=selected)
    metadata = metadata.reset_index(drop=True)
    metadata["log_median_umis"] = np.log1p(metadata["median_umis"])
    metadata["log_n_cells"] = np.log1p(metadata["n_cells"])
    return expression, metadata


def main():
    sea_expr, sea_meta = load_seaad()
    outputs = []
    sea_covariates = ["high_ad", "age_at_death", "sex_male", "rin", "log_median_umis", "log_n_cells"]
    for cell_group in ["microglia", "neurons"]:
        mask = sea_meta["cell_group"] == cell_group
        outputs.append(analyse_dataset(
            "SEA-AD", sea_expr.loc[mask].reset_index(drop=True), sea_meta.loc[mask].reset_index(drop=True),
            sea_covariates, cell_group,
        ))

    gse_expr, gse_meta = load_gse243292()
    outputs.append(analyse_dataset(
        "GSE243292", gse_expr, gse_meta,
        ["pathology_stage", "apoe4_dosage", "trem2_r47h", "log_median_umis", "log_n_cells"],
        "microglia",
    ))

    result = pd.concat(outputs, ignore_index=True)
    result["fdr_global"] = bh_adjust(result["p_value"])
    result.to_csv(OUT / "01_tspan14_targeted_context_correlations.tsv", sep="\t", index=False)

    modules = result[result["feature_type"] == "module"].copy()
    replicated = []
    for feature, frame in modules[modules["cell_context"] == "microglia"].groupby("feature"):
        if set(frame["dataset"]) >= {"SEA-AD", "GSE243292"}:
            same_sign = len(set(np.sign(frame["partial_spearman_rho"]))) == 1
            replicated.append({
                "feature": feature,
                "same_direction_across_microglia_datasets": same_sign,
                "both_nominal_p_lt_0_05": bool((frame["p_value"] < 0.05).all()),
                "both_fdr_lt_0_05": bool((frame["fdr_within_dataset"] < 0.05).all()),
            })
    pd.DataFrame(replicated).to_csv(OUT / "02_module_replication_audit.tsv", sep="\t", index=False)

    summary = [
        "# Targeted TSPAN14 cell-context audit",
        "",
        "This pre-specified analysis tested residual expression covariance, not differential expression and not splice-event mediation.",
        "Associations were adjusted for disease/pathology and technical covariates within each dataset.",
        "Only findings replicated in direction and significance across both independent microglial datasets should be considered supportive.",
        "",
        f"Tests performed: {len(result)}; globally FDR-significant: {(result['fdr_global'] < 0.05).sum()}.",
        f"Module tests replicated at nominal P<0.05 in both microglial datasets: {sum(x['both_nominal_p_lt_0_05'] for x in replicated)}.",
    ]
    (OUT / "CELL_CONTEXT_TARGETED_AUDIT.md").write_text("\n".join(summary) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
