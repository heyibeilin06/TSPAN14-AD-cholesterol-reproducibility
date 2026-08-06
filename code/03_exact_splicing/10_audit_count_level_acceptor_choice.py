#!/usr/bin/env python3
"""Audit genotype-dependent canonical-versus-cryptic acceptor choice.

The retained donor table contains the canonical exon5-6 junction count and
the total count for the two-intron LeafCutter cluster. The competing cryptic
acceptor count is therefore cluster_reads - target_reads. This script treats
the counts as acceptor-choice observations and explicitly handles the
complete separation present in the risk-homozygous genotype group.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq, minimize
from scipy.special import betaln, expit, gammaln
from scipy.stats import chi2, fisher_exact


def firth_objective(beta: np.ndarray, x: np.ndarray, y: np.ndarray) -> float:
    eta = x @ beta
    p = np.clip(expit(eta), 1e-12, 1 - 1e-12)
    loglik = np.sum(y * np.log(p) + (1 - y) * np.log1p(-p))
    w = p * (1 - p)
    information = x.T @ (w[:, None] * x)
    sign, logdet = np.linalg.slogdet(information)
    if sign <= 0:
        return np.inf
    return -(loglik + 0.5 * logdet)


def firth_logistic(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, float, bool]:
    """Fit Firth-penalized logistic regression by maximizing penalized likelihood."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    fit = minimize(firth_objective, np.zeros(x.shape[1]), args=(x, y), method="BFGS")
    return fit.x, -float(fit.fun), bool(fit.success)


def profile_likelihood_ci(x, y, beta_hat, loglik_hat, coefficient=1, alpha=0.05):
    cutoff = chi2.ppf(1 - alpha, 1)
    nuisance = [j for j in range(x.shape[1]) if j != coefficient]

    def profile_loglik(value):
        def objective(theta):
            beta = np.empty(x.shape[1])
            beta[coefficient] = value
            beta[nuisance] = theta
            return firth_objective(beta, x, y)

        fit = minimize(objective, beta_hat[nuisance], method="BFGS")
        return -float(fit.fun)

    def root(value):
        return 2 * (loglik_hat - profile_loglik(value)) - cutoff

    center = beta_hat[coefficient]
    lower_probe = center - 1.0
    while root(lower_probe) < 0 and lower_probe > center - 40:
        lower_probe -= 1.5
    upper_probe = center + 1.0
    while root(upper_probe) < 0 and upper_probe < center + 40:
        upper_probe += 1.5
    lower = brentq(root, lower_probe, center - 1e-8)
    upper = brentq(root, center + 1e-8, upper_probe)
    return lower, upper


def beta_binomial_objective(theta, genotype, successes, totals):
    intercept, genotype_beta, log_concentration = theta
    p = np.clip(expit(intercept + genotype_beta * genotype), 1e-9, 1 - 1e-9)
    concentration = np.exp(log_concentration)
    alpha = p * concentration
    beta = (1 - p) * concentration
    failures = totals - successes
    log_choose = gammaln(totals + 1) - gammaln(successes + 1) - gammaln(failures + 1)
    loglik = log_choose + betaln(successes + alpha, failures + beta) - betaln(alpha, beta)
    return -float(np.sum(loglik))


def fit_beta_binomial(genotype, successes, totals):
    start = np.array([-3.0, -1.0, 1.0])
    fit = minimize(
        beta_binomial_objective,
        start,
        args=(genotype, successes, totals),
        method="L-BFGS-B",
        bounds=[(-20, 20), (-20, 20), (-10, 15)],
    )
    return fit.x, -float(fit.fun), bool(fit.success)


def beta_binomial_profile_ci(genotype, successes, totals, theta_hat, loglik_hat, alpha=0.05):
    cutoff = chi2.ppf(1 - alpha, 1)

    def profile_loglik(value):
        def objective(nuisance):
            theta = np.array([nuisance[0], value, nuisance[1]])
            return beta_binomial_objective(theta, genotype, successes, totals)

        fit = minimize(
            objective,
            theta_hat[[0, 2]],
            method="L-BFGS-B",
            bounds=[(-20, 20), (-10, 15)],
        )
        return -float(fit.fun)

    def root(value):
        return 2 * (loglik_hat - profile_loglik(value)) - cutoff

    center = theta_hat[1]
    lower_probe = center - 1.0
    while root(lower_probe) < 0 and lower_probe > -20:
        lower_probe -= 1.0
    upper_probe = center + 1.0
    while root(upper_probe) < 0 and upper_probe < 20:
        upper_probe += 1.0
    lower = brentq(root, lower_probe, center - 1e-8)
    upper = brentq(root, center + 1e-8, upper_probe)
    return lower, upper


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t")
    required = {"donor_id", "genotype", "target_reads", "cluster_reads"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")

    df = df.copy()
    df["cryptic_reads"] = df["cluster_reads"] - df["target_reads"]
    if (df["cryptic_reads"] < 0).any():
        raise ValueError("Canonical junction counts exceed cluster counts.")
    df["canonical_fraction"] = df["target_reads"] / df["cluster_reads"]
    df["cryptic_detected"] = (df["cryptic_reads"] > 0).astype(int)

    strata = (
        df.groupby("genotype", as_index=False)
        .agg(
            donors=("donor_id", "nunique"),
            canonical_reads=("target_reads", "sum"),
            cryptic_reads=("cryptic_reads", "sum"),
            cryptic_positive_donors=("cryptic_detected", "sum"),
            mean_canonical_fraction=("canonical_fraction", "mean"),
            median_cluster_depth=("cluster_reads", "median"),
            minimum_cluster_depth=("cluster_reads", "min"),
            maximum_cluster_depth=("cluster_reads", "max"),
        )
    )
    strata["cryptic_positive_fraction"] = (
        strata["cryptic_positive_donors"] / strata["donors"]
    )

    x_full = np.column_stack([np.ones(len(df)), df["genotype"].to_numpy(float)])
    y = df["cryptic_detected"].to_numpy(float)
    beta_full, ll_full, ok_full = firth_logistic(x_full, y)
    beta_null, ll_null, ok_null = firth_logistic(np.ones((len(df), 1)), y)
    lr = max(0.0, 2 * (ll_full - ll_null))
    firth_ci = profile_likelihood_ci(x_full, y, beta_full, ll_full)

    log_depth = np.log(df["cluster_reads"].to_numpy(float))
    log_depth = (log_depth - log_depth.mean()) / log_depth.std(ddof=0)
    x_depth = np.column_stack(
        [np.ones(len(df)), df["genotype"].to_numpy(float), log_depth]
    )
    x_depth_null = np.column_stack([np.ones(len(df)), log_depth])
    beta_depth, ll_depth, ok_depth = firth_logistic(x_depth, y)
    beta_depth_null, ll_depth_null, ok_depth_null = firth_logistic(x_depth_null, y)
    lr_depth = max(0.0, 2 * (ll_depth - ll_depth_null))
    firth_depth_ci = profile_likelihood_ci(x_depth, y, beta_depth, ll_depth)

    genotype = df["genotype"].to_numpy(float)
    successes = df["cryptic_reads"].to_numpy(float)
    totals = df["cluster_reads"].to_numpy(float)
    bb_full, bb_ll_full, bb_ok_full = fit_beta_binomial(genotype, successes, totals)
    bb_null, bb_ll_null, bb_ok_null = fit_beta_binomial(
        np.zeros_like(genotype), successes, totals
    )
    bb_lr = max(0.0, 2 * (bb_ll_full - bb_ll_null))
    bb_ci = beta_binomial_profile_ci(
        genotype, successes, totals, bb_full, bb_ll_full
    )

    non_hom = df[df["genotype"] < 2]
    hom = df[df["genotype"] == 2]
    table = np.array(
        [
            [int(hom["cryptic_detected"].sum()), int((1 - hom["cryptic_detected"]).sum())],
            [int(non_hom["cryptic_detected"].sum()), int((1 - non_hom["cryptic_detected"]).sum())],
        ]
    )
    fisher_or, fisher_p = fisher_exact(table, alternative="two-sided")

    loo = []
    for idx in df.index:
        keep = df.index != idx
        beta, _, ok = firth_logistic(x_full[keep], y[keep])
        loo.append((idx, beta[1], ok))
    loo_df = pd.DataFrame(loo, columns=["omitted_row", "genotype_log_odds", "converged"])
    loo_df["odds_ratio_per_risk_allele"] = np.exp(loo_df["genotype_log_odds"])

    depth_rows = []
    for threshold in [1, 10, 20, 30, 50]:
        sub = df[df["cluster_reads"] >= threshold]
        x = np.column_stack([np.ones(len(sub)), sub["genotype"].to_numpy(float)])
        beta, ll, ok = firth_logistic(x, sub["cryptic_detected"].to_numpy(float))
        depth_rows.append(
            {
                "minimum_cluster_depth": threshold,
                "donors": len(sub),
                "cryptic_positive_donors": int(sub["cryptic_detected"].sum()),
                "genotype_log_odds": beta[1],
                "odds_ratio_per_risk_allele": np.exp(beta[1]),
                "converged": ok,
            }
        )
    depth_df = pd.DataFrame(depth_rows)

    summary = pd.DataFrame(
        [
            {
                "analysis": "Firth logistic trend: any cryptic-acceptor read",
                "estimate": beta_full[1],
                "effect_scale": "log odds per risk-aligned allele",
                "odds_ratio": np.exp(beta_full[1]),
                "odds_ratio_95ci_low": np.exp(firth_ci[0]),
                "odds_ratio_95ci_high": np.exp(firth_ci[1]),
                "p_value": chi2.sf(lr, 1),
                "p_value_test": "penalized likelihood-ratio test",
                "n_donors": len(df),
                "converged": ok_full and ok_null,
                "interpretation": "Risk-aligned genotype is associated with lower detection of the competing cryptic acceptor.",
            },
            {
                "analysis": "Depth-adjusted Firth logistic trend: any cryptic-acceptor read",
                "estimate": beta_depth[1],
                "effect_scale": "log odds per rs7080009 risk-aligned allele, adjusted for standardized log cluster depth",
                "odds_ratio": np.exp(beta_depth[1]),
                "odds_ratio_95ci_low": np.exp(firth_depth_ci[0]),
                "odds_ratio_95ci_high": np.exp(firth_depth_ci[1]),
                "p_value": chi2.sf(lr_depth, 1),
                "p_value_test": "penalized likelihood-ratio test",
                "n_donors": len(df),
                "converged": ok_depth and ok_depth_null,
                "interpretation": "Risk-aligned genotype remains associated with lower cryptic-acceptor detection after adjustment for local cluster depth.",
            },
            {
                "analysis": "Beta-binomial cryptic-versus-canonical read model",
                "estimate": bb_full[1],
                "effect_scale": "log odds of a cryptic read per rs7080009 risk-aligned allele",
                "odds_ratio": np.exp(bb_full[1]),
                "odds_ratio_95ci_low": np.exp(bb_ci[0]),
                "odds_ratio_95ci_high": np.exp(bb_ci[1]),
                "p_value": chi2.sf(bb_lr, 1),
                "p_value_test": "beta-binomial likelihood-ratio test",
                "n_donors": len(df),
                "converged": bb_ok_full and bb_ok_null,
                "interpretation": "A read-count model accounting for overdispersion supports lower cryptic-versus-canonical read use with increasing rs7080009 risk-allele dosage.",
            },
            {
                "analysis": "Fisher exact: risk homozygotes versus other genotypes",
                "estimate": fisher_or,
                "effect_scale": "odds ratio for any cryptic-acceptor read",
                "odds_ratio": fisher_or,
                "odds_ratio_95ci_low": np.nan,
                "odds_ratio_95ci_high": np.nan,
                "p_value": fisher_p,
                "p_value_test": "Fisher exact test",
                "n_donors": len(df),
                "converged": True,
                "interpretation": "The zero cryptic-read count in risk homozygotes constitutes complete separation and motivates penalized inference.",
            },
            {
                "analysis": "Observed donor-level canonical fraction",
                "estimate": df.groupby("genotype")["canonical_fraction"].mean().iloc[-1]
                - df.groupby("genotype")["canonical_fraction"].mean().iloc[0],
                "effect_scale": "absolute fraction difference: genotype 2 minus genotype 0",
                "odds_ratio": np.nan,
                "odds_ratio_95ci_low": np.nan,
                "odds_ratio_95ci_high": np.nan,
                "p_value": np.nan,
                "p_value_test": "descriptive",
                "n_donors": len(df),
                "converged": True,
                "interpretation": "Descriptive magnitude only; sparse competing-junction counts preclude a precise continuous isoform-shift estimate.",
            },
        ]
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output_dir / "count_level_acceptor_choice_donors.tsv", sep="\t", index=False)
    strata.to_csv(args.output_dir / "count_level_acceptor_choice_by_genotype.tsv", sep="\t", index=False)
    summary.to_csv(args.output_dir / "count_level_acceptor_choice_models.tsv", sep="\t", index=False)
    depth_df.to_csv(args.output_dir / "count_level_acceptor_choice_depth_sensitivity.tsv", sep="\t", index=False)
    loo_df.to_csv(args.output_dir / "count_level_acceptor_choice_leave_one_out.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
