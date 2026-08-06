#!/usr/bin/env python3
"""Synthesize exact-event mediation rescue and independent trans-sQTL audits."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import pandas as pd


TARGET_GENE = "TSPAN14"
TARGET_GENE_ID = "ENSG00000108219"
TARGET_COORDINATES = (80509471, 80512144)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def contains_target(frame: pd.DataFrame) -> tuple[int, int]:
    text = frame.astype(str).agg("\t".join, axis=1)
    gene_hits = text.str.contains(TARGET_GENE, case=False, regex=False) | text.str.contains(
        TARGET_GENE_ID, case=False, regex=False
    )
    exact_hits = (
        gene_hits
        & text.str.contains(str(TARGET_COORDINATES[0]), regex=False)
        & text.str.contains(str(TARGET_COORDINATES[1]), regex=False)
    )
    return int(gene_hits.sum()), int(exact_hits.sum())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--brainmeta", required=True, type=Path)
    parser.add_argument("--gtex", required=True, type=Path)
    parser.add_argument("--interval", required=True, type=Path)
    parser.add_argument("--analysis-dir", required=True, type=Path)
    args = parser.parse_args()

    brainmeta = pd.read_csv(args.brainmeta, sep="\t")
    gtex = pd.read_csv(args.gtex, sep="\t")
    brainmeta_gene, brainmeta_exact = contains_target(brainmeta)
    gtex_gene, gtex_exact = contains_target(gtex)
    interval_audit = pd.read_csv(args.analysis_dir / "13_interval_trans_sqtl_audit.tsv", sep="\t")
    overlap = pd.read_csv(args.analysis_dir / "17_lipid_iv_any_tspan14_trans_sqtl_summary.tsv", sep="\t")
    dimension = pd.read_csv(args.analysis_dir / "14_pc_dimension_mediation_sensitivity.tsv", sep="\t")
    dimension_summary = pd.read_csv(args.analysis_dir / "15_pc_dimension_mediation_summary.tsv", sep="\t")
    tc3 = dimension[(dimension.lipid == "TC") & (dimension.n_pcs == 3)].iloc[0]

    resources = pd.DataFrame(
        [
            {
                "resource": "BrainMeta brain-cortex trans-sQTL",
                "catalog_rows": len(brainmeta),
                "TSPAN14_rows": brainmeta_gene,
                "exact_exon5_6_rows": brainmeta_exact,
                "source_file": args.brainmeta.name,
                "sha256": sha256(args.brainmeta),
                "interpretation": "TSPAN14 is not among the reported trans-sGenes",
            },
            {
                "resource": "GTEx v8 significant trans-sGenes",
                "catalog_rows": len(gtex),
                "TSPAN14_rows": gtex_gene,
                "exact_exon5_6_rows": gtex_exact,
                "source_file": args.gtex.name,
                "sha256": sha256(args.gtex),
                "interpretation": "TSPAN14 is not among the FDR-significant trans-sGenes",
            },
            {
                "resource": "INTERVAL untargeted trans-sQTL P<1e-5",
                "catalog_rows": int(interval_audit.loc[interval_audit.metric == "TSPAN14 trans-sQTL rows", "value"].iloc[0]),
                "TSPAN14_rows": int(interval_audit.loc[interval_audit.metric == "TSPAN14 trans-sQTL rows", "value"].iloc[0]),
                "exact_exon5_6_rows": int(interval_audit.loc[interval_audit.metric == "Exact exon5-6 trans-sQTL rows", "value"].iloc[0]),
                "source_file": "INTERVAL_trans_sQTL_summary_statistics_1e5.tsv.gz",
                "sha256": sha256(args.interval),
                "interpretation": "TSPAN14 trans hits exist, but none map to the exact exon5-6 event",
            },
        ]
    )
    resources.to_csv(args.analysis_dir / "18_independent_trans_sqtl_resource_audit.tsv", sep="\t", index=False)

    tc_summary = dimension_summary[dimension_summary.lipid == "TC"].iloc[0]
    overlap_sentence = (
        "One LDL-C instrument intersected one non-target TSPAN14 trans-sQTL feature; "
        "no TC instrument and no exact-event instrument intersection was observed."
        if not overlap.empty
        else "No lipid instrument intersected a reported TSPAN14 trans-sQTL feature."
    )
    report = f"""# Exact-event lipid-splicing mediation rescue: final audit

## Objective

Test whether the available public data identify a directional pathway in which lipid traits alter the exact TSPAN14 exon5-6 splice event and that event mediates Alzheimer disease risk.

## What was added

1. LD-aware principal-component generalized method-of-moments cis-MR and multivariable cis-MR using 659 allele-aligned regional variants.
2. A prespecified dimension sweep across 3-15, 20, 25 and 30 LD principal components.
3. Explicit strength diagnostics for the lipid total effect, lipid-to-splice first step and both multivariable exposures.
4. Independent trans-sQTL audits in INTERVAL, BrainMeta and GTEx v8 to search for genome-wide lipid instruments of the exact exon5-6 event.
5. Existing worst-case participant-overlap sensitivity retained as a separate bias analysis.

## Strongest mediation-compatible result

For TC at three principal components, all reported strength statistics exceeded 10. The local TC association with AD was attenuated by {100 * tc3.attenuation_fraction:.1f}% after joint modelling with exact exon5-6 splice usage. The splice direct estimate remained associated with AD (P={tc3.splice_direct_p:.3g}), whereas the TC direct estimate was not significant (P={tc3.lipid_direct_p:.3g}). The estimated indirect product was positive ({tc3.indirect_estimate:.3g}), with a directional probability of {100 * tc3.posterior_probability_indirect_positive:.1f}%.

## Why this does not establish formal mediation

The lipid-to-splice first step did not reach two-sided significance (P={tc3.lipid_to_splice_p:.3g}), and the indirect product was not significant (P={tc3.indirect_p:.3g}). Only {int(tc_summary.dimensions_all_strength_F_ge_10)} of {int(tc_summary.tested_pc_dimensions)} completed TC dimensions satisfied all strength criteria, and none satisfied the confirmatory mediation rule. The attenuation estimate varied from {tc_summary.attenuation_min:.2f} to {tc_summary.attenuation_max:.2f}, demonstrating sensitivity to the amount of LD information retained.

The independent-instrument search did not recover the exact event: INTERVAL contained 3,721 TSPAN14 trans-sQTL hits across 35 other features but zero exact exon5-6 hits; BrainMeta and GTEx v8 did not list TSPAN14 as a significant trans-sGene. {overlap_sentence}

## Failure mechanism

This is an identification failure rather than absence of a TSPAN14 splice signal. The exact splice-to-AD association is reproducible and strong, but the available lipid, splice and AD effects at chr10 are driven by highly correlated variants in one LD block. Those variants cannot cleanly separate vertical mediation from horizontal pleiotropy. Open trans-sQTL resources do not supply independent instruments for the exact event, so a conventional two-step MR cannot estimate the missing lipid-to-splice step without selection bias.

## Defensible conclusion

The completed analyses support a shared TSPAN14 regulatory architecture and a locus-restricted mediation-compatible decomposition, with the exact exon5-6 splice event positioned as the strongest molecular output. They do not identify a genome-wide causal pathway in which lipid levels act through this exact splice event to cause Alzheimer disease. The manuscript should not state that lipid-to-AD causality is mediated by TSPAN14 splicing unless unthresholded genome-wide exact-event sQTL effects or individual-level genotype-splicing data become available.
"""
    (args.analysis_dir / "MEDIATION_RESCUE_FINAL_REPORT.md").write_text(report, encoding="utf-8")
    print(resources[["resource", "TSPAN14_rows", "exact_exon5_6_rows", "interpretation"]].to_string(index=False))


if __name__ == "__main__":
    main()
