"""Synthesize the complete MR module into a claim-level mediation decision."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mr-dir", type=Path, required=True)
    args = parser.parse_args()
    out = args.mr_dir

    bidirectional = pd.read_csv(out / "03_genomewide_bidirectional_mr.tsv", sep="\t")
    apoe = pd.read_csv(out / "04_genomewide_bidirectional_apoe_sensitivity.tsv", sep="\t")
    strict = pd.read_csv(out / "08_local_reciprocal_mr.tsv", sep="\t")
    strength = pd.read_csv(out / "12_local_mvmr_conditional_strength.tsv", sep="\t")
    direct = pd.read_csv(out / "13_local_mvmr_official_estimates.tsv", sep="\t")
    threshold = pd.read_csv(out / "14_local_mvmr_threshold_sensitivity.tsv", sep="\t")
    global_mvmr = pd.read_csv(out / "17_global_joint_lipid_mvmr_estimates.tsv", sep="\t")
    global_strength = pd.read_csv(out / "18_global_joint_lipid_mvmr_strength_sensitivity.tsv", sep="\t")
    presso = pd.read_csv(out / "19_bidirectional_mr_presso.tsv", sep="\t")

    rows = []
    for lipid in ("TC", "LDL", "nonHDL"):
        total = bidirectional.loc[
            bidirectional["exposure"].eq(lipid) & bidirectional["outcome"].eq("AD")
        ].iloc[0]
        total_no_apoe = apoe.loc[
            apoe["exposure"].eq(lipid) & apoe["outcome"].eq("AD")
        ].iloc[0]
        sqtl_ad = strict.loc[
            strict["exposure"].eq("exact_exon5_6_sQTL")
            & strict["outcome"].eq("AD")
        ].iloc[0]
        model = f"exact_exon5_6_sQTL + {lipid} -> AD"
        splice_direct = direct.loc[
            direct["model"].eq(model)
            & direct["exposure"].eq("exact_exon5_6_sQTL")
        ].iloc[0]
        lipid_direct = direct.loc[
            direct["model"].eq(model) & direct["exposure"].eq(lipid)
        ].iloc[0]
        model_strength = strength.loc[strength["model"].eq(model)]
        strict_count = threshold.loc[
            threshold["model"].eq(model) & threshold["sqtl_p_threshold"].eq(5e-8),
            "n_instruments",
        ].iloc[0]
        apparent_attenuation = 1 - splice_direct["estimate"] / sqtl_ad["estimate"]
        rows.append(
            {
                "candidate_path": f"{lipid} -> exact exon5-6 splicing -> AD",
                "genomewide_lipid_to_AD_estimate": total["estimate"],
                "genomewide_lipid_to_AD_p": total["pvalue"],
                "apoe_excluded_lipid_to_AD_estimate": total_no_apoe["estimate"],
                "apoe_excluded_lipid_to_AD_p": total_no_apoe["pvalue"],
                "strict_sqtl_to_AD_estimate": sqtl_ad["estimate"],
                "strict_sqtl_to_AD_p": sqtl_ad["pvalue"],
                "local_mvmr_splice_direct_estimate": splice_direct["estimate"],
                "local_mvmr_splice_direct_p": splice_direct["pvalue"],
                "local_mvmr_lipid_direct_estimate": lipid_direct["estimate"],
                "local_mvmr_lipid_direct_p": lipid_direct["pvalue"],
                "apparent_splice_attenuation_fraction": apparent_attenuation,
                "minimum_conditional_F": model_strength["conditional_F"].min(),
                "strict_sqtl_instrument_count": strict_count,
                "lipid_to_splice_step": "not identifiable with independent genome-wide instruments",
                "indirect_effect": pd.NA,
                "proportion_mediated": pd.NA,
                "mediation_decision": "not_identified",
                "reason": (
                    "The genome-wide lipid total effect is not supported; strict local MVMR has only one "
                    "sQTL instrument, and the six-variant relaxed model has conditional F below 10. "
                    "The apparent coefficient attenuation compares different weak-instrument models and "
                    "cannot be interpreted as mediation."
                ),
            }
        )
    decision = pd.DataFrame(rows)
    decision.to_csv(out / "21_mediation_identification_gate.tsv", sep="\t", index=False)

    report = f"""# Issue 5: complete MR and mediation audit

## Analyses completed

1. Genome-wide bidirectional two-sample MR between AD and TC, LDL-C,
   non-HDL-C, HDL-C and triglycerides using genome-wide-significant instruments
   jointly filtered against a UK Biobank European LD reference at r2 < 0.001.
   Multiplicative random-effects IVW, weighted-median, MR-Egger, Cochran-Q,
   APOE-region exclusion and leave-one-out/MR-PRESSO diagnostics were specified.
2. Strict forward exact exon5-exon6 sQTL MR and locus-restricted reciprocal MR.
3. Global LDL-C + HDL-C + triglyceride MVMR using {int(global_mvmr['n_jointly_clumped_instruments'].iloc[0])}
   jointly clumped instruments.
4. Local exact-splice + lipid MVMR for TC, LDL-C and non-HDL-C, including
   official conditional-F statistics, threshold sensitivity and leave-one-out
   analysis.
5. A formal mediation-identification gate for lipid -> splice -> AD paths.

## Genome-wide bidirectional results

TC, LDL-C and non-HDL-C did not show an IVW total effect on AD (P =
{bidirectional.query("exposure == 'TC' and outcome == 'AD'").iloc[0]['pvalue']:.3g},
{bidirectional.query("exposure == 'LDL' and outcome == 'AD'").iloc[0]['pvalue']:.3g}, and
{bidirectional.query("exposure == 'nonHDL' and outcome == 'AD'").iloc[0]['pvalue']:.3g}).
The reverse AD-to-cholesterol estimates were positive in the unrestricted
analysis but became null after excluding chromosome 19:40-50 Mb, indicating
that the apparent reverse direction was predominantly APOE-region dependent.

The jointly clumped global LDL-C + HDL-C + triglyceride MVMR had conditionally
strong instruments across sampling-correlation sensitivity scenarios
(minimum conditional F = {global_strength['conditional_F'].min():.2f}). None
of the three direct effects was significant (all P >=
{global_mvmr['pvalue'].min():.3g}), and residual heterogeneity was extreme.
This does not support a general circulating-lipid mediation route to AD.

MR-PRESSO global tests were significant across the evaluated directions,
confirming pervasive horizontal heterogeneity. In the prespecified 100-strongest
instrument sensitivity set, outlier-corrected TC-, LDL-C- and non-HDL-C-to-AD
estimates remained nonsignificant (P =
{presso.query("exposure == 'TC' and outcome == 'AD' and estimate_type == 'outlier_corrected'").iloc[0]['pvalue']:.3g},
{presso.query("exposure == 'LDL' and outcome == 'AD' and estimate_type == 'outlier_corrected'").iloc[0]['pvalue']:.3g}, and
{presso.query("exposure == 'nonHDL' and outcome == 'AD' and estimate_type == 'outlier_corrected'").iloc[0]['pvalue']:.3g}).
Individual outlier calls were not used because 1,000 simulations did not provide
sufficient resolution for 64-100 instruments; the global and corrected-model
results are retained as sensitivity evidence only.

## Local splice-lipid MVMR

At the conventional sQTL threshold, only one r2 < 0.1 instrument was available,
so a two-exposure overidentified MVMR could not be fitted. P < 1e-4 yielded only
two instruments. Six instruments were obtained only after relaxing the sQTL
threshold to P < 1e-3. In those models, conditional F statistics ranged from
{strength['conditional_F'].min():.2f} to {strength['conditional_F'].max():.2f},
all below 10. Official MVMR estimates for both splice usage and lipid exposure
were nonsignificant after accounting for residual variance (all P >=
{direct['pvalue'].min():.3g}). Leave-one-out estimates varied substantially and
conditional strength remained weak.

The fitted splice coefficient decreased from approximately 0.065 in strict
single-instrument forward MR to approximately 0.027 in the relaxed MVMR.
However, this is not a valid mediated fraction: the estimators use different
instrument sets, the MVMR instruments are conditionally weak, and the lipid and
splice associations arise from the same cis architecture.

## Causal interpretation

The complete MR module supports locus-level molecular triangulation: the exact
TSPAN14 splice feature, AD and cholesterol traits share a strong local genetic
configuration. It does not identify a directional lipid -> TSPAN14 splicing ->
AD mediation path. A valid two-step mediation analysis would require
genome-wide lipid instruments with measured associations to the exact junction,
or a second set of splice instruments that predicts junction usage independently
of lipid effects. Released BA24 sQTL statistics are cis-only, so the first step
cannot be estimated with independent instruments.

The current evidence is therefore compatible with correlated horizontal
pleiotropy or a shared regulatory program. The exact splice event remains the
best-resolved molecular readout of the locus, but not a proven mediator of a
circulating-lipid effect.

## Figure 5 consequence

Any solid arrow labelled `TSPAN14 splicing -> lipid -> AD` should be removed.
AD, cholesterol traits and the exact sQTL should converge on the regulatory
haplotype as parallel association layers. A possible splice-dependent cellular
lipid phenotype may remain as a dashed, explicitly testable prediction.
"""
    (out / "ISSUE5_COMPLETE_MR_MEDIATION_REPORT.md").write_text(report, encoding="utf-8")

    claims = pd.DataFrame(
        [
            {
                "claim": "Genetically predicted exact exon5-exon6 usage is associated with AD and cholesterol-trait associations at the TSPAN14 locus.",
                "status": "supported_as_locus_level_molecular_triangulation",
                "prohibited_extension": "Do not translate this into proof of a lipid-mediated causal path.",
            },
            {
                "claim": "TC, LDL-C or non-HDL-C has a genome-wide total causal effect on AD in these data.",
                "status": "not_supported",
                "prohibited_extension": "Do not use local colocalization to override the null genome-wide total-effect MR.",
            },
            {
                "claim": "Adding lipid exposure removes the direct effect of TSPAN14 splicing on AD.",
                "status": "not_identified",
                "prohibited_extension": "The relaxed local MVMR has conditional F below 10 and different instruments from the strict cis-MR.",
            },
            {
                "claim": "TSPAN14 splicing mediates a circulating-lipid effect on AD.",
                "status": "not_identified",
                "prohibited_extension": "No indirect effect or proportion mediated should be reported.",
            },
            {
                "claim": "The locus is compatible with shared or correlated pleiotropic regulation.",
                "status": "supported_as_the_most_conservative_model",
                "prohibited_extension": "Compatibility is not proof that all associations arise from one molecular mechanism.",
            },
        ]
    )
    claims.to_csv(out / "22_final_causal_claim_ledger.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
