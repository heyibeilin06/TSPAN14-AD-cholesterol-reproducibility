#!/usr/bin/env python3
"""Build synchronized Supplementary Tables S1-S19 from retained audited data."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs" / "supplement_v19"
SRC = OUT / "source_tables"
RETAINED = ROOT / "tables" / "supplementary" / "source_data"
TITLE = "Local Alzheimer disease-cholesterol colocalization converges on TSPAN14 canonical–cryptic splice choice"


def read(rel: str) -> pd.DataFrame:
    return pd.read_csv(ROOT / rel, sep="\t")


def retained(number: int) -> pd.DataFrame:
    return pd.read_csv(RETAINED / f"Table_S{number:02d}.tsv", sep="\t")


def tagged(rel: str, label: str) -> pd.DataFrame:
    data = read(rel)
    if "analysis_block" in data.columns:
        data = data.rename(columns={"analysis_block": "source_analysis_block"})
    data.insert(0, "analysis_block", label)
    return data


def combine(*frames: pd.DataFrame) -> pd.DataFrame:
    return pd.concat(frames, ignore_index=True, sort=False)


def table4() -> pd.DataFrame:
    old = retained(4)
    old = old[~old.get("analysis_block", pd.Series(index=old.index, dtype=str)).eq("coloc prior sensitivity")].copy()
    result = combine(
        old,
        tagged("audit/reviewer_revision/coloc_prior_sensitivity.tsv", "coloc prior sensitivity"),
        tagged("audit/reviewer_revision/exact_event_coloc_prior_sensitivity.tsv", "exact-event coloc prior sensitivity"),
    )
    return result.drop_duplicates().reset_index(drop=True)


def table6() -> pd.DataFrame:
    return combine(
        tagged("audit/reviewer_revision/count_level_acceptor_choice_donors.tsv", "donor-level canonical and cryptic read counts"),
        tagged("audit/reviewer_revision/count_level_acceptor_choice_by_genotype.tsv", "rs7080009 genotype summary"),
    )


def table7() -> pd.DataFrame:
    old = retained(7)
    old = old[old.get("analysis_block", pd.Series(index=old.index, dtype=str)).isin(["count concordance", "count concordance summary"])].copy()
    return combine(
        tagged("audit/reviewer_revision/count_level_acceptor_choice_models.tsv", "penalized donor-detection and beta-binomial read-count models"),
        tagged("audit/reviewer_revision/count_level_acceptor_choice_depth_sensitivity.tsv", "depth-threshold sensitivity"),
        tagged("audit/reviewer_revision/count_level_acceptor_choice_leave_one_out.tsv", "leave-one-donor-out sensitivity"),
        old,
    )


def table8() -> pd.DataFrame:
    old = retained(8)
    old = old[~old.get("analysis_block", pd.Series(index=old.index, dtype=str)).isin([
        "GTEx focal-tissue donor counts",
        "GTEx focal-tissue donor overlap",
    ])].copy()
    return combine(
        old,
        tagged("audit/reviewer_revision/gtex_focal_tissue_donor_counts.tsv", "GTEx focal-tissue donor counts"),
        tagged("audit/reviewer_revision/gtex_focal_tissue_donor_overlap.tsv", "GTEx focal-tissue donor overlap"),
    )


def table10() -> pd.DataFrame:
    return combine(
        tagged("audit/reviewer_revision/cis_mr_numerical_sensitivity.tsv", "rebuilt-LD cis-MR numerical sensitivity"),
        tagged("audit/reviewer_revision/cis_mr_ld_diagnostics.tsv", "rebuilt-LD diagnostics"),
        tagged("audit/reviewer_revision/cis_mr_ld_perturbation_sensitivity.tsv", "LD perturbation sensitivity"),
        tagged("audit/reviewer_revision/cis_mr_ld_allele_mapping.tsv", "signed allele mapping"),
        tagged("audit/reviewer_revision/cis_mr_overlap_bias_bound_rebuilt_ld.tsv", "conservative sample-overlap bias bound"),
    )


def table13() -> pd.DataFrame:
    data = retained(13)
    data["primary_result"] = data["primary_result"].replace(
        "all_laub_core_ad_effect_alleles_decrease_target_junction_usage",
        "all_laub_core_ad_effect_alleles_decrease_cryptic_acceptor_junction_usage",
    )
    return data


def table17() -> pd.DataFrame:
    data = retained(17)
    mask = data.get("figure_element", pd.Series(index=data.index, dtype=str)).eq("Four neural tissues")
    data.loc[mask, "evidence_class"] = "cross-tissue consistency in partially overlapping GTEx neural tissues"
    data.loc[mask, "supported_claim"] = (
        "The coordinate-identical junction shows directionally consistent associations across partially "
        "overlapping GTEx BA24, hippocampus, putamen and cervical spinal cord donor sets."
    )
    return data


TABLES: dict[int, tuple[str, callable]] = {
    1: ("Data resources and baseline genome-wide genetic correlation", lambda: retained(1)),
    2: ("Extended-APOE conditioning and physical-window sensitivity", lambda: retained(2)),
    3: ("Complete non-APOE regional AD-lipid screen", lambda: retained(3)),
    4: ("Regional, multiple-signal and prior-sensitivity colocalization", table4),
    5: ("TSPAN14 fine-mapping and functional prioritization", lambda: retained(5)),
    6: ("Donor-level canonical–cryptic acceptor read counts", table6),
    7: ("Canonical–cryptic read-count models and robustness", table7),
    8: ("Exact exon5-6 sQTL cross-tissue consistency and GTEx donor overlap", table8),
    9: ("Exon5-6 and exon6-7 neural-tissue co-usage", lambda: retained(9)),
    10: ("Exact-event cis-MR, LD diagnostics and overlap sensitivity", table10),
    11: ("Bidirectional MR, multivariable MR and PC-GMM identification scope", lambda: retained(11)),
    12: ("Neural cell localization and molecular-QTL context", lambda: retained(12)),
    13: ("AD disease-state RNA and differential transcript usage", table13),
    14: ("Cross-study single-nucleus meta-analysis", lambda: retained(14)),
    15: ("Independent pseudobulk disease-state sensitivity", lambda: retained(15)),
    16: ("Transcript, EC2 and AlphaFold reference-model mapping", lambda: retained(16)),
    17: ("Proteomic, perturbation and downstream functional evidence", table17),
    18: ("Signed functional-variant haplotype audit", lambda: tagged(
        "audit/reviewer_revision/signed_functional_haplotype_audit.tsv",
        "1000 Genomes Phase 3 EUR phased haplotypes",
    )),
    19: ("Reproducibility and reviewer-revision audit", lambda: combine(
        tagged("audit/reviewer_revision/coloc_prior_sensitivity.tsv", "coloc prior sensitivity"),
        tagged("audit/reviewer_revision/exact_event_coloc_prior_sensitivity.tsv", "exact-event coloc prior sensitivity"),
        tagged("audit/reviewer_revision/count_level_acceptor_choice_models.tsv", "canonical–cryptic model audit"),
        tagged("audit/reviewer_revision/cis_mr_ld_diagnostics.tsv", "cis-MR LD audit"),
        tagged("audit/reviewer_revision/cis_mr_ld_perturbation_sensitivity.tsv", "cis-MR perturbation audit"),
        tagged("audit/reviewer_revision/gwas_cohort_overlap_audit.tsv", "source-GWAS cohort-overlap audit"),
        tagged("audit/reviewer_revision/cis_mr_overlap_bias_bound_rebuilt_ld.tsv", "sample-overlap bias-bound audit"),
        tagged("audit/claims_to_code.tsv", "claim-to-code map"),
    )),
}


def main() -> None:
    SRC.mkdir(parents=True, exist_ok=True)
    (OUT / "audit").mkdir(parents=True, exist_ok=True)
    manifest = []
    for number, (title, builder) in TABLES.items():
        data = builder()
        path = SRC / f"Table_S{number:02d}.tsv"
        data.to_csv(path, sep="\t", index=False)
        manifest.append({"number": number, "title": title, "file": path.name, "rows": len(data), "columns": len(data.columns)})
    pd.DataFrame(manifest).to_csv(SRC / "Table_S00_Index.tsv", sep="\t", index=False)
    evidence = [
        ("S1", "APOE contribution to AD-HDL-C genetic correlation", "Tables S1-S2"),
        ("S2-S3", "regional screen, fine-mapping and colocalization sensitivity", "Tables S3-S5"),
        ("S4", "canonical–cryptic acceptor counts and models", "Tables S6-S7"),
        ("S5-S6", "exact-event cross-tissue consistency, donor overlap and adjacent-junction co-usage", "Tables S8-S9"),
        ("S7-S8", "cis-MR, overlap, MVMR and mediation scope", "Tables S10-S11"),
        ("S9", "cell context, disease RNA and structural interpretation", "Tables S12-S17"),
        ("Methods and Results", "signed functional haplotype and reproducibility audits", "Tables S18-S19"),
    ]
    pd.DataFrame(evidence, columns=["supplementary_figures", "claim", "supplementary_tables"]).to_csv(
        OUT / "audit" / "manuscript_evidence_map.tsv", sep="\t", index=False
    )
    (OUT / "audit" / "supplement_manifest.json").write_text(
        json.dumps({"manuscript_title": TITLE, "tables": manifest, "figures": 9}, indent=2), encoding="utf-8"
    )
    print(json.dumps({"title": TITLE, "tables": len(manifest), "output": str(OUT)}, indent=2))


if __name__ == "__main__":
    main()
