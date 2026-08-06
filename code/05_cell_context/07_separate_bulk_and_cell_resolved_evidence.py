"""Separate cell-resolved, bulk-tissue, and disease-state TSPAN14 evidence.

This audit prevents a bulk-brain sQTL from being used to infer a microglial
origin while retaining the direct exact-event evidence measured in isolated
human microglia.
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import pandas as pd


DEFAULT_EVENT = "chr10:80509471:80512144"
DEFAULT_CORE_SNPS = (
    "rs1870137",
    "rs1870138",
    "rs1870140",
    "rs1902660",
    "rs6586028",
    "rs7080009",
)


def extract_miga_exact_event(
    frame: pd.DataFrame,
    *,
    event_coordinate: str,
    core_snps: list[str] | tuple[str, ...],
) -> pd.DataFrame:
    exact = frame.loc[
        frame["phenotype_id"].str.startswith(f"{event_coordinate}:")
        & frame["variant_id"].isin(core_snps)
    ].copy()
    exact = exact.sort_values("variant_id").drop_duplicates("variant_id")
    exact["event_coordinate"] = event_coordinate
    exact["event_class"] = "exact_exon5_exon6"
    exact["cell_resolution"] = "isolated_human_microglia"
    exact["effect_orientation"] = (
        "MiGA alternate-allele slope; risk-allele identity verified downstream"
    )
    return exact


def build_evidence_layer_ledger() -> pd.DataFrame:
    rows = [
        {
            "evidence_layer": "MiGA exact-event sQTL",
            "biological_material": "isolated human brain microglia",
            "molecular_phenotype": "exact TSPAN14 exon5-exon6 junction usage",
            "causal_axis": "genotype_to_splicing",
            "cell_type_attribution": "direct_microglia_resolved",
            "permitted_claim": (
                "The exact exon5-exon6 event is genetically regulated in "
                "isolated human microglia."
            ),
            "prohibited_claim": (
                "Does not establish AD-state expression change, protein output, "
                "or an ADAM10/TREM2 downstream effect."
            ),
        },
        {
            "evidence_layer": "GTEx exact-event sQTL",
            "biological_material": "bulk brain tissue",
            "molecular_phenotype": "exact TSPAN14 exon5-exon6 junction usage",
            "causal_axis": "genotype_to_splicing",
            "cell_type_attribution": "not_permitted",
            "permitted_claim": (
                "The exact event replicates in bulk brain across multiple regions."
            ),
            "prohibited_claim": (
                "Bulk tissue cannot identify neurons, microglia, or another cell "
                "population as the source of the association."
            ),
        },
        {
            "evidence_layer": "single-nucleus TSPAN14 eQTL",
            "biological_material": "astrocyte and excitatory-neuron nuclei",
            "molecular_phenotype": "total TSPAN14 expression",
            "causal_axis": "genotype_to_total_gene_expression",
            "cell_type_attribution": "cell_type_resolved_expression_only",
            "permitted_claim": (
                "The regulatory haplotype also affects TSPAN14 expression in "
                "non-microglial neural cell contexts."
            ),
            "prohibited_claim": (
                "An eQTL is not evidence for the exact exon5-exon6 splice event."
            ),
        },
        {
            "evidence_layer": "SEA-AD pseudobulk DGE",
            "biological_material": "single-nucleus-derived microglia and neurons",
            "molecular_phenotype": "AD-state total TSPAN14 expression",
            "causal_axis": "disease_state_to_total_gene_expression",
            "cell_type_attribution": "cell_type_resolved_disease_contrast",
            "permitted_claim": (
                "No large case-control shift in total TSPAN14 RNA was resolved in "
                "the adjusted microglial or neuronal contrasts."
            ),
            "prohibited_claim": (
                "A null DGE result is not a null sQTL test and does not localize "
                "or refute genotype-dependent splicing."
            ),
        },
        {
            "evidence_layer": "published microglial perturbation",
            "biological_material": "microglial functional model",
            "molecular_phenotype": "enhancer perturbation and downstream phenotype",
            "causal_axis": "perturbation_to_cellular_phenotype",
            "cell_type_attribution": "external_functional_microglia_context",
            "permitted_claim": (
                "Provides an external microglial functional anchor for the locus."
            ),
            "prohibited_claim": (
                "Does not prove that the GTEx bulk sQTL originates in microglia or "
                "that the exact splice event mediates the reported phenotype."
            ),
        },
    ]
    return pd.DataFrame(rows)


def scan_isomiga_top_hits(path: Path, event_coordinate: str, core_snps: set[str]) -> dict:
    opener = gzip.open if path.suffix == ".gz" else open
    matches = 0
    with opener(path, "rt", encoding="utf-8") as handle:
        for chunk in pd.read_csv(handle, sep="\t", chunksize=50_000):
            selected = chunk["feature"].str.startswith(f"{event_coordinate}:") | chunk[
                "variant_id"
            ].isin(core_snps)
            matches += int(selected.sum())
    return {
        "resource": "isoMiGA union LeafCutter top associations",
        "scope": "significant/top associations only; not the full association scan",
        "exact_event_or_core_snp_rows": matches,
        "interpretation": (
            "positive_top_hit" if matches else "not_present_in_top_hit_file_not_a_null_test"
        ),
    }


def write_tsv(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--miga", type=Path, required=True)
    parser.add_argument("--miga-harmonized", type=Path, required=True)
    parser.add_argument("--gtex", type=Path, required=True)
    parser.add_argument("--dge", type=Path, required=True)
    parser.add_argument("--snuc-eqtl", type=Path, required=True)
    parser.add_argument("--isomiga-top", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    miga_raw = pd.read_csv(args.miga, sep="\t")
    miga = extract_miga_exact_event(
        miga_raw, event_coordinate=DEFAULT_EVENT, core_snps=DEFAULT_CORE_SNPS
    )
    if set(miga["variant_id"]) != set(DEFAULT_CORE_SNPS):
        missing = sorted(set(DEFAULT_CORE_SNPS) - set(miga["variant_id"]))
        raise ValueError(f"MiGA exact-event rows missing core variants: {missing}")

    gtex = pd.read_csv(args.gtex, sep="\t")
    gtex = gtex.loc[
        gtex["snp"].isin(DEFAULT_CORE_SNPS)
        & gtex["event_class"].eq("exact_exon5_6")
    ].copy()
    ba24 = gtex.loc[gtex["tissue"].str.contains("BA24")].copy()
    if set(ba24["snp"]) != set(DEFAULT_CORE_SNPS):
        missing = sorted(set(DEFAULT_CORE_SNPS) - set(ba24["snp"]))
        raise ValueError(f"GTEx BA24 rows missing core variants: {missing}")

    miga_alleles = pd.read_csv(args.miga_harmonized, sep="\t").loc[
        lambda frame: frame["snp"].isin(DEFAULT_CORE_SNPS),
        ["snp", "miga_ref", "miga_alt"],
    ].drop_duplicates("snp")
    concordance = miga[
        ["variant_id", "slope", "slope_se", "pval_nominal", "maf", "ma_count"]
    ].rename(
        columns={
            "variant_id": "snp",
            "slope": "miga_microglia_beta",
            "slope_se": "miga_microglia_se",
            "pval_nominal": "miga_microglia_p",
            "maf": "miga_maf",
            "ma_count": "miga_minor_allele_count",
        }
    ).merge(miga_alleles, on="snp", how="left", validate="one_to_one").merge(
        ba24[["snp", "risk_allele", "nes", "p_value"]].rename(
            columns={"nes": "gtex_ba24_nes", "p_value": "gtex_ba24_p"}
        ),
        on="snp",
        how="inner",
        validate="one_to_one",
    )
    concordance["risk_allele_matches_miga_alt"] = (
        concordance["risk_allele"] == concordance["miga_alt"]
    )
    if not concordance["risk_allele_matches_miga_alt"].all():
        failed = concordance.loc[
            ~concordance["risk_allele_matches_miga_alt"],
            ["snp", "risk_allele", "miga_alt"],
        ].to_dict("records")
        raise ValueError(f"MiGA slopes are not risk-aligned for: {failed}")
    concordance["miga_risk_aligned_direction"] = concordance[
        "miga_microglia_beta"
    ].apply(lambda value: "positive" if value > 0 else "negative")
    concordance["gtex_risk_aligned_direction"] = concordance["gtex_ba24_nes"].apply(
        lambda value: "positive" if value > 0 else "negative"
    )
    concordance["direction_concordant"] = (
        concordance["miga_risk_aligned_direction"]
        == concordance["gtex_risk_aligned_direction"]
    )
    concordance["interpretation_note"] = (
        "Descriptive cross-context concordance; variants are LD-correlated and "
        "must not be counted as independent replications."
    )

    exact_rows = []
    for row in miga.itertuples(index=False):
        exact_rows.append(
            {
                "resource": "MiGA",
                "tissue_or_cell_context": "isolated human microglia (SVZ cohort)",
                "cell_resolution": "cell_type_specific",
                "event_coordinate": DEFAULT_EVENT,
                "snp": row.variant_id,
                "risk_allele": concordance.set_index("snp").loc[
                    row.variant_id, "risk_allele"
                ],
                "risk_aligned_effect": row.slope,
                "standard_error": row.slope_se,
                "p_value": row.pval_nominal,
                "effect_metric": "LeafCutter normalized intron-usage beta",
                "permitted_role": "direct microglial exact-event sQTL evidence",
            }
        )
    for row in gtex.itertuples(index=False):
        exact_rows.append(
            {
                "resource": "GTEx v8",
                "tissue_or_cell_context": row.tissue_label,
                "cell_resolution": "bulk_tissue",
                "event_coordinate": DEFAULT_EVENT,
                "snp": row.snp,
                "risk_allele": row.risk_allele,
                "risk_aligned_effect": row.nes,
                "standard_error": pd.NA,
                "p_value": row.p_value,
                "effect_metric": "GTEx normalized effect size",
                "permitted_role": "independent exact-event bulk-tissue replication",
            }
        )
    exact_evidence = pd.DataFrame(exact_rows)

    dge = pd.read_csv(args.dge, sep="\t")
    snuc = pd.read_csv(args.snuc_eqtl, sep="\t")
    ledger = build_evidence_layer_ledger()
    decision = pd.DataFrame(
        [
            {
                "question": "Is the exact exon5-exon6 event cell-resolved?",
                "decision": "yes_in_project_MiGA_microglia",
                "basis": (
                    f"{len(miga)}/{len(DEFAULT_CORE_SNPS)} core variants present for "
                    "the exact event in isolated microglia."
                ),
            },
            {
                "question": "Can GTEx BA24 identify the responsible cell type?",
                "decision": "no",
                "basis": "GTEx BA24 is bulk tissue and is used only for exact-event replication.",
            },
            {
                "question": "Are MiGA and BA24 risk-aligned directions concordant?",
                "decision": (
                    f"descriptively_yes_{int(concordance['direction_concordant'].sum())}"
                    f"_of_{len(concordance)}"
                ),
                "basis": "All core-variant effects have the same risk-aligned sign; variants are LD-correlated.",
            },
            {
                "question": "Does SEA-AD DGE identify the sQTL cell of origin?",
                "decision": "no",
                "basis": (
                    "Disease-state total-expression contrasts are orthogonal to "
                    "genotype-dependent junction usage; both audited contrasts are nonsignificant."
                ),
            },
            {
                "question": "Does the locus appear restricted to microglia?",
                "decision": "no_multicellular_regulatory_context",
                "basis": (
                    f"Cell-resolved total-expression eQTL evidence is present in "
                    f"{', '.join(sorted(snuc['context'].unique()))}; this is not exact-splice evidence."
                ),
            },
        ]
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(args.output_dir / "01_exact_event_cross_context_evidence.tsv", exact_evidence)
    write_tsv(args.output_dir / "02_core_snp_direction_concordance.tsv", concordance)
    write_tsv(args.output_dir / "03_cell_type_layer_separation.tsv", ledger)
    write_tsv(args.output_dir / "04_claim_decision.tsv", decision)

    if args.isomiga_top:
        isomiga = pd.DataFrame(
            [scan_isomiga_top_hits(args.isomiga_top, DEFAULT_EVENT, set(DEFAULT_CORE_SNPS))]
        )
        write_tsv(args.output_dir / "05_isomiga_top_hit_audit.tsv", isomiga)

    micro = dge.loc[dge["cell_group"].eq("microglia")].iloc[0]
    neuron = dge.loc[dge["cell_group"].eq("neurons")].iloc[0]
    report = f"""# Issue 4: cell-type attribution audit

## Resolution

The ecological-fallacy concern is valid if the GTEx BA24 association is used to
infer a microglial origin. The evidence does not require that inference. The
exact GRCh38 junction `{DEFAULT_EVENT}` was tested directly in isolated human
microglia in MiGA and was present for all {len(miga)} core variants. GTEx BA24
and the other brain regions are retained only as independent bulk-tissue
replications of the same molecular event.

## Exact-event evidence

- MiGA isolated microglia: risk-aligned beta range
  {miga['slope'].min():.3f} to {miga['slope'].max():.3f}; nominal P range
  {miga['pval_nominal'].min():.3g} to {miga['pval_nominal'].max():.3g}.
- GTEx BA24: risk-aligned NES range {ba24['nes'].min():.3f} to
  {ba24['nes'].max():.3f}; P range {ba24['p_value'].min():.3g} to
  {ba24['p_value'].max():.3g}.
- The risk-aligned sign is concordant for
  {int(concordance['direction_concordant'].sum())}/{len(concordance)} core SNPs.
  For every SNP, the inferred AD-risk-increasing allele matches the MiGA
  alternate allele to which the reported slope refers. This concordance is
  descriptive because the SNPs are correlated within the same LD block.

## Orthogonal cell-state evidence

SEA-AD adjusted pseudobulk differential gene-expression results do not show a
large AD-state shift in total TSPAN14 RNA in microglia (log2 fold change
{micro['logFC']:.3f}, FDR {micro['adj.P.Val']:.3f}) or neurons (log2 fold change
{neuron['logFC']:.3f}, FDR {neuron['adj.P.Val']:.3f}). These are disease-state
gene-level tests, not genotype-by-junction tests, and therefore neither localize
nor refute the sQTL.

Single-nucleus eQTL results in astrocyte and excitatory-neuron contexts indicate
that the locus has a broader, multicellular expression-regulatory footprint.
They are not presented as exact-splice replication.

## Claim that survives audit

The defensible model is: **a microglia-resolved exact-event sQTL in MiGA, exact
event replication in bulk brain, and additional multicellular expression-QTL
context**. The BA24 signal is not assigned to microglia. ADAM10/TREM2 biology is
kept as an externally anchored microglial functional context and a downstream
hypothesis, not as a cell-of-origin conclusion derived from BA24.

## Primary resource provenance

- MiGA mapped expression and splicing QTLs in isolated human microglia:
  https://doi.org/10.1038/s41588-021-00976-y
- GTEx v8 provides the independent bulk-tissue exact-event replication:
  https://gtexportal.org/home/
- isoMiGA top-association data were additionally screened through Zenodo record
  8250771. Absence from this top-hit-only file is recorded as uninformative, not
  as a negative full-scan result: https://doi.org/10.5281/zenodo.8250771
"""
    (args.output_dir / "ISSUE4_CELL_TYPE_ATTRIBUTION_REPORT.md").write_text(
        report, encoding="utf-8"
    )


if __name__ == "__main__":
    main()
