from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path
from statistics import mean


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "outputs" / "targeted_splice_validation"
MIGA_QTL = ROOT / "results" / "tables" / "microglia_full_qtl_target_extracts" / "SVZ_eur_splicing_peer0_gene.cis_qtl_nominal_tabixed.candidate_extract.tsv"
EXON_SOURCE = ROOT / "results" / "figures" / "figure2" / "figure2_panel_b_tspan14_exon_source.tsv"

EVENT_5_6 = "chr10:80509471:80512144:clu_4260_+:ENSG00000108219.15"
EVENT_6_7 = "chr10:80512269:80514019:clu_4261_+:ENSG00000108219.15"
CORE_SNPS = ("rs7080009", "rs1870138", "rs1870137", "rs1902660", "rs6586028", "rs1870140")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def pearson(x: list[float], y: list[float]) -> float | None:
    if len(x) < 3 or len(x) != len(y):
        return None
    x_mean, y_mean = mean(x), mean(y)
    numerator = sum((a - x_mean) * (b - y_mean) for a, b in zip(x, y))
    denominator = math.sqrt(sum((a - x_mean) ** 2 for a in x) * sum((b - y_mean) ** 2 for b in y))
    return numerator / denominator if denominator else None


def rank(values: list[float]) -> list[float]:
    ordered = sorted(enumerate(values), key=lambda item: item[1])
    ranks = [0.0] * len(values)
    index = 0
    while index < len(ordered):
        end = index
        while end + 1 < len(ordered) and ordered[end + 1][1] == ordered[index][1]:
            end += 1
        average_rank = (index + end + 2) / 2
        for pos in range(index, end + 1):
            ranks[ordered[pos][0]] = average_rank
        index = end + 1
    return ranks


def spearman(x: list[float], y: list[float]) -> float | None:
    return pearson(rank(x), rank(y))


def classify_adjacent_feature(
    *,
    psi_available: bool,
    psi_spearman: float | None,
    genotype_effect_relation: str,
    same_transcript_chain: bool,
) -> dict[str, str]:
    if not psi_available or psi_spearman is None:
        return {
            "evidence_class": "unresolved_adjacent_feature",
            "supported_claim": "Adjacent TSPAN14 splice feature provides locus-context evidence only; it is not exact replication of the project junction.",
            "manuscript_action": "Remove target-equivalent and replication wording; do not use this feature as a core direction-consistency layer.",
        }
    if psi_spearman <= -0.2 or genotype_effect_relation == "reciprocal":
        return {
            "evidence_class": "competitive_adjacent_feature",
            "supported_claim": "The two junctions represent distinct splice outputs with evidence compatible with competitive transcript processing.",
            "manuscript_action": "Present separate junction outputs and test exon skipping/isoform consequences; do not aggregate effect directions.",
        }
    if psi_spearman >= 0.2 and genotype_effect_relation == "concordant" and same_transcript_chain:
        return {
            "evidence_class": "co_regulated_adjacent_feature",
            "supported_claim": "The two features support a local TSPAN14 splice module while remaining non-identical junctions.",
            "manuscript_action": "Describe co-regulated adjacent features, never exact replication.",
        }
    return {
        "evidence_class": "unresolved_adjacent_feature",
        "supported_claim": "Adjacent TSPAN14 splice feature has insufficient co-usage evidence for module-level interpretation.",
        "manuscript_action": "Keep as contextual locus evidence only.",
    }


def audit_miga_qtls() -> tuple[list[dict[str, object]], dict[str, object]]:
    by_event: dict[str, dict[str, dict[str, str]]] = defaultdict(dict)
    for row in read_tsv(MIGA_QTL):
        if row["phenotype_id"] in (EVENT_5_6, EVENT_6_7):
            by_event[row["phenotype_id"]][row["variant_id"]] = row

    common = sorted(set(by_event[EVENT_5_6]) & set(by_event[EVENT_6_7]))
    slope_5_6 = [float(by_event[EVENT_5_6][variant]["slope"]) for variant in common]
    slope_6_7 = [float(by_event[EVENT_6_7][variant]["slope"]) for variant in common]
    signs_same = sum((a >= 0) == (b >= 0) for a, b in zip(slope_5_6, slope_6_7))
    correlation = pearson(slope_5_6, slope_6_7)

    rows: list[dict[str, object]] = []
    for variant in CORE_SNPS:
        first, second = by_event[EVENT_5_6].get(variant, {}), by_event[EVENT_6_7].get(variant, {})
        rows.append(
            {
                "variant_id": variant,
                "exon5_6_slope": first.get("slope", ""),
                "exon5_6_pvalue": first.get("pval_nominal", ""),
                "exon6_7_slope": second.get("slope", ""),
                "exon6_7_pvalue": second.get("pval_nominal", ""),
                "interpretation": "Same MiGA SVZ resource; this compares variant-level QTL effects, not individual-level PSI co-usage.",
            }
        )
    summary = {
        "shared_variant_count": len(common),
        "slope_pearson_r": correlation,
        "same_sign_fraction": signs_same / len(common) if common else None,
        "genotype_effect_relation": "weak_or_inconsistent",
        "interpretation": "Weak descriptive slope correlation cannot substitute for a within-individual PSI correlation and does not establish co-regulation.",
    }
    return rows, summary


def audit_transcript_structure() -> tuple[list[dict[str, object]], bool]:
    exons = [row for row in read_tsv(EXON_SOURCE) if row.get("transcript_id") == "ENST00000429989"]
    indexed = {int(row["exon_rank"]): row for row in exons}
    exon6 = indexed[6]
    # GENCODE genomic exon coordinates are 1-based inclusive.
    exon6_length = int(exon6["exon_end"]) - int(exon6["exon_start"]) + 1
    rows = [
        {
            "transcript_id": "ENST00000429989",
            "event": "project_exon5_exon6",
            "junction_coordinates": "chr10:80509471-80512144",
            "structural_status": "canonical_chain_member",
            "protein_transition": "AA150/151",
            "interpretation": "Canonical transcript contains exon5 followed by exon6.",
        },
        {
            "transcript_id": "ENST00000429989",
            "event": "adjacent_exon6_exon7",
            "junction_coordinates": "chr10:80512269-80514018",
            "structural_status": "canonical_chain_member",
            "protein_transition": "AA192/193",
            "interpretation": "Canonical transcript contains exon6 followed by exon7; coordinate end is LeafCutter/BED convention one base before exon7 start.",
        },
        {
            "transcript_id": "ENST00000429989",
            "event": "hypothetical_exon5_exon7_skip",
            "junction_coordinates": "chr10:80509471-80514019",
            "structural_status": "not_observed_in_current_inputs",
            "protein_transition": "exon6_length_nt=%d; mod3=%d" % (exon6_length, exon6_length % 3),
            "interpretation": "Skipping coding exon6 would be in-frame by length alone; this does not demonstrate that exon5-exon7 skipping occurs or excludes other NMD mechanisms.",
        },
    ]
    return rows, True


def read_optional_psi(path: Path | None) -> tuple[bool, float | None, int, str]:
    if path is None or not path.exists():
        return False, None, 0, "No individual-level junction-usage matrix was available locally."
    rows = read_tsv(path)
    lookup = {row.get("feature_id", row.get("phenotype_id", "")): row for row in rows}
    first, second = lookup.get(EVENT_5_6), lookup.get(EVENT_6_7)
    if not first or not second:
        return False, None, 0, "Matrix does not contain both exact MiGA feature identifiers."
    metadata = {"feature_id", "phenotype_id", "chrom", "start", "end", "strand"}
    samples = [key for key in first if key in second and key not in metadata]
    x, y = [], []
    for sample in samples:
        try:
            x.append(float(first[sample]))
            y.append(float(second[sample]))
        except (TypeError, ValueError):
            continue
    return True, spearman(x, y), len(x), "Computed from supplied same-cohort individual-level junction-usage matrix."


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit adjacent TSPAN14 splice-event evidence without conflating QTL slopes and PSI.")
    parser.add_argument("--psi-matrix", type=Path, default=None, help="Optional TSV containing both exact feature IDs as rows and matched sample IDs as columns.")
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    qtl_rows, qtl_summary = audit_miga_qtls()
    transcript_rows, same_chain = audit_transcript_structure()
    psi_available, psi_rho, psi_n, psi_note = read_optional_psi(args.psi_matrix)
    decision = classify_adjacent_feature(
        psi_available=psi_available,
        psi_spearman=psi_rho,
        genotype_effect_relation=str(qtl_summary["genotype_effect_relation"]),
        same_transcript_chain=same_chain,
    )

    write_tsv(
        args.outdir / "01_input_inventory.tsv",
        [
            {"input": str(MIGA_QTL.relative_to(ROOT)), "availability": "available", "role": "MiGA SVZ nominal sQTL summary statistics"},
            {"input": str(EXON_SOURCE.relative_to(ROOT)), "availability": "available", "role": "canonical TSPAN14-207 exon annotation"},
            {"input": str(args.psi_matrix) if args.psi_matrix else "not_supplied", "availability": "available" if psi_available else "not_available", "role": "same-cohort individual-level junction usage required for PSI correlation"},
        ],
        ["input", "availability", "role"],
    )
    write_tsv(args.outdir / "02_miga_adjacent_junction_qtl_audit.tsv", qtl_rows, list(qtl_rows[0]))
    write_tsv(args.outdir / "03_transcript_structure_audit.tsv", transcript_rows, list(transcript_rows[0]))
    write_tsv(
        args.outdir / "04_external_psi_correlation.tsv",
        [{"psi_available": psi_available, "spearman_rho": psi_rho if psi_rho is not None else "", "n_samples": psi_n, "interpretation": psi_note}],
        ["psi_available", "spearman_rho", "n_samples", "interpretation"],
    )
    evidence_row = {**qtl_summary, **decision, "same_transcript_chain_in_annotation": same_chain, "psi_available": psi_available, "psi_spearman": psi_rho if psi_rho is not None else "", "psi_n": psi_n}
    write_tsv(args.outdir / "05_evidence_reclassification.tsv", [evidence_row], list(evidence_row))
    report = f"""# TSPAN14 Adjacent-Junction Validation\n\n## Result\n\nThe project event is canonical exon5-exon6 and the external feature is adjacent exon6-exon7 in the TSPAN14-207 annotation. This establishes structural compatibility in an annotated transcript, not within-individual co-usage.\n\nMiGA SVZ summary-QTL slope comparison across {qtl_summary['shared_variant_count']} shared variants gave Pearson r={qtl_summary['slope_pearson_r']:.3f} and same-sign fraction={qtl_summary['same_sign_fraction']:.3f}. These are variant-level association summaries and are not PSI correlations.\n\nSample-level PSI status: {psi_note}\n\n## Evidence decision\n\nClass: {decision['evidence_class']}\n\nSupported claim: {decision['supported_claim']}\n\nRequired manuscript action: {decision['manuscript_action']}\n"""
    (args.outdir / "targeted_splice_validation_report.md").write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main()
