from __future__ import annotations

import csv
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"
REPORTS = ROOT / "results" / "reports"

TARGET_JUNCTION = "chr10:80509471:80512144:clu_4260:+"
CANONICAL_TRANSCRIPT = "ENST00000429989"
CANONICAL_DISPLAY = "TSPAN14-207"


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def upsert_index(rows: list[dict[str, str]]) -> None:
    index = TABLES / "final_output_index.tsv"
    existing = read_tsv(index)
    by_path = {row.get("path", ""): row for row in existing if row.get("path", "")}
    for row in rows:
        by_path[row["path"]] = row
    write_tsv(index, list(by_path.values()), ["type", "name", "path", "status", "description"])


def as_float(value: object, default: float = 0.0) -> float:
    try:
        if value in (None, "", "NA"):
            return default
        return float(str(value))
    except (TypeError, ValueError):
        return default


def transcript_root(value: str) -> str:
    return value.split(".")[0]


def extract_transcripts(text: str) -> list[str]:
    return sorted(set(re.findall(r"ENST\d+(?:\.\d+)?", text)))


def load_frontiers_neuro_dtu() -> dict[str, str]:
    rows = read_tsv(TABLES / "iteration50_frontiers_tspan14_splicing_hits.tsv")
    for row in rows:
        if row.get("sheet") == "Neuro_DTU" and row.get("support_grade") == "fdr_significant_ad_control_gene_level_dtu":
            return row
    raise RuntimeError("Frontiers Neuro_DTU TSPAN14 row not found")


def build_transcript_bridge() -> list[dict[str, object]]:
    frontiers = load_frontiers_neuro_dtu()
    isoform_rows = read_tsv(TABLES / "iteration17_tspan14_isoform_event_map.tsv")
    exact_rows = [row for row in isoform_rows if row.get("junction_support_class") == "exact_leafcutter_junction_boundaries"]
    canonical = next((row for row in exact_rows if row.get("transcript_id") == CANONICAL_TRANSCRIPT), {})
    coding_exact = [row for row in exact_rows if row.get("transcript_biotype") == "protein_coding" and row.get("junction_in_translation_span") == "yes"]
    nmd_exact = [row for row in exact_rows if row.get("transcript_biotype") == "nonsense_mediated_decay"]
    evidence = {row.get("evidence_domain"): row for row in read_tsv(TABLES / "iteration22_isoform_consequence_master_table.tsv")}
    protein = {row.get("evidence_layer"): row for row in read_tsv(TABLES / "iteration18_integrated_splice_protein_consequence.tsv")}
    iso_miga = {row.get("evidence_layer"): row for row in read_tsv(TABLES / "iteration21_external_isoform_replication_summary.tsv")}
    return [
        {
            "bridge_axis": "frontiers_neuro_gene_level_dtu",
            "source": "Frontiers2023 suptables/tableS3.xlsx Neuro_DTU",
            "published_unit": "gene_level_DTU",
            "gene": frontiers.get("gene", "TSPAN14"),
            "p_value": frontiers.get("p_value", ""),
            "q_or_fdr": frontiers.get("fdr_or_q", ""),
            "effect_or_direction": frontiers.get("effect_or_direction", ""),
            "cell_or_sample_context": frontiers.get("cohort_or_celltype", "Neuro_DTU"),
            "transcript_interpretation": "FDR-significant AD/control neuron DTU indicates altered TSPAN14 transcript usage at gene level.",
            "strength": "disease_state_gene_level_transcript_usage_support",
            "claim_boundary": "No transcript ID, exon structure or target junction is supplied in the public significant-gene DTU row; do not claim exact target-junction validation.",
        },
        {
            "bridge_axis": "canonical_target_junction_interpretation",
            "source": "iteration17_tspan14_isoform_event_map.tsv",
            "published_unit": "project_target_junction_annotation",
            "gene": "TSPAN14",
            "p_value": "",
            "q_or_fdr": "",
            "effect_or_direction": evidence.get("allele_harmonized_splice_direction", {}).get("primary_result", ""),
            "cell_or_sample_context": "MiGA_SVZ_sQTL",
            "transcript_interpretation": (
                f"The project target junction maps exactly to canonical {CANONICAL_TRANSCRIPT}/{CANONICAL_DISPLAY} and "
                f"{len(coding_exact)} translated protein-coding TSPAN14 transcript contexts."
            ),
            "strength": "canonical_transcript_context_for_interpreting_gene_level_DTU",
            "claim_boundary": "This transcript map explains biological plausibility of a DTU signal, but it is not the Frontiers event-level DTU result.",
        },
        {
            "bridge_axis": "splice_to_protein_interpretation",
            "source": "iteration18_integrated_splice_protein_consequence.tsv",
            "published_unit": "protein_feature_annotation",
            "gene": "TSPAN14",
            "p_value": "",
            "q_or_fdr": "",
            "effect_or_direction": "",
            "cell_or_sample_context": "canonical_TSPAN14_protein_Q8NG11",
            "transcript_interpretation": protein.get("ADAM10_interaction_region_overlap", {}).get("observation", ""),
            "strength": protein.get("ADAM10_interaction_region_overlap", {}).get("support_level", ""),
            "claim_boundary": protein.get("ADAM10_interaction_region_overlap", {}).get("claim_boundary", ""),
        },
        {
            "bridge_axis": "external_isoform_qtl_context",
            "source": "iteration21_external_isoform_replication_summary.tsv",
            "published_unit": "isoMiGA_top_association",
            "gene": "TSPAN14",
            "p_value": iso_miga.get("isoMiGA canonical transcript regulation", {}).get("best_Fixed_P", ""),
            "q_or_fdr": iso_miga.get("isoMiGA canonical transcript regulation", {}).get("best_Fixed_FDR", ""),
            "effect_or_direction": iso_miga.get("isoMiGA canonical transcript regulation", {}).get("best_fixed_beta", ""),
            "cell_or_sample_context": "isoMiGA",
            "transcript_interpretation": "External isoMiGA top-association supports genetic regulation of canonical ENST00000429989.",
            "strength": iso_miga.get("isoMiGA canonical transcript regulation", {}).get("support_level", ""),
            "claim_boundary": iso_miga.get("isoMiGA canonical transcript regulation", {}).get("claim_boundary", ""),
        },
        {
            "bridge_axis": "resolved_complexity",
            "source": "iteration17_tspan14_isoform_event_map.tsv",
            "published_unit": "transcript_annotation_count",
            "gene": "TSPAN14",
            "p_value": "",
            "q_or_fdr": "",
            "effect_or_direction": "",
            "cell_or_sample_context": "TSPAN14_transcript_model",
            "transcript_interpretation": f"Exact target-junction annotation includes {len(exact_rows)} transcripts, {len(coding_exact)} coding translated contexts and {len(nmd_exact)} NMD contexts.",
            "strength": "supports_transcript_usage_mechanism_with_explicit_complexity",
            "claim_boundary": "Transcript complexity prevents converting DTU or junction usage directly into net TSPAN14 protein abundance.",
        },
    ]


def build_external_context_rows() -> list[dict[str, object]]:
    rows = read_tsv(TABLES / "iteration49_ad_disease_state_splicing_tspan14_hits.tsv")
    out: list[dict[str, object]] = []
    seen: set[tuple[str, str, str]] = set()
    for row in rows:
        feature = row.get("feature", "")
        transcripts = extract_transcripts(feature)
        relation = row.get("target_relation", "")
        if row.get("resource", "").startswith("ExonSkipAD") and transcripts:
            has_canonical = any(transcript_root(t) == CANONICAL_TRANSCRIPT for t in transcripts)
            key = (row.get("resource", ""), relation, ";".join(transcripts[:8]))
            if key in seen:
                continue
            seen.add(key)
            out.append(
                {
                    "source": row.get("resource", ""),
                    "source_file": row.get("source_file", ""),
                    "evidence_type": row.get("evidence_type", ""),
                    "target_relation": relation,
                    "support_grade": row.get("support_grade", ""),
                    "n_transcript_ids_parsed": len(transcripts),
                    "contains_canonical_ENST00000429989": has_canonical,
                    "parsed_transcript_ids": ";".join(transcripts[:12]),
                    "feature_excerpt": feature[:260],
                    "interpretation": (
                        "External AD exon-skipping resources contain TSPAN14 transcript/exon-skipping context that includes the canonical transcript."
                        if has_canonical
                        else "External AD exon-skipping resources contain same-gene transcript/exon-skipping context."
                    ),
                    "claim_boundary": row.get("claim_boundary", ""),
                }
            )
    return out


def build_summary(bridge_rows: list[dict[str, object]], external_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    frontiers = next(row for row in bridge_rows if row.get("bridge_axis") == "frontiers_neuro_gene_level_dtu")
    canonical_external = [row for row in external_rows if str(row.get("contains_canonical_ENST00000429989")) == "True"]
    target_overlap_external = [row for row in external_rows if row.get("target_relation") == "target_coordinate_overlap"]
    return [
        {
            "summary_item": "gene_level_dtu_transcript_interpretation",
            "frontiers_neuro_dtu_q": frontiers.get("q_or_fdr", ""),
            "n_bridge_rows": len(bridge_rows),
            "n_external_transcript_context_rows": len(external_rows),
            "n_external_canonical_transcript_rows": len(canonical_external),
            "n_external_target_coordinate_overlap_rows": len(target_overlap_external),
            "final_status": "gene_level_dtu_interpretable_as_transcript_usage_support",
            "allowed_claim": "Frontiers Neuro_DTU provides FDR-significant AD/control neuron gene-level TSPAN14 transcript-usage support that is biologically interpretable through the project's canonical TSPAN14 isoform map.",
            "claim_boundary": "This does not close exact target-junction disease-state validation because the public Frontiers significant DTU table is gene-level and lacks event coordinates or transcript-specific abundance.",
        },
        {
            "summary_item": "mechanistic_use",
            "frontiers_neuro_dtu_q": frontiers.get("q_or_fdr", ""),
            "n_bridge_rows": len(bridge_rows),
            "n_external_transcript_context_rows": len(external_rows),
            "n_external_canonical_transcript_rows": len(canonical_external),
            "n_external_target_coordinate_overlap_rows": len(target_overlap_external),
            "final_status": "transcript_story_strengthened_without_target_junction_overclaim",
            "allowed_claim": "Use the DTU result as disease-state transcript-usage context that converges with canonical TSPAN14-207, exon5-exon6, AA150/151 and ADAM10-region isoform biology.",
            "claim_boundary": "Do not present this as direct AD/control validation of MiGA SVZ junction chr10:80509471-80512144 or as protein-level ADAM10/APP/TREM2 evidence.",
        },
    ]


def update_evidence_model(summary_rows: list[dict[str, object]]) -> None:
    path = TABLES / "iteration39_mechanism_figure_edges.tsv"
    rows = read_tsv(path)
    if not rows:
        return
    q_value = summary_rows[0].get("frontiers_neuro_dtu_q", "")
    for row in rows:
        if row.get("edge_id") == "E5":
            label = row.get("figure_label", "")
            if "gene-level DTU transcript interpretation" not in label:
                row["figure_label"] = f"{label}; gene-level DTU transcript interpretation q={q_value}"
            files = row.get("evidence_files", "")
            prefix = (
                "iteration64_gene_level_dtu_transcript_interpretation_summary.tsv; "
                "iteration64_gene_level_dtu_transcript_bridge.tsv; "
                "iteration64_external_ad_transcript_context.tsv"
            )
            if "iteration64_gene_level_dtu_transcript_interpretation_summary.tsv" not in files:
                row["evidence_files"] = f"{prefix}; {files}"
            if "gene-level DTU is transcript-usage context" not in row.get("claim_boundary", ""):
                row["claim_boundary"] = (
                    row.get("claim_boundary", "")
                    + " Iteration64 clarifies that gene-level DTU is transcript-usage context aligned with isoform biology, not direct target-junction validation."
                )
    write_tsv(
        path,
        rows,
        [
            "edge_id",
            "from_node",
            "to_node",
            "main_statement",
            "evidence_strength",
            "figure_line_style",
            "figure_label",
            "evidence_files",
            "must_not_draw",
            "claim_boundary",
        ],
    )


def write_report(bridge_rows: list[dict[str, object]], external_rows: list[dict[str, object]], summary: list[dict[str, object]]) -> None:
    lines = [
        "# Iteration64 gene-level DTU transcript interpretation",
        "",
        "## Bottom line",
        "",
        (
            "Frontiers Neuro_DTU can now be used as disease-state gene-level transcript-usage support for TSPAN14. "
            "The interpretation is mechanistically meaningful because the project target splice event already maps to canonical "
            "TSPAN14-207/ENST00000429989, coding exon5-exon6, AA150/151 and the UniProt ADAM10-interaction region."
        ),
        "",
        "The boundary is unchanged: the Frontiers public significant DTU table is gene-level and does not provide event coordinates or transcript-specific abundances for the target MiGA SVZ junction.",
        "",
        "## Summary",
        "",
        "| Item | Status | Allowed claim | Boundary |",
        "|---|---|---|---|",
    ]
    for row in summary:
        lines.append(f"| {row.get('summary_item')} | {row.get('final_status')} | {row.get('allowed_claim')} | {row.get('claim_boundary')} |")
    lines.extend(["", "## Evidence bridge", "", "| Axis | Interpretation | Strength | Boundary |", "|---|---|---|---|"])
    for row in bridge_rows:
        lines.append(
            f"| {row.get('bridge_axis')} | {row.get('transcript_interpretation')} | {row.get('strength')} | {row.get('claim_boundary')} |"
        )
    lines.extend(
        [
            "",
            "## Outputs",
            "",
            "- results/tables/iteration64_gene_level_dtu_transcript_bridge.tsv",
            "- results/tables/iteration64_external_ad_transcript_context.tsv",
            "- results/tables/iteration64_gene_level_dtu_transcript_interpretation_summary.tsv",
        ]
    )
    REPORTS.mkdir(parents=True, exist_ok=True)
    (REPORTS / "iteration64_gene_level_dtu_transcript_interpretation.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    bridge_rows = build_transcript_bridge()
    external_rows = build_external_context_rows()
    summary = build_summary(bridge_rows, external_rows)
    write_tsv(
        TABLES / "iteration64_gene_level_dtu_transcript_bridge.tsv",
        bridge_rows,
        [
            "bridge_axis",
            "source",
            "published_unit",
            "gene",
            "p_value",
            "q_or_fdr",
            "effect_or_direction",
            "cell_or_sample_context",
            "transcript_interpretation",
            "strength",
            "claim_boundary",
        ],
    )
    write_tsv(
        TABLES / "iteration64_external_ad_transcript_context.tsv",
        external_rows,
        [
            "source",
            "source_file",
            "evidence_type",
            "target_relation",
            "support_grade",
            "n_transcript_ids_parsed",
            "contains_canonical_ENST00000429989",
            "parsed_transcript_ids",
            "feature_excerpt",
            "interpretation",
            "claim_boundary",
        ],
    )
    write_tsv(
        TABLES / "iteration64_gene_level_dtu_transcript_interpretation_summary.tsv",
        summary,
        [
            "summary_item",
            "frontiers_neuro_dtu_q",
            "n_bridge_rows",
            "n_external_transcript_context_rows",
            "n_external_canonical_transcript_rows",
            "n_external_target_coordinate_overlap_rows",
            "final_status",
            "allowed_claim",
            "claim_boundary",
        ],
    )
    update_evidence_model(summary)
    write_report(bridge_rows, external_rows, summary)
    upsert_index(
        [
            {
                "type": "table",
                "name": "iteration64_gene_level_dtu_transcript_bridge",
                "path": "results/tables/iteration64_gene_level_dtu_transcript_bridge.tsv",
                "status": "current",
                "description": "Bridge from Frontiers gene-level Neuro_DTU to the existing TSPAN14 transcript and isoform interpretation.",
            },
            {
                "type": "table",
                "name": "iteration64_external_ad_transcript_context",
                "path": "results/tables/iteration64_external_ad_transcript_context.tsv",
                "status": "current",
                "description": "Parsed external AD exon-skipping transcript context for TSPAN14.",
            },
            {
                "type": "table",
                "name": "iteration64_gene_level_dtu_transcript_interpretation_summary",
                "path": "results/tables/iteration64_gene_level_dtu_transcript_interpretation_summary.tsv",
                "status": "current",
                "description": "Decision summary for using gene-level DTU as transcript-usage support without target-junction overclaim.",
            },
            {
                "type": "report",
                "name": "iteration64_gene_level_dtu_transcript_interpretation",
                "path": "results/reports/iteration64_gene_level_dtu_transcript_interpretation.md",
                "status": "current",
                "description": "Report documenting transcript-level interpretation of the Frontiers gene-level DTU evidence.",
            },
        ]
    )
    print(f"Wrote {len(bridge_rows)} bridge rows, {len(external_rows)} external transcript-context rows and {len(summary)} summary rows.")


if __name__ == "__main__":
    main()
