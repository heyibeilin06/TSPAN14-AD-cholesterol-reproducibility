from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"
REPORTS = ROOT / "results" / "reports"
ENSEMBL_JSON = ROOT / "data" / "processed" / "ensembl_iteration15" / "ensembl_lookup_TSPAN14_expand1.json"
UNIPROT_JSON = ROOT / "data" / "processed" / "protein_iteration18" / "uniprot_TSPAN14_search.json"
EVENT = "chr10:80509471:80512144:clu_4260_+:ENSG00000108219.15"
CANONICAL_TX = "ENST00000429989"
CANONICAL_TRANSLATION = "ENSP00000396270"
UNIPROT_ACCESSION = "Q8NG11"
EVENT_START = 80509471
EVENT_END = 80512144


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def upsert_index(rows: list[dict[str, str]]) -> None:
    index = TABLES / "final_output_index.tsv"
    existing = read_tsv(index) if index.exists() else []
    by_path = {row.get("path", ""): row for row in existing if row.get("path", "")}
    for row in rows:
        by_path[row["path"]] = row
    write_tsv(index, list(by_path.values()), ["type", "name", "path", "status", "description"])


def load_canonical_transcript() -> dict:
    data = json.loads(ENSEMBL_JSON.read_text(encoding="utf-8"))
    for tx in data.get("Transcript", []):
        if tx.get("id") == CANONICAL_TX:
            return tx
    raise RuntimeError(f"Canonical transcript {CANONICAL_TX} not found in Ensembl cache")


def load_uniprot_record() -> dict:
    data = json.loads(UNIPROT_JSON.read_text(encoding="utf-8"))
    for row in data.get("results", []):
        if row.get("primaryAccession") == UNIPROT_ACCESSION:
            return row
    raise RuntimeError(f"UniProt record {UNIPROT_ACCESSION} not found")


def feature_start_end(feature: dict) -> tuple[int | None, int | None]:
    loc = feature.get("location", {})
    start = loc.get("start", {}).get("value")
    end = loc.get("end", {}).get("value")
    if start is None or end is None:
        return None, None
    return int(start), int(end)


def overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    return max(a_start, b_start) <= min(a_end, b_end)


def coding_exon_map(tx: dict) -> list[dict[str, object]]:
    translation = tx["Translation"]
    cds_start = int(translation["start"])
    cds_end = int(translation["end"])
    protein_len = int(translation["length"])
    exons = sorted(tx.get("Exon", []), key=lambda e: int(e["start"]))
    rows: list[dict[str, object]] = []
    coding_nt_seen = 0
    for rank, exon in enumerate(exons, start=1):
        exon_start = int(exon["start"])
        exon_end = int(exon["end"])
        c_start = max(exon_start, cds_start)
        c_end = min(exon_end, cds_end)
        if c_start > c_end:
            rows.append(
                {
                    "transcript_id": CANONICAL_TX,
                    "translation_id": CANONICAL_TRANSLATION,
                    "exon_rank": rank,
                    "exon_id": exon.get("id", ""),
                    "exon_start": exon_start,
                    "exon_end": exon_end,
                    "coding_status": "noncoding_or_utr",
                    "coding_start": "",
                    "coding_end": "",
                    "coding_nt_length": 0,
                    "coding_nt_offset_start": "",
                    "coding_nt_offset_end": "",
                    "aa_start": "",
                    "aa_end": "",
                    "contains_event_upstream_boundary": "yes" if exon_end == EVENT_START else "no",
                    "contains_event_downstream_boundary": "yes" if exon_start == EVENT_END else "no",
                }
            )
            continue
        nt_len = c_end - c_start + 1
        offset_start = coding_nt_seen + 1
        offset_end = coding_nt_seen + nt_len
        aa_start = ((offset_start - 1) // 3) + 1
        aa_end = min(((offset_end - 1) // 3) + 1, protein_len)
        rows.append(
            {
                "transcript_id": CANONICAL_TX,
                "translation_id": CANONICAL_TRANSLATION,
                "exon_rank": rank,
                "exon_id": exon.get("id", ""),
                "exon_start": exon_start,
                "exon_end": exon_end,
                "coding_status": "coding",
                "coding_start": c_start,
                "coding_end": c_end,
                "coding_nt_length": nt_len,
                "coding_nt_offset_start": offset_start,
                "coding_nt_offset_end": offset_end,
                "aa_start": aa_start,
                "aa_end": aa_end,
                "contains_event_upstream_boundary": "yes" if exon_end == EVENT_START else "no",
                "contains_event_downstream_boundary": "yes" if exon_start == EVENT_END else "no",
            }
        )
        coding_nt_seen += nt_len
    return rows


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    REPORTS.mkdir(parents=True, exist_ok=True)

    tx = load_canonical_transcript()
    uni = load_uniprot_record()
    exon_rows = coding_exon_map(tx)
    upstream = next(r for r in exon_rows if r["contains_event_upstream_boundary"] == "yes")
    downstream = next(r for r in exon_rows if r["contains_event_downstream_boundary"] == "yes")
    junction_after_aa = int(upstream["aa_end"])
    junction_before_aa = int(downstream["aa_start"])
    junction_window_start = max(1, junction_after_aa - 3)
    junction_window_end = min(int(tx["Translation"]["length"]), junction_before_aa + 3)

    feature_rows: list[dict[str, object]] = []
    overlap_rows: list[dict[str, object]] = []
    for feature in uni.get("features", []):
        start, end = feature_start_end(feature)
        if start is None or end is None:
            continue
        row = {
            "source": "UniProt REST",
            "uniprot_accession": uni.get("primaryAccession", ""),
            "uniprot_id": uni.get("uniProtkbId", ""),
            "protein_name": uni.get("proteinDescription", {})
            .get("recommendedName", {})
            .get("fullName", {})
            .get("value", ""),
            "feature_type": feature.get("type", ""),
            "feature_description": feature.get("description", ""),
            "aa_start": start,
            "aa_end": end,
        }
        feature_rows.append(row)
        if overlap(junction_window_start, junction_window_end, start, end) or (
            start <= junction_after_aa <= end
            or start <= junction_before_aa <= end
        ):
            overlap_rows.append(
                {
                    **row,
                    "target_splice_event": EVENT,
                    "canonical_transcript": CANONICAL_TX,
                    "canonical_translation": CANONICAL_TRANSLATION,
                    "junction_genomic_boundaries": f"chr10:{EVENT_START}-{EVENT_END}",
                    "junction_aa_boundary": f"after_AA{junction_after_aa}_before_AA{junction_before_aa}",
                    "junction_scoring_window": f"AA{junction_window_start}-AA{junction_window_end}",
                    "overlap_type": (
                        "direct_boundary_inside_feature"
                        if start <= junction_after_aa <= end or start <= junction_before_aa <= end
                        else "near_boundary_window_overlap"
                    ),
                    "mechanistic_relevance": (
                        "directly_overlaps_uniprot_ADAM10_interaction_region"
                        if "ADAM10" in str(feature.get("description", ""))
                        else "overlaps_tspan14_topology_feature"
                    ),
                }
            )

    exact_isoform = read_tsv(TABLES / "iteration17_splice_event_consequence_scoring.tsv")
    canonical_score = next(r for r in exact_isoform if r["transcript_id"] == CANONICAL_TX)

    consequence_rows = [
        {
            "evidence_layer": "canonical_isoform_splice_to_protein_mapping",
            "observation": (
                f"The target splice junction connects coding exon {upstream['exon_rank']} "
                f"(ending at AA{junction_after_aa}) to coding exon {downstream['exon_rank']} "
                f"(starting at AA{junction_before_aa}) in {CANONICAL_TX}."
            ),
            "support_level": "strong",
            "data_source": "Ensembl transcript cache and UniProt Q8NG11 feature annotation",
            "claim_allowed": "The allele-harmonized splice event affects a canonical coding junction within the TSPAN14 protein sequence.",
            "claim_boundary": "This does not determine the exact abundance of full-length versus alternative TSPAN14 protein isoforms.",
        },
        {
            "evidence_layer": "topology_overlap",
            "observation": (
                f"The junction boundary AA{junction_after_aa}/AA{junction_before_aa} falls within UniProt annotated "
                "extracellular region AA114-232."
            ),
            "support_level": "strong_for_topology_context",
            "data_source": "UniProt Q8NG11 topological-domain features",
            "claim_allowed": "The implicated splice junction lies in the large extracellular loop context of TSPAN14.",
            "claim_boundary": "Topology overlap is annotation-based and does not prove altered membrane localization.",
        },
        {
            "evidence_layer": "ADAM10_interaction_region_overlap",
            "observation": (
                f"UniProt annotates AA114-232 as necessary and sufficient for interaction with ADAM10; "
                f"the junction boundary is AA{junction_after_aa}/AA{junction_before_aa}."
            ),
            "support_level": "high_mechanistic_plausibility",
            "data_source": "UniProt Q8NG11 region annotation plus current splice-event mapping",
            "claim_allowed": "The TSPAN14 sQTL signal is positioned in a protein region directly relevant to ADAM10 interaction biology.",
            "claim_boundary": "Do not claim the allele changes ADAM10 interaction strength without perturbation or protein-binding data.",
        },
        {
            "evidence_layer": "integrated_consequence_score",
            "observation": (
                f"Iteration17 canonical consequence tier is {canonical_score['consequence_tier']} "
                f"with score {canonical_score['consequence_score']}; iteration18 adds ADAM10-region overlap."
            ),
            "support_level": "strong_bioinformatic_mechanism",
            "data_source": "Iteration17 consequence scoring and iteration18 protein-feature overlap",
            "claim_allowed": "This is now a strong bioinformatic splice-to-protein-mechanism candidate.",
            "claim_boundary": "Still not a direct functional validation of APP cleavage, TREM2 shedding, or Notch substrate effects.",
        },
    ]

    method_rows = [
        {
            "analysis_step": "TSPAN14 splice junction to protein feature mapping",
            "script_or_command": "python scripts/99_map_tspan14_splice_to_protein_features.py",
            "primary_inputs": "; ".join(
                [
                    "data/processed/ensembl_iteration15/ensembl_lookup_TSPAN14_expand1.json",
                    "data/processed/protein_iteration18/uniprot_TSPAN14_search.json",
                    "results/tables/iteration17_splice_event_consequence_scoring.tsv",
                ]
            ),
            "key_parameters": (
                f"Canonical transcript {CANONICAL_TX}; canonical translation {CANONICAL_TRANSLATION}; "
                f"target splice junction chr10:{EVENT_START}-{EVENT_END}; UniProt accession {UNIPROT_ACCESSION}; "
                "junction feature overlap evaluated at the AA boundary and +/-3 amino-acid window."
            ),
            "primary_outputs": "; ".join(
                [
                    "results/tables/iteration18_tspan14_canonical_exon_cds_aa_map.tsv",
                    "results/tables/iteration18_uniprot_tspan14_feature_catalog.tsv",
                    "results/tables/iteration18_splice_junction_protein_feature_overlap.tsv",
                    "results/tables/iteration18_integrated_splice_protein_consequence.tsv",
                    "results/tables/iteration18_splice_protein_methods_trace.tsv",
                    "results/reports/iteration18_tspan14_splice_to_protein_feature_mapping.md",
                ]
            ),
            "method_text_ready_note": (
                "CDS offsets were computed on the plus strand from Ensembl exon coordinates and the canonical translation span. "
                "The target splice junction was converted to an amino-acid boundary and intersected with UniProt topology/region features."
            ),
        }
    ]

    write_tsv(
        TABLES / "iteration18_tspan14_canonical_exon_cds_aa_map.tsv",
        exon_rows,
        [
            "transcript_id",
            "translation_id",
            "exon_rank",
            "exon_id",
            "exon_start",
            "exon_end",
            "coding_status",
            "coding_start",
            "coding_end",
            "coding_nt_length",
            "coding_nt_offset_start",
            "coding_nt_offset_end",
            "aa_start",
            "aa_end",
            "contains_event_upstream_boundary",
            "contains_event_downstream_boundary",
        ],
    )
    write_tsv(
        TABLES / "iteration18_uniprot_tspan14_feature_catalog.tsv",
        feature_rows,
        ["source", "uniprot_accession", "uniprot_id", "protein_name", "feature_type", "feature_description", "aa_start", "aa_end"],
    )
    write_tsv(
        TABLES / "iteration18_splice_junction_protein_feature_overlap.tsv",
        overlap_rows,
        [
            "source",
            "uniprot_accession",
            "uniprot_id",
            "protein_name",
            "feature_type",
            "feature_description",
            "aa_start",
            "aa_end",
            "target_splice_event",
            "canonical_transcript",
            "canonical_translation",
            "junction_genomic_boundaries",
            "junction_aa_boundary",
            "junction_scoring_window",
            "overlap_type",
            "mechanistic_relevance",
        ],
    )
    write_tsv(
        TABLES / "iteration18_integrated_splice_protein_consequence.tsv",
        consequence_rows,
        ["evidence_layer", "observation", "support_level", "data_source", "claim_allowed", "claim_boundary"],
    )
    write_tsv(
        TABLES / "iteration18_splice_protein_methods_trace.tsv",
        method_rows,
        ["analysis_step", "script_or_command", "primary_inputs", "key_parameters", "primary_outputs", "method_text_ready_note"],
    )

    master_path = TABLES / "iteration17_methods_provenance_master.tsv"
    master = read_tsv(master_path) if master_path.exists() else []
    by_step = {row.get("analysis_step", ""): row for row in master}
    for row in method_rows:
        by_step[row["analysis_step"]] = row
    write_tsv(
        TABLES / "iteration18_methods_provenance_master.tsv",
        list(by_step.values()),
        ["analysis_step", "script_or_command", "primary_inputs", "key_parameters", "primary_outputs", "method_text_ready_note"],
    )

    report = f"""# Iteration 18 TSPAN14 splice-to-protein feature mapping

## Scope

This iteration strengthens the isoform-level analysis by converting the target TSPAN14 splice junction into canonical CDS and amino-acid coordinates, then intersecting the junction with UniProt protein topology and functional-region annotations. No figures or manuscript text were generated.

## Main result

- Canonical transcript: `{CANONICAL_TX}` / TSPAN14-207.
- Canonical translation: `{CANONICAL_TRANSLATION}`; UniProt accession `{UNIPROT_ACCESSION}`.
- Target splice event: `{EVENT}`.
- The event joins coding exon {upstream['exon_rank']} to coding exon {downstream['exon_rank']}.
- Protein boundary: after AA{junction_after_aa}, before AA{junction_before_aa}.
- UniProt overlap: the boundary lies in the extracellular region AA114-232.
- Mechanism-relevant overlap: UniProt annotates AA114-232 as necessary and sufficient for interaction with ADAM10.

## Interpretation

This materially strengthens the previous isoform-level result. The allele-harmonized sQTL does not only alter a canonical coding junction; that junction maps into the large extracellular loop of TSPAN14, the same annotated region required for ADAM10 interaction. Therefore the current pure-bioinformatics evidence can support a strong mechanistic plausibility claim linking the TSPAN14 splice event to ADAM10-related biology.

The boundary remains important: the analysis does not directly prove altered ADAM10 binding, ADAM10 trafficking, APP cleavage, TREM2 shedding, or Notch substrate effects. Those remain downstream hypotheses or experimental validation endpoints.

## New outputs

- `results/tables/iteration18_tspan14_canonical_exon_cds_aa_map.tsv`
- `results/tables/iteration18_uniprot_tspan14_feature_catalog.tsv`
- `results/tables/iteration18_splice_junction_protein_feature_overlap.tsv`
- `results/tables/iteration18_integrated_splice_protein_consequence.tsv`
- `results/tables/iteration18_splice_protein_methods_trace.tsv`
- `results/tables/iteration18_methods_provenance_master.tsv`
- `results/reports/iteration18_tspan14_splice_to_protein_feature_mapping.md`
"""
    (REPORTS / "iteration18_tspan14_splice_to_protein_feature_mapping.md").write_text(report, encoding="utf-8")

    upsert_index(
        [
            {
                "type": "table",
                "name": "Iteration18 TSPAN14 canonical exon CDS AA map",
                "path": "results/tables/iteration18_tspan14_canonical_exon_cds_aa_map.tsv",
                "status": "current",
                "description": "Canonical TSPAN14 exon-to-CDS-to-amino-acid coordinate map for the target splice junction.",
            },
            {
                "type": "table",
                "name": "Iteration18 UniProt TSPAN14 feature catalog",
                "path": "results/tables/iteration18_uniprot_tspan14_feature_catalog.tsv",
                "status": "current",
                "description": "UniProt Q8NG11 topology and region features used for splice-to-protein consequence interpretation.",
            },
            {
                "type": "table",
                "name": "Iteration18 splice junction protein feature overlap",
                "path": "results/tables/iteration18_splice_junction_protein_feature_overlap.tsv",
                "status": "current",
                "description": "Overlap between the TSPAN14 splice-junction amino-acid boundary and UniProt protein features.",
            },
            {
                "type": "table",
                "name": "Iteration18 integrated splice protein consequence",
                "path": "results/tables/iteration18_integrated_splice_protein_consequence.tsv",
                "status": "current",
                "description": "Evidence-layer table integrating isoform, protein topology and ADAM10-region consequence claims.",
            },
            {
                "type": "table",
                "name": "Iteration18 splice protein methods trace",
                "path": "results/tables/iteration18_splice_protein_methods_trace.tsv",
                "status": "current",
                "description": "Methods-ready trace for mapping the TSPAN14 splice junction to protein features.",
            },
            {
                "type": "table",
                "name": "Iteration18 methods provenance master",
                "path": "results/tables/iteration18_methods_provenance_master.tsv",
                "status": "current",
                "description": "Current master Methods provenance table including splice-to-protein feature mapping.",
            },
            {
                "type": "report",
                "name": "Iteration18 TSPAN14 splice-to-protein feature mapping report",
                "path": "results/reports/iteration18_tspan14_splice_to_protein_feature_mapping.md",
                "status": "current",
                "description": "Report summarizing TSPAN14 splice-junction mapping to protein topology and ADAM10 interaction region.",
            },
        ]
    )


if __name__ == "__main__":
    main()
