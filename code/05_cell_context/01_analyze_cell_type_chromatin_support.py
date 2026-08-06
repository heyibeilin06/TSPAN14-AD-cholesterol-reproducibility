from __future__ import annotations

import csv
import json
import urllib.request
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"
REPORTS = ROOT / "results" / "reports"
CACHE = ROOT / "data" / "processed" / "external_metadata_iteration63_chromatin"

CHROM = "chr10"
WINDOW_START = 80480000
WINDOW_END = 80526000
TARGET_SPLICE_START = 80509471
TARGET_SPLICE_END = 80512144
UCSC_CCRE_URL = (
    "https://api.genome.ucsc.edu/getData/track?"
    f"genome=hg38;track=encodeCcreCombined;chrom={CHROM};start={WINDOW_START};end={WINDOW_END}"
)

PRIMARY_CELL_AXES = ["microglia", "neuron", "endothelial"]


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


def as_int(value: object, default: int = 0) -> int:
    try:
        if value in (None, "", "NA"):
            return default
        return int(float(str(value)))
    except (TypeError, ValueError):
        return default


def as_float(value: object, default: float = 0.0) -> float:
    try:
        if value in (None, "", "NA"):
            return default
        return float(str(value))
    except (TypeError, ValueError):
        return default


def fetch_or_load_ucsc_ccre() -> dict[str, Any]:
    CACHE.mkdir(parents=True, exist_ok=True)
    cache_path = CACHE / "ucsc_encode_ccre_tspan14_window.json"
    legacy_path = ROOT / "data" / "processed" / "ucsc_encode_ccre_tspan14_window.json"
    if not cache_path.exists() and legacy_path.exists():
        cache_path.write_text(legacy_path.read_text(encoding="utf-8-sig"), encoding="utf-8")
    if not cache_path.exists():
        req = urllib.request.Request(UCSC_CCRE_URL, headers={"User-Agent": "SLM-TSPAN14-chromatin-audit/1.0"})
        with urllib.request.urlopen(req, timeout=60) as response:
            payload = response.read().decode("utf-8-sig")
        cache_path.write_text(payload, encoding="utf-8")
    return json.loads(cache_path.read_text(encoding="utf-8-sig"))


def is_enhancer_like(row: dict[str, Any]) -> bool:
    ccre = str(row.get("ccre", ""))
    label = str(row.get("ucscLabel", ""))
    return ccre in {"dELS", "pELS", "PLS"} or label.startswith("enh")


def overlap_bp(start1: int, end1: int, start2: int, end2: int) -> int:
    return max(0, min(end1, end2) - max(start1, start2))


def distance_to_interval(pos: int, start: int, end: int) -> int:
    if start <= pos <= end:
        return 0
    return min(abs(pos - start), abs(pos - end))


def distance_between_intervals(start1: int, end1: int, start2: int, end2: int) -> int:
    if overlap_bp(start1, end1, start2, end2) > 0:
        return 0
    if end1 < start2:
        return start2 - end1
    return start1 - end2


def load_snp_rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for row in read_tsv(TABLES / "iteration62_ld_block_functional_annotation.tsv"):
        pos = as_int(row.get("pos_hg38"))
        if pos:
            out: dict[str, object] = dict(row)
            out["pos_hg38"] = pos
            rows.append(out)
    return rows


def build_ccre_rows(ccres: list[dict[str, Any]], snps: list[dict[str, object]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for item in sorted(ccres, key=lambda x: as_int(x.get("chromStart"))):
        start = as_int(item.get("chromStart"))
        end = as_int(item.get("chromEnd"))
        contained = [
            str(snp.get("snp"))
            for snp in snps
            if start <= as_int(snp.get("pos_hg38")) <= end
        ]
        nearest_snp = ""
        nearest_distance = ""
        if snps:
            nearest = min(snps, key=lambda snp: distance_to_interval(as_int(snp.get("pos_hg38")), start, end))
            nearest_snp = str(nearest.get("snp"))
            nearest_distance = distance_to_interval(as_int(nearest.get("pos_hg38")), start, end)
        rows.append(
            {
                "source": "UCSC Genome Browser API encodeCcreCombined hg38",
                "source_url": UCSC_CCRE_URL,
                "chrom": item.get("chrom", CHROM),
                "chromStart_0based": start,
                "chromEnd_0based": end,
                "ccre_id": item.get("name", ""),
                "ccre_class": item.get("ccre", ""),
                "encode_label": item.get("encodeLabel", ""),
                "ucsc_label": item.get("ucscLabel", ""),
                "description": item.get("description", ""),
                "score": item.get("score", ""),
                "zScore": item.get("zScore", ""),
                "is_enhancer_like": is_enhancer_like(item),
                "overlaps_target_splice_interval": overlap_bp(start, end, TARGET_SPLICE_START, TARGET_SPLICE_END) > 0,
                "target_splice_overlap_bp": overlap_bp(start, end, TARGET_SPLICE_START, TARGET_SPLICE_END),
                "distance_to_target_splice_interval_bp": distance_between_intervals(start, end, TARGET_SPLICE_START, TARGET_SPLICE_END),
                "contains_core_or_extended_snp": bool(contained),
                "contained_snps": ";".join(contained),
                "nearest_core_or_extended_snp": nearest_snp,
                "nearest_snp_distance_bp": nearest_distance,
                "evidence_boundary": "ENCODE cCRE supports regulatory-element context in the LD block, not disease-state cell-type specificity.",
            }
        )
    return rows


def build_snp_rows(ccre_rows: list[dict[str, object]], snps: list[dict[str, object]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for snp in snps:
        pos = as_int(snp.get("pos_hg38"))
        containing = [
            c
            for c in ccre_rows
            if as_int(c.get("chromStart_0based")) <= pos <= as_int(c.get("chromEnd_0based"))
        ]
        nearest = min(
            ccre_rows,
            key=lambda c: distance_to_interval(pos, as_int(c.get("chromStart_0based")), as_int(c.get("chromEnd_0based"))),
        )
        enhancer_containing = [c for c in containing if str(c.get("is_enhancer_like")) == "True"]
        rows.append(
            {
                "snp": snp.get("snp", ""),
                "snp_role": snp.get("snp_role", ""),
                "pos_hg38": pos,
                "functional_priority_grade": snp.get("functional_priority_grade", ""),
                "inside_target_splice_interval": snp.get("relation_to_target_splice_interval") == "inside_target_splice_interval",
                "n_overlapping_ccres": len(containing),
                "n_overlapping_enhancer_like_ccres": len(enhancer_containing),
                "overlapping_ccre_ids": ";".join(str(c.get("ccre_id")) for c in containing),
                "overlapping_ccre_classes": ";".join(sorted({str(c.get("ccre_class")) for c in containing if c.get("ccre_class")})),
                "nearest_ccre_id": nearest.get("ccre_id", ""),
                "nearest_ccre_class": nearest.get("ccre_class", ""),
                "nearest_ccre_is_enhancer_like": nearest.get("is_enhancer_like", ""),
                "nearest_ccre_distance_bp": distance_to_interval(
                    pos, as_int(nearest.get("chromStart_0based")), as_int(nearest.get("chromEnd_0based"))
                ),
                "claim_use": "Variant-level cCRE overlap/nearest context only; retain LD-block, not single-SNP causal, interpretation.",
            }
        )
    return rows


def build_summary(ccre_rows: list[dict[str, object]], snp_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    cell61 = {row.get("cell_axis", ""): row for row in read_tsv(TABLES / "iteration61_cell_type_localization_summary.tsv")}
    n_ccre = len(ccre_rows)
    enhancer_like = [row for row in ccre_rows if str(row.get("is_enhancer_like")) == "True"]
    target_overlap = [row for row in ccre_rows if str(row.get("overlaps_target_splice_interval")) == "True"]
    snp_overlap = [row for row in snp_rows if as_int(row.get("n_overlapping_ccres")) > 0]
    rows: list[dict[str, object]] = [
        {
            "axis": "regulatory_block",
            "n_ucsc_ccres_in_window": n_ccre,
            "n_enhancer_like_ccres": len(enhancer_like),
            "n_target_splice_overlapping_ccres": len(target_overlap),
            "n_core_or_extended_snps_with_ccre_overlap": len(snp_overlap),
            "cell_type_chromatin_specificity": "not_cell_type_specific_combined_cCRE_track",
            "paired_cell_expression_context": "",
            "final_status": "ld_block_chromatin_regulatory_context_supported",
            "allowed_claim": "The TSPAN14 LD block lies in a regulatory cCRE-rich region, including enhancer-like cCREs and cCREs near or overlapping the target splice interval.",
            "claim_boundary": "UCSC encodeCcreCombined is a combined regulatory annotation track; do not call it cell-type-specific chromatin evidence.",
        }
    ]
    for axis in PRIMARY_CELL_AXES:
        loc = cell61.get(axis, {})
        if axis == "microglia":
            status = "direct_microglia_crispri_chromatin_anchor_plus_ccre_block"
            specificity = "Laub_microglia_CRISPRi_functional_anchor; UCSC_cCRE_block_context"
            allowed = (
                "Microglial plausibility is strengthened because the LD block has ENCODE cCRE regulatory context and the Laub overlap provides "
                "microglia CRISPRi functional anchoring for the same TSPAN14 enhancer block."
            )
            boundary = (
                "The direct microglia statement comes from Laub CRISPRi/enhancer evidence plus project cell-type expression context; "
                "UCSC cCRE itself is not a microglia-specific ATAC/ChromHMM result."
            )
        elif axis == "neuron":
            status = "ccre_block_plus_public_neuron_disease_context"
            specificity = "UCSC_cCRE_block_context; no_direct_neuron_chromatin_track_in_current_audit"
            allowed = (
                "Neuron relevance can be framed as public neuron AD-state RNA/DTU support occurring on a TSPAN14 locus that is independently annotated as regulatory."
            )
            boundary = "Do not claim neuron-specific chromatin opening or enhancer activity from current cCRE-only evidence."
        else:
            status = "ccre_block_plus_endothelial_expression_context"
            specificity = "UCSC_cCRE_block_context; no_direct_endothelial_chromatin_track_in_current_audit"
            allowed = (
                "Endothelial relevance can be framed as AD427 endothelial detectability occurring on a TSPAN14 locus that is independently annotated as regulatory."
            )
            boundary = "Do not claim endothelial-specific chromatin opening or AD-state endothelial differential expression from current parsed evidence."
        rows.append(
            {
                "axis": axis,
                "n_ucsc_ccres_in_window": n_ccre,
                "n_enhancer_like_ccres": len(enhancer_like),
                "n_target_splice_overlapping_ccres": len(target_overlap),
                "n_core_or_extended_snps_with_ccre_overlap": len(snp_overlap),
                "cell_type_chromatin_specificity": specificity,
                "paired_cell_expression_context": loc.get("final_status", ""),
                "final_status": status,
                "allowed_claim": allowed,
                "claim_boundary": boundary,
            }
        )
    return rows


def update_evidence_model(summary_rows: list[dict[str, object]]) -> None:
    path = TABLES / "iteration39_mechanism_figure_edges.tsv"
    rows = read_tsv(path)
    if not rows:
        return
    block = next((row for row in summary_rows if row.get("axis") == "regulatory_block"), {})
    ccre_n = block.get("n_ucsc_ccres_in_window", "")
    enhancer_n = block.get("n_enhancer_like_ccres", "")
    for row in rows:
        if row.get("edge_id") == "E2":
            label = row.get("figure_label", "")
            if "ENCODE cCREs=" not in label:
                row["figure_label"] = f"{label}; ENCODE cCREs={ccre_n}; enhancer-like={enhancer_n}"
            files = row.get("evidence_files", "")
            prefix = (
                "iteration63_cell_type_chromatin_support_summary.tsv; "
                "iteration63_cell_type_chromatin_ccre_overlap.tsv; "
                "iteration63_cell_type_chromatin_snp_overlap.tsv"
            )
            if "iteration63_cell_type_chromatin_support_summary.tsv" not in files:
                row["evidence_files"] = f"{prefix}; {files}"
            row["claim_boundary"] = (
                "The LD block is functionally annotated, externally CRISPRi-anchored and cCRE-rich; "
                "this remains block-level prioritization rather than single-SNP causal localization, and not cell-type-specific chromatin proof."
            )
        if row.get("edge_id") == "E5":
            files = row.get("evidence_files", "")
            if "iteration63_cell_type_chromatin_support_summary.tsv" not in files:
                row["evidence_files"] = "iteration63_cell_type_chromatin_support_summary.tsv; " + files
            row["claim_boundary"] = (
                row.get("claim_boundary", "")
                + " Iteration63 adds regulatory-block chromatin context; only microglia has a direct external functional chromatin anchor, whereas neuron/endothelial remain RNA/context-linked."
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


def write_report(ccre_rows: list[dict[str, object]], snp_rows: list[dict[str, object]], summary: list[dict[str, object]]) -> None:
    block = next(row for row in summary if row.get("axis") == "regulatory_block")
    lines = [
        "# Iteration63 cell-type chromatin support",
        "",
        "## Bottom line",
        "",
        (
            "The TSPAN14 LD block now has an explicit chromatin-regulatory support layer: the UCSC/ENCODE hg38 "
            f"`encodeCcreCombined` query returned {block.get('n_ucsc_ccres_in_window')} cCREs in chr10:{WINDOW_START}-{WINDOW_END}, "
            f"including {block.get('n_enhancer_like_ccres')} enhancer-like cCREs and "
            f"{block.get('n_target_splice_overlapping_ccres')} cCREs overlapping the target splice interval."
        ),
        "",
        "This strengthens the main mechanism figure as a regulatory-block model. It does not convert the evidence into cell-type-specific ATAC/chromatin proof for neuron or endothelial cells.",
        "",
        "## Cell-axis interpretation",
        "",
        "| Axis | Status | Allowed claim | Boundary |",
        "|---|---|---|---|",
    ]
    for row in summary:
        lines.append(
            f"| {row.get('axis')} | {row.get('final_status')} | {row.get('allowed_claim')} | {row.get('claim_boundary')} |"
        )
    lines.extend(
        [
            "",
            "## Variant-to-cCRE context",
            "",
            "| SNP | Grade | cCRE overlaps | nearest cCRE | distance bp |",
            "|---|---|---:|---|---:|",
        ]
    )
    for row in snp_rows:
        lines.append(
            f"| {row.get('snp')} | {row.get('functional_priority_grade')} | {row.get('n_overlapping_ccres')} | "
            f"{row.get('nearest_ccre_id')} ({row.get('nearest_ccre_class')}) | {row.get('nearest_ccre_distance_bp')} |"
        )
    lines.extend(
        [
            "",
            "## Outputs",
            "",
            "- results/tables/iteration63_cell_type_chromatin_ccre_overlap.tsv",
            "- results/tables/iteration63_cell_type_chromatin_snp_overlap.tsv",
            "- results/tables/iteration63_cell_type_chromatin_support_summary.tsv",
            "- data/processed/external_metadata_iteration63_chromatin/ucsc_encode_ccre_tspan14_window.json",
        ]
    )
    path = REPORTS / "iteration63_cell_type_chromatin_support.md"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    payload = fetch_or_load_ucsc_ccre()
    ccres = payload.get("encodeCcreCombined", [])
    if not isinstance(ccres, list) or not ccres:
        raise RuntimeError("UCSC encodeCcreCombined query returned no cCRE rows")
    snps = load_snp_rows()
    if not snps:
        raise RuntimeError("Missing iteration62 SNP annotation rows")
    ccre_rows = build_ccre_rows(ccres, snps)
    snp_rows = build_snp_rows(ccre_rows, snps)
    summary = build_summary(ccre_rows, snp_rows)
    write_tsv(
        TABLES / "iteration63_cell_type_chromatin_ccre_overlap.tsv",
        ccre_rows,
        [
            "source",
            "source_url",
            "chrom",
            "chromStart_0based",
            "chromEnd_0based",
            "ccre_id",
            "ccre_class",
            "encode_label",
            "ucsc_label",
            "description",
            "score",
            "zScore",
            "is_enhancer_like",
            "overlaps_target_splice_interval",
            "target_splice_overlap_bp",
            "distance_to_target_splice_interval_bp",
            "contains_core_or_extended_snp",
            "contained_snps",
            "nearest_core_or_extended_snp",
            "nearest_snp_distance_bp",
            "evidence_boundary",
        ],
    )
    write_tsv(
        TABLES / "iteration63_cell_type_chromatin_snp_overlap.tsv",
        snp_rows,
        [
            "snp",
            "snp_role",
            "pos_hg38",
            "functional_priority_grade",
            "inside_target_splice_interval",
            "n_overlapping_ccres",
            "n_overlapping_enhancer_like_ccres",
            "overlapping_ccre_ids",
            "overlapping_ccre_classes",
            "nearest_ccre_id",
            "nearest_ccre_class",
            "nearest_ccre_is_enhancer_like",
            "nearest_ccre_distance_bp",
            "claim_use",
        ],
    )
    write_tsv(
        TABLES / "iteration63_cell_type_chromatin_support_summary.tsv",
        summary,
        [
            "axis",
            "n_ucsc_ccres_in_window",
            "n_enhancer_like_ccres",
            "n_target_splice_overlapping_ccres",
            "n_core_or_extended_snps_with_ccre_overlap",
            "cell_type_chromatin_specificity",
            "paired_cell_expression_context",
            "final_status",
            "allowed_claim",
            "claim_boundary",
        ],
    )
    update_evidence_model(summary)
    write_report(ccre_rows, snp_rows, summary)
    upsert_index(
        [
            {
                "type": "table",
                "name": "iteration63_cell_type_chromatin_ccre_overlap",
                "path": "results/tables/iteration63_cell_type_chromatin_ccre_overlap.tsv",
                "status": "current",
                "description": "UCSC/ENCODE cCRE overlap table for the TSPAN14 LD-block window.",
            },
            {
                "type": "table",
                "name": "iteration63_cell_type_chromatin_snp_overlap",
                "path": "results/tables/iteration63_cell_type_chromatin_snp_overlap.tsv",
                "status": "current",
                "description": "Core and extended TSPAN14 SNP overlap/nearest-distance annotation against UCSC/ENCODE cCREs.",
            },
            {
                "type": "table",
                "name": "iteration63_cell_type_chromatin_support_summary",
                "path": "results/tables/iteration63_cell_type_chromatin_support_summary.tsv",
                "status": "current",
                "description": "Cell-axis interpretation of chromatin support with explicit boundaries for microglia, neuron and endothelial evidence.",
            },
            {
                "type": "report",
                "name": "iteration63_cell_type_chromatin_support",
                "path": "results/reports/iteration63_cell_type_chromatin_support.md",
                "status": "current",
                "description": "Report documenting cCRE regulatory-block support and cell-axis chromatin claim boundaries.",
            },
        ]
    )
    print(
        f"Wrote {len(ccre_rows)} cCRE rows, {len(snp_rows)} SNP rows and {len(summary)} chromatin summary rows."
    )


if __name__ == "__main__":
    main()
