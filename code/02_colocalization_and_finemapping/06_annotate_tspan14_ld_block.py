from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"
REPORTS = ROOT / "results" / "reports"

CORE_SNPS = ["rs7080009", "rs1870138", "rs1870137", "rs1902660", "rs6586028", "rs1870140"]
EXTENDED_SNPS = ["rs7922621", "rs7096909"]


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
        return float(value)
    except (TypeError, ValueError):
        return default


def by_key(rows: list[dict[str, str]], key: str) -> dict[str, dict[str, str]]:
    return {row.get(key, ""): row for row in rows if row.get(key, "")}


def load_sources() -> dict[str, object]:
    return {
        "vep": by_key(read_tsv(TABLES / "iteration15_tspan14_variant_functional_annotation.tsv"), "snp"),
        "credible": by_key(read_tsv(TABLES / "iteration13_functional_snp_credible_set_audit.tsv"), "snp"),
        "direction": by_key(read_tsv(TABLES / "iteration44_main_figure_direction_axis.tsv"), "snp"),
        "chain": by_key(read_tsv(TABLES / "iteration35_tspan14_direction_chain.tsv"), "snp"),
        "bridge": by_key(read_tsv(TABLES / "iteration13_enhancer_splice_event_bridge.tsv"), "snp"),
        "niagads": read_tsv(TABLES / "iteration45_niagads_precise_external_sqtl_replication.tsv"),
        "finemap": read_tsv(TABLES / "technical_iteration8_tspan14_causal_snp_variant_ranking.tsv"),
        "laub_ld": read_tsv(TABLES / "laub2026_snp_ld_crosswalk.tsv"),
    }


def niagads_summary(rows: list[dict[str, str]]) -> dict[str, dict[str, object]]:
    out: dict[str, dict[str, object]] = {}
    for row in rows:
        snp = row.get("snp", "")
        if not snp:
            continue
        item = out.setdefault(
            snp,
            {
                "n_niagads_sqtl_regions": 0,
                "niagads_regions": set(),
                "best_niagads_sqtl_fdr": 1.0,
                "n_positive_beta_rows": 0,
                "directionally_upgraded": False,
            },
        )
        item["niagads_regions"].add(row.get("context", ""))
        item["best_niagads_sqtl_fdr"] = min(as_float(item["best_niagads_sqtl_fdr"], 1.0), as_float(row.get("fdr"), 1.0))
        if as_float(row.get("portal_beta")) > 0:
            item["n_positive_beta_rows"] = int(item["n_positive_beta_rows"]) + 1
        if row.get("direction_call") == "risk_allele_matches_variant_alt_and_positive_beta":
            item["directionally_upgraded"] = True
    for item in out.values():
        item["n_niagads_sqtl_regions"] = len(item["niagads_regions"])
        item["niagads_regions"] = ";".join(sorted(str(x) for x in item["niagads_regions"] if x))
    return out


def max_laub_ld(rows: list[dict[str, str]]) -> dict[str, float]:
    out: dict[str, float] = {}
    for row in rows:
        laub = row.get("laub_or_published_snp", "")
        key = row.get("our_key_snp", "")
        r2 = as_float(row.get("r2"))
        if laub:
            out[laub] = max(out.get(laub, 0.0), r2)
        if key:
            out[key] = max(out.get(key, 0.0), r2)
    return out


def score_variant(row: dict[str, object]) -> tuple[int, str, str]:
    score = 0
    reasons: list[str] = []
    role = str(row.get("snp_role", ""))
    consequence = str(row.get("vep_consequence_terms", ""))
    relation = str(row.get("relation_to_target_splice_interval", ""))
    if "Laub_CRISPRi_anchor" in role or "laub_core_functional" in role:
        score += 4
        reasons.append("direct_Laub_CRISPRi_anchor")
    elif "project_AD_lipid_sQTL_proxy" in role or "our_key_coloc_or_sqtl" in role:
        score += 2
        reasons.append("project_AD_lipid_sQTL_proxy")
    if "regulatory_region_variant" in consequence:
        score += 2
        reasons.append("VEP_regulatory_region_variant")
    if "splice_polypyrimidine_tract_variant" in consequence:
        score += 2
        reasons.append("splice_polypyrimidine_tract_annotation")
    if relation == "inside_target_splice_interval":
        score += 2
        reasons.append("inside_target_splice_interval")
    elif relation in {"upstream_of_target_splice_interval", "downstream_of_target_splice_interval"}:
        score += 1
        reasons.append("near_target_splice_interval")
    if as_float(row.get("n_cs_memberships")) >= 4:
        score += 2
        reasons.append("multi_trait_credible_set_member")
    elif as_float(row.get("n_cs_memberships")) >= 2:
        score += 1
        reasons.append("partial_credible_set_member")
    if as_float(row.get("max_bridge_r2")) >= 0.98:
        score += 2
        reasons.append("high_LD_with_functional_block")
    elif as_float(row.get("max_bridge_r2")) >= 0.75:
        score += 1
        reasons.append("moderate_LD_with_functional_block")
    if as_float(row.get("n_niagads_sqtl_regions")) >= 3:
        score += 2
        reasons.append("NIAGADS_three_region_target_equivalent_sQTL")
    if str(row.get("direction_axis_call", "")).startswith("AD_risk_allele_increases"):
        score += 2
        reasons.append("AD_lipid_sQTL_direction_axis")
    if score >= 14:
        grade = "highest_priority_functional_anchor"
    elif score >= 11:
        grade = "high_priority_regulatory_proxy"
    elif score >= 8:
        grade = "moderate_priority_supportive_variant"
    else:
        grade = "context_variant"
    return score, grade, ";".join(reasons)


def build_core_annotation(sources: dict[str, object]) -> list[dict[str, object]]:
    vep: dict[str, dict[str, str]] = sources["vep"]  # type: ignore[assignment]
    credible: dict[str, dict[str, str]] = sources["credible"]  # type: ignore[assignment]
    direction: dict[str, dict[str, str]] = sources["direction"]  # type: ignore[assignment]
    chain: dict[str, dict[str, str]] = sources["chain"]  # type: ignore[assignment]
    bridge: dict[str, dict[str, str]] = sources["bridge"]  # type: ignore[assignment]
    niagads = niagads_summary(sources["niagads"])  # type: ignore[arg-type]
    ldmax = max_laub_ld(sources["laub_ld"])  # type: ignore[arg-type]
    rows: list[dict[str, object]] = []
    for snp in CORE_SNPS + EXTENDED_SNPS:
        v = vep.get(snp, {})
        c = credible.get(snp, {})
        d = direction.get(snp, {})
        ch = chain.get(snp, {})
        b = bridge.get(snp, {})
        n = niagads.get(snp, {})
        role = d.get("snp_axis_role") or c.get("snp_class") or ch.get("snp_role")
        row: dict[str, object] = {
            "snp": snp,
            "snp_role": role,
            "chrom": v.get("chrom", ""),
            "pos_hg38": v.get("pos", b.get("position", "")),
            "relation_to_target_splice_interval": v.get("relation_to_target_splice_interval", ""),
            "nearest_event_boundary_distance_bp": v.get("nearest_event_boundary_distance_bp", ""),
            "vep_consequence_terms": v.get("vep_consequence_terms", ""),
            "vep_impacted_features": v.get("vep_impacted_features", ""),
            "local_fallback_consequence": v.get("local_fallback_consequence", ""),
            "mean_pip": c.get("mean_pip", ""),
            "min_pip": c.get("min_pip", ""),
            "n_cs_memberships": c.get("n_cs_memberships", ""),
            "max_bridge_r2": c.get("max_bridge_r2") or ldmax.get(snp, ""),
            "credible_set_interpretation": c.get("interpretation", ""),
            "laub_functional_anchor": ch.get("laub_functional_anchor", ""),
            "direction_axis_call": d.get("main_figure_axis_call", ""),
            "main_figure_grade": d.get("main_figure_grade", ""),
            "n_niagads_sqtl_regions": n.get("n_niagads_sqtl_regions", 0),
            "niagads_regions": n.get("niagads_regions", ""),
            "best_niagads_sqtl_fdr": n.get("best_niagads_sqtl_fdr", ""),
            "n_positive_beta_rows": n.get("n_positive_beta_rows", 0),
            "directionally_upgraded_external_sqtl": n.get("directionally_upgraded", False),
            "manuscript_use_boundary": "Use as LD-block functional annotation and prioritization, not as single-SNP causal proof.",
        }
        score, grade, reasons = score_variant(row)
        row["functional_prior_score"] = score
        row["functional_priority_grade"] = grade
        row["functional_score_components"] = reasons
        rows.append(row)
    return sorted(rows, key=lambda item: (-as_float(item.get("functional_prior_score")), str(item.get("snp"))))


def build_finemap_top_table(sources: dict[str, object], core_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    annotated = {row["snp"]: row for row in core_rows}
    finemap: list[dict[str, str]] = sources["finemap"]  # type: ignore[assignment]
    top = sorted(finemap, key=lambda row: as_float(row.get("product_pip")), reverse=True)[:20]
    out: list[dict[str, object]] = []
    for rank, row in enumerate(top, start=1):
        snp = row.get("snp", "")
        ann = annotated.get(snp, {})
        out.append(
            {
                "rank_by_product_pip": rank,
                "snp": snp,
                "AD_PIP": row.get("AD", ""),
                "LDL_PIP": row.get("LDL", ""),
                "TC_PIP": row.get("TC", ""),
                "nonHDL_PIP": row.get("nonHDL", ""),
                "mean_pip": row.get("mean_pip", ""),
                "min_pip": row.get("min_pip", ""),
                "product_pip": row.get("product_pip", ""),
                "n_cs_memberships": row.get("n_cs_memberships", ""),
                "functional_priority_grade": ann.get("functional_priority_grade", "not_in_core_annotation_table"),
                "functional_prior_score": ann.get("functional_prior_score", ""),
                "functional_annotation_note": ann.get("functional_score_components", "credible_set_statistical_candidate_without_current_functional_annotation"),
                "claim_boundary": "Top fine-map ranking remains LD-block prioritization; no single causal SNP is resolved.",
            }
        )
    return out


def build_summary(core_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    highest = [row for row in core_rows if row.get("functional_priority_grade") == "highest_priority_functional_anchor"]
    high_or_better = [
        row for row in core_rows
        if row.get("functional_priority_grade") in {"highest_priority_functional_anchor", "high_priority_regulatory_proxy"}
    ]
    niagads_complete = [row for row in core_rows if as_float(row.get("n_niagads_sqtl_regions")) >= 3]
    inside_target = [row for row in core_rows if row.get("relation_to_target_splice_interval") == "inside_target_splice_interval"]
    return [
        {
            "summary_item": "functional_anchor_core",
            "n_variants": len(core_rows),
            "n_highest_priority_functional_anchors": len(highest),
            "n_high_or_better_variants": len(high_or_better),
            "n_inside_target_splice_interval": len(inside_target),
            "n_complete_niagads_three_region_sqtl": len(niagads_complete),
            "top_variants": ";".join(str(row.get("snp")) for row in high_or_better[:6]),
            "final_status": "ld_block_functionally_anchored",
            "allowed_claim": "The TSPAN14 LD block is functionally anchored by Laub CRISPRi variants, regulatory annotations, multi-trait credible-set membership and NIAGADS target-equivalent sQTL direction.",
            "claim_boundary": "Functional annotation strengthens LD-block prioritization; it does not resolve a single causal SNP or prove protein/downstream mechanism.",
        },
        {
            "summary_item": "single_snp_boundary",
            "n_variants": len(core_rows),
            "n_highest_priority_functional_anchors": len(highest),
            "n_high_or_better_variants": len(high_or_better),
            "n_inside_target_splice_interval": len(inside_target),
            "n_complete_niagads_three_region_sqtl": len(niagads_complete),
            "top_variants": ";".join(str(row.get("snp")) for row in highest),
            "final_status": "single_snp_not_resolved",
            "allowed_claim": "Several high-LD variants jointly define the functional block, with rs7080009/rs1870137/rs1870138 as strongest external CRISPRi anchors.",
            "claim_boundary": "Do not claim one SNP is causal; use block-level regulatory mechanism language.",
        },
    ]


def update_evidence_model(summary: list[dict[str, object]]) -> None:
    edge_path = TABLES / "iteration39_mechanism_figure_edges.tsv"
    edges = read_tsv(edge_path)
    anchor = summary[0]
    for row in edges:
        if row.get("edge_id") == "E2":
            row["figure_label"] = (
                "Laub core SNPs rs7080009/rs1870138/rs1870137; "
                f"functional anchors={anchor.get('n_highest_priority_functional_anchors')}; "
                f"high-priority block variants={anchor.get('n_high_or_better_variants')}"
            )
            row["evidence_files"] = (
                "iteration62_ld_block_functional_annotation_summary.tsv; "
                "iteration62_ld_block_functional_annotation.tsv; "
                "iteration62_credible_set_top_functional_candidates.tsv; "
                + row.get("evidence_files", "")
            )
            row["claim_boundary"] = (
                "The LD block is now functionally annotated and externally CRISPRi-anchored, but remains block-level prioritization rather than single-SNP causal localization."
            )
    write_tsv(
        edge_path,
        edges,
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


def write_report(core_rows: list[dict[str, object]], top_rows: list[dict[str, object]], summary: list[dict[str, object]]) -> None:
    lines = [
        "# Iteration62 TSPAN14 LD-block functional annotation",
        "",
        "## Bottom line",
        "",
        "The TSPAN14 LD block is now functionally annotated as a regulatory block rather than only a statistical colocalization interval. The strongest anchors are the Laub CRISPRi SNPs rs7080009, rs1870138 and rs1870137, which sit inside the target splice interval, carry regulatory-region VEP annotations, belong to multi-trait credible sets and are in high LD with project proxy variants.",
        "",
        "This strengthens block-level mechanism confidence but preserves the causal boundary: no single SNP is resolved as causal.",
        "",
        "## Core variant prioritization",
        "",
        "| SNP | Role | Score | Grade | Key components |",
        "|---|---|---:|---|---|",
    ]
    for row in core_rows:
        lines.append(
            f"| {row.get('snp')} | {row.get('snp_role')} | {row.get('functional_prior_score')} | "
            f"{row.get('functional_priority_grade')} | {row.get('functional_score_components')} |"
        )
    lines.extend(
        [
            "",
            "## Summary",
            "",
        ]
    )
    for row in summary:
        lines.append(f"- {row.get('summary_item')}: {row.get('allowed_claim')} Boundary: {row.get('claim_boundary')}")
    lines.extend(
        [
            "",
            "## Outputs",
            "",
            "- results/tables/iteration62_ld_block_functional_annotation.tsv",
            "- results/tables/iteration62_credible_set_top_functional_candidates.tsv",
            "- results/tables/iteration62_ld_block_functional_annotation_summary.tsv",
        ]
    )
    report = REPORTS / "iteration62_ld_block_functional_annotation.md"
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    sources = load_sources()
    core_rows = build_core_annotation(sources)
    top_rows = build_finemap_top_table(sources, core_rows)
    summary = build_summary(core_rows)

    core_fields = [
        "snp",
        "snp_role",
        "chrom",
        "pos_hg38",
        "relation_to_target_splice_interval",
        "nearest_event_boundary_distance_bp",
        "vep_consequence_terms",
        "vep_impacted_features",
        "local_fallback_consequence",
        "mean_pip",
        "min_pip",
        "n_cs_memberships",
        "max_bridge_r2",
        "credible_set_interpretation",
        "laub_functional_anchor",
        "direction_axis_call",
        "main_figure_grade",
        "n_niagads_sqtl_regions",
        "niagads_regions",
        "best_niagads_sqtl_fdr",
        "n_positive_beta_rows",
        "directionally_upgraded_external_sqtl",
        "functional_prior_score",
        "functional_priority_grade",
        "functional_score_components",
        "manuscript_use_boundary",
    ]
    top_fields = [
        "rank_by_product_pip",
        "snp",
        "AD_PIP",
        "LDL_PIP",
        "TC_PIP",
        "nonHDL_PIP",
        "mean_pip",
        "min_pip",
        "product_pip",
        "n_cs_memberships",
        "functional_priority_grade",
        "functional_prior_score",
        "functional_annotation_note",
        "claim_boundary",
    ]
    summary_fields = [
        "summary_item",
        "n_variants",
        "n_highest_priority_functional_anchors",
        "n_high_or_better_variants",
        "n_inside_target_splice_interval",
        "n_complete_niagads_three_region_sqtl",
        "top_variants",
        "final_status",
        "allowed_claim",
        "claim_boundary",
    ]
    write_tsv(TABLES / "iteration62_ld_block_functional_annotation.tsv", core_rows, core_fields)
    write_tsv(TABLES / "iteration62_credible_set_top_functional_candidates.tsv", top_rows, top_fields)
    write_tsv(TABLES / "iteration62_ld_block_functional_annotation_summary.tsv", summary, summary_fields)
    update_evidence_model(summary)
    write_report(core_rows, top_rows, summary)
    upsert_index(
        [
            {
                "type": "table",
                "name": "iteration62_ld_block_functional_annotation",
                "path": "results/tables/iteration62_ld_block_functional_annotation.tsv",
                "status": "current",
                "description": "Functional-prior annotation for TSPAN14 LD-block core and extended SNPs integrating Laub CRISPRi, VEP, credible sets, LD and NIAGADS sQTL support.",
            },
            {
                "type": "table",
                "name": "iteration62_credible_set_top_functional_candidates",
                "path": "results/tables/iteration62_credible_set_top_functional_candidates.tsv",
                "status": "current",
                "description": "Top TSPAN14 credible-set variants by product PIP with functional annotation status.",
            },
            {
                "type": "table",
                "name": "iteration62_ld_block_functional_annotation_summary",
                "path": "results/tables/iteration62_ld_block_functional_annotation_summary.tsv",
                "status": "current",
                "description": "Decision summary for TSPAN14 LD-block functional anchoring and single-SNP boundary.",
            },
            {
                "type": "report",
                "name": "iteration62_ld_block_functional_annotation",
                "path": "results/reports/iteration62_ld_block_functional_annotation.md",
                "status": "current",
                "description": "Report documenting LD-block functional annotation and functional-prior scoring.",
            },
        ]
    )
    print(f"Wrote {len(core_rows)} core annotation rows, {len(top_rows)} top credible-set rows and {len(summary)} summary rows.")


if __name__ == "__main__":
    main()
