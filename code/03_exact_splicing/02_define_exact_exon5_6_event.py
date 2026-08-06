from __future__ import annotations

import csv
import json
import time
import urllib.error
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"
REPORTS = ROOT / "results" / "reports"
CACHE = ROOT / "data" / "processed" / "ensembl_iteration15"

EVENT_ID = "chr10:80509471:80512144:clu_4260_+:ENSG00000108219.15"
EVENT_CHR = "10"
EVENT_START = 80509471
EVENT_END = 80512144
TRACKED_SNPS = ["rs7080009", "rs1870138", "rs1870137", "rs7922621", "rs6586028", "rs1902660", "rs7096909", "rs1870140"]


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
    # Repair older entries whose name/type were accidentally blank in an earlier index rewrite.
    repairs = {
        "results/tables/laub2026_snp_ld_crosswalk.tsv": ("table", "Laub 2026 TSPAN14 SNP LD crosswalk", "Maps Laub 2026 functional SNPs to project TSPAN14 key SNPs using local LD."),
        "results/tables/laub2026_conclusion_conflict_audit.tsv": ("table", "Laub 2026 conclusion conflict audit", "Checks current manuscript claims against Laub 2026."),
        "results/tables/laub2026_evidence_chain_upgrade.tsv": ("table", "Laub 2026 evidence chain upgrade", "Evidence-chain changes after integrating Laub 2026."),
        "results/tables/laub2026_manuscript_rewording_targets.tsv": ("table", "Laub 2026 manuscript rewording targets", "Replacement wording for title, abstract, results, and discussion."),
        "results/reports/laub2026_tspan14_impact_audit.md": ("report", "Laub 2026 TSPAN14 impact audit Markdown", "Markdown audit and evidence-chain report."),
        "manuscript/docx/laub2026_tspan14_impact_audit_report.docx": ("report", "Laub 2026 TSPAN14 impact audit DOCX", "Word report comparing Laub 2026 with project conclusions and evidence chain."),
    }
    for path, (typ, name, desc) in repairs.items():
        if path in by_path:
            by_path[path].update({"type": typ, "name": name, "description": desc, "status": by_path[path].get("status") or "current"})
    for row in rows:
        by_path[row["path"]] = row
    write_tsv(index, list(by_path.values()), ["type", "name", "path", "status", "description"])


def get_json(url: str, cache_name: str, retries: int = 4) -> tuple[object | None, str]:
    CACHE.mkdir(parents=True, exist_ok=True)
    cache_path = CACHE / cache_name
    if cache_path.exists() and cache_path.stat().st_size > 0:
        return json.loads(cache_path.read_text(encoding="utf-8")), "cache"
    headers = {"Content-Type": "application/json", "User-Agent": "SLM-TSPAN14-bioinformatics-audit/1.0"}
    request = urllib.request.Request(url, headers=headers)
    last_error = ""
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                data = response.read().decode("utf-8")
            cache_path.write_text(data, encoding="utf-8")
            return json.loads(data), "downloaded"
        except (urllib.error.URLError, TimeoutError) as exc:
            last_error = str(exc)
            time.sleep(2 + attempt * 2)
    error_path = CACHE / f"{cache_name}.error.txt"
    error_path.write_text(last_error, encoding="utf-8")
    return None, "failed"


def get_tspan14_annotation() -> tuple[dict, str]:
    url = "https://rest.ensembl.org/lookup/symbol/homo_sapiens/TSPAN14?expand=1;content-type=application/json"
    data, status = get_json(url, "ensembl_lookup_TSPAN14_expand1.json")
    return (data or {}), status


def get_vep(rsid: str) -> tuple[list[dict], str]:
    url = f"https://rest.ensembl.org/vep/human/id/{rsid}?content-type=application/json"
    data, status = get_json(url, f"vep_{rsid}.json", retries=2)
    if isinstance(data, list):
        return data, status
    return [], status


def relation_to_event(pos: int) -> str:
    if EVENT_START < pos < EVENT_END:
        return "inside_target_splice_interval"
    if pos == EVENT_START or pos == EVENT_END:
        return "at_target_splice_boundary"
    if pos < EVENT_START:
        return "upstream_of_target_splice_interval"
    return "downstream_of_target_splice_interval"


def nearest_distance(pos: int) -> int:
    if EVENT_START <= pos <= EVENT_END:
        return 0
    return min(abs(pos - EVENT_START), abs(pos - EVENT_END))


def main() -> None:
    annotation, annotation_status = get_tspan14_annotation()
    transcripts = annotation.get("Transcript", []) if isinstance(annotation, dict) else []

    transcript_rows = []
    exon_rows = []
    junction_rows = []
    for tx in transcripts:
        tx_id = tx.get("id", "")
        exons = tx.get("Exon", []) or []
        sorted_exons = sorted(exons, key=lambda x: int(x.get("start", 0)))
        event_overlaps = []
        flanking_left = []
        flanking_right = []
        for exon in sorted_exons:
            start = int(exon.get("start", 0))
            end = int(exon.get("end", 0))
            if start <= EVENT_END and end >= EVENT_START:
                event_overlaps.append(exon.get("id", ""))
                exon_rows.append(
                    {
                        "transcript_id": tx_id,
                        "transcript_biotype": tx.get("biotype", ""),
                        "is_canonical": tx.get("is_canonical", ""),
                        "exon_id": exon.get("id", ""),
                        "exon_rank": exon.get("rank", ""),
                        "exon_start": start,
                        "exon_end": end,
                        "event_overlap": "yes",
                        "boundary_match": "event_end_matches_exon_start" if start == EVENT_END else ("event_start_matches_exon_end" if end == EVENT_START else "overlap_without_exact_boundary"),
                    }
                )
            if end <= EVENT_START:
                flanking_left.append(exon)
            if start >= EVENT_END:
                flanking_right.append(exon)
        left = flanking_left[-1] if flanking_left else {}
        right = flanking_right[0] if flanking_right else {}
        left_end = int(left.get("end", 0) or 0)
        right_start = int(right.get("start", 0) or 0)
        junction_support = "exact_downstream_acceptor_only"
        if left_end == EVENT_START and right_start == EVENT_END:
            junction_support = "exact_leafcutter_junction_boundaries"
        elif right_start == EVENT_END:
            junction_support = "exact_downstream_acceptor_with_nonmatching_upstream_boundary"
        elif left_end == EVENT_START:
            junction_support = "exact_upstream_boundary_with_nonmatching_downstream_acceptor"
        elif left or right:
            junction_support = "nearest_exon_pair_only"
        junction_rows.append(
            {
                "transcript_id": tx_id,
                "transcript_biotype": tx.get("biotype", ""),
                "is_canonical": tx.get("is_canonical", ""),
                "translation_id": (tx.get("Translation") or {}).get("id", ""),
                "transcript_start": tx.get("start", ""),
                "transcript_end": tx.get("end", ""),
                "nearest_upstream_exon_id": left.get("id", ""),
                "nearest_upstream_exon_end": left.get("end", ""),
                "nearest_downstream_exon_id": right.get("id", ""),
                "nearest_downstream_exon_start": right.get("start", ""),
                "junction_support_class": junction_support,
                "event_overlapping_exon_ids": ";".join(event_overlaps),
            }
        )
        transcript_rows.append(
            {
                "gene": "TSPAN14",
                "ensembl_gene_id": annotation.get("id", ""),
                "transcript_id": tx_id,
                "transcript_version": tx.get("version", ""),
                "transcript_biotype": tx.get("biotype", ""),
                "is_canonical": tx.get("is_canonical", ""),
                "gencode_primary": tx.get("gencode_primary", ""),
                "transcript_start": tx.get("start", ""),
                "transcript_end": tx.get("end", ""),
                "n_exons": len(sorted_exons),
                "overlaps_target_splice_interval": "yes" if event_overlaps or left or right else "no",
                "ensembl_lookup_status": annotation_status,
            }
        )

    write_tsv(
        TABLES / "iteration15_tspan14_transcript_context.tsv",
        transcript_rows,
        ["gene", "ensembl_gene_id", "transcript_id", "transcript_version", "transcript_biotype", "is_canonical", "gencode_primary", "transcript_start", "transcript_end", "n_exons", "overlaps_target_splice_interval", "ensembl_lookup_status"],
    )
    write_tsv(
        TABLES / "iteration15_tspan14_splice_event_exon_overlap.tsv",
        exon_rows,
        ["transcript_id", "transcript_biotype", "is_canonical", "exon_id", "exon_rank", "exon_start", "exon_end", "event_overlap", "boundary_match"],
    )
    write_tsv(
        TABLES / "iteration15_tspan14_splice_junction_transcript_support.tsv",
        junction_rows,
        ["transcript_id", "transcript_biotype", "is_canonical", "translation_id", "transcript_start", "transcript_end", "nearest_upstream_exon_id", "nearest_upstream_exon_end", "nearest_downstream_exon_id", "nearest_downstream_exon_start", "junction_support_class", "event_overlapping_exon_ids"],
    )

    allele_rows = {row["snp"]: row for row in read_tsv(TABLES / "iteration14_miga_snp_allele_lookup.tsv")}
    variant_rows = []
    for snp in TRACKED_SNPS:
        allele = allele_rows.get(snp, {})
        pos = int(allele.get("pos") or 0)
        vep, vep_status = get_vep(snp)
        consequences = set()
        impacted_genes = set()
        for item in vep:
            for block_name in ["transcript_consequences", "regulatory_feature_consequences", "intergenic_consequences"]:
                for cons in item.get(block_name, []) or []:
                    impacted_genes.add(cons.get("gene_symbol", "") or cons.get("gene_id", "") or cons.get("impact", ""))
                    for term in cons.get("consequence_terms", []) or []:
                        consequences.add(term)
        variant_rows.append(
            {
                "snp": snp,
                "chrom": allele.get("chrom", ""),
                "pos": pos,
                "ref_grch38": allele.get("ref_grch38", ""),
                "alt_grch38": allele.get("alt_grch38_effect_allele", ""),
                "relation_to_target_splice_interval": relation_to_event(pos),
                "nearest_event_boundary_distance_bp": nearest_distance(pos),
                "ensembl_vep_status": vep_status,
                "vep_consequence_terms": ";".join(sorted(consequences)) if consequences else "not_available_or_not_returned",
                "vep_impacted_features": ";".join(sorted(x for x in impacted_genes if x)) if impacted_genes else "not_available_or_not_returned",
                "local_fallback_consequence": "intronic_or_junction_interval_variant" if EVENT_START < pos < EVENT_END else "flanking_regulatory_LD_variant",
            }
        )
    write_tsv(
        TABLES / "iteration15_tspan14_variant_functional_annotation.tsv",
        variant_rows,
        ["snp", "chrom", "pos", "ref_grch38", "alt_grch38", "relation_to_target_splice_interval", "nearest_event_boundary_distance_bp", "ensembl_vep_status", "vep_consequence_terms", "vep_impacted_features", "local_fallback_consequence"],
    )

    methods_rows = [
        {
            "analysis_step": "Project data manifest and public-resource audit",
            "script_or_command": "scripts/01_build_data_manifest.py; scripts/03_build_processed_manifest.py; scripts/04_build_study_summary_tables.py",
            "primary_inputs": "project.yaml; local raw GWAS/QTL/single-cell/proteomic resource paths; public source metadata",
            "key_parameters": "record source name, accession/URL, sample ancestry when available, trait label, and processing status",
            "primary_outputs": "results/tables/table1_gwas_datasets.tsv; results/reports/data_processing_status.md; results/tables/final_output_index.tsv",
            "method_text_ready_note": "Use this as the data-source paragraph and Supplementary Table provenance base.",
        },
        {
            "analysis_step": "GWAS harmonization and APOE-aware candidate extraction",
            "script_or_command": "scripts/11_harmonize_gwas.py; scripts/21_extract_shared_significant_snps.py; scripts/23_build_non_apoe_candidate_regions.py; scripts/24_build_apoe_candidate_regions.py",
            "primary_inputs": "AD GWAS summary statistics; GLGC lipid GWAS summary statistics for HDL, LDL, TC, TG and nonHDL",
            "key_parameters": "rsID-based join; genome-wide significant shared SNP screen; APOE-region separated from non-APOE loci",
            "primary_outputs": "results/tables/ad_*_shared_significant_snps.tsv; results/tables/ad_lipid_non_apoe_candidate_regions.tsv; results/tables/ad_lipid_apoe_candidate_regions.tsv",
            "method_text_ready_note": "Emphasize APOE-aware design because raw overlap alone was not interpreted as shared causality.",
        },
        {
            "analysis_step": "Initial local colocalization and prior sensitivity",
            "script_or_command": "scripts/31_extract_coloc_loci.py; scripts/32_run_coloc_abf.R; scripts/34_run_coloc_susie.R; scripts/38_run_coloc_prior_sensitivity.R",
            "primary_inputs": "lead-centered AD-lipid locus windows; harmonized AD/lipid summary statistics",
            "key_parameters": "coloc.abf and coloc.susie; prior sensitivity grid including p12 values from 1e-08 to 1e-04",
            "primary_outputs": "results/tables/coloc_loci/coloc_abf_results.tsv; results/tables/coloc_susie_results.tsv; results/tables/coloc_prior_sensitivity.tsv",
            "method_text_ready_note": "Use as the statistical screen before LD-refined TSPAN14 analyses.",
        },
        {
            "analysis_step": "AD-lipid LDSC",
            "script_or_command": "scripts/12_run_iteration6_ldsc_pair_harmonized.py; scripts/13_parse_iteration6_ldsc_logs.py",
            "primary_inputs": "harmonized AD and lipid GWAS summary statistics; HapMap3/LDSC-compatible SNPs",
            "key_parameters": "pair-harmonized LDSC; AD compared with HDL, LDL, TC, TG, nonHDL",
            "primary_outputs": "results/tables/technical_iteration6_ldsc_pair_rg_results.tsv",
            "method_text_ready_note": "Report rg, SE, z and P; only AD-HDL reached nominal significance in this run.",
        },
        {
            "analysis_step": "Local LD matrix construction and allele alignment",
            "script_or_command": "scripts/17_build_local_1000g_ld_matrices.py; scripts/36_diagnose_local_ld_mismatch.py",
            "primary_inputs": "candidate SNP lists; 1000 Genomes EUR LD reference; AD/lipid allele and frequency columns",
            "key_parameters": "signed LD matrix; panel coverage and MAF concordance checks; lead-SNP LD diagnostics",
            "primary_outputs": "results/tables/technical_iteration7_ld_matrix_status.tsv; results/tables/local_ld_mismatch_diagnostics.tsv; results/tables/technical_iteration7_ld_allele_alignment.tsv",
            "method_text_ready_note": "State that all LD-aware analyses used signed EUR LD matrices and passed local coverage checks.",
        },
        {
            "analysis_step": "LD-aware coloc.susie and SuSiE fine-mapping",
            "script_or_command": "scripts/18_run_iteration7_ldaware_coloc_susie.R; scripts/22_run_iteration7_ldaware_susie_finemapping.R",
            "primary_inputs": "AD/lipid local GWAS windows; local 1000G EUR signed LD matrices",
            "key_parameters": "TSPAN14 non-APOE loci TC/LDL/nonHDL; SuSiE credible sets retained as LD-block evidence",
            "primary_outputs": "results/tables/technical_iteration7_ldaware_coloc_susie_results.tsv; results/tables/technical_iteration8_tspan14_causal_snp_variant_ranking.tsv",
            "method_text_ready_note": "Do not claim single-SNP causal resolution; report LD-block and credible-set support.",
        },
        {
            "analysis_step": "LD-aware HyPrColoc sensitivity",
            "script_or_command": "scripts/19_run_iteration7_ldaware_hyprcoloc.R; scripts/24_run_iteration8_tspan14_hyprcoloc_strengthening.R",
            "primary_inputs": "TSPAN14 AD, TC, LDL and nonHDL local summary-statistic matrices; signed LD-aligned SNPs",
            "key_parameters": "default-prior subset models; full multitrait and pair/subset posterior comparisons",
            "primary_outputs": "results/tables/technical_iteration8_tspan14_hyprcoloc_default_prior_subset.tsv; results/tables/technical_iteration8_tspan14_hyprcoloc_strengthening_results.tsv",
            "method_text_ready_note": "Report strong pair/subset support and moderate full multitrait support rather than oversimplifying into one shared variant claim.",
        },
        {
            "analysis_step": "Microglia/QTL resource inventory and TSPAN14 sQTL colocalization",
            "script_or_command": "scripts/47_build_microglia_full_qtl_resource_inventory.py; scripts/48_stream_extract_microglia_full_qtl.py; scripts/37_run_miga_splicing_full_coloc.R; scripts/25_build_iteration9_tspan14_qtl_strengthening.py",
            "primary_inputs": "MiGA/isoMiGA Zenodo resources; streamed SVZ and THA splicing-QTL extracts; AD/lipid TSPAN14 locus summary statistics",
            "key_parameters": "targeted streaming rather than full blind download; three-way min PP.H4 across AD-QTL and lipid-QTL coloc",
            "primary_outputs": "results/tables/microglia_full_qtl_resource_inventory.tsv; results/tables/miga_full_splicing_qtl_coloc_results.tsv; results/tables/technical_iteration9_tspan14_qtl_layer_summary.tsv",
            "method_text_ready_note": "The strongest molecular layer is SVZ TSPAN14 splicing-QTL colocalization; eQTL is weak/negative.",
        },
        {
            "analysis_step": "Laub functional SNP integration",
            "script_or_command": "scripts/93_laub2026_impact_audit.py; scripts/94_pure_bioinfo_strengthening_after_laub.py",
            "primary_inputs": "Laub 2026 PDF/text extraction; local LD crosswalk; TSPAN14 SuSiE ranking and QTL tables",
            "key_parameters": "tracked SNPs rs7080009, rs1870138, rs1870137, rs7922621 plus project key SNPs",
            "primary_outputs": "results/tables/iteration13_functional_snp_credible_set_audit.tsv; results/tables/iteration13_enhancer_splice_event_bridge.tsv",
            "method_text_ready_note": "Use Laub SNPs as external functional anchors, not as proof that our data alone resolved a causal SNP.",
        },
        {
            "analysis_step": "MiGA SVZ sQTL allele harmonization",
            "script_or_command": "scripts/95_harmonize_miga_sqtl_alleles.py",
            "primary_inputs": "SVZ MiGA nominal splicing extract; MiGA snp_alleles.bed.gz from Zenodo 4301005; AD/lipid allele-aligned variant table",
            "key_parameters": "MiGA slope interpreted as ALT relative to REF; slopes flipped when AD effect allele equals REF",
            "primary_outputs": "results/tables/iteration14_sqtl_allele_harmonized_directionality.tsv; results/tables/iteration14_mechanism_direction_resolution.tsv",
            "method_text_ready_note": "Report that AD/lipid effect alleles decrease the target SVZ TSPAN14 splice event for Laub core SNPs.",
        },
        {
            "analysis_step": "TSPAN14 splice event transcript annotation",
            "script_or_command": "scripts/96_annotate_splice_event_and_methods_trace.py",
            "primary_inputs": "Ensembl REST TSPAN14 transcript/exon annotation; MiGA allele lookup; target LeafCutter event coordinates",
            "key_parameters": "GRCh38 chr10:80509471-80512144; Ensembl lookup/symbol expand=1; VEP rsID annotation with local fallback",
            "primary_outputs": "results/tables/iteration15_tspan14_transcript_context.tsv; results/tables/iteration15_tspan14_variant_functional_annotation.tsv",
            "method_text_ready_note": "Use as Supplementary Methods provenance for transcript-level interpretation.",
        },
        {
            "analysis_step": "pQTL/protein-layer audit",
            "script_or_command": "scripts/25_build_iteration9_tspan14_qtl_strengthening.py; scripts/28_build_iteration12_tspan14_target_prioritization.py",
            "primary_inputs": "NIAGADS NG00102 SOMAscan1.3k analyte annotations for brain, CSF and plasma",
            "key_parameters": "TSPAN14 analyte presence/absence checked before pQTL colocalization attempt",
            "primary_outputs": "results/tables/technical_iteration9_tspan14_pqtl_panel_audit.tsv; results/tables/technical_iteration7_pqtl_resource_audit.tsv",
            "method_text_ready_note": "Protein-layer evidence is unavailable in the audited panel because TSPAN14 analyte was absent; do not report as negative pQTL association.",
        },
        {
            "analysis_step": "AD427 single-cell/single-nucleus expression validation",
            "script_or_command": "scripts/20_run_iteration7_ad427_immune_target_expression.R; scripts/26_run_iteration10_ad427_singlecell_strengthening.R; scripts/27_run_iteration11_tspan14_ad_pathology_significance_sensitivity.R",
            "primary_inputs": "AD427 immune and vasculature single-cell/single-nucleus objects; target gene list including TSPAN14",
            "key_parameters": "cell-type detection fraction and pseudobulk/pathology-associated tests; FDR-controlled interpretation",
            "primary_outputs": "results/tables/technical_iteration10_ad427_target_expression_by_celltype.tsv; results/tables/technical_iteration11_tspan14_ad_pathology_significance_decision.tsv",
            "method_text_ready_note": "Use as expression plausibility, not genotype-aware single-cell QTL evidence; AD-pathology differential expression is not FDR-significant.",
        },
        {
            "analysis_step": "Target prioritization and risk register",
            "script_or_command": "scripts/28_build_iteration12_tspan14_target_prioritization.py",
            "primary_inputs": "genetic, QTL, single-cell, pQTL availability, druggability and mechanism-risk evidence layers",
            "key_parameters": "separate discovery priority from translational maturity; explicit ADAM10/Notch risk boundary",
            "primary_outputs": "results/tables/technical_iteration12_tspan14_target_priority_summary.tsv; results/tables/technical_iteration12_tspan14_target_risk_register.tsv",
            "method_text_ready_note": "After Laub integration, target confidence is stronger, but drug-development readiness remains limited by protein and functional-direction gaps.",
        },
        {
            "analysis_step": "Run-level output validation",
            "script_or_command": "scripts/90_validate_final_outputs.py",
            "primary_inputs": "results/tables/final_output_index.tsv and indexed local outputs",
            "key_parameters": "non-empty file checks plus selected biological/statistical sanity checks",
            "primary_outputs": "terminal validation status; current indexed outputs count",
            "method_text_ready_note": "Use this only as reproducibility QA, not as a biological result.",
        },
    ]
    write_tsv(
        TABLES / "iteration15_methods_provenance_trace.tsv",
        methods_rows,
        ["analysis_step", "script_or_command", "primary_inputs", "key_parameters", "primary_outputs", "method_text_ready_note"],
    )

    exact_junction_count = sum(1 for row in junction_rows if row["junction_support_class"] == "exact_leafcutter_junction_boundaries")
    acceptor_count = sum(
        1
        for row in junction_rows
        if row["nearest_downstream_exon_start"] == str(EVENT_END)
        or row["nearest_downstream_exon_start"] == EVENT_END
    )
    inside_variants = sum(1 for row in variant_rows if row["relation_to_target_splice_interval"] == "inside_target_splice_interval")
    report = f"""# Iteration 15 TSPAN14 splice-event annotation and methods provenance

## Scope

This iteration adds transcript-level annotation and explicit Methods provenance for the pure-bioinformatics workflow. It does not generate figures or manuscript text.

## Main result

- Target event: `{EVENT_ID}`.
- Ensembl transcript lookup status: `{annotation_status}`.
- Transcripts with exact LeafCutter boundary support: {exact_junction_count}.
- Transcripts with exact downstream acceptor support at `{EVENT_END}`: {acceptor_count}.
- Tracked variants inside the target splice interval: {inside_variants}/{len(TRACKED_SNPS)}.

## Interpretation

The Laub core variants remain positioned inside the target TSPAN14 splice interval, supporting a direct spatial bridge between the external CRISPRi enhancer evidence and our strongest MiGA SVZ splicing-QTL event. Transcript annotation provides biological context but does not by itself prove isoform-level protein consequences.

## Methods trace

`results/tables/iteration15_methods_provenance_trace.tsv` records scripts, inputs, parameters, outputs and manuscript-ready wording boundaries for the main pure-bioinformatics analyses. This table should be used when drafting the final Methods and Supplementary Methods sections.

## New outputs

- `results/tables/iteration15_tspan14_transcript_context.tsv`
- `results/tables/iteration15_tspan14_splice_event_exon_overlap.tsv`
- `results/tables/iteration15_tspan14_splice_junction_transcript_support.tsv`
- `results/tables/iteration15_tspan14_variant_functional_annotation.tsv`
- `results/tables/iteration15_methods_provenance_trace.tsv`
"""
    (REPORTS / "iteration15_tspan14_splice_event_annotation_and_methods_trace.md").write_text(report, encoding="utf-8")

    upsert_index(
        [
            {
                "type": "table",
                "name": "Iteration15 TSPAN14 transcript context",
                "path": "results/tables/iteration15_tspan14_transcript_context.tsv",
                "status": "current",
                "description": "Ensembl transcript context for TSPAN14 around the target SVZ splicing event.",
            },
            {
                "type": "table",
                "name": "Iteration15 TSPAN14 splice event exon overlap",
                "path": "results/tables/iteration15_tspan14_splice_event_exon_overlap.tsv",
                "status": "current",
                "description": "Exon overlaps and boundary matches for the target TSPAN14 splicing event.",
            },
            {
                "type": "table",
                "name": "Iteration15 TSPAN14 splice junction transcript support",
                "path": "results/tables/iteration15_tspan14_splice_junction_transcript_support.tsv",
                "status": "current",
                "description": "Transcript-level nearest exon and junction-boundary support for the target TSPAN14 event.",
            },
            {
                "type": "table",
                "name": "Iteration15 TSPAN14 variant functional annotation",
                "path": "results/tables/iteration15_tspan14_variant_functional_annotation.tsv",
                "status": "current",
                "description": "Ensembl VEP/local fallback functional annotation for Laub and project key TSPAN14 variants.",
            },
            {
                "type": "table",
                "name": "Iteration15 methods provenance trace",
                "path": "results/tables/iteration15_methods_provenance_trace.tsv",
                "status": "current",
                "description": "Methods-ready provenance table linking scripts, inputs, parameters, outputs and claim boundaries.",
            },
            {
                "type": "report",
                "name": "Iteration15 splice-event annotation and methods trace report",
                "path": "results/reports/iteration15_tspan14_splice_event_annotation_and_methods_trace.md",
                "status": "current",
                "description": "Report summarizing transcript-level annotation and methods provenance outputs.",
            },
        ]
    )


if __name__ == "__main__":
    main()
