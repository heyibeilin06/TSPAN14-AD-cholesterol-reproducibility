"""Retrieve and classify GTEx v8 TSPAN14 brain sQTLs for the exact project junction."""

from __future__ import annotations

import re
import os
from pathlib import Path

import pandas as pd
import requests


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs" / "p1_cell_context"
DATA = Path(os.environ.get("SLM_DATA_ROOT", ROOT / "data" / "raw"))
LD_FILE = Path(os.environ.get(
    "TSPAN14_LD_FILE",
    DATA / "p0_ld_1000g_eur" / "matrices" / "TC_nonAPOE_004.ld.tsv",
))
API = "https://gtexportal.org/api/v2/association/singleTissueSqtl"
GENCODE_ID = "ENSG00000108219.14"
EXACT_INTERVAL = (80509471, 80512144)
ADJACENT_INTERVAL = (80512269, 80514018)
P_THRESHOLD = 5e-8
R2_THRESHOLD = 0.1


def event_relation(phenotype_id: str) -> str:
    match = re.search(r"(?:chr)?10:(\d+):(\d+):", str(phenotype_id))
    if not match:
        return "unmapped_splice_feature"
    interval = (int(match.group(1)), int(match.group(2)))
    if interval == EXACT_INTERVAL:
        return "exact_project_exon5_exon6"
    if interval == ADJACENT_INTERVAL:
        return "adjacent_exon6_exon7_not_exact"
    return "other_tspan14_splice_feature"


def is_brain_tissue(frame: pd.DataFrame) -> pd.Series:
    return frame["tissueSiteDetailId"].astype(str).str.startswith("Brain_")


def ld_clump(frame: pd.DataFrame, ld: pd.DataFrame, *, r2_threshold: float) -> list[str]:
    candidates = frame[frame["snpId"].isin(ld.index)].sort_values("pValue")
    selected: list[str] = []
    for snp in candidates["snpId"].drop_duplicates():
        if all(float(ld.loc[snp, prior]) ** 2 < r2_threshold for prior in selected):
            selected.append(snp)
    return selected


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    response = requests.get(
        API,
        params={"gencodeId": GENCODE_ID, "datasetId": "gtex_v8", "itemsPerPage": 100000, "page": 0},
        timeout=180,
    )
    response.raise_for_status()
    payload = response.json()
    frame = pd.DataFrame(payload["data"])
    if frame.empty or payload["paging_info"]["totalNumberOfItems"] != len(frame):
        raise RuntimeError("GTEx API did not return the complete requested TSPAN14 sQTL result set")
    frame["event_relation"] = frame["phenotypeId"].map(event_relation)
    exact_brain = frame[is_brain_tissue(frame) & (frame["event_relation"] == "exact_project_exon5_exon6")].copy()
    exact_brain.to_csv(OUT / "p1_gtex_v8_tspan14_exact_brain_sqtl.tsv", sep="\t", index=False)

    summary = (
        exact_brain.groupby("tissueSiteDetailId", dropna=False)
        .agg(
            n_significant_pairs=("snpId", "size"),
            n_unique_variants=("snpId", "nunique"),
            min_p=("pValue", "min"),
            n_genome_wide_pairs=("pValue", lambda values: int((values <= P_THRESHOLD).sum())),
        )
        .reset_index()
        .sort_values("min_p")
    )
    summary["source"] = "GTEx Portal API v2 singleTissueSqtl; gtex_v8; ENSG00000108219.14"
    summary.to_csv(OUT / "p1_gtex_v8_tspan14_exact_brain_summary.tsv", sep="\t", index=False)

    top_tissue = summary.iloc[0]["tissueSiteDetailId"]
    top = exact_brain[(exact_brain["tissueSiteDetailId"] == top_tissue) & (exact_brain["pValue"] <= P_THRESHOLD)].copy()
    ld = pd.read_csv(LD_FILE, sep="\t", index_col=0)
    selected = ld_clump(top, ld, r2_threshold=R2_THRESHOLD)
    selected_rows = top[top["snpId"].isin(selected)].copy()
    selected_rows["selection_reason"] = "LD-independent lead under r2<0.1 for the exact GTEx brain exon5-exon6 event"
    selected_rows.to_csv(OUT / "p1_gtex_exact_brain_sqtl_instrument_candidates.tsv", sep="\t", index=False)
    gate = {
        "exposure_feature": "TSPAN14 exact exon5-exon6 junction",
        "event_relation_to_project": "exact_project_exon5_exon6",
        "external_resource": "GTEx Portal v8 single-tissue sQTL",
        "top_brain_tissue": top_tissue,
        "exact_external_brain_splice_replication": True,
        "n_brain_tissues_with_exact_significant_pairs": int(summary.shape[0]),
        "n_ld_independent_instruments_top_tissue": len(selected),
        "selected_instruments_top_tissue": ";".join(selected),
        "full_feature_summary_available": False,
        "feature_matched_AD_colocalization_passed": False,
        "eligible_for_cis_MR": False,
        "eligible_for_SMR_HEIDI": False,
        "eligible_for_Steiger": False,
        "eligible_for_coloc_constrained_mediation_sensitivity": False,
        "decision": "not_initiated",
        "allowed_interpretation": "GTEx provides cross-tissue consistency for the exact TSPAN14 exon5-exon6 splice feature across partially overlapping donor sets.",
        "prohibited_interpretation": "The significant-pair API response and one independent top-tissue instrument do not identify lipid-to-AD mediation or allow HEIDI diagnostics.",
        "reason": "Only one LD-independent strong instrument was observed in the best-supported brain tissue, and the endpoint is a significant-association resource rather than all tested feature statistics for AD colocalization.",
    }
    pd.DataFrame([gate]).to_csv(OUT / "p1_gtex_exact_brain_sqtl_causal_gate.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
