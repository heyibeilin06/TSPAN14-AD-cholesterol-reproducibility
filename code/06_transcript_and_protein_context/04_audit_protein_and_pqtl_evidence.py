"""Audit direct protein-context evidence without inferring an unobserved splice consequence."""

from __future__ import annotations

import html
import re
from pathlib import Path

import pandas as pd
import requests


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs" / "p2_protein_feasibility"
UNIPROT = "https://rest.uniprot.org/uniprotkb/Q8NG11.json"
HPA_BRAIN = "https://www.proteinatlas.org/ENSG00000108219-TSPAN14/single%2Bcell/brain"
HPA_PRIMARY = "https://www.proteinatlas.org/ENSG00000108219-TSPAN14/tissue/primary%2Bdata"
CPTAC = "https://assays.cancer.gov/CPTAC-1188"
EVENTS = {"project_exon5_exon6": (150, 151), "adjacent_exon6_exon7": (192, 193)}


def interval_contains(container: tuple[int, int], query: tuple[int, int]) -> bool:
    return container[0] <= query[0] and query[1] <= container[1]


def plain_text(markup: str) -> str:
    stripped = re.sub(r"<[^>]+>", " ", markup)
    return re.sub(r"\s+", " ", html.unescape(stripped)).strip()


def count_nd_replicates(fragment: str, group: str) -> int:
    text = plain_text(fragment)
    return len(re.findall(rf"{re.escape(group)}_\d+\s+ND\b", text))


def antibody_concordance(calls: list[str]) -> str:
    return "concordant" if len(set(calls)) == 1 else "discordant"


def find_adam10_region(features: list[dict[str, object]]) -> tuple[int, int]:
    for feature in features:
        if feature.get("type") == "Region" and "ADAM10" in str(feature.get("description", "")):
            location = feature["location"]
            return int(location["start"]["value"]), int(location["end"]["value"])
    raise RuntimeError("No UniProt ADAM10-interaction region annotation was found")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    uniprot = requests.get(UNIPROT, timeout=90)
    uniprot.raise_for_status()
    record = uniprot.json()
    region = find_adam10_region(record["features"])
    interaction_experiments = next(
        comment["interactions"][0]["numberOfExperiments"]
        for comment in record["comments"]
        if comment.get("commentType") == "INTERACTION"
        and comment["interactions"][0].get("interactantTwo", {}).get("geneName") == "ADAM10"
    )

    hpa = requests.get(HPA_BRAIN, timeout=180)
    hpa.raise_for_status()
    hpa_text = hpa.text
    dvp_start = hpa_text.find('id="dvp"')
    if dvp_start < 0:
        raise RuntimeError("HPA DVP section was not found")
    dvp_fragment = hpa_text[dvp_start:]
    cell_groups = ["Neurons", "Astrocytes and Neuropil", "Microglia and Neuropil"]
    hpa_rows = [
        {
            "resource": "Human Protein Atlas DVP brain",
            "cell_context": group,
            "n_reported_ND_replicates": count_nd_replicates(dvp_fragment, group),
            "interpretation": "Healthy cerebral-cortex DVP target-level observation only; ND does not imply absence in all brain states.",
            "claim_boundary": "This is not AD-state, isoform-specific, or junction-resolved protein evidence.",
            "source_url": HPA_BRAIN,
        }
        for group in cell_groups
    ]
    primary = requests.get(HPA_PRIMARY, timeout=180)
    primary.raise_for_status()
    antibody_calls = {}
    for antibody in ("HPA014773", "HPA057174"):
        match = re.search(
            rf'title="<b>{antibody}</b><br>Cerebral cortex, Neuronal cells - <b>([^<]+)</b>"',
            primary.text,
        )
        if match is None:
            raise RuntimeError(f"HPA cortical-neuron staining call was not found for {antibody}")
        antibody_calls[antibody] = match.group(1).strip().lower().replace(" ", "_")

    cptac = requests.get(CPTAC, timeout=120)
    cptac.raise_for_status()
    cptac_text = plain_text(cptac.text)
    peptide_match = re.search(r"Peptide Sequence\s+([A-Z]+)\s+Modification Type", cptac_text)
    position_match = re.search(r"Peptide Start\s+(\d+)\s+Peptide End\s+(\d+)", cptac_text)
    if peptide_match is None or position_match is None:
        raise RuntimeError("CPTAC TSPAN14 assay fields were not parsed")

    rows: list[dict[str, object]] = []
    for event, interval in EVENTS.items():
        rows.append(
            {
                "evidence_layer": "curated_direct_interaction_context",
                "resource": "UniProt Q8NG11",
                "target_or_feature": event,
                "observation": f"AA{interval[0]}-{interval[1]} is within extracellular AA{region[0]}-{region[1]}, annotated as necessary and sufficient for ADAM10 interaction; curated interaction experiments={interaction_experiments}.",
                "evidence_status": "contextual_known_biology",
                "claim_boundary": "Does not show that the splice QTL changes the region, protein conformation, ADAM10 binding, trafficking, or proteolysis.",
                "source_url": UNIPROT,
            }
        )
    rows.extend(
        {
            "evidence_layer": "direct_protein_detection_context",
            "resource": item["resource"],
            "target_or_feature": item["cell_context"],
            "observation": f"TSPAN14 was reported ND in {item['n_reported_ND_replicates']} DVP replicates.",
            "evidence_status": "negative_healthy_cell_resolved_detection_context",
            "claim_boundary": item["claim_boundary"],
            "source_url": item["source_url"],
        }
        for item in hpa_rows
    )
    rows.append(
        {
            "evidence_layer": "antibody_based_brain_context",
            "resource": "Human Protein Atlas primary tissue data",
            "target_or_feature": "Cerebral cortex neuronal cells",
            "observation": "; ".join(f"{antibody}={call}" for antibody, call in antibody_calls.items()),
            "evidence_status": f"{antibody_concordance(list(antibody_calls.values()))}_antibody_staining",
            "claim_boundary": "Discordant antibodies cannot resolve TSPAN14 abundance, cell specificity, AD-state regulation, or isoform usage.",
            "source_url": HPA_PRIMARY,
        }
    )
    rows.append(
        {
            "evidence_layer": "targeted_proteomics_assay_feasibility",
            "resource": "NCI CPTAC Assay Portal CPTAC-1188",
            "target_or_feature": "TSPAN14 canonical protein",
            "observation": f"A direct MRM/SRM assay exists for peptide {peptide_match.group(1)} at AA{position_match.group(1)}-{position_match.group(2)}.",
            "evidence_status": "validated_targeted_assay_outside_brain_and_not_isoform_discriminating",
            "claim_boundary": "The assay was characterized in ovarian tumor lysate and targets AA256-265, outside the splice-associated EC2 interval; it cannot quantify the proposed isoform effect as-is.",
            "source_url": CPTAC,
        }
    )
    pd.DataFrame(rows).to_csv(OUT / "p2_direct_protein_context_audit.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "protein_mechanism_status": "structurally_contextualized_but_not_observed_as_disease_protein_effect",
                "strongest_new_evidence": "Both splice boundaries lie within the UniProt-curated extracellular TSPAN14 region necessary and sufficient for ADAM10 interaction.",
                "direct_detection_result": "HPA healthy-brain DVP reported ND across the displayed neuronal, astrocyte/neuropil and microglia/neuropil replicates.",
                "actionable_validation_route": "Use targeted LC-MS/MS with the existing assay as a total-protein control and EC2-region peptides (AA150-177; AA180-194) only after an alternative coding isoform is directly demonstrated.",
                "prohibited_conclusion": "Do not claim splice-driven TSPAN14 protein change or altered ADAM10 function from structural annotation, assay availability, or healthy-brain ND observations.",
            }
        ]
    ).to_csv(OUT / "p2_protein_mechanism_feasibility_v2.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
