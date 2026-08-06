#!/usr/bin/env python3
"""Audit shared-prior sensitivity for exact exon5-6 colocalization."""

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
COLOC = ROOT / "figures" / "Figure_02" / "data" / "Figure_2_exact_event_coloc.tsv"
SCATTER = ROOT / "figures" / "Figure_02" / "data" / "Figure_2_colocalization_scatter.tsv"
OUTPUT = ROOT / "audit" / "reviewer_revision" / "exact_event_coloc_prior_sensitivity.tsv"


def rescale_h4(pp_h4: float, multiplier: float) -> float:
    odds = pp_h4 / (1.0 - pp_h4)
    scaled_odds = odds * multiplier
    return scaled_odds / (1.0 + scaled_odds)


def main() -> None:
    coloc = pd.read_csv(COLOC, sep="\t")
    scatter = pd.read_csv(SCATTER, sep="\t")
    counts = scatter.groupby("comparison")["snp"].nunique().to_dict()
    count_labels = {"AD": "AD", "TC": "TC", "LDL-C": "LDL", "non-HDL-C": "nonHDL"}

    rows = []
    for row in coloc.itertuples(index=False):
        pp_h4 = float(row.pph4)
        rows.append(
            {
                "trait_pair": f"canonical exon5-6 sQTL-{row.trait}",
                "region_grch38": "chr10:80450000-80550000",
                "allele_compatible_overlapping_snps": int(counts[count_labels[row.trait]]),
                "p1": 1e-4,
                "p2": 1e-4,
                "p12_default": 1e-5,
                "pp_h4_p12_div_10": rescale_h4(pp_h4, 0.1),
                "pp_h4_default": pp_h4,
                "pp_h4_p12_times_10": rescale_h4(pp_h4, 10.0),
                "method": "coloc.abf posterior-odds rescaling with variant ABFs fixed",
                "interpretation": "exact coordinate-matched canonical-branch colocalization prior sensitivity",
            }
        )

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(OUTPUT, sep="\t", index=False)
    print(OUTPUT)


if __name__ == "__main__":
    main()
