#!/usr/bin/env python3
"""Report a conservative overlap-bias bound on the rebuilt signed-LD baseline."""

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
SENSITIVITY = ROOT / "audit" / "reviewer_revision" / "cis_mr_numerical_sensitivity.tsv"
ATLAS = ROOT / "figures" / "Figure_04" / "data" / "Figure_4_instrument_effect_atlas.tsv"
OUTPUT = ROOT / "audit" / "reviewer_revision" / "cis_mr_overlap_bias_bound_rebuilt_ld.tsv"


def main() -> None:
    sensitivity = pd.read_csv(SENSITIVITY, sep="\t")
    atlas = pd.read_csv(ATLAS, sep="\t")

    parameter = pd.to_numeric(sensitivity["parameter"], errors="coerce")
    baseline = sensitivity[sensitivity["analysis"].eq("ridge") & parameter.eq(0)].copy()
    exposure = atlas[atlas["layer"].eq("Exact exon5-6 sQTL")]
    minimum_f = float(exposure["f_statistic"].min())

    output = baseline[["outcome", "n_instruments", "estimate", "se", "p_value"]].copy()
    output.insert(1, "baseline_model", "rebuilt signed 1000 Genomes EUR LD; ridge = 0")
    output["molecular_qtl_donors"] = 147
    output["maximum_assumed_overlap_fraction"] = 1.0
    output["minimum_instrument_f"] = minimum_f
    output["relative_weak_instrument_bias_bound_percent"] = 100.0 / minimum_f
    output["interpretation"] = (
        "complete-overlap upper bound; no overlap-adjusted estimate or P value is inferred "
        "without cohort-resolved sampling covariance"
    )

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(OUTPUT, sep="\t", index=False)
    print(OUTPUT)


if __name__ == "__main__":
    main()
