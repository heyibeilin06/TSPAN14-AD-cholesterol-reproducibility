#!/usr/bin/env python3
"""Calculate signed r and D-prime for phased TSPAN14 functional variants."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


VARIANTS = {
    82251544: {"rsid": "rs7922621", "risk": "A", "grch38": 80491788},
    82269461: {"rsid": "rs7080009", "risk": "C", "grch38": 80509705},
    82269611: {"rsid": "rs1870138", "risk": "G", "grch38": 80509855},
    82269848: {"rsid": "rs1870137", "risk": "G", "grch38": 80510092},
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    raw = pd.read_csv(args.input, sep="\t", header=None, dtype=str)
    rows = {}
    for _, row in raw.iterrows():
        pos = int(row.iloc[1])
        if pos not in VARIANTS:
            continue
        alleles = [row.iloc[3]] + row.iloc[4].split(",")
        risk = VARIANTS[pos]["risk"]
        if risk not in alleles:
            raise ValueError(f"Risk allele {risk} absent at {pos}: {alleles}")
        risk_index = alleles.index(risk)
        haps = []
        for gt in row.iloc[5:]:
            if gt in {".", "./.", ".|."}:
                haps.extend([np.nan, np.nan])
                continue
            phased = gt.replace("/", "|").split("|")
            haps.extend([int(x) == risk_index for x in phased])
        rows[pos] = {"alleles": "/".join(alleles), "risk_haplotype": np.asarray(haps, float)}

    sentinel = rows[82269461]["risk_haplotype"]
    output = []
    for pos, meta in VARIANTS.items():
        x = sentinel
        y = rows[pos]["risk_haplotype"]
        keep = np.isfinite(x) & np.isfinite(y)
        x = x[keep]
        y = y[keep]
        p_a = x.mean()
        p_b = y.mean()
        p_ab = np.mean((x == 1) & (y == 1))
        d = p_ab - p_a * p_b
        if d >= 0:
            d_max = min(p_a * (1 - p_b), (1 - p_a) * p_b)
        else:
            d_max = min(p_a * p_b, (1 - p_a) * (1 - p_b))
        d_prime = d / d_max if d_max > 0 else np.nan
        signed_r = np.corrcoef(x, y)[0, 1]
        output.append(
            {
                "variant": meta["rsid"],
                "position_grch37": pos,
                "position_grch38": meta["grch38"],
                "1000g_alleles": rows[pos]["alleles"],
                "ad_risk_allele_grch38_forward": meta["risk"],
                "risk_allele_frequency_1000g_eur": p_b,
                "sentinel": "rs7080009",
                "signed_r_risk_allele_dosage": signed_r,
                "r_squared": signed_r**2,
                "d_prime_risk_phase": d_prime,
                "n_phased_chromosomes": len(x),
                "phase_interpretation": "same risk phase" if signed_r > 0 else "opposite risk phase",
            }
        )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(output).to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
