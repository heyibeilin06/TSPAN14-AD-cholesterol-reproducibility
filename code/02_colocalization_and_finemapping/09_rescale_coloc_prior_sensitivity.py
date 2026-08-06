#!/usr/bin/env python3
"""Rescale retained coloc posterior odds across a prespecified p12 grid.

With variant-level approximate Bayes factors held fixed, changing only p12
multiplies the H4 posterior odds by p12_new / p12_default. This exactly
reproduces the posterior rescaling for the retained ABF model while avoiding
reconstruction of deleted intermediate files. Regional raw GWAS files are
streamed only to audit the allele-compatible overlapping SNP count.
"""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

import pandas as pd


def read_region(path: Path, start: int, end: int):
    rows = {}
    with gzip.open(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            chrom = row["CHR"].strip()
            if chrom not in {"10", "chr10"}:
                continue
            position = int(row["BP"].strip())
            if start <= position <= end:
                rows[row["SNP"]] = (row["A1"].upper(), row["A2"].upper())
            elif position > end:
                break
    return rows


def allele_compatible(a, b):
    complement = str.maketrans("ACGT", "TGCA")
    a_set = set(a)
    b_set = set(b)
    if a_set == b_set:
        return True
    b_comp = {x.translate(complement) for x in b_set}
    return a_set == b_comp


def rescale_h4(pp_h4, multiplier):
    odds = pp_h4 / (1 - pp_h4)
    new_odds = odds * multiplier
    return new_odds / (1 + new_odds)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--regional-results", type=Path, required=True)
    parser.add_argument("--ad", type=Path, required=True)
    parser.add_argument("--trait", action="append", nargs=2, metavar=("LABEL", "PATH"), required=True)
    parser.add_argument("--start", type=int, default=80244228)
    parser.add_argument("--end", type=int, default=80744228)
    parser.add_argument("--trait-start", type=int, default=82003984)
    parser.add_argument("--trait-end", type=int, default=82503984)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    retained = pd.read_csv(args.regional_results, sep="\t")
    retained = retained[retained["locus"].eq("chr10 TSPAN14")].copy()
    requested = [("AD", args.ad, args.start, args.end)] + [
        (label, Path(path), args.trait_start, args.trait_end) for label, path in args.trait
    ]
    with ProcessPoolExecutor(max_workers=len(requested)) as pool:
        futures = {
            label: pool.submit(read_region, path, start, end)
            for label, path, start, end in requested
        }
        regions = {label: future.result() for label, future in futures.items()}
    ad = regions["AD"]
    outputs = []
    p1 = 1e-4
    p2 = 1e-4
    default_p12 = 1e-5
    for label, path in args.trait:
        trait = regions[label]
        overlap = [s for s in ad.keys() & trait.keys() if allele_compatible(ad[s], trait[s])]
        source = retained[retained["trait"].eq(label)].iloc[0]
        pp_h4 = float(source["PP.H4"])
        outputs.append(
            {
                "trait_pair": f"AD-{source['trait_label']}",
                "region_grch38": f"chr10:{args.start}-{args.end}",
                "allele_compatible_overlapping_snps": len(overlap),
                "p1": p1,
                "p2": p2,
                "p12_default": default_p12,
                "pp_h4_p12_div_10": rescale_h4(pp_h4, 0.1),
                "pp_h4_default": pp_h4,
                "pp_h4_p12_times_10": rescale_h4(pp_h4, 10.0),
                "default_pp_h3": source["PP.H3"],
                "method": "coloc.abf posterior-odds rescaling with variant ABFs fixed",
                "ld_panel": "not required by coloc.abf; multiple-signal sensitivity used rebuilt 1000 Genomes Phase 3 EUR LD",
            }
        )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(outputs).to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
