"""Jointly clump LDL/HDL/TG instruments and prepare a global lipid MVMR input."""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

import pandas as pd


def load_core(path: Path):
    spec = importlib.util.spec_from_file_location("complete_mr_core", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--core-script", type=Path, required=True)
    parser.add_argument("--bidirectional", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--ld-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    core = load_core(args.core_script)

    bidirectional = pd.read_csv(args.bidirectional, sep="\t")
    selected_traits = ("LDL", "HDL", "TG")
    candidates = bidirectional.loc[
        bidirectional["exposure"].isin(selected_traits)
        & bidirectional["outcome"].eq("AD")
    ].copy()
    candidates = (
        candidates.sort_values("p_x")
        .drop_duplicates("SNP")
        .rename(columns={"p_x": "P"})
    )
    joint = core.ld_clump_by_blocks(
        candidates[["SNP", "CHR", "BP", "P"]], args.ld_root, 0.001
    )
    snps = set(joint["SNP"])

    common = (
        candidates.loc[candidates["SNP"].isin(snps)]
        .drop_duplicates("SNP")
        .rename(
            columns={
                "effect_allele": "A1",
                "other_allele": "A2",
                "beta_y": "BETA",
                "se_y": "SE",
                "p_y": "P_OUTCOME",
                "eaf_y": "EAF",
                "n_y": "N",
            }
        )
    )
    common["P"] = common["P_OUTCOME"]
    common = common[["SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "EAF", "N"]]

    result = common.rename(columns={"BETA": "beta_AD", "SE": "se_AD"})[
        ["SNP", "CHR", "BP", "A1", "A2", "beta_AD", "se_AD"]
    ]
    file_names = {
        "LDL": "glgc_ldl_eur_harmonized.tsv.gz",
        "HDL": "glgc_hdl_eur_harmonized.tsv.gz",
        "TG": "glgc_tg_eur_harmonized.tsv.gz",
    }
    for trait, name in file_names.items():
        raw = core.read_selected(args.data_root / name, snps)
        aligned = core.harmonize(common, raw)[["SNP", "beta_y", "se_y"]].rename(
            columns={"beta_y": f"beta_{trait}", "se_y": f"se_{trait}"}
        )
        result = result.merge(aligned, on="SNP", how="inner")
    result["joint_instrument_rule"] = (
        "union of LDL/HDL/TG genome-wide-significant instruments, jointly UKB-EUR clumped at r2<0.001"
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
