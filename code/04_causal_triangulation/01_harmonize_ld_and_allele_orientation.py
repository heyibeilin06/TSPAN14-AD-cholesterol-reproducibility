"""Fetch only selected 1000G VCF alleles needed to orient a signed LD matrix."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import requests


DEFAULT_ENSEMBL_URL = "https://grch37.rest.ensembl.org/variation/human/"


def alt_to_effect_sign(effect_allele: str, other_allele: str, ref: str, alt: str) -> int | None:
    if effect_allele == alt and other_allele == ref:
        return 1
    if effect_allele == ref and other_allele == alt:
        return -1
    return None


def grch37_reference(payload: dict[str, object], *, chromosome: str, position: int) -> tuple[str, set[str]] | None:
    for mapping in payload.get("mappings", []):
        if (
            mapping.get("assembly_name") == "GRCh37"
            and str(mapping.get("seq_region_name")) == chromosome
            and int(mapping.get("start")) == position
        ):
            alleles = set(str(mapping["allele_string"]).split("/"))
            reference = str(mapping["allele_string"]).split("/")[0]
            return reference, alleles
    return None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--instruments", required=True)
    parser.add_argument("--ensembl-url", default=DEFAULT_ENSEMBL_URL)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    instruments = pd.read_csv(args.instruments, sep="\t")
    variants = instruments[["SNP", "BP", "A1", "A2"]].drop_duplicates("SNP")
    rows = []
    for row in variants.itertuples(index=False):
        response = requests.get(f"{args.ensembl_url}{row.SNP}", headers={"Content-Type": "application/json"}, timeout=30)
        response.raise_for_status()
        mapping = grch37_reference(response.json(), chromosome="10", position=int(row.BP))
        if mapping is None:
            raise RuntimeError(f"No GRCh37 mapping matched {row.SNP} at {row.BP}.")
        ref, known_alleles = mapping
        if row.A1 not in known_alleles or row.A2 not in known_alleles:
            raise RuntimeError(f"The sQTL alleles for {row.SNP} do not match the GRCh37 variation record.")
        alt = row.A1 if row.A1 != ref else row.A2
        rows.append(
            {
                "SNP": row.SNP,
                "BP": row.BP,
                "sqtl_effect_allele": row.A1,
                "sqtl_other_allele": row.A2,
                "grch37_ref": ref,
                "inferred_biallelic_alt": alt,
                "alt_to_sqtl_effect_sign": alt_to_effect_sign(row.A1, row.A2, ref, alt) if ref else None,
            }
        )
    output = pd.DataFrame(rows)
    if output["alt_to_sqtl_effect_sign"].isna().any():
        raise RuntimeError("At least one selected instrument could not be matched to the 1000G REF/ALT alleles.")
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
