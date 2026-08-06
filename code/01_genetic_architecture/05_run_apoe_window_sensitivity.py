"""Recompute AD-lipid LDSC after excluding graded APOE lead-anchor windows.

The external data root is supplied through ``SLM_DATA_ROOT``. Each trait is
munged once, then small munged summary statistics are filtered for APOE
windows before pairwise LDSC.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DATA = Path(os.environ.get("SLM_DATA_ROOT", ROOT / "data" / "raw"))
LDSC_ROOT = DATA / "p0_ldsc"
LDSC_CODE = LDSC_ROOT / "cbii_ldsc"
LDSC_PYTHON_DEPS = LDSC_ROOT / "python_minimal"
REFERENCE = LDSC_ROOT / "eur_w_ld_chr"
WORK = LDSC_ROOT / "apoe_window_reanalysis"
OUT = ROOT / "outputs" / "p0_reanalysis"

APOE_LEAD_HG38 = 44_890_259
WINDOWS = {
    "w250kb": 250_000,
    "w500kb": 500_000,
    "w1Mb": 1_000_000,
    "w2Mb": 2_000_000,
    "w5Mb": 5_000_000,
}
TRAIT_FILES = {
    "AD": DATA / "processed" / "ad_bellenguez_2022_harmonized.tsv.gz",
    "HDL": DATA / "processed" / "glgc_hdl_eur_harmonized.tsv.gz",
    "LDL": DATA / "processed" / "glgc_ldl_eur_harmonized.tsv.gz",
    "TG": DATA / "processed" / "glgc_tg_eur_harmonized.tsv.gz",
    "TC": DATA / "processed" / "glgc_tc_eur_harmonized.tsv.gz",
    "nonHDL": DATA / "processed" / "glgc_nonhdl_eur_harmonized.tsv.gz",
}
PAIRS = ("HDL", "LDL", "TG", "TC", "nonHDL")


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def append_tsv(path: Path, rows: list[dict[str, object]], key_fields: tuple[str, ...]) -> None:
    old: list[dict[str, str]] = []
    if path.exists():
        with path.open(encoding="utf-8", newline="") as handle:
            old = list(csv.DictReader(handle, delimiter="\t"))
    merged = {tuple(str(row[field]) for field in key_fields): row for row in old}
    for row in rows:
        merged[tuple(str(row[field]) for field in key_fields)] = row
    write_tsv(path, list(merged.values()))


def run(command: list[str], log_path: Path) -> tuple[int, str]:
    env = os.environ.copy()
    env["PYTHONPATH"] = str(LDSC_PYTHON_DEPS)
    completed = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=env)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(completed.stdout, encoding="utf-8", errors="replace")
    return completed.returncode, completed.stdout


def munged_path(trait: str) -> Path:
    return WORK / "munged" / f"{trait}.sumstats.gz"


def munge_trait(trait: str) -> dict[str, object]:
    output = munged_path(trait)
    source = TRAIT_FILES[trait]
    if output.exists():
        return {"trait": trait, "status": "reused", "source": str(source), "munged": str(output), "log": ""}
    output.parent.mkdir(parents=True, exist_ok=True)
    prefix = output.parent / trait
    command = [
        sys.executable,
        str(LDSC_CODE / "munge_sumstats.py"),
        "--sumstats", str(source), "--out", str(prefix),
        "--snp", "SNP", "--a1", "A1", "--a2", "A2", "--p", "P",
        "--signed-sumstats", "BETA,0", "--N-col", "N", "--chunksize", "500000",
        "--merge-alleles", str(REFERENCE / "w_hm3.snplist"),
    ]
    code, log = run(command, WORK / "logs" / f"munge_{trait}.log")
    return {"trait": trait, "status": "completed" if code == 0 and output.exists() else "failed", "source": str(source), "munged": str(output) if output.exists() else "", "log": str(WORK / "logs" / f"munge_{trait}.log")}


def build_exclusion_sets() -> list[dict[str, object]]:
    sets = {label: set() for label in WINDOWS}
    for source in TRAIT_FILES.values():
        with gzip.open(source, "rt", encoding="utf-8", errors="replace", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                if row.get("CHR", "").strip() != "19":
                    continue
                try:
                    position = int(row["BP"])
                except (KeyError, TypeError, ValueError):
                    continue
                for label, radius in WINDOWS.items():
                    if abs(position - APOE_LEAD_HG38) <= radius and row.get("SNP"):
                        sets[label].add(row["SNP"])
    rows = []
    exclusion_dir = WORK / "exclude"
    exclusion_dir.mkdir(parents=True, exist_ok=True)
    for label, rsids in sets.items():
        path = exclusion_dir / f"{label}_rsids.txt"
        path.write_text("\n".join(sorted(rsids)) + "\n", encoding="utf-8")
        rows.append({
            "window": label, "radius_bp": WINDOWS[label], "anchor_hg38": APOE_LEAD_HG38,
            "excluded_rsids_union": len(rsids), "rsid_file": str(path),
            "scope": "Union of rsIDs within the stated hg38 chr19 lead-anchor window across retained AD and lipid processed GWAS inputs.",
        })
    write_tsv(OUT / "p0_apoe_window_exclusion_manifest.tsv", rows)
    return rows


def filter_munged_sumstats(trait: str, window: str) -> tuple[Path, int, int]:
    exclusion = set((WORK / "exclude" / f"{window}_rsids.txt").read_text(encoding="utf-8").splitlines())
    source = munged_path(trait)
    output = WORK / "filtered" / window / f"{trait}.sumstats.gz"
    if output.exists():
        with gzip.open(output, "rt", encoding="utf-8") as handle:
            retained = sum(1 for _ in handle) - 1
        with gzip.open(source, "rt", encoding="utf-8") as handle:
            base_total = sum(1 for _ in handle) - 1
        return output, retained, base_total - retained
    output.parent.mkdir(parents=True, exist_ok=True)
    retained = removed = 0
    with gzip.open(source, "rt", encoding="utf-8") as reader, gzip.open(output, "wt", encoding="utf-8") as writer:
        header = reader.readline()
        writer.write(header)
        snp_index = header.rstrip("\n").split("\t").index("SNP")
        for line in reader:
            fields = line.rstrip("\n").split("\t")
            if fields[snp_index] in exclusion:
                removed += 1
            else:
                writer.write(line)
                retained += 1
    return output, retained, removed


def parse_rg(log: str) -> tuple[str, str, str]:
    match = re.search(r"^Genetic Correlation:\s*([^\s]+)\s*\(([^)]+)\)", log, flags=re.MULTILINE)
    p_match = re.search(r"^P:\s*([^\s]+)", log, flags=re.MULTILINE)
    if match and p_match:
        return match.group(1), match.group(2), p_match.group(1)
    return "", "", ""


def run_window(window: str) -> None:
    if window != "baseline" and not (WORK / "exclude" / f"{window}_rsids.txt").exists():
        raise RuntimeError("Build APOE exclusion sets before LDSC window runs")
    results: list[dict[str, object]] = []
    for trait in PAIRS:
        if window == "baseline":
            ad, lipid = munged_path("AD"), munged_path(trait)
            ad_retained = lipid_retained = "base_munged"
            ad_removed = lipid_removed = 0
        else:
            ad, ad_retained, ad_removed = filter_munged_sumstats("AD", window)
            lipid, lipid_retained, lipid_removed = filter_munged_sumstats(trait, window)
        prefix = WORK / "rg" / window / f"AD_{trait}"
        prefix.parent.mkdir(parents=True, exist_ok=True)
        command = [
            sys.executable, str(LDSC_CODE / "ldsc.py"), "--rg", f"{ad},{lipid}",
            "--ref-ld-chr", str(REFERENCE) + os.sep, "--w-ld-chr", str(REFERENCE) + os.sep,
            "--out", str(prefix),
        ]
        code, log = run(command, WORK / "logs" / f"rg_{window}_AD_{trait}.log")
        rg, se, p = parse_rg(log)
        results.append({
            "window": window, "comparison": f"AD-{trait}", "status": "completed" if code == 0 and rg else "failed",
            "rg": rg, "se": se, "p": p, "ad_retained": ad_retained, "lipid_retained": lipid_retained,
            "ad_removed": ad_removed, "lipid_removed": lipid_removed,
            "log": str(WORK / "logs" / f"rg_{window}_AD_{trait}.log"),
            "interpretation_boundary": "Chromosome-19 lead-anchor sensitivity for genome-wide LDSC; does not condition individual variants or establish mediation.",
        })
    append_tsv(OUT / "p0_apoe_window_ldsc_results.tsv", results, ("window", "comparison"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("prepare", "run"))
    parser.add_argument("--traits", default="AD,HDL,LDL,TG,TC,nonHDL")
    parser.add_argument("--window", choices=("baseline", *WINDOWS), default="baseline")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    if args.stage == "prepare":
        traits = tuple(item for item in args.traits.split(",") if item)
        status = [munge_trait(trait) for trait in traits]
        write_tsv(OUT / "p0_ldsc_munge_status.tsv", status)
        if set(traits) == set(TRAIT_FILES):
            build_exclusion_sets()
    else:
        run_window(args.window)


if __name__ == "__main__":
    main()
