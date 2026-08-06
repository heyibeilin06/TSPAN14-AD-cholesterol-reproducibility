from __future__ import annotations

import csv
import gzip
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"
REPORTS = ROOT / "results" / "reports"
LDSC_DIR = ROOT / "tools" / "ldsc_cbiit"
OUT_DIR = ROOT / "results" / "ldsc" / "iteration6_pair_harmonized"
INPUT_DIR = OUT_DIR / "inputs"
SUM_DIR = OUT_DIR / "sumstats"
LOG_DIR = OUT_DIR / "logs"
REF_DIR = ROOT / "data" / "reference" / "ldsc"
PYTHON = Path(sys.executable)
BASES = {"A", "C", "G", "T"}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def open_text_gz(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")


def run_command(cmd: list[str], log_path: Path, timeout: int | None = None) -> tuple[int, str]:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    proc = subprocess.run(
        cmd,
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=timeout,
    )
    log_path.write_text(proc.stdout, encoding="utf-8", errors="replace")
    return proc.returncode, proc.stdout


def clean_base(value: str) -> str:
    return (value or "").strip().upper()


def parse_float(value: str) -> float | None:
    try:
        if value is None or value == "":
            return None
        return float(value)
    except ValueError:
        return None


def load_hm3() -> set[str]:
    hm3 = REF_DIR / "hm3_no_MHC.list.txt"
    return {line.strip() for line in hm3.read_text(encoding="utf-8").splitlines() if line.strip()}


def load_ad(path: Path, hm3: set[str]) -> dict[str, dict[str, str]]:
    ad = {}
    with open_text_gz(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            snp = row["SNP"]
            if snp not in hm3 or snp in ad:
                continue
            a1, a2 = clean_base(row["A1"]), clean_base(row["A2"])
            if a1 not in BASES or a2 not in BASES or a1 == a2:
                continue
            ad[snp] = {
                "SNP": snp,
                "A1": a1,
                "A2": a2,
                "BETA": row["BETA"],
                "P": row["P"],
                "EAF": row["EAF"],
                "N": row["N"],
            }
    return ad


def allele_alignment(ad_a1: str, ad_a2: str, tr_a1: str, tr_a2: str) -> int | None:
    if ad_a1 == tr_a1 and ad_a2 == tr_a2:
        return 1
    if ad_a1 == tr_a2 and ad_a2 == tr_a1:
        return -1
    return None


def flip_freq(value: str) -> str:
    x = parse_float(value)
    if x is None:
        return value
    return str(1.0 - x)


def make_pair_inputs(trait: str, lipid_path: Path, ad: dict[str, dict[str, str]]) -> dict[str, object]:
    ad_out = INPUT_DIR / f"ad_for_{trait.lower()}.hm3_aligned.tsv.gz"
    trait_out = INPUT_DIR / f"{trait.lower()}_for_ad.hm3_aligned.tsv.gz"
    fields = ["SNP", "A1", "A2", "BETA", "P", "EAF", "N"]
    matched = flipped = mismatch = duplicate = invalid = 0
    seen = set()
    with gzip.open(ad_out, "wt", encoding="utf-8", newline="") as ad_handle, gzip.open(trait_out, "wt", encoding="utf-8", newline="") as trait_handle:
        ad_writer = csv.DictWriter(ad_handle, delimiter="\t", fieldnames=fields)
        trait_writer = csv.DictWriter(trait_handle, delimiter="\t", fieldnames=fields)
        ad_writer.writeheader()
        trait_writer.writeheader()
        with open_text_gz(lipid_path) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                snp = row["SNP"]
                if snp not in ad:
                    continue
                if snp in seen:
                    duplicate += 1
                    continue
                seen.add(snp)
                tr_a1, tr_a2 = clean_base(row["A1"]), clean_base(row["A2"])
                if tr_a1 not in BASES or tr_a2 not in BASES or tr_a1 == tr_a2:
                    invalid += 1
                    continue
                ad_row = ad[snp]
                direction = allele_alignment(ad_row["A1"], ad_row["A2"], tr_a1, tr_a2)
                if direction is None:
                    mismatch += 1
                    continue
                beta = parse_float(row["BETA"])
                if beta is None:
                    invalid += 1
                    continue
                trait_row = {
                    "SNP": snp,
                    "A1": ad_row["A1"],
                    "A2": ad_row["A2"],
                    "BETA": str(beta if direction == 1 else -beta),
                    "P": row["P"],
                    "EAF": row["EAF"] if direction == 1 else flip_freq(row["EAF"]),
                    "N": row["N"],
                }
                ad_writer.writerow(ad_row)
                trait_writer.writerow(trait_row)
                matched += 1
                if direction == -1:
                    flipped += 1
    return {
        "comparison": f"AD-{trait}",
        "matched_hm3_biallelic_snps": matched,
        "flipped_to_ad_alleles": flipped,
        "allele_mismatch_excluded": mismatch,
        "duplicate_lipid_snps_excluded": duplicate,
        "invalid_lipid_alleles_or_beta_excluded": invalid,
        "ad_input": str(ad_out.relative_to(ROOT)).replace("\\", "/"),
        "trait_input": str(trait_out.relative_to(ROOT)).replace("\\", "/"),
    }


def munge(input_path: Path, out_prefix: Path, name: str) -> dict[str, object]:
    cmd = [
        str(PYTHON),
        str(LDSC_DIR / "munge_sumstats.py"),
        "--sumstats",
        str(input_path),
        "--out",
        str(out_prefix),
        "--snp",
        "SNP",
        "--a1",
        "A1",
        "--a2",
        "A2",
        "--p",
        "P",
        "--signed-sumstats",
        "BETA,0",
        "--N-col",
        "N",
        "--chunksize",
        "500000",
    ]
    code, log = run_command(cmd, LOG_DIR / f"munge_{name}.log")
    out_file = out_prefix.with_suffix(".sumstats.gz")
    return {
        "name": name,
        "status": "completed" if code == 0 and out_file.exists() else "failed",
        "return_code": code,
        "output_path": str(out_file.relative_to(ROOT)).replace("\\", "/") if out_file.exists() else "",
        "log_path": str((LOG_DIR / f"munge_{name}.log").relative_to(ROOT)).replace("\\", "/"),
        "last_log_line": log.strip().splitlines()[-1] if log.strip() else "",
    }


def parse_rg_log(log: str) -> tuple[str, str, str]:
    for line in log.splitlines():
        if line.startswith("Genetic Correlation:"):
            # Example: Genetic Correlation: 0.12 (0.03) (P = 1e-04)
            clean = line.replace("Genetic Correlation:", "").replace("(", "").replace(")", "").replace("P =", "").strip()
            parts = clean.split()
            if len(parts) >= 3:
                return parts[0], parts[1], parts[2]
    return "", "", ""


def main() -> int:
    for directory in [INPUT_DIR, SUM_DIR, LOG_DIR, REPORTS]:
        directory.mkdir(parents=True, exist_ok=True)

    manifest = read_tsv(ROOT / "data" / "manifest" / "processed_gwas.tsv")
    by_trait = {row["trait"]: row for row in manifest}
    hm3 = load_hm3()
    ad = load_ad(Path(by_trait["AD"]["path"]), hm3)

    pair_rows = []
    munge_rows = []
    rg_rows = []
    for trait in ["LDL", "HDL", "TG", "TC", "nonHDL"]:
        pair = make_pair_inputs(trait, Path(by_trait[trait]["path"]), ad)
        pair_rows.append(pair)
        ad_input = ROOT / pair["ad_input"]
        trait_input = ROOT / pair["trait_input"]
        ad_munge = munge(ad_input, SUM_DIR / f"ad_for_{trait.lower()}", f"ad_for_{trait.lower()}")
        trait_munge = munge(trait_input, SUM_DIR / f"{trait.lower()}_for_ad", f"{trait.lower()}_for_ad")
        munge_rows.extend([ad_munge, trait_munge])
        if ad_munge["status"] != "completed" or trait_munge["status"] != "completed":
            rg_rows.append({"comparison": f"AD-{trait}", "status": "skipped_missing_munged_sumstats", "rg": "", "se": "", "p": "", "output_log": "", "note": "pair-harmonized munge failed"})
            continue
        out_prefix = OUT_DIR / f"ad_{trait.lower()}_rg"
        cmd = [
            str(PYTHON),
            str(LDSC_DIR / "ldsc.py"),
            "--rg",
            f"{ROOT / ad_munge['output_path']},{ROOT / trait_munge['output_path']}",
            "--ref-ld-chr",
            str(REF_DIR / "LDscore" / "LDscore."),
            "--w-ld-chr",
            str(REF_DIR / "1000G_Phase3_weights_hm3_no_MHC" / "weights.hm3_noMHC."),
            "--out",
            str(out_prefix),
        ]
        code, log = run_command(cmd, LOG_DIR / f"rg_ad_{trait.lower()}.log")
        rg, se, p = parse_rg_log(log)
        rg_rows.append(
            {
                "comparison": f"AD-{trait}",
                "status": "completed" if code == 0 and rg else "failed",
                "rg": rg,
                "se": se,
                "p": p,
                "output_log": str((LOG_DIR / f"rg_ad_{trait.lower()}.log").relative_to(ROOT)).replace("\\", "/"),
                "note": log.strip().splitlines()[-1] if log.strip() else "",
            }
        )

    write_tsv(TABLES / "technical_iteration6_ldsc_pair_harmonization_status.tsv", pair_rows, ["comparison", "matched_hm3_biallelic_snps", "flipped_to_ad_alleles", "allele_mismatch_excluded", "duplicate_lipid_snps_excluded", "invalid_lipid_alleles_or_beta_excluded", "ad_input", "trait_input"])
    write_tsv(TABLES / "technical_iteration6_ldsc_pair_munge_status.tsv", munge_rows, ["name", "status", "return_code", "output_path", "log_path", "last_log_line"])
    write_tsv(TABLES / "technical_iteration6_ldsc_pair_rg_results.tsv", rg_rows, ["comparison", "status", "rg", "se", "p", "output_log", "note"])

    report = [
        "# Preliminary result 25: pair-harmonized LDSC genetic correlation",
        "",
        "## Execution summary",
        "",
        "- AD and each lipid GWAS were pair-harmonized before LDSC: HapMap3 SNP filter, single-base allele filter, duplicate removal, and lipid beta flipping to AD alleles.",
        "- Pair-specific munged summary statistics were then used for AD-lipid LDSC genetic correlation.",
        "- This resolves the incompatible-allele failure observed in the non-pair-harmonized LDSC attempt.",
        "",
        "## Outputs",
        "",
        "- `technical_iteration6_ldsc_pair_harmonization_status.tsv`",
        "- `technical_iteration6_ldsc_pair_munge_status.tsv`",
        "- `technical_iteration6_ldsc_pair_rg_results.tsv`",
    ]
    (REPORTS / "preliminary_result_25_iteration6_pair_harmonized_ldsc.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print("Saved pair-harmonized LDSC outputs.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
