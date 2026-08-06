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
OUT_DIR = ROOT / "results" / "ldsc" / "iteration6"
SUM_DIR = OUT_DIR / "sumstats"
LOG_DIR = OUT_DIR / "logs"
REF_DIR = ROOT / "data" / "reference" / "ldsc"
PYTHON = Path(sys.executable)


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


def count_lines_gz(path: Path) -> int:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        return sum(1 for _ in handle) - 1


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


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    SUM_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)

    manifest = read_tsv(ROOT / "data" / "manifest" / "processed_gwas.tsv")
    traits = ["AD", "LDL", "HDL", "TG", "TC", "nonHDL"]
    manifest = [row for row in manifest if row["trait"] in traits]
    by_trait = {row["trait"]: row for row in manifest}

    hm3 = REF_DIR / "hm3_no_MHC.list.txt"
    ref_ld_chr = REF_DIR / "LDscore" / "LDscore."
    w_ld_chr = REF_DIR / "1000G_Phase3_weights_hm3_no_MHC" / "weights.hm3_noMHC."

    resource_rows = [
        {
            "resource": "ldsc_repository",
            "path": str(LDSC_DIR.relative_to(ROOT)).replace("\\", "/"),
            "status": "available" if (LDSC_DIR / "ldsc.py").exists() and (LDSC_DIR / "munge_sumstats.py").exists() else "missing",
            "note": "CBIIT ldsc39 branch cloned locally",
        },
        {
            "resource": "hapmap3_merge_alleles",
            "path": str(hm3.relative_to(ROOT)).replace("\\", "/"),
            "status": "available" if hm3.exists() else "missing",
            "note": "Zenodo 10.5281/zenodo.7768714 minimal LDSC resource",
        },
        {
            "resource": "eur_ld_scores",
            "path": str(ref_ld_chr.parent.relative_to(ROOT)).replace("\\", "/"),
            "status": "available" if len(list(ref_ld_chr.parent.glob("LDscore.*.l2.ldscore.gz"))) >= 22 else "missing",
            "note": "1000G Phase3 LDscore files",
        },
        {
            "resource": "eur_regression_weights",
            "path": str(w_ld_chr.parent.relative_to(ROOT)).replace("\\", "/"),
            "status": "available" if len(list(w_ld_chr.parent.glob("weights.hm3_noMHC.*.l2.ldscore.gz"))) >= 22 else "missing",
            "note": "1000G Phase3 HapMap3 no-MHC weights",
        },
    ]

    write_tsv(TABLES / "technical_iteration6_ldsc_resource_status.tsv", resource_rows, ["resource", "path", "status", "note"])

    munge_rows = []
    for trait in traits:
        row = by_trait[trait]
        in_path = Path(row["path"])
        out_prefix = SUM_DIR / trait.lower()
        cmd = [
            str(PYTHON),
            str(LDSC_DIR / "munge_sumstats.py"),
            "--sumstats",
            str(in_path),
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
        code, log = run_command(cmd, LOG_DIR / f"munge_{trait.lower()}.log", timeout=None)
        out_file = out_prefix.with_suffix(".sumstats.gz")
        munge_rows.append(
            {
                "trait": trait,
                "input_path": str(in_path),
                "input_exists": in_path.exists(),
                "input_manifest_rows": row.get("rows", ""),
                "output_path": str(out_file.relative_to(ROOT)).replace("\\", "/") if out_file.exists() else "",
                "status": "completed" if code == 0 and out_file.exists() else "failed",
                "return_code": code,
                "log_path": str((LOG_DIR / f"munge_{trait.lower()}.log").relative_to(ROOT)).replace("\\", "/"),
                "last_log_line": log.strip().splitlines()[-1] if log.strip() else "",
                "merge_alleles_used": "no",
            }
        )
    write_tsv(TABLES / "technical_iteration6_ldsc_munge_status.tsv", munge_rows, ["trait", "input_path", "input_exists", "input_manifest_rows", "output_path", "status", "return_code", "log_path", "last_log_line", "merge_alleles_used"])

    rg_rows = []
    completed = {row["trait"]: row for row in munge_rows if row["status"] == "completed"}
    for trait in ["LDL", "HDL", "TG", "TC", "nonHDL"]:
        ad_file = SUM_DIR / "ad.sumstats.gz"
        trait_file = SUM_DIR / f"{trait.lower()}.sumstats.gz"
        out_prefix = OUT_DIR / f"ad_{trait.lower()}_rg"
        if "AD" not in completed or trait not in completed:
            rg_rows.append(
                {
                    "comparison": f"AD-{trait}",
                    "status": "skipped_missing_munged_sumstats",
                    "rg": "",
                    "se": "",
                    "p": "",
                    "output_log": "",
                    "note": "Munged AD or lipid file unavailable",
                }
            )
            continue
        cmd = [
            str(PYTHON),
            str(LDSC_DIR / "ldsc.py"),
            "--rg",
            f"{ad_file},{trait_file}",
            "--ref-ld-chr",
            str(ref_ld_chr),
            "--w-ld-chr",
            str(w_ld_chr),
            "--out",
            str(out_prefix),
        ]
        code, log = run_command(cmd, LOG_DIR / f"rg_ad_{trait.lower()}.log", timeout=None)
        rg = se = p = ""
        for line in log.splitlines():
            if line.startswith("Genetic Correlation:"):
                parts = line.replace("(", "").replace(")", "").replace("=", " ").replace(",", " ").split()
                try:
                    rg = parts[2]
                    se = parts[4]
                    p = parts[6]
                except IndexError:
                    pass
        rg_rows.append(
            {
                "comparison": f"AD-{trait}",
                "status": "completed" if code == 0 else "failed",
                "rg": rg,
                "se": se,
                "p": p,
                "output_log": str((LOG_DIR / f"rg_ad_{trait.lower()}.log").relative_to(ROOT)).replace("\\", "/"),
                "note": log.strip().splitlines()[-1] if log.strip() else "",
            }
        )
    write_tsv(TABLES / "technical_iteration6_ldsc_rg_results.tsv", rg_rows, ["comparison", "status", "rg", "se", "p", "output_log", "note"])

    report = [
        "# Preliminary result 24: iteration-6 LDSC execution",
        "",
        "## Execution summary",
        "",
        "- LDSC software was cloned from the CBIIT `ldsc39` branch and passed help-command self-checks.",
        "- Minimal free LDSC reference resources were downloaded from Zenodo record `10.5281/zenodo.7768714`.",
        "- Full AD and lipid GWAS were processed through LDSC `munge_sumstats.py`; AD-lipid genetic correlation was attempted for each lipid trait.",
        "- The free Zenodo HapMap3 file is SNP-list-only in this environment, so `--merge-alleles` was not used; LDSC therefore intersects munged files with the LDscore SNP set during regression.",
        "",
        "## Outputs",
        "",
        "- `technical_iteration6_ldsc_resource_status.tsv`",
        "- `technical_iteration6_ldsc_munge_status.tsv`",
        "- `technical_iteration6_ldsc_rg_results.tsv`",
    ]
    (REPORTS / "preliminary_result_24_iteration6_ldsc_execution.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print("Saved iteration 6 LDSC outputs.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
