from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"
LOG_DIR = ROOT / "results" / "ldsc" / "iteration6_pair_harmonized" / "logs"


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def parse_log(path: Path) -> dict[str, str]:
    rg = se = z = p = ""
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    for line in lines:
        if line.startswith("Genetic Correlation:"):
            clean = line.replace("Genetic Correlation:", "").replace("(", "").replace(")", "").strip()
            parts = clean.split()
            if len(parts) >= 2:
                rg, se = parts[0], parts[1]
        elif line.startswith("Z-score:"):
            z = line.split(":", 1)[1].strip()
        elif line.startswith("P:"):
            p = line.split(":", 1)[1].strip()
    status = "completed" if rg and se and p else "failed"
    return {"status": status, "rg": rg, "se": se, "z": z, "p": p}


def main() -> int:
    rows = []
    for trait in ["ldl", "hdl", "tg", "tc", "nonhdl"]:
        log_path = LOG_DIR / f"rg_ad_{trait}.log"
        parsed = parse_log(log_path)
        rows.append(
            {
                "comparison": f"AD-{trait.upper() if trait != 'nonhdl' else 'nonHDL'}",
                **parsed,
                "output_log": str(log_path.relative_to(ROOT)).replace("\\", "/"),
                "note": "pair-harmonized LDSC rg parsed from completed log" if parsed["status"] == "completed" else "rg not found in LDSC log",
            }
        )
    write_tsv(TABLES / "technical_iteration6_ldsc_pair_rg_results.tsv", rows, ["comparison", "status", "rg", "se", "z", "p", "output_log", "note"])
    print("Parsed pair-harmonized LDSC logs.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
