from __future__ import annotations

import csv
import hashlib
import re
import sys
from pathlib import Path


PACKAGE = Path(__file__).resolve().parents[1]


def read_tsv(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    errors: list[str] = []
    required = [
        "README.md",
        "LICENSE",
        "CITATION.cff",
        "config/data_sources.tsv",
        "config/script_migration.tsv",
        "audit/claims_to_code.tsv",
        "audit/figure_panel_lineage.tsv",
        "audit/file_checksums.tsv",
    ]
    for relative in required:
        if not (PACKAGE / relative).is_file():
            errors.append(f"missing required file: {relative}")

    claim_path = PACKAGE / "audit" / "claims_to_code.tsv"
    if claim_path.exists():
        for row in read_tsv(claim_path):
            if row["reproducibility_status"] != "PASS":
                errors.append(f"claim not reproducible: {row['claim_id']}")
            for script in [item.strip() for item in row["public_scripts"].split(";")]:
                if not (PACKAGE / script).is_file():
                    errors.append(f"claim {row['claim_id']} missing script: {script}")

    lineage_path = PACKAGE / "audit" / "figure_panel_lineage.tsv"
    if lineage_path.exists():
        seen = set()
        for row in read_tsv(lineage_path):
            key = (row["figure"], row["panel"], row["public_data"])
            if key in seen:
                errors.append(f"duplicate figure lineage: {key}")
            seen.add(key)
            for field in ("public_data", "compatibility_data", "plot_script"):
                if not (PACKAGE / row[field]).is_file():
                    errors.append(f"missing figure lineage file: {row[field]}")

    absolute_path_pattern = re.compile(r"[A-Za-z]:[/\\](?:Users|weixin|Desktop)[/\\]", re.I)
    for directory in (PACKAGE / "code", PACKAGE / "figures"):
        for path in directory.rglob("*"):
            if path.suffix.lower() not in {".py", ".r", ".ps1", ".mjs", ".md"} or not path.is_file():
                continue
            text = path.read_text(encoding="utf-8", errors="ignore")
            if absolute_path_pattern.search(text):
                errors.append(f"machine-specific absolute path in {path.relative_to(PACKAGE).as_posix()}")

    for number in range(1, 7):
        directory = PACKAGE / "figures" / f"Figure_{number:02d}"
        for name in ("README.md", "legend.md"):
            if not (directory / name).is_file():
                errors.append(f"Figure {number} missing {name}")
        for extension in ("png", "pdf", "svg", "tiff"):
            if not (directory / "output" / f"Figure_{number}.{extension}").is_file():
                errors.append(f"Figure {number} missing {extension} output")

    checksum_path = PACKAGE / "audit" / "file_checksums.tsv"
    if checksum_path.exists():
        for row in read_tsv(checksum_path):
            path = PACKAGE / row["relative_path"]
            if not path.is_file():
                errors.append(f"checksum target missing: {row['relative_path']}")
            elif sha256(path) != row["sha256"]:
                errors.append(f"checksum mismatch: {row['relative_path']}")

    if errors:
        print("FAIL")
        for error in errors:
            print(f"- {error}")
        return 1
    print("PASS")
    print("All manuscript claims, figure panels, public scripts, source tables, and checksums validated.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
