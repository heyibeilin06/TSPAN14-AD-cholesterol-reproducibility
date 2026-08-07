#!/usr/bin/env python3
"""Refresh the portable SHA-256 ledger for repository release files."""

from __future__ import annotations

import hashlib
import csv
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "audit" / "file_checksums.tsv"
TEXT_SUFFIXES = {
    "", ".cff", ".csv", ".json", ".md", ".mjs", ".py", ".r", ".sh",
    ".svg", ".tsv", ".txt", ".yaml", ".yml", ".ps1",
}


def digest(path: Path) -> str:
    sha = hashlib.sha256()
    if path.suffix.lower() in TEXT_SUFFIXES:
        sha.update(path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n"))
    else:
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                sha.update(chunk)
    return sha.hexdigest()


def main() -> None:
    lineage_path = ROOT / "audit" / "figure_panel_lineage.tsv"
    with lineage_path.open("r", encoding="utf-8", newline="") as handle:
        lineage = list(csv.DictReader(handle, delimiter="\t"))
    for row in lineage:
        row["sha256"] = digest(ROOT / row["compatibility_data"])
    with lineage_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=lineage[0].keys(), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(lineage)

    tracked = subprocess.run(
        ["git", "ls-files", "--cached"], cwd=ROOT, check=True, capture_output=True, text=True
    ).stdout.splitlines()
    rows = []
    for relative in sorted(tracked):
        path = ROOT / relative
        if not path.is_file() or path == OUTPUT:
            continue
        rows.append((relative, path.stat().st_size, digest(path)))
    OUTPUT.write_text(
        "relative_path\tsize_bytes\tsha256\n" +
        "\n".join(f"{relative}\t{size}\t{sha}" for relative, size, sha in rows) + "\n",
        encoding="utf-8",
    )
    print(f"Wrote {len(rows)} checksum records to {OUTPUT}")


if __name__ == "__main__":
    main()
