#!/usr/bin/env python3
"""Export manuscript tables as machine-readable TSV source files."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from docx import Document


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manuscript", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    document = Document(args.manuscript)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for number, table in enumerate(document.tables, start=1):
        path = args.output_dir / f"Table_{number}_source.tsv"
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            for row in table.rows:
                writer.writerow([cell.text.strip() for cell in row.cells])
        print(path)


if __name__ == "__main__":
    main()
