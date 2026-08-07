#!/usr/bin/env python3
"""Verify manuscript legends and canonical main-figure files."""

from __future__ import annotations

import argparse
import hashlib
import re
from pathlib import Path

from docx import Document


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--submission-root", required=True, type=Path)
    parser.add_argument("--package-root", required=True, type=Path)
    parser.add_argument("--manuscript", required=True, type=Path)
    args = parser.parse_args()

    document = Document(args.manuscript)
    legends = [
        p for p in document.paragraphs
        if re.match(r"^Figure [1-6] \|", p.text)
    ]
    if len(legends) != 6:
        raise SystemExit(f"Expected 6 main-figure legends; found {len(legends)}")

    for number, paragraph in enumerate(legends, 1):
        expected = (
            args.package_root
            / f"figures/Figure_{number:02d}/legend.md"
        ).read_text(encoding="utf-8").strip()
        if paragraph.text != expected:
            raise SystemExit(f"Figure {number} legend differs from legend.md")
        for run in (run for run in paragraph.runs if run.text):
            if run.font.name != "Times New Roman":
                raise SystemExit(f"Figure {number} legend has a non-Times run")
            if run.font.size is None or abs(run.font.size.pt - 12) > 0.01:
                raise SystemExit(f"Figure {number} legend has a non-12-pt run")
        if paragraph.paragraph_format.line_spacing != 2.0:
            raise SystemExit(f"Figure {number} legend is not double-spaced")

        submission_pdf = args.submission_root / f"Figures/Figure_{number}.pdf"
        package_pdf = (
            args.package_root
            / f"figures/Figure_{number:02d}/output/Figure_{number}.pdf"
        )
        if sha256(submission_pdf) != sha256(package_pdf):
            raise SystemExit(f"Figure {number} submission/package PDFs differ")

    print("Main-figure synchronization passed: 6 legends and 6 PDFs")


if __name__ == "__main__":
    main()
