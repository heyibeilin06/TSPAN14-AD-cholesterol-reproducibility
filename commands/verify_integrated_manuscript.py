#!/usr/bin/env python3
"""Structural checks for the reviewer-readable integrated manuscript."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

from docx import Document


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manuscript", required=True, type=Path)
    args = parser.parse_args()

    document = Document(args.manuscript)
    legends = [
        p for p in document.paragraphs
        if re.match(r"^Figure [1-6] \|", p.text)
    ]
    if len(legends) != 6:
        raise SystemExit(f"Expected 6 legends; found {len(legends)}")
    if len(document.inline_shapes) != 6:
        raise SystemExit(f"Expected 6 inline figures; found {len(document.inline_shapes)}")

    for number, legend in enumerate(legends, 1):
        previous = legend._p.getprevious()
        if previous is None or not previous.xpath(".//w:drawing"):
            raise SystemExit(f"Figure {number} is not immediately before its legend")
        if not legend.paragraph_format.page_break_before:
            raise SystemExit(f"Figure {number} legend does not start on a new page")

    for number, shape in enumerate(document.inline_shapes, 1):
        doc_pr = shape._inline.docPr
        if doc_pr.get("title") != f"Figure {number}" or not doc_pr.get("descr"):
            raise SystemExit(f"Figure {number} alternative text is incomplete")

    print("Integrated manuscript passed: 6 figures, 6 legends, alt text and ordering")


if __name__ == "__main__":
    main()
