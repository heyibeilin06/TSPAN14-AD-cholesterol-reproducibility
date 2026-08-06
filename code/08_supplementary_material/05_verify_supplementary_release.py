#!/usr/bin/env python3
"""Release gates for the independently rebuilt supplement v19."""
from __future__ import annotations
import hashlib, json, re
from pathlib import Path
from zipfile import ZipFile
from docx import Document
from pypdf import PdfReader

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/"outputs"/"supplement_v19"
MANUSCRIPT=ROOT/"outputs"/"final_manuscript"/"手稿_主表重制版.docx"
TITLE="A non-APOE Alzheimer disease-cholesterol locus converges on exact TSPAN14 splice regulation"
OLD_TITLE="An APOE-aware Alzheimer disease-lipid screen identifies exact TSPAN14 splice regulation"

def text(path):
    d=Document(path); parts=[p.text for p in d.paragraphs]
    for t in d.tables:
        for r in t.rows: parts.extend(c.text for c in r.cells)
    return "\n".join(parts)

def refs(body, kind):
    label=r"Fig(?:ure)?s?\.?" if kind=="fig" else r"Tables?"
    pattern=rf"Supplementary\s+{label}\s+((?:S?\d+)(?:\s*(?:-|–|—|,|and)\s*S?\d+)*)"
    out=set()
    for m in re.finditer(pattern,body,re.I):
        s=m.group(1)
        for a,b in re.findall(r"S?(\d+)\s*(?:-|–|—)\s*S?(\d+)",s,re.I): out.update(range(int(a),int(b)+1))
        s=re.sub(r"S?\d+\s*(?:-|–|—)\s*S?\d+","",s,flags=re.I)
        out.update(int(x) for x in re.findall(r"S?(\d+)",s,re.I))
    return out

def digest(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for c in iter(lambda:f.read(1<<20),b""):h.update(c)
    return h.hexdigest()

def main():
    m=text(MANUSCRIPT); s=text(BASE/"Supplementary_Information.docx")
    expected_f=set(range(1,10)); expected_t=set(range(1,18))
    required=[BASE/"Supplementary_Information.docx",BASE/"Supplementary_Information.pdf",BASE/"Supplementary_Tables.xlsx",BASE/"source_tables"/"Table_S00_Index.tsv"]
    required += [BASE/"source_tables"/f"Table_S{i:02d}.tsv" for i in expected_t]
    required += [BASE/"figures"/f"Supplementary_Figure_S{i}.{e}" for i in expected_f for e in ("png","pdf","svg","tiff")]
    zip_ok={}
    for p in required[:3]:
        if p.suffix in {".docx",".xlsx"}:
            with ZipFile(p) as z: zip_ok[p.name]=z.testzip() is None
    report={
        "status":"PASS",
        "supplement_title_exact": TITLE in s.splitlines()[:5],
        "obsolete_title_absent": OLD_TITLE not in s,
        "manuscript_figure_references":sorted(refs(m,"fig")),
        "manuscript_table_references":sorted(refs(m,"table")),
        "figure_legends":sorted(int(x) for x in re.findall(r"Supplementary Figure S(\d+)",s)),
        "table_legends":sorted(int(x) for x in re.findall(r"Supplementary Table S(\d+)",s)),
        "missing_or_empty_files":[str(p) for p in required if not p.exists() or p.stat().st_size==0],
        "office_zip_integrity":zip_ok,
        "pdf_pages":len(PdfReader(BASE/"Supplementary_Information.pdf").pages),
        "workbook_qa":json.loads((BASE/"qa"/"workbook"/"workbook_qa.json").read_text())["status"],
        "docx_rendered_pages":len(list((BASE/"qa"/"docx").glob("page-*.png"))),
    }
    report["missing_manuscript_figures"]=sorted(expected_f-set(report["manuscript_figure_references"]))
    report["missing_manuscript_tables"]=sorted(expected_t-set(report["manuscript_table_references"]))
    report["missing_figure_legends"]=sorted(expected_f-set(report["figure_legends"]))
    report["missing_table_legends"]=sorted(expected_t-set(report["table_legends"]))
    failures=[not report["supplement_title_exact"],not report["obsolete_title_absent"],report["missing_or_empty_files"],report["missing_manuscript_figures"],report["missing_manuscript_tables"],report["missing_figure_legends"],report["missing_table_legends"],not all(zip_ok.values()),report["pdf_pages"]!=14,report["workbook_qa"]!="PASS",report["docx_rendered_pages"]!=14]
    if any(failures):report["status"]="FAIL"
    (BASE/"qa_report.json").write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding="utf-8")
    files=[]
    for p in sorted(BASE.rglob("*")):
        if p.is_file() and "qa" not in p.parts and p.name not in {"release_manifest.tsv"} and not p.name.endswith(".inspect.ndjson"):
            files.append((p.relative_to(BASE),p.stat().st_size,digest(p)))
    (BASE/"release_manifest.tsv").write_text("file\tsize_bytes\tsha256\n"+"\n".join(f"{p}\t{n}\t{h}" for p,n,h in files)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2,ensure_ascii=False)); raise SystemExit(0 if report["status"]=="PASS" else 1)

if __name__=="__main__":main()
