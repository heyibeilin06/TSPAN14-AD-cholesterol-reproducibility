#!/usr/bin/env node
import fs from "node:fs/promises";
import path from "node:path";
import { pathToFileURL } from "node:url";

const artifactToolModule = process.env.ARTIFACT_TOOL_MODULE
  ? pathToFileURL(process.env.ARTIFACT_TOOL_MODULE).href
  : "@oai/artifact-tool";
const { SpreadsheetFile, Workbook } = await import(artifactToolModule);

const root = path.resolve(process.argv[2] ?? ".");
const base = path.join(root, "outputs", "supplement_v19");
const sourceDir = path.join(base, "source_tables");
const qaDir = path.join(base, "qa", "workbook");
const title = "A non-APOE Alzheimer disease-cholesterol locus converges on TSPAN14 splice choice";

function colName(index) {
  let n = index + 1, out = "";
  while (n) { n--; out = String.fromCharCode(65 + (n % 26)) + out; n = Math.floor(n / 26); }
  return out;
}

function parseCell(value) {
  const s = value.trim();
  if (!s || ["NA", "NaN", "nan", "None"].includes(s)) return null;
  if (/^(TRUE|FALSE)$/i.test(s)) return /^TRUE$/i.test(s);
  if (/^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?$/i.test(s)) return Number(s);
  return /^[=+@]/.test(s) ? `'${s}` : s;
}

function parseTsv(text) {
  return text.replace(/^\uFEFF/, "").trimEnd().split(/\r?\n/).map(line => line.split("\t").map(parseCell));
}

const indexRows = parseTsv(await fs.readFile(path.join(sourceDir, "Table_S00_Index.tsv"), "utf8"));
const indexHeader = indexRows[0];
const records = indexRows.slice(1).map(row => Object.fromEntries(indexHeader.map((key, i) => [key, row[i]])));
const workbook = Workbook.create();
const previewRanges = new Map();
const readme = workbook.worksheets.add("README");
readme.showGridLines = false;
readme.getRange("A1:D1").merge();
readme.getRange("A1").values = [["Supplementary Tables"]];
readme.getRange("A2:D2").merge();
readme.getRange("A2").values = [[title]];
readme.getRange("A4:D4").values = [["Table", "Content", "Rows", "Source file"]];
readme.getRange(`A5:D${records.length + 4}`).values = records.map(r => [`Table S${r.number}`, r.title, r.rows, r.file]);
readme.getRange(`A${records.length + 6}:D${records.length + 7}`).merge();
readme.getRange(`A${records.length + 6}`).values = [["Each sheet is rebuilt from the current audited analysis outputs. Blank cells denote non-applicable fields; numerical values retain source precision."]];
readme.getRange("A1:D1").format = { font: { name: "Arial", size: 16, bold: true, color: "#123B5D" }, rowHeight: 28 };
readme.getRange("A2:D2").format = { font: { name: "Arial", size: 10, italic: true, color: "#465A6B" }, wrapText: true };
readme.getRange("A4:D4").format = { font: { name: "Arial", size: 10, bold: true }, borders: { top: { style: "medium", color: "#111111" }, bottom: { style: "thin", color: "#111111" } } };
readme.getRange(`A${records.length + 4}:D${records.length + 4}`).format.borders = { bottom: { style: "medium", color: "#111111" } };
readme.getRange(`A5:D${records.length + 7}`).format = { font: { name: "Arial", size: 9 }, wrapText: true, verticalAlignment: "top" };
readme.getRange("A:A").format.columnWidth = 15;
readme.getRange("B:B").format.columnWidth = 62;
readme.getRange("C:C").format.columnWidth = 10;
readme.getRange("D:D").format.columnWidth = 24;
readme.freezePanes.freezeRows(4);
previewRanges.set("README", `A1:D${records.length + 7}`);

for (const rec of records) {
  const number = Number(rec.number);
  const matrix = parseTsv(await fs.readFile(path.join(sourceDir, rec.file), "utf8"));
  const cols = matrix[0].length;
  const rows = matrix.length;
  const last = colName(cols - 1);
  const sheet = workbook.worksheets.add(`Table S${number}`);
  sheet.showGridLines = false;
  sheet.getRange(`A1:${last}1`).merge();
  sheet.getRange("A1").values = [[`Supplementary Table S${number} | ${rec.title}`]];
  sheet.getRange(`A2:${last}2`).merge();
  sheet.getRange("A2").values = [["Current-analysis source data; values are reported at source precision."]];
  sheet.getRange(`A4:${last}${rows + 3}`).values = matrix;
  sheet.getRange(`A1:${last}1`).format = { font: { name: "Arial", size: 12, bold: true, color: "#123B5D" }, wrapText: true, rowHeight: 28 };
  sheet.getRange(`A2:${last}2`).format = { font: { name: "Arial", size: 8, italic: true, color: "#5F6B73" }, wrapText: true };
  sheet.getRange(`A4:${last}4`).format = { font: { name: "Arial", size: 9, bold: true }, wrapText: true, verticalAlignment: "bottom", borders: { top: { style: "medium", color: "#111111" }, bottom: { style: "thin", color: "#111111" } } };
  if (rows > 1) {
    sheet.getRange(`A5:${last}${rows + 3}`).format = { font: { name: "Arial", size: 8 }, wrapText: true, verticalAlignment: "top" };
    sheet.getRange(`A${rows + 3}:${last}${rows + 3}`).format.borders = { bottom: { style: "medium", color: "#111111" } };
  }
  for (let c = 0; c < cols; c++) {
    const values = matrix.slice(0, Math.min(rows, 150)).map(row => String(row[c] ?? ""));
    const maxLen = Math.max(8, ...values.map(v => v.length));
    const header = String(matrix[0][c] ?? "").toLowerCase();
    let width = Math.min(30, Math.max(10, maxLen * 0.82 + 2));
    if (/interpret|note|reason|description|url|scope|boundary/.test(header)) width = Math.min(42, Math.max(width, 24));
    sheet.getRange(`${colName(c)}:${colName(c)}`).format.columnWidth = width;
  }
  sheet.freezePanes.freezeRows(4);
  previewRanges.set(`Table S${number}`, `A1:${last}${Math.min(rows + 3, 60)}`);
}

await fs.mkdir(qaDir, { recursive: true });
const output = await SpreadsheetFile.exportXlsx(workbook);
const workbookPath = path.join(base, "Supplementary_Tables.xlsx");
await output.save(workbookPath);

const previews = [];
for (const name of ["README", ...records.map(r => `Table S${r.number}`)]) {
  const image = await workbook.render({ sheetName: name, range: previewRanges.get(name), scale: 0.65, format: "png" });
  const target = path.join(qaDir, `${name.replaceAll(" ", "_")}.png`);
  await fs.writeFile(target, new Uint8Array(await image.arrayBuffer()));
  previews.push({ sheet: name, bytes: (await fs.stat(target)).size });
}
const errors = await workbook.inspect({ kind: "match", searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A", options: { useRegex: true, maxResults: 200 }, maxChars: 8000 });
const report = { status: (errors.ndjson ?? "").includes("#REF!") ? "FAIL" : "PASS", sheets: previews.length, previews, formulaErrors: errors.ndjson ?? "" };
await fs.writeFile(path.join(qaDir, "workbook_qa.json"), JSON.stringify(report, null, 2));
console.log(JSON.stringify({ workbookPath, status: report.status, sheets: previews.length }));
