"""Condition APOE-region summary statistics with a large UKB-European LD matrix.

This workflow replaces physical-window deletion with summary-statistic
conditional residualization inside the LDetect block containing APOE. It then
re-runs genome-wide LDSC using the residualized chr19 statistics.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import os
import re
import subprocess
import sys
from pathlib import Path
from statistics import NormalDist

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
DATA = Path(os.environ.get("SLM_DATA_ROOT", ROOT / "data" / "raw"))
LDSC_ROOT = DATA / "p0_ldsc"
LDSC_CODE = LDSC_ROOT / "cbii_ldsc"
LDSC_DEPS = DATA / "tools" / "ldsc_py312"
LDSC_REFERENCE = LDSC_ROOT / "eur_w_ld_chr"
MUNGED = LDSC_ROOT / "apoe_window_reanalysis" / "munged"
WORK = LDSC_ROOT / "apoe_conditional_reanalysis"
OUT = ROOT / "outputs" / "mentor_revision" / "apoe_conditional"

MAGENPY_PACKAGES = DATA / "tools" / "python_packages"
UKB_LD_ARCHIVE = DATA / "reference" / "ukb_eur_ld" / "EUR.tar.gz"
UKB_LD_CHR19 = DATA / "reference" / "ukb_eur_ld" / "EUR" / "chr_19"
LDSC_PYTHON = Path(os.environ.get("LDSC_PYTHON", sys.executable))

TRAITS = ("AD", "HDL", "LDL", "TG", "TC", "nonHDL")
LIPIDS = TRAITS[1:]
MODELS = ("own_lead", "pair_union_leads")
RAW_FILES = {
    "AD": DATA / "processed" / "ad_bellenguez_2022_harmonized.tsv.gz",
    "HDL": DATA / "processed" / "glgc_hdl_eur_harmonized.tsv.gz",
    "LDL": DATA / "processed" / "glgc_ldl_eur_harmonized.tsv.gz",
    "TG": DATA / "processed" / "glgc_tg_eur_harmonized.tsv.gz",
    "TC": DATA / "processed" / "glgc_tc_eur_harmonized.tsv.gz",
    "nonHDL": DATA / "processed" / "glgc_nonhdl_eur_harmonized.tsv.gz",
}
GWAS_P_THRESHOLD = 5e-8
MIN_CONDITIONAL_VARIANCE = 0.10
MAX_SIGNALS = 20


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def allele_sign(gwas_a1: str, gwas_a2: str, ld_a1: str, ld_a2: str) -> int | None:
    gwas = (gwas_a1.upper(), gwas_a2.upper())
    ld = (ld_a1.upper(), ld_a2.upper())
    if gwas == ld:
        return 1
    if gwas == ld[::-1]:
        return -1
    return None


def conditional_z(
    z: np.ndarray,
    n: np.ndarray,
    ld: np.ndarray,
    conditioning_indices: list[int],
    ridge: float = 1e-6,
) -> tuple[np.ndarray, np.ndarray]:
    """Residualize marginal standardized effects against conditioning SNPs."""
    if not conditioning_indices:
        return z.copy(), np.ones_like(z, dtype=float)
    cond = np.asarray(conditioning_indices, dtype=int)
    r_xy = z / np.sqrt(n)
    r_cc = ld[np.ix_(cond, cond)].astype(float)
    r_cc = r_cc + np.eye(len(cond)) * ridge
    solved_effect = np.linalg.solve(r_cc, r_xy[cond])
    solved_ld = np.linalg.solve(r_cc, ld[cond, :])
    residual_effect = r_xy - ld[:, cond] @ solved_effect
    variance_fraction = 1.0 - np.einsum("ij,ji->i", ld[:, cond], solved_ld)
    variance_fraction = np.clip(variance_fraction, 1e-8, 1.0)
    result = np.sqrt(n) * residual_effect / np.sqrt(variance_fraction)
    result[cond] = 0.0
    return result, variance_fraction


def joint_selected_z(z: np.ndarray, n: np.ndarray, ld: np.ndarray, selected: list[int]) -> np.ndarray:
    if not selected:
        return np.array([], dtype=float)
    idx = np.asarray(selected, dtype=int)
    r_ss = ld[np.ix_(idx, idx)].astype(float) + np.eye(len(idx)) * 1e-6
    inverse = np.linalg.inv(r_ss)
    joint_effect = inverse @ (z[idx] / np.sqrt(n[idx]))
    return joint_effect * np.sqrt(n[idx]) / np.sqrt(np.diag(inverse))


def z_threshold(p_threshold: float) -> float:
    return NormalDist().inv_cdf(1.0 - p_threshold / 2.0)


def stepwise_select(
    z: np.ndarray,
    n: np.ndarray,
    ld: np.ndarray,
    p_threshold: float = GWAS_P_THRESHOLD,
    min_conditional_variance: float = MIN_CONDITIONAL_VARIANCE,
    max_signals: int = MAX_SIGNALS,
) -> list[int]:
    """Approximate COJO forward selection with backward joint-significance checks."""
    cutoff = z_threshold(p_threshold)
    selected: list[int] = []
    while len(selected) < max_signals:
        current_z, variance_fraction = conditional_z(z, n, ld, selected)
        eligible = [
            i for i in range(len(z))
            if i not in selected and variance_fraction[i] >= min_conditional_variance
        ]
        if not eligible:
            break
        best = max(eligible, key=lambda i: abs(current_z[i]))
        if abs(current_z[best]) < cutoff:
            break
        selected.append(best)
        while len(selected) > 1:
            joint_z = joint_selected_z(z, n, ld, selected)
            weakest = int(np.argmin(np.abs(joint_z)))
            if abs(joint_z[weakest]) >= cutoff:
                break
            selected.pop(weakest)
    return selected


def p_from_z(z: float) -> float:
    return math.erfc(abs(z) / math.sqrt(2.0))


def load_ld_block() -> dict[str, object]:
    sys.path.insert(0, str(MAGENPY_PACKAGES))
    from magenpy import LDMatrix  # type: ignore

    ld_store = LDMatrix.from_path(str(UKB_LD_CHR19))
    blocks = ld_store.zarr_group.attrs["Estimator properties"]["LD blocks"]
    apoe_position_grch37 = 45_411_941
    block = next(pair for pair in blocks if pair[0] <= apoe_position_grch37 < pair[1])
    bp = np.asarray(ld_store.bp_position)
    global_indices = np.where((bp >= block[0]) & (bp < block[1]))[0]
    start, end = int(global_indices[0]), int(global_indices[-1]) + 1
    matrix = ld_store.load_data(
        start_row=start,
        end_row=end,
        dtype="float64",
        return_square=True,
        return_symmetric=True,
        return_as_csr=True,
    ).toarray()
    return {
        "matrix": matrix,
        "snps": np.asarray(ld_store.snps)[global_indices].astype(str),
        "a1": np.asarray(ld_store.a1)[global_indices].astype(str),
        "a2": np.asarray(ld_store.a2)[global_indices].astype(str),
        "maf": np.asarray(ld_store.maf)[global_indices].astype(float),
        "bp": bp[global_indices].astype(int),
        "block_start": int(block[0]),
        "block_end": int(block[1]),
        "sample_size": int(ld_store.sample_size),
        "genome_build": str(ld_store.genome_build),
    }


def ensure_raw_block_caches(ld_block: dict[str, object]) -> None:
    block_snps = set(np.asarray(ld_block["snps"]).tolist())
    cache_dir = WORK / "raw_block"
    cache_dir.mkdir(parents=True, exist_ok=True)
    manifest_rows: list[dict[str, object]] = []
    for trait, source in RAW_FILES.items():
        output = cache_dir / f"{trait}.tsv.gz"
        if not output.exists():
            rows: list[dict[str, object]] = []
            with gzip.open(source, "rt", encoding="utf-8", errors="replace") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                for row in reader:
                    if row.get("SNP") not in block_snps:
                        continue
                    try:
                        beta = float(row["BETA"])
                        se = float(row["SE"])
                        n = float(row["N"])
                    except (KeyError, TypeError, ValueError):
                        continue
                    if se <= 0 or not all((row.get("A1"), row.get("A2"))):
                        continue
                    rows.append({
                        "SNP": row["SNP"], "A1": row["A1"], "A2": row["A2"],
                        "Z": beta / se, "N": n,
                    })
            with gzip.open(output, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=("SNP", "A1", "A2", "Z", "N"), delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)
        with gzip.open(output, "rt", encoding="utf-8") as handle:
            retained = sum(1 for _ in handle) - 1
        manifest_rows.append({
            "trait": trait, "source": str(source), "source_size_bytes": source.stat().st_size,
            "cache": str(output), "retained_ld_block_rows": retained,
            "selection": "rsID intersection with UKB EUR LDetect APOE block",
        })
    write_tsv(OUT / "raw_block_cache_manifest.tsv", manifest_rows)


def load_trait_block(trait: str, ld_block: dict[str, object]) -> dict[str, object]:
    snps = np.asarray(ld_block["snps"])
    lookup = {snp: i for i, snp in enumerate(snps)}
    records: list[dict[str, object]] = []
    blank = mismatch = 0
    with gzip.open(WORK / "raw_block" / f"{trait}.tsv.gz", "rt", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            snp = row["SNP"]
            if snp not in lookup:
                continue
            if not row.get("Z") or not row.get("N") or not row.get("A1") or not row.get("A2"):
                blank += 1
                continue
            ld_i = lookup[snp]
            sign = allele_sign(row["A1"], row["A2"], str(ld_block["a1"][ld_i]), str(ld_block["a2"][ld_i]))
            if sign is None:
                mismatch += 1
                continue
            records.append({
                "snp": snp,
                "ld_index": ld_i,
                "gwas_a1": row["A1"],
                "gwas_a2": row["A2"],
                "sign": sign,
                "z_ld": float(row["Z"]) * sign,
                "n": float(row["N"]),
            })
    indices = np.asarray([int(row["ld_index"]) for row in records], dtype=int)
    return {
        "records": records,
        "z": np.asarray([float(row["z_ld"]) for row in records]),
        "n": np.asarray([float(row["n"]) for row in records]),
        "matrix": np.asarray(ld_block["matrix"])[np.ix_(indices, indices)],
        "blank": blank,
        "mismatch": mismatch,
    }


def marginal_lead(dataset: dict[str, object]) -> str:
    best = int(np.argmax(np.abs(np.asarray(dataset["z"]))))
    return str(dataset["records"][best]["snp"])


def map_anchors_to_dataset(
    dataset: dict[str, object],
    ld_block: dict[str, object],
    anchors: list[tuple[str, str]],
    minimum_proxy_r2: float = 0.80,
) -> tuple[list[int], list[dict[str, object]]]:
    records = dataset["records"]
    record_by_ld = {int(row["ld_index"]): i for i, row in enumerate(records)}
    snps = np.asarray(ld_block["snps"])
    full_ld = np.asarray(ld_block["matrix"])
    selected: list[int] = []
    rows: list[dict[str, object]] = []
    for anchor_trait, anchor in anchors:
        anchor_match = np.where(snps == anchor)[0]
        if not len(anchor_match):
            rows.append({"anchor_trait": anchor_trait, "anchor_snp": anchor, "selected_snp": "", "r2_to_anchor": "", "status": "anchor_not_in_ld"})
            continue
        anchor_ld = int(anchor_match[0])
        candidates = list(record_by_ld)
        best_ld = max(candidates, key=lambda i: abs(full_ld[anchor_ld, i]))
        r2 = float(full_ld[anchor_ld, best_ld] ** 2)
        local_i = record_by_ld[best_ld]
        if r2 < minimum_proxy_r2:
            rows.append({"anchor_trait": anchor_trait, "anchor_snp": anchor, "selected_snp": records[local_i]["snp"], "r2_to_anchor": r2, "status": "below_r2_threshold"})
            continue
        if any(dataset["matrix"][local_i, prior] ** 2 >= 0.90 for prior in selected):
            rows.append({"anchor_trait": anchor_trait, "anchor_snp": anchor, "selected_snp": records[local_i]["snp"], "r2_to_anchor": r2, "status": "collinear_with_selected"})
            continue
        selected.append(local_i)
        rows.append({
            "anchor_trait": anchor_trait,
            "anchor_snp": anchor,
            "selected_snp": records[local_i]["snp"],
            "r2_to_anchor": r2,
            "status": "exact" if records[local_i]["snp"] == anchor else "proxy",
        })
    return selected, rows


def write_conditioned_sumstats(
    trait: str,
    model: str,
    comparison: str,
    dataset: dict[str, object],
    conditioned_z_ld: np.ndarray,
    variance_fraction: np.ndarray,
    selected: list[int],
    block_snps: set[str],
) -> tuple[Path, int, int]:
    output = WORK / "conditional" / model / comparison / f"{trait}.sumstats.gz"
    output.parent.mkdir(parents=True, exist_ok=True)
    selected_set = set(selected)
    replacement = {
        str(row["snp"]): (
            str(row["gwas_a1"]), str(row["gwas_a2"]),
            float(conditioned_z_ld[i]) * int(row["sign"]), float(row["n"]),
        )
        for i, row in enumerate(dataset["records"])
        if i in selected_set or variance_fraction[i] >= MIN_CONDITIONAL_VARIANCE
    }
    retained = removed = 0
    with gzip.open(MUNGED / f"{trait}.sumstats.gz", "rt", encoding="utf-8") as source, gzip.open(
        output, "wt", encoding="utf-8"
    ) as target:
        header = source.readline()
        target.write(header)
        fields = header.rstrip("\n").split("\t")
        snp_i, a1_i, a2_i = fields.index("SNP"), fields.index("A1"), fields.index("A2")
        z_i, n_i = fields.index("Z"), fields.index("N")
        for line in source:
            values = line.rstrip("\n").split("\t")
            snp = values[snp_i]
            if snp not in block_snps:
                target.write(line)
                retained += 1
            elif snp in replacement and len(values) == len(fields):
                a1, a2, z_value, n_value = replacement[snp]
                values[a1_i], values[a2_i] = a1, a2
                values[z_i], values[n_i] = f"{z_value:.8g}", f"{n_value:.8g}"
                target.write("\t".join(values) + "\n")
                retained += 1
            else:
                removed += 1
    return output, retained, removed


def run_command(command: list[str], log_path: Path) -> tuple[int, str]:
    env = os.environ.copy()
    env["PYTHONPATH"] = str(LDSC_DEPS)
    complete = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=env)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(complete.stdout, encoding="utf-8", errors="replace")
    return complete.returncode, complete.stdout


def parse_rg(log: str) -> tuple[str, str, str]:
    estimate = re.search(r"^Genetic Correlation:\s*([^\s]+)\s*\(([^)]+)\)", log, flags=re.MULTILINE)
    p_value = re.search(r"^P:\s*([^\s]+)", log, flags=re.MULTILINE)
    if estimate and p_value:
        return estimate.group(1), estimate.group(2), p_value.group(1)
    return "", "", ""


def prepare_conditioned_files() -> None:
    ld_block = load_ld_block()
    ensure_raw_block_caches(ld_block)
    block_snps = set(np.asarray(ld_block["snps"]).tolist())
    manifest = [{
        "resource": "UK Biobank EUR LD matrix",
        "doi": "10.5281/zenodo.14614207",
        "ancestry": "EUR",
        "reference_n": ld_block["sample_size"],
        "genome_build": ld_block["genome_build"],
        "chromosome": 19,
        "ld_block_start": ld_block["block_start"],
        "ld_block_end": ld_block["block_end"],
        "ld_block_variants": len(block_snps),
        "archive_md5": "41826edf74f9cc14b3e97024119ad2e6",
        "archive": str(UKB_LD_ARCHIVE),
    }]
    write_tsv(OUT / "ld_reference_manifest.tsv", manifest)

    datasets = {trait: load_trait_block(trait, ld_block) for trait in TRAITS}
    leads = {trait: marginal_lead(datasets[trait]) for trait in TRAITS}
    signal_rows: list[dict[str, object]] = []
    qc_rows: list[dict[str, object]] = []
    residual_rows: list[dict[str, object]] = []
    for lipid in LIPIDS:
        comparison = f"AD_{lipid}"
        for model in MODELS:
            for trait in ("AD", lipid):
                dataset = datasets[trait]
                if model == "own_lead":
                    anchors = [(trait, leads[trait])]
                else:
                    anchors = [("AD", leads["AD"]), (lipid, leads[lipid])]
                selected, mapping_rows = map_anchors_to_dataset(dataset, ld_block, anchors)
                for row in mapping_rows:
                    signal_rows.append({
                        "comparison": f"AD-{lipid}",
                        "conditioned_trait": trait,
                        "model": model,
                        **row,
                    })
                z = np.asarray(dataset["z"])
                n = np.asarray(dataset["n"])
                matrix = np.asarray(dataset["matrix"])
                conditioned, variance_fraction = conditional_z(z, n, matrix, selected)
                output, retained, removed = write_conditioned_sumstats(
                    trait, model, comparison, dataset, conditioned, variance_fraction, selected, block_snps
                )
                stable = variance_fraction >= MIN_CONDITIONAL_VARIANCE
                selected_mask = np.zeros(len(z), dtype=bool)
                selected_mask[selected] = True
                analyzed = stable & ~selected_mask
                residual_rows.append({
                    "comparison": f"AD-{lipid}",
                    "trait": trait,
                    "model": model,
                    "conditioning_signals": len(selected),
                    "max_abs_z_before": float(np.max(np.abs(z))),
                    "max_abs_z_after_stable_variants": float(np.max(np.abs(conditioned[analyzed]))) if np.any(analyzed) else "",
                    "significant_before": int(np.sum(np.abs(z) >= z_threshold(GWAS_P_THRESHOLD))),
                    "significant_after_stable_variants": int(np.sum(np.abs(conditioned[analyzed]) >= z_threshold(GWAS_P_THRESHOLD))),
                    "low_conditional_variance_removed": int(np.sum(~stable & ~selected_mask)),
                    "conditioned_sumstats": str(output),
                })
                qc_rows.append({
                    "comparison": f"AD-{lipid}",
                    "trait": trait,
                    "model": model,
                    "valid_aligned_block_variants": len(dataset["records"]),
                    "blank_block_rows_removed": dataset["blank"],
                    "allele_mismatch_rows_removed": dataset["mismatch"],
                    "genomewide_rows_retained": retained,
                    "block_rows_removed_from_ldsc": removed,
                })
    write_tsv(OUT / "conditioning_signals.tsv", signal_rows)
    write_tsv(OUT / "block_harmonization_qc.tsv", qc_rows)
    write_tsv(OUT / "conditional_residual_qc.tsv", residual_rows)


def run_ldsc() -> None:
    if not LDSC_PYTHON.exists():
        raise FileNotFoundError(f"LDSC Python not found: {LDSC_PYTHON}")
    rows: list[dict[str, object]] = []
    for model in MODELS:
        for lipid in LIPIDS:
            comparison = f"AD_{lipid}"
            ad = WORK / "conditional" / model / comparison / "AD.sumstats.gz"
            other = WORK / "conditional" / model / comparison / f"{lipid}.sumstats.gz"
            prefix = WORK / "rg" / model / f"AD_{lipid}"
            log_path = WORK / "logs" / f"rg_{model}_AD_{lipid}.log"
            prefix.parent.mkdir(parents=True, exist_ok=True)
            command = [
                str(LDSC_PYTHON), str(LDSC_CODE / "ldsc.py"),
                "--rg", f"{ad},{other}",
                "--ref-ld-chr", str(LDSC_REFERENCE) + os.sep,
                "--w-ld-chr", str(LDSC_REFERENCE) + os.sep,
                "--out", str(prefix),
            ]
            code, log = run_command(command, log_path)
            rg, se, p = parse_rg(log)
            rows.append({
                "model": model,
                "comparison": f"AD-{lipid}",
                "status": "completed" if code == 0 and rg else "failed",
                "rg": rg,
                "se": se,
                "p": p,
                "log": str(log_path),
            })
    write_tsv(OUT / "conditional_ldsc_results.tsv", rows)
    write_integrated_summary(rows)


def write_integrated_summary(conditional_rows: list[dict[str, object]]) -> None:
    rows: list[dict[str, object]] = []
    window_results = ROOT / "outputs" / "p0_reanalysis" / "p0_apoe_window_ldsc_results.tsv"
    if window_results.exists():
        with window_results.open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row["comparison"] == "AD-HDL" and row["window"] in {"baseline", "w1Mb", "w2Mb", "w5Mb"}:
                    rows.append({
                        "analysis_class": "baseline" if row["window"] == "baseline" else "physical-window sensitivity",
                        "model": row["window"], "comparison": row["comparison"],
                        "rg": row["rg"], "se": row["se"], "p": row["p"],
                        "scope": "unmodified genome-wide statistics" if row["window"] == "baseline" else "APOE lead-anchor window exclusion",
                    })
    for row in conditional_rows:
        if row["comparison"] == "AD-HDL" and row["status"] == "completed":
            rows.append({
                "analysis_class": "LD-based summary conditional analysis",
                "model": row["model"], "comparison": row["comparison"],
                "rg": row["rg"], "se": row["se"], "p": row["p"],
                "scope": "GRCh37 LDetect APOE block conditioned on pre-specified observed lead signal(s)",
            })
    write_tsv(OUT / "apoe_sensitivity_summary.tsv", rows)


def write_method_note() -> None:
    text = "# APOE summary-statistic conditional analysis\n\n"
    text += "Two pre-specified lead-signal models were evaluated within the GRCh37 LDetect block chr19:44,744,108-46,102,697. The own-lead model conditions each GWAS on its strongest available association in the block. The pair-union model conditions both members of each AD-lipid pair on the union of the AD and lipid lead signals. If an anchor is absent from a GWAS input, an r2>=0.80 proxy from the same UKB matrix is used; r2>=0.90 redundant anchors are represented once. This design corresponds to fixed-signal summary conditional analysis rather than automated signal discovery.\n\n"
    text += "Marginal standardized effects were aligned to UK Biobank LD alleles, residualized against the pre-specified signals and transformed to conditional Z statistics. Conditioning variants were assigned Z=0. Variants retaining less than 10% conditional genotype variance were excluded as unstable, matching the collinearity principle used by COJO. The residualized block replaced the original block in each genome-wide LDSC input.\n\n"
    text += "The APOE coding variants rs429358 and rs7412 are absent from the 1.4-million-variant UKB LD release and from the AD LDSC input. The analysis therefore conditions the observed extended APOE-region association signals and is not described as direct APOE-genotype conditioning.\n\n"
    text += "The LD reference comprises 362,063 unrelated European-ancestry UK Biobank participants (Zenodo DOI 10.5281/zenodo.14614207). Physical-window exclusions remain separate sensitivity analyses and are not interpreted as complete APOE removal.\n"
    (OUT / "METHOD_NOTE.md").write_text(text, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("prepare", "ldsc", "all"), default="all", nargs="?")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    if args.stage in ("prepare", "all"):
        prepare_conditioned_files()
        write_method_note()
    if args.stage in ("ldsc", "all"):
        run_ldsc()


if __name__ == "__main__":
    main()
