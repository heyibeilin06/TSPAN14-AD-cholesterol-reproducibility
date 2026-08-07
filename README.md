# TSPAN14 AD-cholesterol reproducibility package

This release accompanies the manuscript *A non-APOE Alzheimer disease-cholesterol locus converges on TSPAN14 splice choice*.

## Scope

The repository reproduces the reported derived-data figures and tables and maps every code-derived manuscript claim to a purpose-named analysis script. Release 1.2.1 synchronizes the Round 10 manuscript, main figures, Supplementary Figures S1-S9 and Supplementary Tables S1-S19. It supports a shared local regulatory configuration and a candidate canonical-versus-cryptic TSPAN14 splice readout. It does not encode a serial lipid-to-splicing-to-AD mediation claim or splice-specific ADAM10/TREM2 causality.

## Contents

- `code/`: upstream analysis scripts organized by scientific purpose;
- `figures/`: panel-ordered R code, derived source tables, legends and publication outputs for Figures 1-6;
- `tables/`: source data for Tables 1-2 and Supplementary Tables S1-S19;
- `config/`: analysis manifest, data-source registry and script migration map;
- `commands/`: cross-platform validation and figure-reproduction commands;
- `audit/`: claim-to-code, panel-lineage, asset-decision and checksum records.

Raw GWAS, genotype, individual-level transcriptomic and controlled-access data are not redistributed. Accessions and URLs are listed in `config/data_sources.tsv`.

## Quick start

```bash
python tests/validate_release.py
bash commands/rebuild_figures.sh
```

Supplementary source tables, the Excel workbook and revised Supplementary Figures S2, S3 and S9 can be rebuilt with the scripts in `code/08_supplementary_material/`. The publication-ready supplementary outputs are stored in `tables/supplementary/` and `figures/Supplementary_Figures/`.

On Windows, use `powershell -ExecutionPolicy Bypass -File commands/rebuild_figures.ps1`. See `commands/README.md` and `data/README.md` for the two-level reproduction model.
