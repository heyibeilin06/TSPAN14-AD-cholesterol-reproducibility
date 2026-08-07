# Reproduction commands

The package has two reproducibility levels.

1. **Figure/table reproduction from redistributed derived data** requires no controlled-access input. Run the validation command first, then rebuild the figures.
2. **Full upstream analysis** requires the public or registered-access resources listed in `config/data_sources.tsv`. Set `SLM_DATA_ROOT` to a local directory arranged as described in `data/README.md`, then run the purpose-specific script recorded in `config/analysis_manifest.tsv`.

Windows PowerShell:

```powershell
python tests/validate_release.py
powershell -ExecutionPolicy Bypass -File commands/rebuild_figures.ps1
```

Linux/macOS:

```bash
python tests/validate_release.py
bash commands/rebuild_figures.sh
```

All plotting is performed in R. Outputs are written to `reproduced/Figure_01` through `reproduced/Figure_06` and do not overwrite the publication files.

The synchronized supplementary release is rebuilt in this order:

```powershell
python code/08_supplementary_material/01_build_supplementary_source_tables.py
node code/08_supplementary_material/02_build_supplementary_workbook.mjs .
Rscript code/08_supplementary_material/08_regenerate_supplementary_figure_s2.R .
Rscript code/08_supplementary_material/06_regenerate_supplementary_figure_s3.R .
Rscript code/08_supplementary_material/07_regenerate_supplementary_figure_s9.R .
python code/08_supplementary_material/04_build_supplementary_information.py
```

`02_build_supplementary_workbook.mjs` uses `@oai/artifact-tool`; install that package in the active Node environment or set `ARTIFACT_TOOL_MODULE` to its module file.
