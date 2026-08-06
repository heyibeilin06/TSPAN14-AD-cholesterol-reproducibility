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
