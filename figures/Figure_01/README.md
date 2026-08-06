# Figure 1

## Panel order and source data

- **A:** `Figure_1_ldsc_apoe_conditioning.tsv`
- **B:** `Figure_1_regional_screen.tsv`
- **C:** `Figure_1_evidence_matrix.tsv`
- **D:** `Figure_1_variant_fingerprint.tsv`

## Code order

1. `01_prepare_panel_source_data.R` records how final panel source tables were assembled.
2. `02_plot_panels_A_B_C_D.R` draws panels in manuscript order and assembles the figure.
3. `03_quality_control.R` performs figure-specific quality control.

## Reproduction command

```bash
Rscript figures/Figure_01/code/02_plot_panels_A_B_C_D.R --project-root . --source-dir figures/Figure_01/data --output-dir reproduced/Figure_01
```

The panel-prefixed files provide explicit public lineage. Compatibility copies retain the filenames expected by the plotting script.
