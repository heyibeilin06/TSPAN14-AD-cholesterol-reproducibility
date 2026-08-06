# Figure 4

## Panel order and source data

- **A:** `Figure_4_instrument_effect_atlas.tsv`
- **B:** `Figure_4_ld_aware_cis_mr.tsv`
- **C:** `Figure_4_cis_mr_diagnostics.tsv`
- **D:** `Figure_4_genomewide_lipid_to_ad.tsv`
- **E:** `Figure_4_global_joint_mvmr.tsv`, `Figure_4_global_joint_mvmr_strength.tsv`
- **F:** `Figure_4_pc_gmm_dimension_sensitivity.tsv`

## Code order

1. `01_prepare_panel_source_data.R` records how final panel source tables were assembled.
2. `02_plot_panels_A_B_C_D_E_F.R` draws panels in manuscript order and assembles the figure.
3. `03_quality_control.R` performs figure-specific quality control.

## Reproduction command

```bash
Rscript figures/Figure_04/code/02_plot_panels_A_B_C_D_E_F.R --project-root . --source-dir figures/Figure_04/data --output-dir reproduced/Figure_04
```

The panel-prefixed files provide explicit public lineage. Compatibility copies retain the filenames expected by the plotting script.
