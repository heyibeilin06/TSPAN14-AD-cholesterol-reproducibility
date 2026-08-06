# Figure 5

## Panel order and source data

- **A:** `Figure_5_cell_context_atlas.tsv`, `Figure_5_exact_event_cross_context.tsv`
- **B:** `Figure_5_single_nucleus_eqtl.tsv`
- **C:** `Figure_5_disease_state_rna.tsv`
- **D:** `Figure_5_ec2_structure.tsv`, `Figure_5_transcript_events.tsv`, `Figure_5_structure_summary.tsv`
- **E:** `Figure_5_structure_summary.tsv`

## Code order

1. `01_prepare_panel_source_data.R` records how final panel source tables were assembled.
2. `02_plot_panels_A_B_C_D_E.R` draws panels in manuscript order and assembles the figure.
3. `03_quality_control.R` performs figure-specific quality control.

## Reproduction command

```bash
Rscript figures/Figure_05/code/02_plot_panels_A_B_C_D_E.R --project-root . --source-dir figures/Figure_05/data --output-dir reproduced/Figure_05
```

The panel-prefixed files provide explicit public lineage. Compatibility copies retain the filenames expected by the plotting script.
