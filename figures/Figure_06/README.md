# Figure 6

## Panel order and source data

- **MODEL:** `Figure_5_exact_event_cross_context.tsv`, `Figure_6_structure_summary.tsv`, `Figure_6_edges.tsv`, `Figure_5_exact_event_triangulation.tsv`

## Code order

1. `01_prepare_panel_source_data.R` records how final panel source tables were assembled.
2. `02_plot_panels_MODEL.R` draws panels in manuscript order and assembles the figure.
3. `03_quality_control.R` performs figure-specific quality control.

## Reproduction command

```bash
Rscript figures/Figure_06/code/02_plot_panels_MODEL.R --project-root . --source-dir figures/Figure_06/data --output-dir reproduced/Figure_06
```

The panel-prefixed files provide explicit public lineage. Compatibility copies retain the filenames expected by the plotting script.
