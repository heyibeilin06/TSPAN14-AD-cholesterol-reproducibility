# Figure 3

## Panel order and source data

- **A:** `Figure_3_exact_replication_matrix.tsv`
- **B:** `Figure_3_ba24_donor_ier.tsv`, `Figure_3_depth_sensitivity_delta_ier.tsv`
- **C:** `Figure_3_brain_junction_counts.tsv`, `Figure_3_brain_cousage_summary.tsv`
- **D:** `Figure_3_primary_exons.tsv`, `Figure_3_splice_events.tsv`, `Figure_3_alignment_metric_audit.tsv`

## Code order

1. `01_prepare_panel_source_data.R` records how final panel source tables were assembled.
2. `02_plot_panels_A_B_C_D.R` draws panels in manuscript order and assembles the figure.
3. `03_quality_control.R` performs figure-specific quality control.

## Reproduction command

```bash
Rscript figures/Figure_03/code/02_plot_panels_A_B_C_D.R --project-root . --source-dir figures/Figure_03/data --output-dir reproduced/Figure_03
```

The panel-prefixed files provide explicit public lineage. Compatibility copies retain the filenames expected by the plotting script.
