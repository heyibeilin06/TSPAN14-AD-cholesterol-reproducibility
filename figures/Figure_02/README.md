# Figure 2

## Panel order and source data

- **A:** `Figure_2_locus_tracks_grch38.tsv`, `Figure_2_annotation_markers_grch38.tsv`
- **B:** `Figure_2_gencode_v38_transcripts.tsv`, `Figure_2_gencode_v38_exons.tsv`, `Figure_2_regulatory_elements.tsv`, `Figure_2_annotation_markers_grch38.tsv`
- **C:** `Figure_2_colocalization_scatter.tsv`, `Figure_2_exact_event_coloc.tsv`
- **D:** `Figure_2_variant_annotation_matrix.tsv`

## Code order

1. `01_prepare_panel_source_data.R` records how final panel source tables were assembled.
2. `02_plot_panels_A_B_C_D.R` draws panels in manuscript order and assembles the figure.
3. `03_quality_control.R` performs figure-specific quality control.

## Reproduction command

```bash
Rscript figures/Figure_02/code/02_plot_panels_A_B_C_D.R --project-root . --source-dir figures/Figure_02/data --output-dir reproduced/Figure_02
```

The panel-prefixed files provide explicit public lineage. Compatibility copies retain the filenames expected by the plotting script.
