# External data layout

Raw GWAS, molecular-QTL, genotype, junction-count and controlled-access files are not redistributed. Obtain them from the resources and releases listed in `config/data_sources.tsv` and set the environment variable `SLM_DATA_ROOT` to their local root.

The scripts use descriptive subdirectories under that root, including `processed/`, `p0_ldsc/`, `p0_ld_1000g_eur/`, `reference/` and `tools/`. Dataset-specific accessions, coordinate systems and access conditions are recorded in the data-source registry. `TSPAN14_LD_FILE` can override the default local LD-matrix path for the GTEx exact-event workflow.

The `figures/*/data` and `tables/*/source_data` directories contain only derived, manuscript-level source tables that may be redistributed.
