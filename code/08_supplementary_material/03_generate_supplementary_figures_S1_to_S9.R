#!/usr/bin/env Rscript
# Portable orchestration entry point for all supplementary figures.

args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(ifelse(length(args), args[1], "."), winslash = "/")
script_dir <- file.path(root, "code", "08_supplementary_material")
rscript <- file.path(R.home("bin"), "Rscript")
scripts <- c(
  "09_regenerate_supplementary_figures_s1_s5_s6_s8.R",
  "08_regenerate_supplementary_figure_s2.R",
  "06_regenerate_supplementary_figure_s3.R",
  "05_update_supplementary_figures_S4_S7.R",
  "07_regenerate_supplementary_figure_s9.R"
)

for (script in scripts) {
  status <- system2(rscript, c(file.path(script_dir, script), root))
  if (!identical(status, 0L)) {
    stop("Supplementary figure build failed in ", script, " (exit status ", status, ")")
  }
}

cat("Generated Supplementary Figures S1-S9 from package-local data.\n")
