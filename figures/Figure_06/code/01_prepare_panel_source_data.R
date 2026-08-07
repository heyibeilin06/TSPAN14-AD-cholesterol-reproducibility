#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop("Missing value for ", flag, call. = FALSE)
  args[[i + 1L]]
}

root <- normalizePath(value_after("--project-root", "."), winslash = "/", mustWork = TRUE)
out_dir <- value_after("--output-dir", file.path(root, "figures", "Figure_06", "data"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_checked <- function(path, required) {
  if (!file.exists(path)) stop("Missing audited source: ", path, call. = FALSE)
  x <- fread(path, sep = "\t", na.strings = c("", "NA"))
  missing <- setdiff(required, names(x))
  if (length(missing)) stop("Missing columns in ", basename(path), ": ", paste(missing, collapse = ", "), call. = FALSE)
  x
}

context_source <- read_checked(
  file.path(root, "figures", "Figure_05", "data", "Figure_5_exact_event_cross_context.tsv"),
  c("resource", "context_label", "cell_resolution", "event_coordinate", "snp", "risk_allele",
    "risk_aligned_effect", "p_value")
)
context <- context_source[, .(
  resource, context_label, cell_resolution, event_coordinate, snp, risk_allele,
  risk_aligned_effect, p_value,
  permitted_role = fcase(
    cell_resolution == "cell_type_specific", "Direct cell-resolved exact-event sQTL",
    resource == "GTEx v8", "Cross-tissue exact-event consistency",
    default = "Contextual exact-event evidence"
  )
)]

coloc_source <- read_checked(
  file.path(root, "figures", "Figure_02", "data", "Figure_2_exact_event_coloc.tsv"),
  c("trait", "pph4", "pph4_conservative", "prior_interpretation")
)
triangulation <- coloc_source[, .(
  trait,
  pph4_default = pph4,
  pph4_conservative,
  prior_interpretation
)]

structure <- read_checked(
  file.path(root, "figures", "Figure_05", "data", "Figure_5_structure_summary.tsv"),
  c("interval", "mean_pLDDT", "local_PAE")
)[1]

edges <- data.table(
  source = c(
    "ad", "cholesterol", "ld_block", "exact_splice", "ld_block", "ld_block", "ld_block",
    "exact_splice", "ld_block", "published_editing", "published_editing",
    "ec2_boundary", "ec2_boundary", "ec2_boundary"
  ),
  target = c(
    "ld_block", "ld_block", "exact_splice", "cross_tissue_context", "microglia", "neural_eqtl",
    "ad_state_rna", "ec2_boundary", "published_editing", "adam10", "trem2",
    "isoform_prediction", "substrate_prediction", "lipid_state_prediction"
  ),
  edge_class = c(
    rep("present_study", 3), "context_support", "present_study", "context_support", "context_support",
    "reference_annotation", rep("published_perturbation", 3), rep("prediction", 3)
  ),
  interpretation = c(
    "AD association maps to the local configuration",
    "Cholesterol associations map to the same local configuration",
    "Exact exon5-6 processing is the leading candidate molecular readout",
    "Identical-coordinate cross-tissue consistency in partially overlapping GTEx tissues",
    "Cell-resolved exact-event sQTL in isolated microglia",
    "Single-nucleus total-expression eQTL",
    "Disease-state abundance assessed separately",
    "Reference localization of AA150/151 within EC2",
    "Published rs7922621 editing anchor",
    "Published cell-surface ADAM10 result",
    "Published soluble TREM2 result",
    "Protein-isoform consequence to be tested",
    "Trafficking and substrate processing to be tested",
    "Lipid-state cellular phenotype to be tested"
  )
)

stopifnot(
  nrow(triangulation) == 4L,
  all(triangulation$pph4_default > 0.95),
  all(triangulation$pph4_conservative > 0.69),
  uniqueN(context[resource == "GTEx v8", context_label]) == 4L,
  structure$mean_pLDDT > 90,
  nrow(edges[source == "cholesterol" & target %chin% c("ad", "exact_splice")]) == 0L
)

tables <- list(
  Figure_6_cross_tissue_context.tsv = context,
  Figure_6_exact_event_triangulation.tsv = triangulation,
  Figure_6_structure_summary.tsv = structure,
  Figure_6_edges.tsv = edges
)
for (filename in names(tables)) {
  fwrite(tables[[filename]], file.path(out_dir, filename), sep = "\t", quote = FALSE, na = "NA")
}

manifest <- data.table(
  artifact = names(tables),
  role = c(
    "Cell-resolved and cross-tissue exact-event context",
    "Default and conservative exact-event colocalization posteriors",
    "EC2 structural-confidence summary",
    "Evidence-class and prediction-edge ledger"
  )
)
fwrite(manifest, file.path(out_dir, "source_manifest.tsv"), sep = "\t", quote = FALSE)
message("Prepared Figure 6 source tables in: ", out_dir)
