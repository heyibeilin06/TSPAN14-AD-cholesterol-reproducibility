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
out_dir <- file.path(root, "outputs", "main_figures_v9", "source_data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_checked <- function(path, required) {
  if (!file.exists(path)) stop("Missing audited source: ", path, call. = FALSE)
  x <- fread(path, sep = "\t", na.strings = c("", "NA"))
  missing <- setdiff(required, names(x))
  if (length(missing)) stop("Missing columns in ", basename(path), ": ", paste(missing, collapse = ", "), call. = FALSE)
  x
}

claim_map <- read_checked(
  file.path(root, "outputs", "main_figures_v8", "source_data", "Figure_5_evidence_claim_map.tsv"),
  c("figure_element", "evidence_class", "supported_claim", "claim_boundary")
)
causal <- read_checked(
  file.path(root, "outputs", "mentor_revision", "complete_mr", "22_final_causal_claim_ledger.tsv"),
  c("claim", "status", "prohibited_extension")
)
cell <- read_checked(
  file.path(root, "outputs", "mentor_revision", "cell_type_attribution", "03_cell_type_layer_separation.tsv"),
  c("evidence_layer", "biological_material", "molecular_phenotype", "causal_axis", "cell_type_attribution", "permitted_claim", "prohibited_claim")
)
replication <- read_checked(
  file.path(out_dir, "Figure_3_exact_replication_matrix.tsv"),
  c("snp", "tissue", "p_value", "nes", "event_class")
)
structure <- read_checked(
  file.path(out_dir, "Figure_5_structure_summary.tsv"),
  c("interval", "mean_pLDDT", "local_PAE_AA147_154")
)

nodes <- data.table(
  node_id = c(
    "ad", "cholesterol", "ld_block", "exact_splice", "replication",
    "microglia", "neural_eqtl", "ad_state_rna", "published_editing",
    "ec2_boundary", "adam10", "trem2",
    "isoform_prediction", "substrate_prediction", "lipid_state_prediction"
  ),
  display_label = c(
    "Alzheimer disease", "TC | LDL-C | non-HDL-C", "TSPAN14 regulatory LD block",
    "Exact exon5-6 processing", "Identical junction in four neural tissues",
    "Isolated microglia", "Excitatory neurons / astrocytes", "AD-state total RNA",
    "rs7922621 editing", "EC2 boundary AA150/151", "ADAM10", "TREM2",
    "TSPAN14 protein-isoform composition", "ADAM10 trafficking / substrate processing",
    "Lipid-state-dependent cellular phenotypes"
  ),
  evidence_class = c(
    "present_study", "present_study", "present_study", "present_study", "independent_context",
    "present_study", "independent_context", "independent_context", "published_perturbation",
    "reference_annotation", "published_perturbation", "published_perturbation",
    "prediction", "prediction", "prediction"
  ),
  supported_scope = c(
    "local association", "local association", "shared regulatory configuration", "genotype-regulated exact junction",
    "bulk neural replication", "cell-resolved exact sQTL", "cell-resolved total-expression eQTL",
    "no uniform FDR-significant shift", "variant-to-gene functional anchor", "reference structural localization",
    "published cell-surface phenotype", "published soluble-shedding phenotype",
    "testable exact-splice consequence", "testable exact-splice consequence", "testable cellular consequence"
  )
)

edges <- data.table(
  source = c(
    "ad", "cholesterol", "ld_block", "exact_splice", "exact_splice",
    "ld_block", "ld_block", "ld_block", "ld_block",
    "published_editing", "published_editing", "exact_splice",
    "ec2_boundary", "ec2_boundary", "ec2_boundary"
  ),
  target = c(
    "ld_block", "ld_block", "exact_splice", "replication", "ec2_boundary",
    "microglia", "neural_eqtl", "ad_state_rna", "published_editing",
    "adam10", "trem2", "ec2_boundary",
    "isoform_prediction", "substrate_prediction", "lipid_state_prediction"
  ),
  edge_class = c(
    "present_study", "present_study", "present_study", "independent_context", "reference_annotation",
    "present_study", "independent_context", "independent_context", "published_perturbation",
    "published_perturbation", "published_perturbation", "reference_annotation",
    "prediction", "prediction", "prediction"
  ),
  interpretation = c(
    "AD association maps to the local configuration", "cholesterol associations map to the same local configuration",
    "the locus resolves to exact transcript processing", "identical-coordinate replication",
    "coding boundary localization", "isolated-microglia exact sQTL", "single-nucleus total-expression eQTL",
    "disease-state abundance assessed separately", "published editing anchor",
    "published cell-surface ADAM10 result", "published soluble TREM2 result", "reference localization",
    "protein isoform prediction", "trafficking and cleavage prediction", "lipid-state phenotype prediction"
  )
)

replication_summary <- replication[snp == "rs6586028" & event_class == "exact_exon5_6", .(
  n_tissues = uniqueN(tissue),
  minimum_nes = min(nes), maximum_nes = max(nes),
  maximum_p_value = max(p_value)
)]
structure_summary <- structure[, .(
  interval = interval[1], mean_pLDDT = mean_pLDDT[1], local_PAE = local_PAE_AA147_154[1]
)]

# Guard against drawing a serial lipid-mediated causal path that the MR analyses did not identify.
forbidden_pairs <- edges[source == "cholesterol" & target %chin% c("ad", "exact_splice")]
stopifnot(
  nrow(forbidden_pairs) == 0L,
  nrow(edges[source == "exact_splice" & target %chin% c("adam10", "trem2") & edge_class != "prediction"]) == 0L,
  replication_summary$n_tissues == 4L,
  structure_summary$mean_pLDDT > 90,
  any(causal$status == "not_identified"),
  nrow(claim_map) >= 7L,
  nrow(cell) >= 5L
)

tables <- list(
  Figure_6_nodes.tsv = nodes,
  Figure_6_edges.tsv = edges,
  Figure_6_replication_summary.tsv = replication_summary,
  Figure_6_structure_summary.tsv = structure_summary,
  Figure_6_claim_boundaries.tsv = claim_map,
  Figure_6_causal_scope.tsv = causal,
  Figure_6_cell_context_scope.tsv = cell
)
for (filename in names(tables)) fwrite(tables[[filename]], file.path(out_dir, filename), sep = "\t", na = "NA")

manifest <- data.table(
  artifact = names(tables),
  role = c(
    "Mechanism-model nodes and evidence classes", "Directed evidence and prediction edges",
    "Four-tissue exact-event replication summary", "EC2 structural-confidence summary",
    "Supported-claim and boundary ledger", "MR causal-scope guardrail", "Cell-attribution guardrail"
  )
)
fwrite(manifest, file.path(out_dir, "Figure_6_source_manifest.tsv"), sep = "\t")
message("Prepared Figure 6 v9 source tables in: ", out_dir)
