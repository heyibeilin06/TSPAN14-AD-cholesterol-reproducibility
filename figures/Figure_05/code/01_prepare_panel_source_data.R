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

cell_dir <- file.path(root, "outputs", "mentor_revision", "cell_type_attribution")
p1_dir <- file.path(root, "outputs", "p1_cell_context")
p1_v16_dir <- file.path(root, "outputs", "p1_cell_context_v16")
splice_dir <- file.path(root, "outputs", "targeted_splice_validation")

cross_context <- read_checked(
  file.path(cell_dir, "01_exact_event_cross_context_evidence.tsv"),
  c("resource", "tissue_or_cell_context", "cell_resolution", "event_coordinate", "snp", "risk_allele",
    "risk_aligned_effect", "standard_error", "p_value", "effect_metric", "permitted_role")
)
snuc <- read_checked(
  file.path(p1_dir, "p1_niagads_laub_core_snuc_eqtl_direction.tsv"),
  c("rsID", "xQTLtype", "context", "study", "pvalue", "FDR", "beta", "se", "risk_allele_direction")
)
disease <- read_checked(
  file.path(p1_dir, "p1_seaad_pseudobulk_dge_tspan14.tsv"),
  c("logFC", "t", "P.Value", "adj.P.Val", "gene", "cell_group", "contrast", "n_comparator", "n_high_ad")
)
meta <- read_checked(
  file.path(p1_v16_dir, "nakatsuka2025_tspan14_meta_results.tsv"),
  c("source", "method", "tested_direction", "cell_type", "gene", "effect", "effect_metric", "p_value", "adjusted_p_value")
)
structure <- read_checked(
  file.path(splice_dir, "08_tspan14_ec2_secondary_structure.tsv"),
  c("residue", "region", "sse_psea", "pLDDT", "splice_context")
)
transcript <- read_checked(
  file.path(splice_dir, "03_transcript_structure_audit.tsv"),
  c("transcript_id", "event", "junction_coordinates", "structural_status", "protein_transition", "interpretation")
)
structure_summary <- read_checked(
  file.path(splice_dir, "09_tspan14_ec2_structural_context_summary.tsv"),
  c("interval", "mean_pLDDT", "local_PAE_AA147_154", "secondary_structure_counts", "interpretation", "claim_boundary")
)

# One lead marker is used for cross-context display to avoid counting LD-correlated variants as replications.
exact_context <- cross_context[snp == "rs1870137"]
exact_context[, context_label := fcase(
  resource == "MiGA", "Isolated microglia",
  grepl("Anterior cingulate", tissue_or_cell_context), "Anterior cingulate BA24",
  grepl("Hippocampus", tissue_or_cell_context), "Hippocampus",
  grepl("Putamen", tissue_or_cell_context), "Putamen",
  grepl("Spinal", tissue_or_cell_context), "Cervical spinal cord",
  default = tissue_or_cell_context
)]
exact_context[, evidence_scope := ifelse(resource == "MiGA", "Cell-resolved exact sQTL", "Bulk-tissue exact sQTL")]

# The same lead risk allele is displayed across available single-nucleus contexts and studies.
snuc_lead <- snuc[rsID == "rs7080009"]
snuc_lead[, `:=`(
  cell_context = fcase(context == "ast", "Astrocytes", context == "exc", "Excitatory neurons", default = context),
  study_label = fcase(
    study == "ROSMAP_CUIMC1_2_MIT", "ROSMAP combined",
    study == "ROSMAP_CUIMC1", "ROSMAP CUIMC",
    study == "ROSMAP_MIT", "ROSMAP MIT",
    default = study
  ),
  lo = beta - 1.96 * se,
  hi = beta + 1.96 * se
)]

disease[, `:=`(
  se = abs(logFC / t),
  lo = logFC - 1.96 * abs(logFC / t),
  hi = logFC + 1.96 * abs(logFC / t),
  source_label = "SEA-AD adjusted pseudobulk",
  cell_label = fcase(cell_group == "microglia", "Microglia", cell_group == "neurons", "Neurons", default = cell_group),
  estimate = logFC,
  fdr = adj.P.Val,
  interval_available = TRUE
)]

meta_selected <- meta[
  gene == "TSPAN14" & method == "SumRank" & tested_direction == "upregulated" &
    cell_type %chin% c("Microglia", "Oligodendrocytes")
]
meta_selected[, `:=`(
  source_label = "17-study single-cell meta-analysis",
  cell_label = cell_type,
  estimate = effect,
  fdr = adjusted_p_value,
  lo = NA_real_, hi = NA_real_, se = NA_real_,
  interval_available = FALSE,
  n_comparator = NA_real_, n_high_ad = NA_real_
)]

disease_display <- rbindlist(list(
  disease[, .(source_label, cell_label, estimate, se, lo, hi, P.Value, fdr,
              interval_available, n_comparator, n_high_ad)],
  meta_selected[, .(source_label, cell_label, estimate, se, lo, hi, P.Value = p_value, fdr,
                    interval_available, n_comparator, n_high_ad)]
), use.names = TRUE, fill = TRUE)

atlas <- data.table(
  evidence_layer = c(
    "Exact exon5-6 sQTL", "Exact exon5-6 sQTL",
    "Total-expression eQTL", "Total-expression eQTL",
    "Adjusted AD-state RNA", "Adjusted AD-state RNA",
    "Published perturbation"
  ),
  context = c(
    "Microglia", "Bulk neural tissue",
    "Excitatory neurons", "Astrocytes",
    "Microglia", "Excitatory neurons",
    "Microglia"
  ),
  evidence_class = c(
    "Exact splice QTL", "Exact splice QTL",
    "Gene-expression QTL", "Gene-expression QTL",
    "Disease-state test", "Disease-state test",
    "Perturbation"
  ),
  evidence_strength = c(4, 3, 3, 3, 2, 2, 3),
  interpretation = c(
    "Cell-resolved genetic regulation", "Cross-tissue consistency",
    "Risk allele increases total expression", "Risk allele increases total expression",
    "No FDR-significant abundance shift", "No FDR-significant abundance shift",
    "External regulatory and cellular anchor"
  )
)

stopifnot(
  nrow(exact_context) == 5L,
  nrow(snuc_lead) == 4L,
  nrow(disease_display) == 4L,
  all(disease_display$fdr > 0.05),
  nrow(structure[residue %between% c(114, 232)]) == 119L,
  transcript[event == "project_exon5_exon6", protein_transition] == "AA150/151",
  abs(structure[residue == 150, pLDDT] - 93.19) < 0.01,
  abs(structure[residue == 151, pLDDT] - 93.38) < 0.01
)

tables <- list(
  Figure_5_cell_context_atlas.tsv = atlas,
  Figure_5_exact_event_cross_context.tsv = exact_context,
  Figure_5_single_nucleus_eqtl.tsv = snuc_lead,
  Figure_5_disease_state_rna.tsv = disease_display,
  Figure_5_ec2_structure.tsv = structure,
  Figure_5_transcript_events.tsv = transcript,
  Figure_5_structure_summary.tsv = structure_summary
)
for (filename in names(tables)) fwrite(tables[[filename]], file.path(out_dir, filename), sep = "\t", na = "NA")

manifest <- data.table(
  artifact = names(tables),
  role = c(
    "Evidence-layer by cell-context atlas",
    "Lead-variant exact exon5-6 sQTL across isolated microglia and bulk neural tissues",
    "Lead-risk-allele single-nucleus TSPAN14 eQTL estimates",
    "Adjusted disease-state RNA and cross-study single-cell meta-analysis",
    "Residue-level EC2 secondary structure and AlphaFold confidence",
    "Exact and adjacent transcript-event coordinate audit",
    "AA150-193 structural-context summary"
  )
)
fwrite(manifest, file.path(out_dir, "Figure_5_source_manifest.tsv"), sep = "\t")
message("Prepared Figure 5 v9 source tables in: ", out_dir)
