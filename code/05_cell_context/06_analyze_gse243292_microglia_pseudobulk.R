#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(edgeR))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop(paste(
    "Usage: 231_gse243292_microglia_edger.R",
    "donor_counts.tsv donor_metadata.tsv state_counts.tsv state_metadata.tsv output_dir"
  ))
}

donor_counts <- read.delim(args[[1L]], row.names = 1L, check.names = FALSE)
donor_metadata <- read.delim(args[[2L]], stringsAsFactors = FALSE, check.names = FALSE)
state_counts <- read.delim(args[[3L]], row.names = 1L, check.names = FALSE)
state_metadata <- read.delim(args[[4L]], stringsAsFactors = FALSE, check.names = FALSE)
output_dir <- args[[5L]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
donor_metadata$sample <- as.character(donor_metadata$sample)
state_metadata$sample <- as.character(state_metadata$sample)
state_metadata$donor <- as.character(state_metadata$donor)
stopifnot(setequal(colnames(donor_counts), donor_metadata$sample))
stopifnot(setequal(colnames(state_counts), state_metadata$sample))
donor_metadata <- donor_metadata[match(colnames(donor_counts), donor_metadata$sample), , drop = FALSE]
state_metadata <- state_metadata[match(colnames(state_counts), state_metadata$sample), , drop = FALSE]

donor_metadata$log_median_umis <- log2(donor_metadata$median_umis + 1)
donor_metadata$log_median_umis_centered <- as.numeric(scale(
  donor_metadata$log_median_umis, center = TRUE, scale = FALSE
))
donor_metadata$apoe4_dosage_centered <- as.numeric(scale(
  donor_metadata$apoe4_dosage, center = TRUE, scale = FALSE
))
donor_metadata$reactive_fraction_centered <- as.numeric(scale(
  donor_metadata$reactive_associated_fraction, center = TRUE, scale = FALSE
))

models <- list(
  pathology_stage_adjusted = list(
    formula = ~ apoe4_dosage_centered + trem2_r47h +
      log_median_umis_centered + pathology_stage,
    coefficient = "pathology_stage"
  ),
  pathology_stage_composition_adjusted = list(
    formula = ~ apoe4_dosage_centered + trem2_r47h +
      log_median_umis_centered + reactive_fraction_centered + pathology_stage,
    coefficient = "pathology_stage"
  ),
  ad_binary_adjusted = list(
    formula = ~ apoe4_dosage_centered + trem2_r47h +
      log_median_umis_centered + ad_binary,
    coefficient = "ad_binary"
  )
)

run_edger <- function(counts, metadata, design, coefficient, model_name) {
  if (qr(design)$rank != ncol(design)) {
    stop(sprintf("Rank-deficient design for %s", model_name))
  }
  y <- DGEList(counts = counts)
  keep <- filterByExpr(y, design = design, min.count = 5L)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- normLibSizes(y, method = "TMM")
  y <- estimateDisp(y, design = design, robust = TRUE)
  fit <- glmQLFit(y, design = design, robust = TRUE)
  test <- glmQLFTest(fit, coef = which(colnames(design) == coefficient))
  result <- topTags(test, n = Inf, sort.by = "none")$table
  result$gene <- rownames(result)
  result$model <- model_name
  result$coefficient <- coefficient
  list(
    result = result,
    audit = data.frame(
      model = model_name,
      coefficient = coefficient,
      n_samples = nrow(metadata),
      n_genes_tested = nrow(y),
      design_rank = qr(design)$rank,
      design_columns = paste(colnames(design), collapse = ";"),
      condition_number = kappa(design),
      residual_df = min(fit$df.residual),
      stringsAsFactors = FALSE
    )
  )
}

donor_results <- list()
donor_audits <- list()
for (model_name in names(models)) {
  specification <- models[[model_name]]
  design <- model.matrix(specification$formula, data = donor_metadata)
  fit <- run_edger(
    donor_counts, donor_metadata, design, specification$coefficient, model_name
  )
  donor_results[[model_name]] <- fit$result
  donor_audits[[model_name]] <- fit$audit
}

paired_donors <- intersect(
  state_metadata$donor[state_metadata$state_group == "Reactive_associated"],
  state_metadata$donor[state_metadata$state_group == "Other_microglia"]
)
state_metadata <- state_metadata[state_metadata$donor %in% paired_donors, , drop = FALSE]
state_counts <- state_counts[, state_metadata$sample, drop = FALSE]
state_metadata$donor <- factor(state_metadata$donor)
state_metadata$state_group <- relevel(factor(state_metadata$state_group), ref = "Other_microglia")
state_design <- model.matrix(~ donor + state_group, data = state_metadata)
state_fit <- run_edger(
  state_counts,
  state_metadata,
  state_design,
  "state_groupReactive_associated",
  "within_donor_reactive_state"
)

all_results <- rbind(do.call(rbind, donor_results), state_fit$result)
tspan14 <- all_results[all_results$gene == "TSPAN14", , drop = FALSE]
audits <- rbind(do.call(rbind, donor_audits), state_fit$audit)

write.table(
  all_results,
  file.path(output_dir, "gse243292_microglia_edger_all_genes.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  tspan14,
  file.path(output_dir, "gse243292_microglia_edger_tspan14.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  audits,
  file.path(output_dir, "gse243292_microglia_edger_model_audit.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
