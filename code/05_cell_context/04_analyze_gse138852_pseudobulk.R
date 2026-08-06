#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: 229_gse138852_pool_pseudobulk_edger.R counts.tsv metadata.tsv output_dir")
}

count_path <- args[[1L]]
metadata_path <- args[[2L]]
output_dir <- args[[3L]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

counts <- read.delim(count_path, row.names = 1L, check.names = FALSE)
metadata <- read.delim(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(identical(colnames(counts), metadata$sample))

all_results <- list()
tspan14_results <- list()
model_audit <- list()

for (cell_type in unique(metadata$cell_type)) {
  keep_samples <- metadata$cell_type == cell_type
  sample_data <- droplevels(transform(
    metadata[keep_samples, , drop = FALSE],
    condition = relevel(factor(condition), ref = "Control"),
    pair = factor(pair)
  ))
  cell_counts <- counts[, keep_samples, drop = FALSE]
  design <- model.matrix(~ pair + condition, data = sample_data)
  stopifnot(qr(design)$rank == ncol(design))

  y <- DGEList(counts = cell_counts, samples = sample_data)
  keep_genes <- filterByExpr(y, design = design, min.count = 5L)
  y <- y[keep_genes, , keep.lib.sizes = FALSE]
  y <- normLibSizes(y, method = "TMM")
  y <- estimateDisp(y, design = design, robust = TRUE)
  fit <- glmQLFit(y, design = design, robust = TRUE)
  coefficient <- which(colnames(design) == "conditionAD")
  test <- glmQLFTest(fit, coef = coefficient)
  result <- topTags(test, n = Inf, sort.by = "none")$table
  result$gene <- rownames(result)
  result$cell_type <- cell_type
  result$contrast <- "AD_vs_Control"
  result$n_ad_pools <- sum(sample_data$condition == "AD")
  result$n_control_pools <- sum(sample_data$condition == "Control")
  result$minimum_cells_per_pool <- min(sample_data$n_cells)
  all_results[[cell_type]] <- result
  tspan14_results[[cell_type]] <- result[result$gene == "TSPAN14", , drop = FALSE]
  model_audit[[cell_type]] <- data.frame(
    cell_type = cell_type,
    n_genes_input = nrow(cell_counts),
    n_genes_tested = nrow(y),
    design_rank = qr(design)$rank,
    residual_df = fit$df.residual[1L],
    n_ad_pools = sum(sample_data$condition == "AD"),
    n_control_pools = sum(sample_data$condition == "Control"),
    minimum_cells_per_pool = min(sample_data$n_cells),
    inference_level = "pooled library; paired-block sensitivity analysis"
  )
}

write.table(
  do.call(rbind, all_results),
  file.path(output_dir, "gse138852_pool_pseudobulk_edger_all_genes.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  do.call(rbind, tspan14_results),
  file.path(output_dir, "gse138852_pool_pseudobulk_edger_tspan14.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  do.call(rbind, model_audit),
  file.path(output_dir, "gse138852_pool_pseudobulk_edger_model_audit.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
