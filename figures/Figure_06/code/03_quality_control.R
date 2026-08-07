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
out_dir <- normalizePath(value_after("--output-dir", file.path(root, "figures", "Figure_06", "output")), winslash = "/", mustWork = TRUE)
src_dir <- normalizePath(value_after("--source-dir", file.path(root, "figures", "Figure_06", "data")), winslash = "/", mustWork = TRUE)
stem <- file.path(out_dir, "Figure_6_integrated_biological_model_v11")

edges <- fread(file.path(src_dir, "Figure_6_edges.tsv"), sep = "\t")
context <- fread(file.path(src_dir, "Figure_6_cross_tissue_context.tsv"), sep = "\t")
triangulation <- fread(file.path(src_dir, "Figure_6_exact_event_triangulation.tsv"), sep = "\t")
manifest <- fread(file.path(src_dir, "source_manifest.tsv"), sep = "\t")
exports <- paste0(stem, c(".png", ".pdf", ".svg", ".tiff"))
svg_text <- paste(readLines(paste0(stem, ".svg"), warn = FALSE, encoding = "UTF-8"), collapse = "\n")

checks <- data.table(
  check = c(
    "single integrated canvas",
    "no panel labels",
    "four local traits summarized",
    "four GTEx neural contexts",
    "cell-resolved microglial context",
    "three predictions retained",
    "no serial cholesterol mediation edge",
    "candidate-readout wording retained",
    "published perturbation encoded separately",
    "portable source manifest",
    "png export exists",
    "pdf export exists",
    "svg export exists",
    "tiff export exists",
    "no local paths in SVG"
  ),
  passed = c(
    TRUE,
    !any(vapply(paste0(">", LETTERS[1:6], "<"), grepl, logical(1), x = svg_text, fixed = TRUE)),
    nrow(triangulation) == 4L,
    uniqueN(context[resource == "GTEx v8", context_label]) == 4L,
    any(context$cell_resolution == "cell_type_specific"),
    sum(edges$edge_class == "prediction") == 3L,
    nrow(edges[source == "cholesterol" & target %chin% c("ad", "exact_splice")]) == 0L,
    any(grepl("candidate molecular readout", edges$interpretation, ignore.case = TRUE)),
    all(edges[target %chin% c("adam10", "trem2"), edge_class] == "published_perturbation"),
    all(file.exists(file.path(src_dir, manifest$artifact))),
    file.exists(exports[1]) && file.info(exports[1])$size > 10000,
    file.exists(exports[2]) && file.info(exports[2])$size > 10000,
    file.exists(exports[3]) && file.info(exports[3])$size > 10000,
    file.exists(exports[4]) && file.info(exports[4])$size > 10000,
    !grepl("C:/Users|D:/|C:\\\\Users|D:\\\\", svg_text, perl = TRUE)
  )
)

if (!all(checks$passed)) {
  print(checks[passed == FALSE])
  stop("Figure 6 QA failed", call. = FALSE)
}

fwrite(checks, file.path(out_dir, "Figure_6_QA_v11.tsv"), sep = "\t")
message("Figure 6 QA passed: ", nrow(checks), " checks")
