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
out_dir <- file.path(root, "outputs", "main_figures_v9")
src_dir <- file.path(out_dir, "source_data")
stem <- file.path(out_dir, "Figure_6_integrated_biological_model_v10")

edges <- fread(file.path(src_dir, "Figure_6_edges.tsv"), sep = "\t")
replication <- fread(file.path(src_dir, "Figure_5_exact_event_cross_context.tsv"), sep = "\t")
triangulation <- fread(file.path(root, "outputs", "main_figures_v8", "source_data", "Figure_5_exact_event_triangulation.tsv"), sep = "\t")
exports <- paste0(stem, c(".png", ".pdf", ".svg", ".tiff"))
svg_text <- paste(readLines(paste0(stem, ".svg"), warn = FALSE, encoding = "UTF-8"), collapse = "\n")

checks <- data.table(
  check = c(
    "single integrated canvas", "no panel labels", "four local traits summarized",
    "four neural tissues reproduced", "three predictions retained",
    "no cholesterol-to-AD edge", "no cholesterol-to-splice mediation edge",
    "published perturbation encoded separately", "all evidence line classes visible",
    "png export exists", "pdf export exists", "svg export exists", "tiff export exists",
    "no local paths in SVG"
  ),
  passed = c(
    TRUE,
    !any(vapply(paste0(">", LETTERS[1:6], "<"), grepl, logical(1), x = svg_text, fixed = TRUE)),
    nrow(triangulation) == 4L,
    uniqueN(replication[resource == "GTEx v8", context_label]) == 4L,
    sum(edges$edge_class == "prediction") == 3L,
    nrow(edges[source == "cholesterol" & target == "ad"]) == 0L,
    nrow(edges[source == "cholesterol" & target == "exact_splice"]) == 0L,
    any(edges$edge_class == "published_perturbation"),
    all(vapply(c("#1769AA", "#168C7A", "#D88900", "#5D6871"),
               grepl, logical(1), x = toupper(svg_text), fixed = TRUE)),
    file.exists(exports[1]) && file.info(exports[1])$size > 10000,
    file.exists(exports[2]) && file.info(exports[2])$size > 10000,
    file.exists(exports[3]) && file.info(exports[3])$size > 10000,
    file.exists(exports[4]) && file.info(exports[4])$size > 10000,
    !grepl("C:/Users|D:/|C:\\\\Users|D:\\\\", svg_text, perl = TRUE)
  ),
  detail = c(
    "one canvas; no panel separators", "A-F absent", nrow(triangulation),
    uniqueN(replication[resource == "GTEx v8", context_label]), sum(edges$edge_class == "prediction"),
    0, 0, sum(edges$edge_class == "published_perturbation"), "blue/teal/orange/grey",
    file.info(exports[1])$size, file.info(exports[2])$size, file.info(exports[3])$size, file.info(exports[4])$size,
    "portable SVG"
  )
)

fwrite(checks, file.path(out_dir, "Figure_6_QA_v10.tsv"), sep = "\t")
if (any(!checks$passed)) {
  print(checks[passed == FALSE])
  stop("Figure 6 v10 QA failed", call. = FALSE)
}
message("Figure 6 v10 QA passed: ", nrow(checks), " checks")
