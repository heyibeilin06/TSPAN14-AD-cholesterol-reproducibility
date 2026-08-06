#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args) >= 2 && args[1] == "--project-root") args[2] else ".",
                      winslash = "/", mustWork = TRUE)
out_dir <- file.path(root, "outputs", "main_figures_v9")
source_dir <- file.path(out_dir, "source_data")

checks <- list()
add_check <- function(name, passed, detail) {
  checks[[length(checks) + 1L]] <<- data.table(check = name, passed = isTRUE(passed), detail = as.character(detail))
}

tracks <- fread(file.path(source_dir, "Figure_2_locus_tracks_grch38.tsv"), sep = "\t")
scatter <- fread(file.path(source_dir, "Figure_2_colocalization_scatter.tsv"), sep = "\t")
transcripts <- fread(file.path(source_dir, "Figure_2_gencode_v38_transcripts.tsv"), sep = "\t")
ccre <- fread(file.path(source_dir, "Figure_2_regulatory_elements.tsv"), sep = "\t")
coloc <- fread(file.path(source_dir, "Figure_2_exact_event_coloc.tsv"), sep = "\t")
annotation <- fread(file.path(source_dir, "Figure_2_variant_annotation_matrix.tsv"), sep = "\t")

add_check("five aligned tracks", setequal(unique(tracks$track), c("AD", "TC", "LDL", "nonHDL", "BA24 exact exon5-6 sQTL")), paste(unique(tracks$track), collapse = ", "))
add_check("regional variants present", nrow(tracks) > 10000, nrow(tracks))
add_check("four locus comparisons", setequal(unique(scatter$comparison), c("AD", "TC", "LDL", "nonHDL")), paste(unique(scatter$comparison), collapse = ", "))
add_check("exact-event posteriors complete", nrow(coloc) == 4 && all(coloc$pph4 > 0.95), paste(round(coloc$pph4, 4), collapse = ", "))
add_check("canonical transcript present", any(grepl("ENST00000429989", transcripts$display_name)), "ENST00000429989")
add_check("47 cCREs retained", nrow(ccre) == 47, nrow(ccre))
add_check("annotation matrix complete", nrow(annotation) == length(unique(annotation$snp)) * 7, nrow(annotation))
add_check("CRISPRi anchors retained", sum(annotation$metric == "CRISPRi" & annotation$value == 1) == 3, sum(annotation$metric == "CRISPRi" & annotation$value == 1))
add_check("prime-editing anchor retained", sum(annotation$metric == "Prime editing" & annotation$value == 1) == 1, sum(annotation$metric == "Prime editing" & annotation$value == 1))
add_check("lead SNP present in every track", tracks[snp == "rs7080009", uniqueN(track)] == 5, tracks[snp == "rs7080009", uniqueN(track)])

stem <- file.path(out_dir, "Figure_2_variant_level_TSPAN14_locus_v9")
for (extension in c("png", "pdf", "svg", "tiff")) {
  path <- paste0(stem, ".", extension)
  add_check(paste(extension, "export exists"), file.exists(path) && file.info(path)$size > 10000, if (file.exists(path)) file.info(path)$size else 0)
}

svg_path <- paste0(stem, ".svg")
if (file.exists(svg_path)) {
  svg_text <- paste(readLines(svg_path, warn = FALSE), collapse = "\n")
  add_check("uppercase panel labels", all(vapply(c("A", "B", "C", "D"), function(x) grepl(paste0(">", x, "<"), svg_text, fixed = TRUE), logical(1))), "A-D")
  add_check("no local paths in SVG", !grepl("C:/Users|D:/|C:\\\\Users|D:\\\\", svg_text, perl = TRUE), "portable SVG")
}

qa <- rbindlist(checks)
fwrite(qa, file.path(out_dir, "Figure_2_QA.tsv"), sep = "\t")
if (any(!qa$passed)) {
  print(qa[passed == FALSE])
  stop("Figure 2 QA failed", call. = FALSE)
}
message("Figure 2 QA passed: ", nrow(qa), " checks")
