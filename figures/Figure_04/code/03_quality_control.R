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
out_dir <- normalizePath(value_after("--output-dir", file.path(root, "figures", "Figure_04", "output")), winslash = "/", mustWork = TRUE)
source_dir <- normalizePath(value_after("--source-dir", file.path(root, "figures", "Figure_04", "data")), winslash = "/", mustWork = TRUE)

checks <- list()
add_check <- function(name, passed, detail) {
  checks[[length(checks) + 1L]] <<- data.table(check = name, passed = isTRUE(passed), detail = as.character(detail))
}

instrument <- fread(file.path(source_dir, "Figure_4_instrument_effect_atlas.tsv"))
cis <- fread(file.path(source_dir, "Figure_4_ld_aware_cis_mr.tsv"))
diag <- fread(file.path(source_dir, "Figure_4_cis_mr_diagnostics.tsv"))
global <- fread(file.path(source_dir, "Figure_4_genomewide_lipid_to_ad.tsv"))
joint <- fread(file.path(source_dir, "Figure_4_global_joint_mvmr.tsv"))
strength <- fread(file.path(source_dir, "Figure_4_global_joint_mvmr_strength.tsv"))
pc <- fread(file.path(source_dir, "Figure_4_pc_gmm_dimension_sensitivity.tsv"))

add_check("five cis instruments", uniqueN(instrument$SNP) == 5, uniqueN(instrument$SNP))
add_check("five instrument layers", uniqueN(instrument$display_layer) == 5, uniqueN(instrument$display_layer))
add_check("eight cis-MR estimates", nrow(cis) == 8, nrow(cis))
add_check("all cis-MR estimates positive", all(cis$estimate > 0), range(cis$estimate))
add_check("diagnostic matrix complete", nrow(diag) == 20, nrow(diag))
add_check("two non-HDL diagnostic flags", sum(!diag$pass) == 2 && all(diag[pass == FALSE, outcome] == "nonHDL"), sum(!diag$pass))
add_check("two global models by five lipids", nrow(global) == 10 && uniqueN(global$exposure) == 5, nrow(global))
add_check("three global MVMR estimates", nrow(joint) == 3, nrow(joint))
add_check("global MVMR estimates non-significant", all(joint$pvalue > 0.05), range(joint$pvalue))
add_check("global MVMR conditional strength passes", all(strength$conditional_F > 10), range(strength$conditional_F))
add_check("three PC-GMM lipid paths", uniqueN(pc$lipid) == 3, uniqueN(pc$lipid))
add_check("no confirmatory mediation model", !any(pc$confirmatory_mediation), sum(pc$confirmatory_mediation))
add_check("PC-GMM strength instability visible", any(pc$all_strength_F_ge_10) && any(!pc$all_strength_F_ge_10), paste(table(pc$all_strength_F_ge_10), collapse = ","))

stem <- file.path(out_dir, "Figure_4_causal_scope_v9")
for (extension in c("png", "pdf", "svg", "tiff")) {
  path <- paste0(stem, ".", extension)
  add_check(paste(extension, "export exists"), file.exists(path) && file.info(path)$size > 10000, if (file.exists(path)) file.info(path)$size else 0)
}

svg_path <- paste0(stem, ".svg")
if (file.exists(svg_path)) {
  svg_text <- paste(readLines(svg_path, warn = FALSE), collapse = "\n")
  add_check("uppercase panel labels", all(vapply(LETTERS[1:6], function(x) grepl(paste0(">", x, "<"), svg_text, fixed = TRUE), logical(1))), "A-F")
  add_check("lipid attenuation label", grepl("Lipid-coefficient attenuation", svg_text, fixed = TRUE), "correct Panel F interpretation")
  add_check("obsolete splice attenuation label absent", !grepl("Splice-coefficient attenuation", svg_text, fixed = TRUE), "obsolete label removed")
  add_check("no local paths in SVG", !grepl("C:/Users|D:/|C:\\\\Users|D:\\\\", svg_text, perl = TRUE), "portable SVG")
}

qa <- rbindlist(checks)
fwrite(qa, file.path(out_dir, "Figure_4_QA.tsv"), sep = "\t")
if (any(!qa$passed)) {
  print(qa[passed == FALSE])
  stop("Figure 4 QA failed", call. = FALSE)
}
message("Figure 4 QA passed: ", nrow(qa), " checks")
