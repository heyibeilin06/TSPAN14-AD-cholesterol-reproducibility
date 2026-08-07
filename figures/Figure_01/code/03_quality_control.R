suppressPackageStartupMessages({
  library(data.table)
  library(magick)
})

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) return(default)
  if (idx == length(args)) stop("Missing value for ", flag, call. = FALSE)
  args[[idx + 1L]]
}

root <- normalizePath(value_after("--project-root", "."), winslash = "/", mustWork = TRUE)
out_dir <- normalizePath(value_after("--output-dir", file.path(root, "figures", "Figure_01", "output")), winslash = "/", mustWork = TRUE)
source_dir <- normalizePath(value_after("--source-dir", file.path(root, "figures", "Figure_01", "data")), winslash = "/", mustWork = TRUE)
stem <- "Figure_1_APOE_aware_TSPAN14_prioritization_v9"

checks <- list()
record_check <- function(name, passed, detail) {
  checks[[length(checks) + 1L]] <<- data.table(check = name, passed = isTRUE(passed), detail = as.character(detail))
  if (!isTRUE(passed)) stop("QA failed: ", name, " — ", detail, call. = FALSE)
}

ldsc <- fread(file.path(source_dir, "Figure_1_ldsc_apoe_conditioning.tsv"), sep = "\t")
regional <- fread(file.path(source_dir, "Figure_1_regional_screen.tsv"), sep = "\t")
prior <- fread(file.path(source_dir, "Figure_1_regional_prior_sensitivity.tsv"), sep = "\t")

get_one <- function(x, condition, column) {
  value <- x[eval(condition), get(column)]
  if (length(value) != 1L) stop("Expected one value for ", column, call. = FALSE)
  value
}

expected_ldsc <- data.table(
  model = c("baseline", "own_lead", "pair_union_leads", "w5Mb"),
  expected = c(0.1394, 0.1155, 0.1181, 0.1059)
)
for (i in seq_len(nrow(expected_ldsc))) {
  model_i <- expected_ldsc$model[[i]]
  observed <- ldsc[model == model_i & comparison == "AD-HDL", rg]
  record_check(
    paste0("LDSC value: ", model_i),
    length(observed) == 1L && abs(observed - expected_ldsc$expected[[i]]) < 1e-10,
    paste("observed", paste(observed, collapse = ","), "expected", expected_ldsc$expected[[i]])
  )
}

expected_coloc <- data.table(
  trait_label = c("TC", "LDL-C", "non-HDL-C"),
  expected = c(0.975651086233109, 0.961885946151768, 0.958600618864801)
)
for (i in seq_len(nrow(expected_coloc))) {
  trait_i <- expected_coloc$trait_label[[i]]
  observed <- regional[focus_class == "TSPAN14" & trait_label == trait_i, PP.H4]
  record_check(
    paste0("TSPAN14 regional PP.H4: ", trait_i),
    length(observed) == 1L && abs(observed - expected_coloc$expected[[i]]) < 1e-12,
    paste("observed", paste(observed, collapse = ","), "expected", expected_coloc$expected[[i]])
  )
}

extensions <- c("svg", "pdf", "tiff", "png")
for (extension in extensions) {
  path <- file.path(out_dir, paste0(stem, ".", extension))
  record_check(
    paste0("Export exists: ", extension),
    file.exists(path) && file.info(path)$size > 10000,
    if (file.exists(path)) paste("bytes", file.info(path)$size) else "missing"
  )
}

png_path <- file.path(out_dir, paste0(stem, ".png"))
tiff_path <- file.path(out_dir, paste0(stem, ".tiff"))
png_info <- image_info(image_read(png_path))[1, ]
tiff_info <- image_info(image_read(tiff_path))[1, ]

expected_png <- c(width = round(183 / 25.4 * 300), height = round(168 / 25.4 * 300))
expected_tiff <- c(width = round(183 / 25.4 * 600), height = round(168 / 25.4 * 600))
record_check(
  "PNG final-size geometry",
  abs(png_info$width - expected_png[["width"]]) <= 2 && abs(png_info$height - expected_png[["height"]]) <= 2,
  paste(png_info$width, "x", png_info$height)
)
record_check(
  "TIFF 600-dpi geometry",
  abs(tiff_info$width - expected_tiff[["width"]]) <= 2 && abs(tiff_info$height - expected_tiff[["height"]]) <= 2,
  paste(tiff_info$width, "x", tiff_info$height)
)

svg_path <- file.path(out_dir, paste0(stem, ".svg"))
svg_text <- paste(readLines(svg_path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
record_check(
  "SVG has no local paths",
  !grepl("C:/Users|D:/|C:\\\\Users|D:\\\\", svg_text, perl = TRUE),
  "path scan"
)
for (label in LETTERS[1:4]) {
  record_check(
    paste0("SVG panel label ", label),
    grepl(paste0(">", label, "</text>"), svg_text, fixed = TRUE),
    "text-node scan"
  )
}

ledger <- file.path(out_dir, "Figure_1_manuscript_citation_ledger.tsv")
ledger_data <- fread(ledger, sep = "\t")
record_check(
  "Citation ledger order",
  identical(ledger_data$figure_citation, paste0("Fig. 1", LETTERS[1:4])),
  paste(ledger_data$figure_citation, collapse = ", ")
)

qa <- rbindlist(checks)
fwrite(qa, file.path(out_dir, "Figure_1_QA.tsv"), sep = "\t")
cat("Figure 1 QA passed:", nrow(qa), "checks\n")
