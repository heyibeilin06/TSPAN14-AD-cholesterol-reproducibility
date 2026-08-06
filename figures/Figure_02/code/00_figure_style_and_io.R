suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(scales)
})

v9_palette <- c(
  ink = "#20262B",
  graphite = "#5D6871",
  grey = "#9CA6AC",
  pale = "#EDF1F3",
  grid = "#D7DEE2",
  blue = "#1769AA",
  blue_dark = "#13496F",
  blue_soft = "#BFD8E8",
  teal = "#168C7A",
  teal_dark = "#0B665D",
  teal_soft = "#C7E2DC",
  purple = "#75529B",
  purple_soft = "#D8CAE4",
  amber = "#D88900",
  amber_soft = "#F2D69C",
  red = "#B24A43",
  red_soft = "#E9C5C1",
  white = "#FFFFFF"
)

v9_trait_colours <- c(
  "HDL-C" = v9_palette[["amber"]],
  "LDL-C" = v9_palette[["blue"]],
  "TG" = "#6E7D86",
  "TC" = v9_palette[["teal"]],
  "non-HDL-C" = v9_palette[["purple"]]
)

theme_v9 <- function(base_size = 6.6) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.32, colour = v9_palette[["ink"]]),
      axis.ticks = element_line(linewidth = 0.30, colour = v9_palette[["ink"]]),
      axis.ticks.length = grid::unit(1.1, "mm"),
      axis.title = element_text(size = base_size, colour = v9_palette[["ink"]]),
      axis.text = element_text(size = base_size - 0.2, colour = v9_palette[["ink"]]),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold", colour = v9_palette[["ink"]]),
      panel.grid = element_blank(),
      legend.title = element_text(size = base_size - 0.2),
      legend.text = element_text(size = base_size - 0.5),
      legend.key.height = grid::unit(2.7, "mm"),
      legend.key.width = grid::unit(3.5, "mm"),
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      plot.margin = margin(2, 2, 2, 2)
    )
}

theme_set(theme_v9())

panel_tag_theme_v9 <- theme(
  plot.tag = element_text(
    size = 9.0, face = "bold", family = "Arial", colour = v9_palette[["ink"]],
    margin = margin(0, 1.3, 0, 0)
  )
)

parse_v9_args <- function(default_output = "outputs/main_figures_v9") {
  args <- commandArgs(trailingOnly = TRUE)
  value_after <- function(flag, default = NULL) {
    idx <- match(flag, args)
    if (is.na(idx)) return(default)
    if (idx == length(args)) stop("Missing value for ", flag, call. = FALSE)
    args[[idx + 1L]]
  }
  root <- normalizePath(value_after("--project-root", "."), winslash = "/", mustWork = TRUE)
  output_arg <- value_after("--output-dir", file.path(root, default_output))
  source_arg <- value_after("--source-dir", file.path(output_arg, "source_data"))
  output_dir <- normalizePath(output_arg, winslash = "/", mustWork = FALSE)
  source_dir <- normalizePath(source_arg, winslash = "/", mustWork = TRUE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  list(root = root, output_dir = output_dir, source_dir = source_dir)
}

read_v9_source <- function(source_dir, filename, required = character()) {
  path <- file.path(source_dir, filename)
  if (!file.exists(path)) stop("Missing source data: ", path, call. = FALSE)
  x <- fread(path, sep = "\t", na.strings = c("NA", ""))
  missing <- setdiff(required, names(x))
  if (length(missing)) stop("Missing columns in ", filename, ": ", paste(missing, collapse = ", "), call. = FALSE)
  numeric_columns <- names(x)[vapply(x, is.numeric, logical(1))]
  for (column in numeric_columns) {
    if (any(!is.finite(x[[column]]) & !is.na(x[[column]]))) {
      stop("Non-finite value in ", filename, " column ", column, call. = FALSE)
    }
  }
  x
}

assert_unique_v9 <- function(x, keys, label) {
  if (x[, anyDuplicated(.SD), .SDcols = keys] > 0) {
    stop("Duplicate keys in ", label, ": ", paste(keys, collapse = ", "), call. = FALSE)
  }
}

format_p_v9 <- function(p) {
  fifelse(
    is.na(p), "",
    fifelse(p < 0.001, format(p, scientific = TRUE, digits = 2), sprintf("%.3f", p))
  )
}

save_pub_v9 <- function(plot, output_dir, stem, width_mm = 183, height_mm = 165) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4

  svglite::svglite(file.path(output_dir, paste0(stem, ".svg")), width = width_in, height = height_in)
  print(plot)
  grDevices::dev.off()

  grDevices::cairo_pdf(file.path(output_dir, paste0(stem, ".pdf")), width = width_in, height = height_in, family = "Arial")
  print(plot)
  grDevices::dev.off()

  ragg::agg_tiff(
    file.path(output_dir, paste0(stem, ".tiff")), width = width_in, height = height_in,
    units = "in", res = 600, compression = "lzw", background = "white"
  )
  print(plot)
  grDevices::dev.off()

  ragg::agg_png(
    file.path(output_dir, paste0(stem, ".png")), width = width_in, height = height_in,
    units = "in", res = 300, background = "white"
  )
  print(plot)
  grDevices::dev.off()
}

assert_clean_svg_v9 <- function(path) {
  text <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  if (grepl("C:/Users|D:/|C:\\\\Users|D:\\\\", text, perl = TRUE)) {
    stop("Local absolute path found in SVG", call. = FALSE)
  }
  invisible(TRUE)
}
