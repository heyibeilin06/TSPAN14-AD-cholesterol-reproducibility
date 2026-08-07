#!/usr/bin/env Rscript

script_arg <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)])
script_dir <- dirname(normalizePath(script_arg, winslash = "/", mustWork = TRUE))
source(file.path(script_dir, "00_figure_style_and_io.R"))

paths <- parse_v9_args()
context <- read_v9_source(
  paths$source_dir, "Figure_6_cross_tissue_context.tsv",
  c("resource", "context_label", "cell_resolution", "event_coordinate", "risk_aligned_effect", "p_value")
)
triangulation <- read_v9_source(
  paths$source_dir, "Figure_6_exact_event_triangulation.tsv",
  c("trait", "pph4_default", "pph4_conservative", "prior_interpretation")
)
structure <- read_v9_source(
  paths$source_dir, "Figure_6_structure_summary.tsv",
  c("interval", "mean_pLDDT", "local_PAE")
)
edges <- read_v9_source(
  paths$source_dir, "Figure_6_edges.tsv",
  c("source", "target", "edge_class", "interpretation")
)

stopifnot(
  nrow(triangulation) == 4L,
  all(triangulation$pph4_default > 0.95),
  all(triangulation$pph4_conservative > 0.69),
  uniqueN(context[resource == "GTEx v8", context_label]) == 4L,
  any(context$cell_resolution == "cell_type_specific"),
  structure$mean_pLDDT[1] > 90,
  nrow(edges[source == "cholesterol" & target %chin% c("ad", "exact_splice")]) == 0L
)

blue <- v9_palette[["blue"]]
blue_dark <- v9_palette[["blue_dark"]]
blue_soft <- v9_palette[["blue_soft"]]
teal <- v9_palette[["teal"]]
teal_dark <- v9_palette[["teal_dark"]]
teal_soft <- v9_palette[["teal_soft"]]
amber <- v9_palette[["amber"]]
amber_soft <- v9_palette[["amber_soft"]]
purple <- v9_palette[["purple"]]
purple_soft <- v9_palette[["purple_soft"]]
red <- v9_palette[["red"]]
grey <- v9_palette[["grey"]]
graphite <- v9_palette[["graphite"]]
ink <- v9_palette[["ink"]]
white <- v9_palette[["white"]]

arrow_closed <- grid::arrow(length = grid::unit(1.35, "mm"), type = "closed")

add_rect <- function(p, xmin, xmax, ymin, ymax, fill, colour = NA, linewidth = 0.35,
                     linetype = "solid") {
  p + annotate("rect", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
               fill = fill, colour = colour, linewidth = linewidth, linetype = linetype)
}

add_text <- function(p, x, y, label, size = 1.65, colour = ink, fontface = "plain",
                     hjust = 0.5, vjust = 0.5, lineheight = 0.92) {
  p + annotate("text", x = x, y = y, label = label, size = size, colour = colour,
               family = "Arial", fontface = fontface, hjust = hjust, vjust = vjust,
               lineheight = lineheight)
}

add_arrow <- function(p, x, y, xend, yend, colour = blue, linewidth = 0.65,
                      linetype = "solid") {
  p + annotate("segment", x = x, y = y, xend = xend, yend = yend,
               colour = colour, linewidth = linewidth, linetype = linetype,
               lineend = "round", arrow = arrow_closed)
}

add_curve <- function(p, x, y, xend, yend, curvature = 0.25, colour = blue,
                      linewidth = 0.65, linetype = "solid") {
  p + annotate("curve", x = x, y = y, xend = xend, yend = yend,
               curvature = curvature, colour = colour, linewidth = linewidth,
               linetype = linetype, lineend = "round", arrow = arrow_closed)
}

p <- ggplot() +
  coord_cartesian(xlim = c(0, 120), ylim = c(0, 78), clip = "off", expand = FALSE) +
  theme_void(base_family = "Arial", base_size = 6) +
  theme(plot.margin = margin(2, 2, 2, 2))

# Module backgrounds establish whitespace and prevent connector/text collisions.
p <- add_rect(p, 1.2, 32.0, 45.0, 75.8, "#F3F7FA", "#C7D7E2", 0.42)
p <- add_rect(p, 34.0, 77.0, 45.0, 75.8, "#F7F4FA", "#D8CCE3", 0.42)
p <- add_rect(p, 79.0, 118.8, 45.0, 75.8, "#F2F8F6", "#C8DFD9", 0.42)
p <- add_rect(p, 1.2, 32.0, 18.0, 42.2, "#FAFBFB", "#DCE2E5", 0.42)
p <- add_rect(p, 34.0, 77.0, 18.0, 42.2, "#F3F7FA", "#C7D7E2", 0.42)
p <- add_rect(p, 79.0, 118.8, 18.0, 42.2, "#FFF8EF", "#E9D0A5", 0.42)

# -----------------------------------------------------------------------------
# Regulatory locus: parallel phenotype inputs converge on one local configuration.
# -----------------------------------------------------------------------------
p <- add_text(p, 3.0, 72.7, "REGULATORY LD BLOCK", 2.15, blue_dark, "bold", hjust = 0)
p <- add_rect(p, 3.5, 10.5, 66.0, 69.5, "#F5E7E5", red, 0.42)
p <- add_text(p, 7.0, 67.75, "AD", 1.75, red, "bold")
p <- add_rect(p, 12.0, 29.5, 66.0, 69.5, "#E8F3F1", teal, 0.42)
p <- add_text(p, 20.75, 67.75, "TC | LDL-C | non-HDL-C", 1.43, teal_dark, "bold")
p <- add_arrow(p, 7.0, 65.7, 11.0, 62.0, blue, 0.55)
p <- add_arrow(p, 20.7, 65.7, 18.5, 62.0, blue, 0.55)

p <- p + annotate("segment", x = 4.0, xend = 29.2, y = 59.0, yend = 59.0,
                  colour = graphite, linewidth = 0.72)
variant_x <- c(5.2, 8.6, 12.0, 15.5, 19.0, 22.7, 26.2, 28.5)
p <- p + geom_point(
  data = data.table(x = variant_x, y = 59.0), aes(x, y), inherit.aes = FALSE,
  shape = 21, size = 2.2, fill = white, colour = blue, stroke = 0.5
)
p <- add_rect(p, 9.2, 22.5, 53.2, 56.3, amber_soft, amber, 0.5)
p <- add_text(p, 15.85, 54.75, "Published regulatory\nperturbation anchors", 1.35, amber, "bold")
p <- p + annotate("segment", x = 15.85, xend = 15.85, y = 56.3, yend = 58.7,
                  colour = amber, linewidth = 0.48, linetype = "dashed")
p <- add_text(p, 16.5, 49.0, "Shared local configuration", 1.62, blue_dark, "bold")
p <- add_text(p, 16.5, 46.8, "Parallel associations; no serial mediation inferred", 1.25, graphite)

# Present-study connector occupies the gutter, not either module.
p <- add_arrow(p, 32.25, 60.0, 33.7, 60.0, blue, 0.72)

# -----------------------------------------------------------------------------
# Candidate molecular readout: exact canonical-versus-cryptic acceptor balance.
# -----------------------------------------------------------------------------
p <- add_text(p, 36.0, 72.7, "CANDIDATE MOLECULAR READOUT", 2.15, blue_dark, "bold", hjust = 0)

exons <- data.table(
  xmin = c(37.0, 41.5, 48.0, 54.0), xmax = c(39.5, 45.5, 52.0, 56.5),
  label = c("4", "5", "6", "7"), fill = c("#E1E7EA", blue_dark, "#55A3D2", "#E1E7EA")
)
p <- p + geom_segment(
  data = data.table(x = c(39.5, 45.5, 52.0), xend = c(41.5, 48.0, 54.0), y = 62.5),
  aes(x = x, xend = xend, y = y, yend = y), inherit.aes = FALSE,
  colour = graphite, linewidth = 0.58
)
p <- p + geom_rect(
  data = exons, aes(xmin = xmin, xmax = xmax, ymin = 60.6, ymax = 64.4),
  inherit.aes = FALSE, fill = exons$fill, colour = white, linewidth = 0.35
)
p <- p + geom_text(
  data = exons, aes(x = (xmin + xmax) / 2, y = 62.5, label = label),
  inherit.aes = FALSE, family = "Arial", size = 1.85, fontface = "bold",
  colour = c(ink, white, white, ink)
)
p <- add_curve(p, 43.5, 64.8, 50.0, 64.8, -0.45, blue, 0.78)
p <- add_curve(p, 43.5, 60.2, 46.6, 60.2, 0.58, amber, 0.72)
p <- add_text(p, 47.0, 68.2, "Canonical exon5-6 | AA150/151", 1.42, blue_dark, "bold")
p <- add_text(p, 43.9, 57.7, "Competing cryptic acceptor", 1.30, amber, "bold")

risk_x <- c(38.0, 41.0, 44.0, 47.0, 50.0, 53.0)
p <- p + geom_point(
  data = data.table(x = risk_x, y = 52.8), aes(x, y), inherit.aes = FALSE,
  shape = 21, size = 1.85, fill = white, colour = blue, stroke = 0.45
)
p <- add_text(p, 45.5, 49.2, "Risk haplotype shifts local acceptor balance", 1.38, blue, "bold")
p <- add_text(p, 45.5, 46.9, "Count models resolve canonical versus cryptic reads", 1.20, graphite)

triangulation[, y := c(66.8, 62.8, 58.8, 54.8)]
p <- add_text(p, 59.0, 70.0, "Trait", 1.30, graphite, "bold", hjust = 0)
p <- add_text(p, 68.7, 70.0, "Default", 1.25, graphite, "bold")
p <- add_text(p, 74.0, 70.0, "p12/10", 1.25, graphite, "bold")
p <- p + geom_text(
  data = triangulation, aes(x = 59.0, y = y, label = trait), inherit.aes = FALSE,
  family = "Arial", size = 1.45, hjust = 0, colour = ink
)
p <- p + geom_text(
  data = triangulation, aes(x = 68.7, y = y, label = sprintf("%.3f", pph4_default)),
  inherit.aes = FALSE, family = "Arial", size = 1.35, colour = graphite
)
p <- p + geom_text(
  data = triangulation, aes(x = 74.0, y = y, label = sprintf("%.3f", pph4_conservative)),
  inherit.aes = FALSE, family = "Arial", size = 1.35, fontface = "bold",
  colour = ifelse(triangulation$pph4_conservative >= 0.80, blue, graphite)
)

# -----------------------------------------------------------------------------
# Contextual evidence: cell-resolved molecular phenotypes remain separated.
# -----------------------------------------------------------------------------
p <- add_text(p, 81.0, 72.7, "CELL AND TISSUE CONTEXT", 2.15, teal_dark, "bold", hjust = 0)

angles <- seq(0, 2 * pi, length.out = 9)[-9]
micro <- data.table(
  x = 88.0, y = 63.0,
  xend = 88.0 + 3.0 * cos(angles), yend = 63.0 + 3.0 * sin(angles)
)
p <- p + geom_segment(
  data = micro, aes(x, y, xend = xend, yend = yend), inherit.aes = FALSE,
  colour = "#5B8E52", linewidth = 0.62, lineend = "round"
)
p <- p + annotate("point", x = 88.0, y = 63.0, shape = 21, size = 5.5,
                  fill = "#E8F2E5", colour = "#5B8E52", stroke = 0.65)
p <- add_text(p, 88.0, 57.6, "Isolated microglia", 1.50, "#4E7E49", "bold")
p <- add_text(p, 88.0, 55.2, "exact-event sQTL", 1.28, blue)

p <- p + annotate("point", x = 105.0, y = 63.0, shape = 21, size = 5.4,
                  fill = "#E8F2F8", colour = blue, stroke = 0.65)
p <- p + annotate("segment", x = 105.0, y = 63.0, xend = 111.0, yend = 59.0,
                  colour = blue, linewidth = 0.65)
p <- p + annotate("segment", x = 104.5, y = 64.0, xend = 101.0, yend = 67.0,
                  colour = blue, linewidth = 0.55)
p <- p + annotate("segment", x = 105.5, y = 64.0, xend = 108.7, yend = 67.2,
                  colour = blue, linewidth = 0.55)
p <- add_text(p, 105.0, 57.6, "Neurons / astrocytes", 1.50, blue_dark, "bold")
p <- add_text(p, 105.0, 55.2, "total-expression eQTL", 1.28, teal)

gtex <- context[resource == "GTEx v8"]
gtex[, `:=`(
  x = c(83.5, 91.5, 99.5, 107.5),
  short_label = c("BA24", "Hippocampus", "Putamen", "Spinal cord")
)]
p <- p + geom_point(
  data = gtex, aes(x = x, y = 51.6, size = -log10(p_value), fill = risk_aligned_effect),
  inherit.aes = FALSE, shape = 21, colour = white, stroke = 0.35
) +
  scale_size_continuous(range = c(2.0, 3.5), guide = "none") +
  scale_fill_gradient(low = teal_soft, high = teal_dark, guide = "none")
p <- p + geom_text(
  data = gtex, aes(x = x, y = 49.5, label = short_label), inherit.aes = FALSE,
  family = "Arial", size = 1.08, colour = graphite
)
p <- add_rect(p, 82.0, 116.0, 45.8, 48.1, "#F1F3F4", "#D0D6D9", 0.32)
p <- add_text(p, 99.0, 46.95, "AD-state RNA: no uniform FDR-significant shift", 1.18, graphite, "bold")

# -----------------------------------------------------------------------------
# Evidence grammar remains explicit and separate from the biological drawing.
# -----------------------------------------------------------------------------
p <- add_text(p, 3.0, 39.2, "EVIDENCE GRAMMAR", 1.85, ink, "bold", hjust = 0)
legend_y <- c(34.8, 31.3, 27.8, 24.3, 20.8)
legend_col <- c(blue, teal, purple, amber, graphite)
legend_lty <- c("solid", "solid", "dotted", "dashed", "dashed")
legend_lab <- c(
  "Present-study association",
  "Cross-tissue / cell context",
  "Reference localization",
  "Published perturbation",
  "Prediction requiring validation"
)
for (i in seq_along(legend_y)) {
  p <- p + annotate("segment", x = 3.0, xend = 7.0, y = legend_y[i], yend = legend_y[i],
                    colour = legend_col[i], linewidth = 0.76, linetype = legend_lty[i])
  p <- add_text(p, 8.2, legend_y[i], legend_lab[i], 1.30, legend_col[i],
                ifelse(i == 1, "bold", "plain"), hjust = 0)
}

# -----------------------------------------------------------------------------
# Reference EC2 localization. No protein consequence is asserted.
# -----------------------------------------------------------------------------
p <- add_text(p, 36.0, 39.2, "REFERENCE EC2 LOCALIZATION", 1.85, blue_dark, "bold", hjust = 0)

lipid_x <- seq(37.0, 74.0, by = 2.0)
p <- p + geom_point(
  data = data.table(x = lipid_x, y = 28.2), aes(x, y), inherit.aes = FALSE,
  shape = 21, size = 2.0, fill = "#D7D0AE", colour = "#B9AE7B", stroke = 0.3
)
p <- p + geom_point(
  data = data.table(x = lipid_x, y = 25.7), aes(x, y), inherit.aes = FALSE,
  shape = 21, size = 2.0, fill = "#D7D0AE", colour = "#B9AE7B", stroke = 0.3
)

helix_x <- c(40.0, 44.5, 51.5, 57.0)
p <- p + geom_rect(
  data = data.table(xmin = helix_x, xmax = helix_x + 2.1),
  aes(xmin = xmin, xmax = xmax, ymin = 23.5, ymax = 30.4), inherit.aes = FALSE,
  fill = blue_dark, colour = white, linewidth = 0.3
)
ec2 <- data.table(
  x = seq(53.6, 57.0, length.out = 100),
  y = 30.4 + 5.0 * sin(seq(0, pi, length.out = 100))
)
p <- p + geom_line(data = ec2, aes(x, y), inherit.aes = FALSE, colour = purple, linewidth = 0.95)
p <- p + annotate("point", x = 55.1, y = 34.9, shape = 21, size = 2.8,
                  fill = white, colour = purple, stroke = 0.65)
p <- add_text(p, 55.1, 37.0, "AA150/151 within EC2", 1.43, purple, "bold")
p <- add_text(p, 48.7, 21.3, "TSPAN14 four-pass topology", 1.35, blue_dark, "bold")
p <- add_rect(p, 61.0, 75.0, 32.0, 36.0, white, purple_soft, 0.45, "dotted")
p <- add_text(p, 68.0, 34.0, "Reference boundary\nmean interval pLDDT > 90", 1.25, purple, "bold")
p <- add_arrow(p, 55.5, 44.7, 55.5, 42.4, purple, 0.58, "dotted")

# -----------------------------------------------------------------------------
# Published perturbation context is spatially separated from splice predictions.
# -----------------------------------------------------------------------------
p <- add_text(p, 81.0, 39.2, "PUBLISHED FUNCTIONAL CONTEXT", 1.85, amber, "bold", hjust = 0)
p <- add_rect(p, 82.0, 116.0, 33.0, 36.6, "#FFF3DF", amber, 0.46)
p <- add_text(p, 99.0, 34.8, "rs7922621 editing: TSPAN14 regulation", 1.38, amber, "bold")

p <- p + geom_point(
  data = data.table(x = seq(82.0, 116.0, by = 2.0), y = 26.8), aes(x, y),
  inherit.aes = FALSE, shape = 21, size = 1.8, fill = "#D7D0AE", colour = "#B9AE7B", stroke = 0.25
)
p <- add_rect(p, 89.0, 91.2, 24.0, 30.0, amber, white, 0.3)
p <- p + annotate("point", x = 90.1, y = 31.3, shape = 21, size = 4.3,
                  fill = amber_soft, colour = amber, stroke = 0.65)
p <- add_text(p, 90.1, 21.6, "Cell-surface ADAM10", 1.30, amber, "bold")

p <- add_rect(p, 106.0, 108.0, 24.0, 30.0, "#5B8E52", white, 0.3)
p <- p + annotate("point", x = 107.0, y = 31.3, shape = 21, size = 4.0,
                  fill = "#E7F1E3", colour = "#5B8E52", stroke = 0.65)
p <- add_text(p, 107.0, 21.6, "Soluble TREM2", 1.30, "#4E7E49", "bold")
p <- add_arrow(p, 90.1, 32.8, 90.1, 31.8, amber, 0.60, "dashed")
p <- add_arrow(p, 107.0, 32.8, 107.0, 31.8, amber, 0.60, "dashed")
p <- add_text(p, 99.0, 19.3, "Published editing did not target the splice acceptors", 1.15, graphite)

# -----------------------------------------------------------------------------
# Bottom validation roadmap: predictions are grey dashed and non-overlapping.
# -----------------------------------------------------------------------------
p <- add_text(p, 2.0, 11.0, "EXACT-SPLICE-SPECIFIC\nPREDICTIONS", 1.65, graphite, "bold", hjust = 0)

prediction <- data.table(
  xmin = c(34.0, 62.0, 91.0), xmax = c(59.0, 88.0, 118.8),
  xmid = c(46.5, 75.0, 104.9),
  label = c(
    "TSPAN14 protein-isoform\ncomposition",
    "ADAM10 trafficking and\nsubstrate processing",
    "Lipid-state-dependent\ncellular phenotypes"
  )
)
p <- p + geom_rect(
  data = prediction, aes(xmin = xmin, xmax = xmax, ymin = 3.0, ymax = 11.5),
  inherit.aes = FALSE, fill = "#F2F4F5", colour = grey, linewidth = 0.45,
  linetype = "dashed"
)
p <- p + geom_text(
  data = prediction, aes(x = xmid, y = 7.25, label = label), inherit.aes = FALSE,
  family = "Arial", size = 1.42, fontface = "bold", colour = graphite, lineheight = 0.9
)
p <- add_curve(p, 48.0, 17.8, 46.5, 11.7, 0.08, graphite, 0.58, "dashed")
p <- add_curve(p, 55.5, 17.8, 75.0, 11.7, -0.12, graphite, 0.58, "dashed")
p <- add_curve(p, 63.0, 17.8, 104.9, 11.7, -0.18, graphite, 0.58, "dashed")

save_pub_v9(p, paths$output_dir, "Figure_6_integrated_biological_model_v11", 183, 122)
assert_clean_svg_v9(file.path(paths$output_dir, "Figure_6_integrated_biological_model_v11.svg"))
message("Exported integrated Figure 6 v11 to: ", paths$output_dir)
