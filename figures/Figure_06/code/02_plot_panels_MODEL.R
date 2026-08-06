#!/usr/bin/env Rscript

script_arg <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)])
script_dir <- dirname(normalizePath(script_arg, winslash = "/", mustWork = TRUE))
source(file.path(script_dir, "00_figure_style_and_io.R"))

paths <- parse_v9_args()
replication <- read_v9_source(paths$source_dir, "Figure_5_exact_event_cross_context.tsv")
triangulation <- read_v9_source(
  paths$source_dir,
  "Figure_5_exact_event_triangulation.tsv"
)
structure <- read_v9_source(paths$source_dir, "Figure_6_structure_summary.tsv")
edges <- read_v9_source(paths$source_dir, "Figure_6_edges.tsv")

stopifnot(
  uniqueN(replication[resource == "GTEx v8", context_label]) == 4L,
  nrow(triangulation) == 4L,
  min(triangulation$pph4) > 0.95,
  min(triangulation$pph4_conservative) > 0.69,
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
pale <- v9_palette[["pale"]]
ink <- v9_palette[["ink"]]
white <- v9_palette[["white"]]

arrow_closed <- grid::arrow(length = grid::unit(1.55, "mm"), type = "closed")

add_rect <- function(p, xmin, xmax, ymin, ymax, fill, colour = NA, linewidth = 0.4,
                     linetype = "solid") {
  p + annotate("rect", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
               fill = fill, colour = colour, linewidth = linewidth, linetype = linetype)
}
add_text <- function(p, x, y, label, size = 1.9, colour = ink, fontface = "plain",
                     hjust = 0.5, vjust = 0.5, lineheight = 0.92) {
  p + annotate("text", x = x, y = y, label = label, size = size, colour = colour,
               fontface = fontface, family = "Arial", hjust = hjust, vjust = vjust,
               lineheight = lineheight)
}
add_label <- function(p, x, y, label, size = 1.8, colour = ink, fill = white,
                      fontface = "plain", hjust = 0.5, lineheight = 0.92) {
  p + annotate("label", x = x, y = y, label = label, size = size, colour = colour,
               fill = fill, linewidth = 0, label.padding = grid::unit(0.55, "mm"),
               fontface = fontface, family = "Arial", hjust = hjust, lineheight = lineheight)
}
add_arrow <- function(p, x, y, xend, yend, colour = blue, linewidth = 0.7,
                      linetype = "solid") {
  p + annotate("segment", x = x, y = y, xend = xend, yend = yend,
               colour = colour, linewidth = linewidth, linetype = linetype,
               lineend = "round", arrow = arrow_closed)
}
add_curve <- function(p, x, y, xend, yend, curvature = 0.2, colour = blue,
                      linewidth = 0.7, linetype = "solid") {
  p + annotate("curve", x = x, y = y, xend = xend, yend = yend,
               curvature = curvature, colour = colour, linewidth = linewidth,
               linetype = linetype, arrow = arrow_closed)
}

p <- ggplot() +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 70), clip = "off", expand = FALSE) +
  theme_void(base_family = "Arial", base_size = 6) +
  theme(plot.margin = margin(2, 2, 2, 2))

# Quiet bands organize a single continuous figure without panel boundaries.
p <- add_rect(p, 1.2, 98.8, 46.0, 68.5, "#FAFBFB", "#E0E5E7", 0.35)
p <- add_rect(p, 1.2, 98.8, 14.0, 44.3, white, NA)
p <- add_rect(p, 1.2, 98.8, 1.2, 12.6, "#FBFCFC", "#E0E5E7", 0.35)

# --- Evidence foundation: locus, exact event and replication. ---
p <- add_text(p, 3.0, 66.7, "REGULATORY LOCUS", 2.25, blue_dark, "bold", hjust = 0)
p <- add_text(p, 34.5, 66.7, "EXACT-EVENT TRIANGULATION", 2.25, blue_dark, "bold", hjust = 0)
p <- add_text(p, 70.0, 66.7, "CROSS-TISSUE NEURAL CONSISTENCY", 2.25, blue_dark, "bold", hjust = 0)

# Compact chromosome and local regulatory block.
p <- add_rect(p, 3.2, 4.8, 52.0, 62.5, "#B9C1C5", NA)
p <- add_rect(p, 3.2, 4.8, 56.5, 58.2, blue, NA)
p <- add_text(p, 5.5, 63.2, "Chr10 q23.1", 1.48, graphite, hjust = 0)
p <- p + annotate("segment", x = 7.0, xend = 30.5, y = 58.0, yend = 58.0,
                  colour = grey, linewidth = 0.75)
variant_x <- c(8.0, 11.0, 14.2, 18.0, 22.0, 26.0, 29.0)
p <- p + geom_point(data = data.table(x = variant_x, y = 58.0), aes(x, y),
                    shape = 21, size = 2.45, fill = white, colour = blue, stroke = 0.52)
p <- add_rect(p, 14.5, 23.5, 52.7, 54.8, amber_soft, amber, 0.5)
p <- add_text(p, 19.0, 53.75, "CRISPRi enhancer", 1.43, amber, "bold")
p <- p + annotate("segment", x = 19.0, xend = 19.0, y = 54.8, yend = 57.7,
                  colour = amber, linewidth = 0.48)
p <- add_text(p, 10.2, 60.7, "high-LD regulatory variants", 1.45, graphite, hjust = 0)

# Parallel input tags have enough width to contain their labels.
p <- add_rect(p, 5.7, 14.0, 47.6, 50.8, "#F5E7E5", red, 0.4)
p <- add_text(p, 9.85, 49.2, "AD", 1.75, red, "bold")
p <- add_rect(p, 15.0, 30.8, 47.6, 50.8, "#E8F3F1", teal, 0.4)
p <- add_text(p, 22.9, 49.2, "TC · LDL-C · non-HDL-C", 1.55, teal_dark, "bold")
p <- add_text(p, 18.2, 51.7, "parallel local associations", 1.42, blue, "bold")

# Exact transcript event and compact triangulation matrix.
exon_top <- data.table(
  xmin = c(35.2, 39.6, 44.8, 50.2), xmax = c(38.0, 43.2, 48.6, 53.0),
  label = c("4", "5", "6", "7"), fill = c("#E3E8EA", blue_dark, "#55A3D2", "#E3E8EA")
)
p <- p + geom_segment(data = data.table(x = c(38.0, 43.2, 48.6), xend = c(39.6, 44.8, 50.2), y = 58.0),
                      aes(x = x, xend = xend, y = y, yend = y), colour = graphite, linewidth = 0.62)
p <- p + geom_rect(data = exon_top, aes(xmin = xmin, xmax = xmax, ymin = 56.0, ymax = 60.0, fill = fill),
                   inherit.aes = FALSE, colour = white, linewidth = 0.4, show.legend = FALSE) +
  scale_fill_identity()
p <- p + geom_text(data = exon_top, aes(x = (xmin + xmax) / 2, y = 58.0, label = label),
                   inherit.aes = FALSE, family = "Arial", size = 2.0, fontface = "bold",
                   colour = c(ink, white, white, ink))
p <- p + annotate("curve", x = 41.4, y = 60.5, xend = 46.7, yend = 60.5,
                  curvature = -0.42, colour = blue, linewidth = 0.78, arrow = arrow_closed)
p <- add_label(p, 44.0, 63.2, "Exact exon5-6 | AA150/151", 1.62, blue, white, "bold")

triangulation[, y := c(60.8, 57.8, 54.8, 51.8)]
p <- add_text(p, 54.5, 63.3, "Trait", 1.5, graphite, "bold", hjust = 0)
p <- add_text(p, 62.0, 63.3, "Default", 1.42, graphite, "bold")
p <- add_text(p, 67.0, 63.3, "p12/10", 1.42, graphite, "bold")
p <- p + geom_text(data = triangulation, aes(x = 54.5, y = y, label = trait),
                   inherit.aes = FALSE, family = "Arial", size = 1.7, colour = ink, hjust = 0)
p <- p + geom_text(data = triangulation, aes(x = 62.0, y = y, label = sprintf("%.3f", pph4)),
                   inherit.aes = FALSE, family = "Arial", size = 1.55, colour = graphite)
p <- p + geom_text(data = triangulation, aes(x = 67.0, y = y, label = sprintf("%.3f", pph4_conservative)),
                   inherit.aes = FALSE, family = "Arial", size = 1.55,
                   fontface = "bold", colour = ifelse(triangulation$pph4_conservative >= 0.80, blue, graphite))

# Four neural tissues with mini NES tracks.
rep_plot <- replication[resource == "GTEx v8"]
rep_plot[, y := c(61.0, 57.8, 54.6, 51.4)[match(context_label,
  c("Anterior cingulate BA24", "Hippocampus", "Putamen", "Cervical spinal cord"))]]
rep_plot[, xend := 86.0 + 8.0 * risk_aligned_effect / 1.25]
p <- p + geom_point(data = rep_plot, aes(x = 72.0, y = y), shape = 21, size = 2.2,
                    fill = teal, colour = white, stroke = 0.35, inherit.aes = FALSE)
p <- p + geom_text(data = rep_plot, aes(x = 73.5, y = y, label = context_label),
                   family = "Arial", size = 1.62, colour = ink, hjust = 0, inherit.aes = FALSE)
p <- p + geom_segment(data = rep_plot, aes(x = 86.0, xend = xend, y = y, yend = y),
                      colour = teal, linewidth = 1.05, lineend = "round", inherit.aes = FALSE)
p <- p + geom_point(data = rep_plot, aes(x = xend, y = y), colour = teal, size = 1.8, inherit.aes = FALSE)
p <- add_text(p, 95.5, 48.7, "risk-aligned NES", 1.45, graphite, hjust = 1)

# --- Continuous molecular and cell-context mechanism. ---
p <- add_text(p, 3.0, 42.5, "STUDY-SUPPORTED MOLECULAR ROUTE", 2.25, blue_dark, "bold", hjust = 0)

# Transcript module: pre-mRNA, risk variants and altered processing.
p <- add_rect(p, 3.0, 34.5, 16.0, 40.5, "#F5F1F8", "#D4C7E2", 0.45)
p <- add_text(p, 5.0, 38.6, "TSPAN14 pre-mRNA", 1.95, purple, "bold", hjust = 0)
exon_mid <- data.table(
  xmin = c(6.0, 11.0, 17.0, 23.5), xmax = c(9.0, 15.0, 21.2, 27.0),
  label = c("4", "5", "6", "7"), fill = c("#E3E8EA", blue_dark, "#55A3D2", "#E3E8EA")
)
p <- p + geom_segment(data = data.table(x = c(9.0, 15.0, 21.2), xend = c(11.0, 17.0, 23.5), y = 34.0),
                      aes(x = x, xend = xend, y = y, yend = y), colour = graphite, linewidth = 0.62)
p <- p + geom_rect(data = exon_mid, aes(xmin = xmin, xmax = xmax, ymin = 31.8, ymax = 36.2, fill = fill),
                   inherit.aes = FALSE, colour = white, linewidth = 0.4, show.legend = FALSE)
p <- p + geom_text(data = exon_mid, aes(x = (xmin + xmax) / 2, y = 34.0, label = label),
                   inherit.aes = FALSE, family = "Arial", size = 2.05, fontface = "bold",
                   colour = c(ink, white, white, ink))
p <- p + annotate("curve", x = 13.0, y = 36.7, xend = 19.1, yend = 36.7,
                  curvature = -0.42, colour = blue, linewidth = 0.82, arrow = arrow_closed)
p <- add_label(p, 30.0, 38.7, "exact exon5-6 sQTL", 1.50, blue, white, "bold", hjust = 1)

risk_x <- c(10.2, 12.8, 15.4, 18.0, 20.6, 23.2)
p <- p + geom_point(data = data.table(x = risk_x, y = 28.7), aes(x, y), shape = 21,
                    size = 2.25, fill = white, colour = blue, stroke = 0.5)
p <- p + geom_segment(data = data.table(x = risk_x, xend = 16.1, y = 29.8, yend = 31.5),
                      aes(x = x, y = y, xend = xend, yend = yend), colour = blue,
                      linewidth = 0.42, linetype = "dashed")
p <- add_text(p, 16.7, 27.1, "risk-aligned linked variants", 1.52, graphite)

# Two compact transcript strips show ratio remodeling without claiming a protein product.
for (yy in c(22.7, 19.2)) {
  p <- p + annotate("segment", x = 8.0, xend = 27.0, y = yy, yend = yy,
                    colour = grey, linewidth = 0.55)
  p <- p + annotate("rect", xmin = 8.0, xmax = 10.0, ymin = yy - 1.15, ymax = yy + 1.15,
                    fill = "#E3E8EA", colour = white, linewidth = 0.35)
  p <- p + annotate("rect", xmin = 12.0, xmax = 15.0, ymin = yy - 1.15, ymax = yy + 1.15,
                    fill = blue_dark, colour = white, linewidth = 0.35)
  p <- p + annotate("rect", xmin = 17.0, xmax = ifelse(yy > 21, 21.0, 19.5),
                    ymin = yy - 1.15, ymax = yy + 1.15,
                    fill = "#55A3D2", colour = white, linewidth = 0.35)
  p <- p + annotate("rect", xmin = 24.0, xmax = 27.0, ymin = yy - 1.15, ymax = yy + 1.15,
                    fill = "#E3E8EA", colour = white, linewidth = 0.35)
}
p <- add_text(p, 16.7, 16.9, "altered relative transcript processing", 1.42, blue, "bold")

# Cell-context module. Direct labels replace ambiguous connector arrows.
p <- add_text(p, 37.0, 39.2, "Cell-context-dependent regulation", 1.95, blue_dark, "bold", hjust = 0)
p <- add_text(p, 37.0, 36.8, "Genotype-resolved molecular contexts", 1.45, graphite, hjust = 0)

angles <- seq(0, 2 * pi, length.out = 9)[-9]
micro <- data.table(x = 41.5, y = 29.2,
                    xend = 41.5 + 2.8 * cos(angles), yend = 29.2 + 2.8 * sin(angles))
p <- p + geom_segment(data = micro, aes(x, y, xend = xend, yend = yend),
                      colour = "#5B8E52", linewidth = 0.65, lineend = "round") +
  annotate("point", x = 41.5, y = 29.2, shape = 21, size = 5.7,
           fill = "#E8F2E5", colour = "#5B8E52", stroke = 0.7) +
  annotate("point", x = 41.5, y = 29.2, size = 1.5, colour = "#5B8E52")
p <- add_text(p, 41.5, 24.8, "Microglia", 1.72, "#4E7E49", "bold")
p <- add_text(p, 41.5, 22.9, "exact sQTL", 1.45, blue)

p <- p + annotate("point", x = 54.0, y = 29.2, shape = 21, size = 5.7,
                  fill = "#E8F2F8", colour = blue, stroke = 0.7) +
  annotate("segment", x = 54.0, y = 29.2, xend = 59.2, yend = 25.5,
           colour = blue, linewidth = 0.7) +
  annotate("segment", x = 53.5, y = 30.2, xend = 50.7, yend = 33.2,
           colour = blue, linewidth = 0.58) +
  annotate("segment", x = 54.5, y = 30.1, xend = 57.0, yend = 33.3,
           colour = blue, linewidth = 0.58)
p <- add_text(p, 53.5, 22.8, "Neurons / astrocytes", 1.65, blue_dark, "bold")
p <- add_text(p, 53.5, 20.9, "total-expression eQTL", 1.42, teal)

p <- add_rect(p, 37.0, 59.5, 16.2, 19.0, "#F1F3F4", "#C8CED1", 0.38)
p <- add_text(p, 48.25, 17.6, "AD-state RNA: no uniform FDR-significant shift", 1.38, graphite, "bold")

# Protein / membrane scene.
p <- add_text(p, 63.0, 42.4, "EC2-CENTRED WORKING MODEL", 2.2, blue_dark, "bold", hjust = 0)
p <- add_rect(p, 62.0, 97.5, 25.3, 28.2, "#E8E1C7", "#CFC49B", 0.4)

helix_x <- c(65.0, 68.0, 72.0, 76.0)
p <- p + geom_rect(data = data.table(xmin = helix_x, xmax = helix_x + 1.45),
                   aes(xmin = xmin, xmax = xmax, ymin = 23.5, ymax = 30.0),
                   inherit.aes = FALSE, fill = blue_dark, colour = white, linewidth = 0.3)
ec2 <- data.table(x = seq(72.72, 76.72, length.out = 100),
                  y = 30.0 + 4.7 * sin(seq(0, pi, length.out = 100)))
p <- p + geom_line(data = ec2, aes(x, y), colour = purple, linewidth = 1.0)
p <- p + annotate("point", x = 74.0, y = 34.0, shape = 21, size = 3.0,
                  fill = white, colour = purple, stroke = 0.7)
p <- add_label(p, 74.0, 38.0, "AA150/151 in EC2", 1.62, purple, white, "bold")
p <- add_text(p, 70.7, 21.8, "TSPAN14", 1.85, blue_dark, "bold")

p <- add_rect(p, 83.0, 85.0, 23.4, 30.1, amber, white, 0.35)
p <- p + annotate("point", x = 84.0, y = 32.3, shape = 21, size = 5.5,
                  fill = amber_soft, colour = amber, stroke = 0.7)
p <- add_text(p, 82.5, 35.9, "ADAM10", 1.72, amber, "bold", hjust = 1)
p <- add_rect(p, 91.5, 93.3, 23.5, 30.0, "#5B8E52", white, 0.35)
p <- p + annotate("point", x = 92.4, y = 32.2, shape = 21, size = 5.0,
                  fill = "#E7F1E3", colour = "#5B8E52", stroke = 0.7)
p <- add_text(p, 94.0, 35.9, "TREM2", 1.72, "#4E7E49", "bold", hjust = 0)

# Published perturbation evidence forms a local module above ADAM10 and TREM2.
# Three short lines keep the annotation inside the box at final print size.
p <- add_rect(p, 77.8, 98.0, 37.6, 41.2, "#FFF8EF", amber, 0.42)
p <- add_text(p, 87.9, 39.4,
              "Published microglial rs7922621 editing\nTSPAN14 regulation\ncell-surface ADAM10 | soluble TREM2",
              1.32, amber, "bold", lineheight = 0.86)
p <- add_arrow(p, 84.0, 37.6, 84.0, 34.9, amber, 0.64, "dashed")
p <- add_arrow(p, 92.4, 37.6, 92.4, 34.7, amber, 0.64, "dashed")

# --- Bottom evidence grammar and validation roadmap. ---
p <- add_text(p, 3.0, 11.0, "EVIDENCE GRAMMAR", 1.8, ink, "bold", hjust = 0)
key_y <- c(8.7, 6.5, 4.3, 2.3)
key_col <- c(blue, teal, amber, graphite)
key_lty <- c("solid", "solid", "dashed", "dashed")
key_lab <- c("Present study", "Cross-tissue / cell context",
             "Published perturbation", "Exact-splice-specific prediction")
for (i in seq_along(key_y)) {
  p <- p + annotate("segment", x = 3.0, xend = 7.0, y = key_y[i], yend = key_y[i],
                    colour = key_col[i], linewidth = 0.78, linetype = key_lty[i])
  p <- add_text(p, 8.0, key_y[i], key_lab[i], 1.42, key_col[i],
                ifelse(i == 1, "bold", "plain"), hjust = 0)
}

pred <- data.table(
  xmin = c(43.0, 61.5, 80.0), xmax = c(59.0, 78.5, 97.0),
  label = c("TSPAN14 protein-isoform\ncomposition",
            "ADAM10 trafficking /\nAPP, TREM2, Notch processing",
            "Lipid-state-dependent\ncellular phenotypes")
)
p <- p + geom_rect(data = pred, aes(xmin = xmin, xmax = xmax, ymin = 2.4, ymax = 10.7),
                   inherit.aes = FALSE, fill = "#F2F4F5", colour = grey,
                   linewidth = 0.45, linetype = "dashed")
p <- p + geom_text(data = pred, aes(x = (xmin + xmax) / 2, y = 6.4, label = label),
                   inherit.aes = FALSE, family = "Arial", size = 1.52,
                   colour = graphite, fontface = "bold", lineheight = 0.9)
p <- add_arrow(p, 70.8, 21.5, 51.0, 10.7, graphite, 0.62, "dashed")
p <- add_arrow(p, 84.0, 23.0, 70.0, 10.7, graphite, 0.62, "dashed")
p <- add_arrow(p, 92.4, 23.0, 88.5, 10.7, graphite, 0.62, "dashed")

save_pub_v9(p, paths$output_dir, "Figure_6_integrated_biological_model_v10", 183, 128)
message("Exported integrated Figure 6 v10 to: ", paths$output_dir)
