#!/usr/bin/env Rscript

script_arg <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)])
script_dir <- dirname(normalizePath(script_arg, winslash = "/", mustWork = TRUE))
source(file.path(script_dir, "00_figure_style_and_io.R"))

paths <- parse_v9_args()
atlas <- read_v9_source(paths$source_dir, "Figure_5_cell_context_atlas.tsv")
exact <- read_v9_source(paths$source_dir, "Figure_5_exact_event_cross_context.tsv")
snuc <- read_v9_source(paths$source_dir, "Figure_5_single_nucleus_eqtl.tsv")
disease <- read_v9_source(paths$source_dir, "Figure_5_disease_state_rna.tsv")
structure <- read_v9_source(paths$source_dir, "Figure_5_ec2_structure.tsv")

context_levels <- c("Microglia", "Excitatory neurons", "Astrocytes", "Bulk neural tissue")
layer_levels <- c("Published perturbation", "Adjusted AD-state RNA", "Total-expression eQTL", "Exact exon5-6 sQTL")
class_colours <- c(
  "Exact splice QTL" = v9_palette[["blue"]],
  "Gene-expression QTL" = v9_palette[["teal"]],
  "Disease-state test" = v9_palette[["grey"]],
  "Perturbation" = v9_palette[["amber"]]
)

# A. Evidence atlas: blank cells are deliberately retained to prevent ecological over-attribution.
atlas_grid <- CJ(evidence_layer = layer_levels, context = context_levels, unique = TRUE)
atlas_grid <- merge(atlas_grid, atlas, by = c("evidence_layer", "context"), all.x = TRUE)
atlas_grid[, `:=`(
  evidence_layer = factor(evidence_layer, levels = layer_levels),
  context = factor(context, levels = context_levels)
)]

p_a <- ggplot(atlas_grid, aes(context, evidence_layer)) +
  geom_tile(fill = v9_palette[["pale"]], colour = "white", linewidth = 0.7) +
  geom_point(
    data = atlas_grid[!is.na(evidence_class)],
    aes(size = evidence_strength, fill = evidence_class),
    shape = 21, colour = "white", stroke = 0.5
  ) +
  scale_size_continuous(range = c(3.0, 6.2), guide = "none") +
  scale_fill_manual(values = class_colours, name = NULL) +
  labs(tag = "A", x = NULL, y = NULL) +
  theme_v9(5.35) +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_text(size = 4.75, face = "bold", angle = 25, hjust = 1),
    axis.text.y = element_text(size = 4.9),
    legend.position = "bottom", legend.direction = "horizontal",
    legend.text = element_text(size = 4.25), legend.key.width = grid::unit(2.6, "mm"),
    plot.margin = margin(2, 2, 2, 2)
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 3.2)))

# B. Single-nucleus total-expression eQTL: distinct from the exact splice event.
snuc[, row_label := factor(
  paste(cell_context, study_label, sep = " | "),
  levels = rev(c(
    "Astrocytes | ROSMAP combined",
    "Excitatory neurons | ROSMAP CUIMC",
    "Excitatory neurons | ROSMAP combined",
    "Excitatory neurons | ROSMAP MIT"
  ))
)]
cell_colours <- c("Astrocytes" = v9_palette[["teal"]], "Excitatory neurons" = v9_palette[["blue"]])

p_b <- ggplot(snuc, aes(beta, row_label, colour = cell_context, fill = cell_context)) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.16, linewidth = 0.72) +
  geom_point(shape = 21, size = 3.1, colour = "white", stroke = 0.42) +
  scale_colour_manual(values = cell_colours, guide = "none") +
  scale_fill_manual(values = cell_colours, name = NULL) +
  scale_x_continuous(limits = c(0, 0.5), breaks = c(0, 0.2, 0.4)) +
  labs(tag = "B", x = "Risk-allele effect on TSPAN14 expression", y = NULL) +
  theme_v9(5.25) +
  theme(
    legend.position = "bottom", legend.direction = "horizontal",
    legend.text = element_text(size = 4.45), legend.key.width = grid::unit(3, "mm"),
    axis.text.y = element_text(size = 4.65), plot.margin = margin(2, 2, 2, 2)
  )

# C. Disease-state RNA: adjusted pseudobulk and independent cross-study synthesis.
disease[, row_label := factor(
  paste(cell_label, source_label, sep = " | "),
  levels = rev(c(
    "Microglia | SEA-AD adjusted pseudobulk",
    "Neurons | SEA-AD adjusted pseudobulk",
    "Microglia | 17-study single-cell meta-analysis",
    "Oligodendrocytes | 17-study single-cell meta-analysis"
  ))
)]
disease[, source_short := fifelse(grepl("SEA-AD", source_label), "SEA-AD pseudobulk", "17-study meta-analysis")]
disease_colours <- c("SEA-AD pseudobulk" = v9_palette[["blue"]], "17-study meta-analysis" = v9_palette[["purple"]])

p_c <- ggplot(disease, aes(estimate, row_label, colour = source_short, fill = source_short)) +
  geom_vline(xintercept = 0, linewidth = 0.32, colour = v9_palette[["graphite"]]) +
  geom_errorbar(
    data = disease[interval_available == TRUE],
    aes(xmin = lo, xmax = hi), orientation = "y", width = 0.16, linewidth = 0.72
  ) +
  geom_point(shape = 21, size = 3.15, colour = "white", stroke = 0.48) +
  scale_colour_manual(values = disease_colours, guide = "none") +
  scale_fill_manual(values = disease_colours, name = NULL) +
  scale_x_continuous(limits = c(-0.42, 0.72), breaks = c(-0.4, 0, 0.4)) +
  labs(tag = "C", subtitle = "No estimate passed FDR < 0.05", x = "AD-state RNA effect", y = NULL) +
  theme_v9(5.15) +
  theme(
    plot.subtitle = element_text(size = 4.7, face = "bold", colour = v9_palette[["graphite"]]),
    legend.position = "bottom", legend.box = "vertical", legend.direction = "horizontal",
    legend.text = element_text(size = 4.25), legend.title = element_text(size = 4.45),
    legend.key.width = grid::unit(2.7, "mm"), axis.text.y = element_text(size = 4.45),
    plot.margin = margin(2, 2, 2, 2)
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 3.0)))

# D. Transcript-to-protein coordinate track and residue-level structural confidence.
structure <- structure[residue %between% c(114, 232)]
structure[, sse_short := fifelse(sse_psea == "alpha_helix", "Alpha helix", "Coil")]
exons <- data.table(
  exon = c("Exon 5", "Exon 6", "Exon 7"),
  xmin = c(114, 151, 193), xmax = c(150, 192, 232),
  fill = c(v9_palette[["blue_soft"]], "#72ACCE", v9_palette[["purple_soft"]])
)
sse_colours <- c("Alpha helix" = v9_palette[["amber"]], "Coil" = "#AEB8BD")

p_d <- ggplot(structure, aes(residue, pLDDT)) +
  annotate("rect", xmin = 114, xmax = 232, ymin = 55, ymax = 98.5,
           fill = v9_palette[["teal_soft"]], alpha = 0.35) +
  geom_ribbon(aes(ymin = 55, ymax = pLDDT), fill = v9_palette[["blue_soft"]], alpha = 0.55) +
  geom_line(colour = v9_palette[["blue_dark"]], linewidth = 0.72) +
  geom_rect(data = exons, aes(xmin = xmin, xmax = xmax, ymin = 100, ymax = 104, fill = exon),
            inherit.aes = FALSE, colour = "white", linewidth = 0.45) +
  geom_text(data = exons, aes(x = (xmin + xmax) / 2, y = 102, label = exon),
            inherit.aes = FALSE, size = 2.0, fontface = "bold", colour = v9_palette[["ink"]]) +
  geom_tile(aes(x = residue, y = 57.1, fill = sse_short), width = 1.02, height = 3.3) +
  geom_hline(yintercept = 90, linetype = "dashed", linewidth = 0.34, colour = v9_palette[["graphite"]]) +
  geom_vline(xintercept = 150.5, colour = v9_palette[["blue"]], linewidth = 1.0) +
  geom_point(data = structure[residue == 150], shape = 21, size = 3.2,
             fill = "white", colour = v9_palette[["blue"]], stroke = 0.7) +
  annotate("label", x = 152.2, y = 106.4, label = "Exact exon5-6 | AA150/151",
           hjust = 0, size = 2.1, fontface = "bold", colour = v9_palette[["blue"]],
           fill = "white", linewidth = 0, label.padding = grid::unit(0.7, "mm")) +
  annotate("label", x = 115.5, y = 97.0,
           label = "EC2 / ADAM10-interaction\nannotation (AA114-232)",
           hjust = 0, size = 1.85, lineheight = 0.9, fontface = "bold",
           colour = v9_palette[["teal_dark"]], fill = "#EFF7F5", linewidth = 0,
           label.padding = grid::unit(0.65, "mm")) +
  annotate("label", x = 153.0, y = 92.4, label = "Boundary pLDDT 93.3",
           hjust = 0, size = 1.85, colour = v9_palette[["graphite"]],
           fill = "#EFF7F5", linewidth = 0, label.padding = grid::unit(0.55, "mm")) +
  scale_fill_manual(
    values = c(setNames(exons$fill, exons$exon), sse_colours),
    breaks = names(sse_colours), name = "Secondary structure"
  ) +
  scale_x_continuous(breaks = c(120, 150, 180, 210, 232), expand = expansion(mult = c(0.015, 0.02))) +
  scale_y_continuous(limits = c(54.5, 107.5), breaks = c(60, 75, 90), expand = c(0, 0)) +
  labs(tag = "D", x = "TSPAN14-207 residue", y = "AlphaFold pLDDT") +
  theme_v9(5.4) +
  theme(
    legend.position = "bottom", legend.direction = "horizontal",
    legend.text = element_text(size = 4.6), legend.title = element_text(size = 4.8),
    legend.key.width = grid::unit(4.5, "mm"), plot.margin = margin(2, 2, 2, 2)
  )

# E. Restrained topology localization inset; schematic positions are not residue-scaled.
helix <- data.table(xmin = c(0.9, 1.8, 3.2, 4.4), xmax = c(1.25, 2.15, 3.55, 4.75))
ec2_curve <- data.table(
  x = seq(3.375, 4.575, length.out = 120),
  y = 0.72 + 0.46 * sin(seq(0, pi, length.out = 120))
)
p_e <- ggplot() +
  annotate("rect", xmin = 0.35, xmax = 5.3, ymin = 0.34, ymax = 0.66,
           fill = "#E7E1C9", colour = "#CFC5A0", linewidth = 0.4) +
  geom_rect(data = helix, aes(xmin = xmin, xmax = xmax, ymin = 0.22, ymax = 0.78),
            fill = v9_palette[["blue"]], colour = v9_palette[["blue_dark"]], linewidth = 0.45) +
  geom_line(data = ec2_curve, aes(x, y), colour = v9_palette[["purple"]], linewidth = 1.05) +
  annotate("segment", x = 1.25, xend = 1.8, y = 0.22, yend = 0.22,
           colour = v9_palette[["blue_dark"]], linewidth = 0.65) +
  annotate("segment", x = 2.15, xend = 3.2, y = 0.78, yend = 0.78,
           colour = v9_palette[["blue_dark"]], linewidth = 0.65) +
  annotate("point", x = 3.75, y = 1.08, shape = 21, size = 3.2,
           fill = "white", colour = v9_palette[["purple"]], stroke = 0.75) +
  annotate("segment", x = 3.75, xend = 3.75, y = 1.08, yend = 1.35,
           colour = v9_palette[["purple"]], linewidth = 0.55) +
  annotate("label", x = 3.87, y = 1.47, label = "AA150/151", size = 2.15,
           hjust = 0, fontface = "bold", colour = v9_palette[["purple"]],
           fill = "white", linewidth = 0, label.padding = grid::unit(0.55, "mm")) +
  annotate("label", x = 3.62, y = 1.22, label = "EC2 (AA114-232)", size = 2.0,
           hjust = 1, fontface = "bold", colour = v9_palette[["purple"]],
           fill = "white", linewidth = 0, label.padding = grid::unit(0.55, "mm")) +
  annotate("text", x = c(1.075, 1.975, 3.375, 4.575), y = 0.49,
           label = c("1", "2", "3", "4"), size = 1.9, fontface = "bold", colour = "white") +
  annotate("label", x = 2.55, y = -0.02,
           label = "Reference localization\npLDDT 93.3 | local PAE 0.875",
           size = 1.85, linewidth = 0.25, label.padding = grid::unit(1.2, "mm"),
           fill = v9_palette[["pale"]], colour = v9_palette[["ink"]]) +
  annotate("text", x = 5.2, y = -0.27, label = "Schematic topology; not to scale",
           hjust = 1, size = 1.55, colour = v9_palette[["graphite"]]) +
  coord_cartesian(xlim = c(0.25, 5.35), ylim = c(-0.34, 1.58), clip = "off") +
  labs(tag = "E") +
  theme_void(base_family = "Arial", base_size = 6) +
  theme(plot.margin = margin(2, 2, 2, 2), plot.tag = panel_tag_theme_v9$plot.tag)

top_row <- p_a | p_b | p_c
top_row <- top_row + plot_layout(widths = c(1.05, 1.00, 1.08))
bottom_row <- p_d | p_e
bottom_row <- bottom_row + plot_layout(widths = c(1.55, 0.72))
figure_5 <- top_row / bottom_row + plot_layout(heights = c(0.95, 1.08)) & panel_tag_theme_v9

save_pub_v9(figure_5, paths$output_dir, "Figure_5_cell_context_structure_v9", 183, 175)
message("Exported Figure 5 v9 to: ", paths$output_dir)
