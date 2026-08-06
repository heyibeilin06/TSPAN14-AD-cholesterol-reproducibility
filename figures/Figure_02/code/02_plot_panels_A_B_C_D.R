#!/usr/bin/env Rscript

script_arg <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)])
script_dir <- dirname(normalizePath(script_arg, winslash = "/", mustWork = TRUE))
source(file.path(script_dir, "00_figure_style_and_io.R"))

paths <- parse_v9_args()

tracks <- read_v9_source(paths$source_dir, "Figure_2_locus_tracks_grch38.tsv")
markers <- read_v9_source(paths$source_dir, "Figure_2_annotation_markers_grch38.tsv")
scatter <- read_v9_source(paths$source_dir, "Figure_2_colocalization_scatter.tsv")
transcripts <- read_v9_source(paths$source_dir, "Figure_2_gencode_v38_transcripts.tsv")
exons <- read_v9_source(paths$source_dir, "Figure_2_gencode_v38_exons.tsv")
ccre <- read_v9_source(paths$source_dir, "Figure_2_regulatory_elements.tsv")
coloc <- read_v9_source(paths$source_dir, "Figure_2_exact_event_coloc.tsv")
annotation <- read_v9_source(paths$source_dir, "Figure_2_variant_annotation_matrix.tsv")

region <- c(80.45, 80.53)
lead_position <- markers[snp == "rs7080009", position_grch38][1] / 1e6
ld_breaks <- c(-Inf, 0.2, 0.4, 0.6, 0.8, 0.95, Inf)
ld_labels <- c("<0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-0.95", ">=0.95")
ld_colours <- c(
  "<0.2" = "#D6E0E6", "0.2-0.4" = "#8FC2D9", "0.4-0.6" = "#55A982",
  "0.6-0.8" = "#E3A72F", "0.8-0.95" = "#D55B4D", ">=0.95" = "#7650A1"
)

# A. Five aligned regional association tracks.
hero <- tracks[
  track %chin% c("AD", "TC", "LDL", "nonHDL", "BA24 exact exon5-6 sQTL") &
    position_grch38 / 1e6 >= region[1] & position_grch38 / 1e6 <= region[2]
]
hero[, track_label := factor(
  track,
  levels = c("AD", "TC", "LDL", "nonHDL", "BA24 exact exon5-6 sQTL"),
  labels = c("AD GWAS", "Total cholesterol", "LDL cholesterol", "Non-HDL cholesterol", "Exact exon5-6 sQTL")
)]
hero[, ld_bin := cut(r2_to_rs7080009, breaks = ld_breaks, labels = ld_labels, right = FALSE)]
lead <- hero[snp == "rs7080009"]
lead_label <- lead[track == "AD"]
lead_label[, label_y := neg_log10_p - 0.45]

p_a <- ggplot(hero, aes(position_grch38 / 1e6, neg_log10_p)) +
  geom_vline(xintercept = lead_position, linewidth = 0.32, linetype = "dotted", colour = v9_palette[["purple"]]) +
  geom_point(aes(fill = ld_bin), shape = 21, size = 0.84, stroke = 0, alpha = 0.92) +
  geom_point(data = lead, shape = 21, size = 2.25, fill = "white", colour = v9_palette[["purple"]], stroke = 0.68) +
  geom_text(data = lead_label, aes(y = label_y, label = snp), hjust = -0.08, vjust = 1,
            size = 1.55, fontface = "bold", colour = v9_palette[["purple"]]) +
  facet_grid(track_label ~ ., scales = "free_y", switch = "y") +
  scale_fill_manual(values = ld_colours, drop = FALSE, name = expression(EUR~LD~r^2)) +
  scale_x_continuous(
    position = "top", limits = region, breaks = seq(80.46, 80.52, 0.02),
    labels = scales::number_format(accuracy = 0.01), expand = c(0, 0)
  ) +
  labs(tag = "A", x = "Chr10 position (Mb, GRCh38)", y = expression(-log[10](P))) +
  theme_v9(6.15) +
  theme(
    strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1, size = 5.65, face = "bold"),
    panel.spacing.y = grid::unit(0.85, "mm"), legend.position = "bottom", legend.direction = "horizontal",
    legend.text = element_text(size = 4.9), legend.key.width = grid::unit(2.5, "mm"),
    axis.title.x = element_text(margin = margin(b = 1.5)), plot.margin = margin(2, 4, 0, 10)
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 2.2)))

# B. Transcript and regulatory tracks on the same physical coordinate.
transcripts[, `:=`(display_name = sub("\\..*$", "", display_name), is_primary = as.logical(is_primary))]
exons[, is_primary := as.logical(is_primary)]
ccre[, `:=`(
  display_class = "ENCODE cCRE",
  functional_anchor = as.logical(functional_anchor)
)]
selected_transcripts <- transcripts[order(-as.integer(is_primary), -track_y)][1:min(7, .N)]
selected_names <- selected_transcripts$display_name
selected_exons <- exons[transcript %chin% selected_names]
selected_transcripts[, plot_y := seq(.N, 1)]
selected_exons <- merge(selected_exons, selected_transcripts[, .(display_name, plot_y)], by.x = "transcript", by.y = "display_name")
primary_y <- selected_transcripts[is_primary == TRUE, plot_y][1]
transcript_label_x <- max(selected_transcripts$txEnd, na.rm = TRUE) / 1e6 + 0.0022
track_right_limit <- transcript_label_x + 0.011

marker_tracks <- copy(markers)
marker_tracks[, marker_class := fcase(
  grepl("Laub_CRISPRi", snp_role), "CRISPRi anchor",
  snp == "rs7922621", "Prime-editing anchor",
  default = "Project proxy"
)]
marker_tracks[, marker_y := -2.25]

p_b <- ggplot() +
  geom_segment(
    data = selected_transcripts,
    aes(txStart / 1e6, plot_y, xend = txEnd / 1e6, yend = plot_y, colour = is_primary), linewidth = 0.42
  ) +
  geom_rect(
    data = selected_exons,
    aes(xmin = start_mb, xmax = end_mb, ymin = plot_y - 0.19, ymax = plot_y + 0.19, fill = is_primary), colour = NA
  ) +
  geom_text(
    data = selected_transcripts,
    aes(x = transcript_label_x, y = plot_y, label = display_name, colour = is_primary), hjust = 0, size = 1.62,
    fontface = ifelse(selected_transcripts$is_primary, "bold", "plain")
  ) +
  geom_curve(
    aes(x = 80.509471, xend = 80.512144, y = primary_y + 0.26, yend = primary_y + 0.26),
    curvature = -0.58, linewidth = 0.72, colour = v9_palette[["blue"]]
  ) +
  annotate("text", x = 80.5108, y = primary_y + 0.86, label = "exact exon5-6", size = 1.68,
           fontface = "bold", colour = v9_palette[["blue_dark"]]) +
  geom_rect(
    data = ccre,
    aes(xmin = start_mb, xmax = end_mb, ymin = -1.28, ymax = -0.82, fill = display_class), colour = NA
  ) +
  geom_rect(
    data = ccre[functional_anchor == TRUE],
    aes(xmin = start_mb, xmax = end_mb, ymin = -1.46, ymax = -0.64),
    fill = NA, colour = v9_palette[["amber"]], linewidth = 0.50
  ) +
  geom_segment(
    data = marker_tracks,
    aes(x = position_grch38 / 1e6, xend = position_grch38 / 1e6, y = -2.08, yend = -1.60, colour = marker_class),
    linewidth = 0.42
  ) +
  geom_point(
    data = marker_tracks,
    aes(position_grch38 / 1e6, marker_y, shape = marker_class, colour = marker_class), size = 2.0, stroke = 0.52
  ) +
  annotate("text", x = region[1] + 0.001, y = -1.05, label = "ENCODE cCRE", hjust = 0,
           size = 1.62, colour = v9_palette[["graphite"]]) +
  annotate("text", x = region[1] + 0.001, y = -2.25, label = "Functional anchors", hjust = 0,
           size = 1.62, colour = v9_palette[["graphite"]]) +
  scale_colour_manual(values = c(
    `TRUE` = v9_palette[["blue_dark"]], `FALSE` = "#8B989F",
    `CRISPRi anchor` = v9_palette[["amber"]], `Prime-editing anchor` = v9_palette[["red"]],
    `Project proxy` = v9_palette[["blue"]]
  ), guide = "none") +
  scale_fill_manual(values = c(
    `TRUE` = v9_palette[["blue_dark"]], `FALSE` = "#AAB5BA",
    `ENCODE cCRE` = "#B9DCCF"
  ), guide = "none") +
  scale_shape_manual(
    values = c(`CRISPRi anchor` = 8, `Prime-editing anchor` = 23, `Project proxy` = 21),
    name = NULL
  ) +
  scale_x_continuous(limits = c(region[1], track_right_limit), breaks = seq(80.46, 80.52, 0.02),
                     labels = scales::number_format(accuracy = 0.01), expand = c(0, 0)) +
  coord_cartesian(ylim = c(-2.75, max(selected_transcripts$plot_y) + 1.05), clip = "off") +
  labs(tag = "B", x = "Chr10 position (Mb, GRCh38)", y = NULL) +
  theme_v9(5.9) +
  theme(
    axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
    legend.position = "bottom", legend.direction = "horizontal",
    legend.text = element_text(size = 4.8), legend.key.width = grid::unit(2.6, "mm"),
    legend.margin = margin(t = -2, b = -1), plot.margin = margin(1, 3, 1, 5)
  ) +
  guides(shape = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 2.1)))

# C. Compact locus-comparison atlas with high-LD variants emphasized.
comparison_levels <- c("AD", "TC", "LDL", "nonHDL")
comparison_labels <- c("AD", "Total cholesterol", "LDL cholesterol", "Non-HDL cholesterol")
scatter <- scatter[comparison %chin% comparison_levels]
scatter[, comparison_label := factor(comparison, levels = comparison_levels, labels = comparison_labels)]
scatter[, ld_bin := cut(r2_to_rs7080009, breaks = ld_breaks, labels = ld_labels, right = FALSE)]
scatter[, high_ld := r2_to_rs7080009 >= 0.6]
pp_map <- setNames(coloc$pph4, coloc$trait)
pp_labels <- data.table(
  comparison_label = factor(comparison_labels, levels = comparison_labels),
  label = sprintf("Default PP.H4 %.3f", pp_map[c("AD", "TC", "LDL-C", "non-HDL-C")])
)

p_c <- ggplot(scatter, aes(trait_neg_log10_p, sqtl_neg_log10_p)) +
  geom_point(data = scatter[high_ld == FALSE], colour = "#CCD6DC", size = 0.42, alpha = 0.55) +
  geom_point(data = scatter[high_ld == TRUE], aes(fill = ld_bin), shape = 21, size = 0.82, stroke = 0, alpha = 0.94) +
  geom_point(data = scatter[snp == "rs7080009"], shape = 21, size = 2.0, fill = "white",
             colour = v9_palette[["purple"]], stroke = 0.62) +
  geom_text(data = pp_labels, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
            hjust = -0.08, vjust = 1.25, size = 1.55, fontface = "bold", colour = v9_palette[["ink"]]) +
  facet_wrap(~comparison_label, nrow = 2, scales = "free_x") +
  scale_fill_manual(values = ld_colours, guide = "none") +
  labs(tag = "C", x = "Trait -log10(P)", y = "Exact exon5-6 sQTL -log10(P)") +
  theme_v9(5.45) +
  theme(
    strip.text = element_text(size = 5.15, face = "bold"), panel.spacing = grid::unit(2.2, "mm"),
    axis.text = element_text(size = 4.9), plot.margin = margin(2, 2, 2, 2)
  )

# D. Functional evidence matrix for sentinel and proxy variants.
metric_levels <- c("EUR LD r2", "Exact interval", "Within 20 kb", "4-trait CS", "sQTL tissues", "CRISPRi", "Prime editing")
metric_labels <- c("EUR LD\nr2", "Exact\ninterval", "Within\n20 kb", "4-trait\nCS", "sQTL\ntissues", "CRISPRi", "Prime\nediting")
annotation[, metric_factor := factor(metric, levels = metric_levels, labels = metric_labels)]
snp_order <- unique(annotation[order(-functional_prior_score, position_grch38), snp])
annotation[, snp_factor := factor(snp, levels = rev(snp_order))]
annotation[, evidence_class := fcase(
  metric == "CRISPRi", "Published perturbation",
  metric == "Prime editing", "Published editing",
  metric == "EUR LD r2", "LD",
  metric == "sQTL tissues", "Replication",
  default = "Project annotation"
)]

matrix_colours <- c(
  LD = v9_palette[["purple"]], Replication = v9_palette[["teal"]],
  `Project annotation` = v9_palette[["blue"]],
  `Published perturbation` = v9_palette[["amber"]], `Published editing` = v9_palette[["red"]]
)

p_d <- ggplot(annotation, aes(metric_factor, snp_factor)) +
  geom_tile(fill = v9_palette[["pale"]], colour = "white", linewidth = 0.55) +
  geom_point(aes(size = value, fill = evidence_class, alpha = value), shape = 21, colour = "white", stroke = 0.35) +
  scale_size_continuous(range = c(0, 4.4), limits = c(0, 1), guide = "none") +
  scale_alpha_continuous(range = c(0, 1), limits = c(0, 1), guide = "none") +
  scale_fill_manual(values = matrix_colours, guide = "none") +
  labs(tag = "D", x = NULL, y = NULL) +
  theme_v9(5.5) +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_text(size = 4.7, face = "bold", lineheight = 0.88, margin = margin(t = 1.2)),
    axis.text.y = element_text(size = 5.0), plot.margin = margin(2, 1, 2, 2)
  )

figure_2 <- p_a / p_b / (p_c | p_d) +
  plot_layout(heights = c(1.10, 0.96, 1.03), widths = c(1.45, 0.85)) & panel_tag_theme_v9

save_pub_v9(figure_2, paths$output_dir, "Figure_2_variant_level_TSPAN14_locus_v9", 183, 190)
message("Exported Figure 2 v9 to: ", paths$output_dir)
