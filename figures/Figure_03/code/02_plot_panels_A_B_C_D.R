#!/usr/bin/env Rscript

script_arg <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)])
script_dir <- dirname(normalizePath(script_arg, winslash = "/", mustWork = TRUE))
source(file.path(script_dir, "00_figure_style_and_io.R"))

paths <- parse_v9_args()
set.seed(20260723)
replication <- read_v9_source(paths$source_dir, "Figure_3_exact_replication_matrix.tsv")
donor <- read_v9_source(paths$source_dir, "Figure_3_ba24_canonical_cryptic_counts.tsv")
models <- read_v9_source(paths$source_dir, "Figure_3_canonical_cryptic_models.tsv")
alignment <- read_v9_source(paths$source_dir, "Figure_3_alignment_metric_audit.tsv")
counts <- read_v9_source(paths$source_dir, "Figure_3_brain_junction_counts.tsv")
cousage <- read_v9_source(paths$source_dir, "Figure_3_brain_cousage_summary.tsv")
exons <- read_v9_source(paths$source_dir, "Figure_3_primary_exons.tsv")
events <- read_v9_source(paths$source_dir, "Figure_3_splice_events.tsv")

# A. Cross-tissue exact-event consistency matrix.
tissue_order <- c("BA24", "Hippocampus", "Putamen", "Spinal cord C1")
snp_order <- unique(replication[order(position), snp_label])
replication[, `:=`(
  tissue_factor = factor(tissue_label, levels = tissue_order),
  snp_factor = factor(snp_label, levels = rev(snp_order))
)]

p_a <- ggplot(replication, aes(tissue_factor, snp_factor)) +
  geom_tile(fill = v9_palette[["pale"]], colour = "white", linewidth = 0.7) +
  geom_point(aes(size = neg_log10_p, fill = nes), shape = 21, colour = "white", stroke = 0.42) +
  scale_size_continuous(range = c(2.2, 7.1), name = expression(-log[10](P))) +
  scale_fill_gradient(low = v9_palette[["blue_soft"]], high = v9_palette[["blue_dark"]], name = "NES") +
  labs(tag = "A", x = NULL, y = NULL) +
  theme_v9(6.0) +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_text(size = 5.2, face = "bold", margin = margin(t = 1.5)),
    axis.text.y = element_text(size = 5.1), legend.position = "bottom",
    legend.box = "horizontal", legend.key.width = grid::unit(3.0, "mm"),
    legend.text = element_text(size = 4.7), legend.title = element_text(size = 4.9),
    plot.margin = margin(2, 2, 2, 2)
  ) +
  guides(
    size = guide_legend(order = 1, title.position = "left", override.aes = list(fill = v9_palette[["blue"]])),
    fill = guide_colourbar(order = 2, title.position = "left", barwidth = grid::unit(16, "mm"), barheight = grid::unit(2.3, "mm"))
  )

# B. Donor-level local acceptor balance and complementary count models.
genotype_colours <- c("0 copies" = v9_palette[["blue_soft"]], "1 copy" = v9_palette[["blue"]], "2 copies" = v9_palette[["blue_dark"]])
donor[, `:=`(
  genotype_label = factor(genotype, levels = 0:2, labels = c("0 copies", "1 copy", "2 copies")),
  genotype_factor = factor(genotype, levels = 0:2, labels = c("0 copies", "1 copy", "2 copies")),
  canonical_percent = 100 * canonical_fraction,
  x_genotype = genotype + 1
)]
group_n <- donor[, .N, by = genotype_label]
group_n[, x := match(genotype_label, c("0 copies", "1 copy", "2 copies"))]

p_b1 <- ggplot(donor, aes(x_genotype, canonical_percent, colour = genotype_factor)) +
  geom_jitter(width = 0.10, height = 0, size = 0.82, alpha = 0.55) +
  stat_summary(aes(group = genotype_factor), fun = median, geom = "crossbar", width = 0.42,
               linewidth = 0.52, colour = v9_palette[["ink"]]) +
  geom_text(data = group_n, aes(x = x, y = 82.4, label = paste0("n=", N)), inherit.aes = FALSE,
            size = 1.65, colour = v9_palette[["graphite"]]) +
  annotate("text", x = 0.92, y = 101.0, label = "Local acceptor balance", hjust = 0,
           size = 1.68, fontface = "bold", colour = v9_palette[["ink"]]) +
  annotate("text", x = 2.55, y = 100.45, label = "106 C/C donors: 100%", hjust = 0.5,
           size = 1.35, colour = v9_palette[["blue_dark"]]) +
  scale_colour_manual(values = genotype_colours, guide = "none") +
  scale_x_continuous(breaks = 1:3, labels = c("0 copies", "1 copy", "2 copies"), limits = c(0.55, 3.45)) +
  scale_y_continuous(limits = c(82, 101.2), breaks = c(85, 90, 95, 100), labels = function(x) paste0(x, "%")) +
  labs(tag = "B", x = "rs7080009 AD-risk C-allele dosage", y = "Canonical read fraction") +
  theme_v9(5.7) +
  theme(axis.text.x = element_text(size = 5.0), plot.margin = margin(2, 2, 2, 2))

models <- models[grepl("Firth|Beta-binomial", analysis)]
models[, model_label := factor(fcase(
  grepl("Depth-adjusted", analysis), "Depth-adjusted Firth",
  grepl("Beta-binomial", analysis), "Beta-binomial counts",
  default = "Firth donor detection"
), levels = c("Beta-binomial counts", "Depth-adjusted Firth", "Firth donor detection"))]
models[, model_colour := fcase(grepl("Beta-binomial", analysis), "Beta-binomial", grepl("Depth-adjusted", analysis), "Depth-adjusted", default = "Firth")]

p_b2 <- ggplot(models, aes(odds_ratio, model_label, colour = model_colour)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.4, colour = v9_palette[["grey"]]) +
  geom_errorbar(aes(xmin = odds_ratio_95ci_low, xmax = odds_ratio_95ci_high), orientation = "y", width = 0.18, linewidth = 0.64) +
  geom_point(size = 2.6) +
  annotate("text", x = 0.00105, y = 3.48, label = "Count and detection models", hjust = 0,
           size = 1.82, fontface = "bold", colour = v9_palette[["ink"]]) +
  scale_colour_manual(values = c("Firth" = v9_palette[["blue"]], "Depth-adjusted" = v9_palette[["teal"]], "Beta-binomial" = v9_palette[["amber"]]), guide = "none") +
  scale_x_log10(limits = c(0.001, 1.2), breaks = c(0.001, 0.01, 0.1, 1)) +
  labs(x = "OR per rs7080009 risk allele (95% CI)", y = NULL) +
  theme_v9(5.45) +
  theme(axis.text.y = element_text(size = 4.55), plot.margin = margin(2, 2, 2, 2))

p_b <- p_b1 | p_b2

# C. Junction co-usage in 2,642 GTEx brain samples and across brain regions.
global_rho <- cor(counts$log_ex5_6, counts$log_ex6_7, method = "spearman")
p_c1 <- ggplot(counts, aes(log_ex5_6, log_ex6_7)) +
  geom_hex(bins = 36, aes(fill = after_stat(count)), colour = NA) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.42, linetype = "dashed", colour = v9_palette[["graphite"]]) +
  annotate("text", x = 0.25, y = 5.55, label = sprintf("Spearman rho = %.3f\nn = %s", global_rho, scales::comma(nrow(counts))),
           hjust = 0, vjust = 1, size = 1.75, fontface = "bold", colour = v9_palette[["ink"]]) +
  scale_fill_gradientn(colours = c(v9_palette[["pale"]], v9_palette[["teal_soft"]], v9_palette[["teal"]], v9_palette[["teal_dark"]]),
                       trans = "sqrt", name = "Samples") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_fixed(xlim = c(0, 5.9), ylim = c(0, 5.9), ratio = 1) +
  labs(tag = "C", x = "log(1 + exon5-6 reads)", y = "log(1 + exon6-7 reads)") +
  theme_v9(5.8) +
  theme(
    legend.position = "bottom", plot.margin = margin(2, 2, 2, 2)
  ) +
  guides(fill = guide_colourbar(
    title.position = "left", barwidth = grid::unit(24, "mm"),
    barheight = grid::unit(2.3, "mm"), ticks.colour = "white", frame.colour = NA
  ))

cousage[, tissue_factor := factor(tissue_label, levels = tissue_label[order(spearman_rho)])]
p_c2 <- ggplot(cousage, aes(spearman_rho, tissue_factor)) +
  geom_segment(aes(x = 0.70, xend = spearman_rho, yend = tissue_factor), linewidth = 0.52, colour = v9_palette[["teal_soft"]]) +
  geom_point(aes(size = n_samples, fill = spearman_rho), shape = 21, colour = "white", stroke = 0.4) +
  scale_size_continuous(range = c(2.1, 4.4), guide = "none") +
  scale_fill_gradient(low = v9_palette[["teal_soft"]], high = v9_palette[["teal_dark"]], guide = "none") +
  scale_x_continuous(limits = c(0.70, 0.94), breaks = c(0.7, 0.8, 0.9), expand = c(0, 0)) +
  labs(x = "Spearman rho", y = NULL) +
  theme_v9(5.35) +
  theme(axis.text.y = element_text(size = 4.55), plot.margin = margin(2, 2, 2, 2))

p_c <- p_c1 | p_c2

# D. Coordinate-accurate splice graph.
exons[, `:=`(xmin = start_grch38 / 1e6, xmax = end_grch38 / 1e6)]
event_colours <- c(exact = v9_palette[["blue"]], competing = v9_palette[["amber"]], adjacent = v9_palette[["teal"]])
variant_x <- 80509705 / 1e6

p_d <- ggplot() +
  geom_segment(aes(x = min(exons$xmin), xend = max(exons$xmax), y = 0, yend = 0),
               linewidth = 0.55, colour = v9_palette[["graphite"]]) +
  geom_rect(data = exons, aes(xmin = xmin, xmax = xmax, ymin = -0.22, ymax = 0.22),
            fill = v9_palette[["blue_dark"]], colour = "white", linewidth = 0.45) +
  geom_text(data = exons, aes(x = (xmin + xmax) / 2, y = -0.39, label = paste0("Exon ", exon_rank)),
            colour = v9_palette[["ink"]], size = 1.65, fontface = "bold") +
  geom_curve(data = events[event_class == "exact"],
             aes(x = start_grch38 / 1e6, xend = end_grch38 / 1e6, y = 0.26, yend = 0.26, colour = event_class),
             curvature = -0.46, linewidth = 1.0) +
  geom_curve(data = events[event_class == "competing"],
             aes(x = start_grch38 / 1e6, xend = end_grch38 / 1e6, y = -0.25, yend = -0.25, colour = event_class),
             curvature = 0.55, linewidth = 0.9) +
  geom_curve(data = events[event_class == "adjacent"],
             aes(x = start_grch38 / 1e6, xend = end_grch38 / 1e6, y = 0.25, yend = 0.25, colour = event_class),
             curvature = -0.46, linewidth = 0.9) +
  geom_vline(xintercept = variant_x, linetype = "dotted", linewidth = 0.48, colour = v9_palette[["purple"]]) +
  annotate("label", x = variant_x + 0.00006, y = 0.42, label = "rs7080009", hjust = 0,
           size = 1.55, fontface = "bold", colour = v9_palette[["purple"]], fill = "white",
           linewidth = 0, label.padding = grid::unit(0.35, "mm")) +
  annotate("text", x = mean(events[event_class == "exact", c(start_grch38, end_grch38)]) / 1e6,
           y = 1.15, label = "Canonical exon5-6 | AA150/151\ncross-tissue-consistent sQTL; BA24 counts modelled",
           size = 1.8, fontface = "bold", colour = v9_palette[["blue_dark"]]) +
  annotate("text", x = mean(events[event_class == "competing", c(start_grch38, end_grch38)]) / 1e6,
           y = -1.02, label = "Competing cryptic acceptor 1", size = 1.65,
           fontface = "bold", colour = v9_palette[["amber"]]) +
  annotate("text", x = mean(events[event_class == "adjacent", c(start_grch38, end_grch38)]) / 1e6,
           y = 0.92, label = "Adjacent exon6-7 | AA192/193\nco-used structural context",
           size = 1.65, colour = v9_palette[["teal_dark"]]) +
  scale_colour_manual(values = event_colours, guide = "none") +
  scale_x_continuous(breaks = c(80.5095, 80.5105, 80.5115, 80.5125, 80.5135),
                     labels = scales::number_format(accuracy = 0.0001), expand = expansion(mult = c(0.03, 0.03))) +
  coord_cartesian(ylim = c(-1.2, 1.32), clip = "off") +
  labs(tag = "D", x = "Chr10 position (Mb, GRCh38)", y = NULL) +
  theme_v9(5.8) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
        plot.margin = margin(4, 4, 2, 4))

figure_3 <- (p_a | p_b) / p_c / p_d +
  plot_layout(heights = c(1.02, 1.00, 0.78), widths = c(0.82, 1.45)) & panel_tag_theme_v9

save_pub_v9(figure_3, paths$output_dir, "Figure_3_exact_splicing_consistency_v10", 183, 190)
message("Exported Figure 3 v9 to: ", paths$output_dir)
