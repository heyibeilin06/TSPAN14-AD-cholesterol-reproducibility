#!/usr/bin/env Rscript

script_arg <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)])
script_dir <- dirname(normalizePath(script_arg, winslash = "/", mustWork = TRUE))
source(file.path(script_dir, "00_figure_style_and_io.R"))

paths <- parse_v9_args()
instrument <- read_v9_source(paths$source_dir, "Figure_4_instrument_effect_atlas.tsv")
cis_mr <- read_v9_source(paths$source_dir, "Figure_4_ld_aware_cis_mr.tsv")
diagnostics <- read_v9_source(paths$source_dir, "Figure_4_cis_mr_diagnostics.tsv")
global <- read_v9_source(paths$source_dir, "Figure_4_genomewide_lipid_to_ad.tsv")
joint_mvmr <- read_v9_source(paths$source_dir, "Figure_4_global_joint_mvmr.tsv")
joint_strength <- read_v9_source(paths$source_dir, "Figure_4_global_joint_mvmr_strength.tsv")
pc <- read_v9_source(paths$source_dir, "Figure_4_pc_gmm_dimension_sensitivity.tsv")

lipid_colours <- c(
  "TC" = v9_palette[["teal"]], "LDL-C" = v9_palette[["blue"]],
  "non-HDL-C" = v9_palette[["purple"]], "HDL-C" = v9_palette[["amber"]],
  "TG" = v9_palette[["red"]]
)

# A. Risk-aligned cis-instrument effect atlas.
layer_levels <- c("Exact sQTL", "AD", "TC", "LDL-C", "non-HDL-C")
snp_levels <- unique(instrument[order(BP), SNP])
instrument[, `:=`(
  layer_factor = factor(display_layer, levels = layer_levels),
  snp_factor = factor(SNP, levels = rev(snp_levels)),
  plot_z = pmax(-10, pmin(10, risk_aligned_z))
)]

p_a <- ggplot(instrument, aes(layer_factor, snp_factor)) +
  geom_tile(fill = v9_palette[["pale"]], colour = "white", linewidth = 0.65) +
  geom_point(aes(size = abs(plot_z), fill = plot_z), shape = 21, colour = "white", stroke = 0.38) +
  scale_size_continuous(range = c(1.7, 6.0), guide = "none") +
  scale_fill_gradient(low = v9_palette[["blue_soft"]], high = v9_palette[["blue_dark"]],
                      limits = c(0, 10), name = "Risk-aligned Z") +
  labs(tag = "A", x = NULL, y = NULL) +
  theme_v9(5.7) +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_text(size = 4.9, face = "bold", angle = 25, hjust = 1),
    axis.text.y = element_text(size = 5.0), legend.position = "bottom",
    legend.key.width = grid::unit(17, "mm"), legend.key.height = grid::unit(2.2, "mm"),
    plot.margin = margin(2, 2, 2, 2)
  ) +
  guides(fill = guide_colourbar(title.position = "left", barwidth = grid::unit(17, "mm"), barheight = grid::unit(2.2, "mm")))

# B. LD-aware cis-MR and lead-instrument estimates.
cis_mr[, `:=`(
  outcome_factor = factor(outcome_label, levels = c("AD", "TC", "LDL-C", "non-HDL-C")),
  method_factor = factor(method_label, levels = c("LD-aware generalized IVW", "Lead-instrument Wald")),
  panel_group = factor(ifelse(outcome_label == "AD", "AD", "Lipids"), levels = c("AD", "Lipids"))
)]
method_colours <- c("LD-aware generalized IVW" = v9_palette[["blue"]], "Lead-instrument Wald" = v9_palette[["grey"]])
format_cis_effect <- function(x) {
  if (max(abs(x), na.rm = TRUE) < 0.02) {
    out <- formatC(x, format = "e", digits = 0)
    out <- gsub("e-0", "e-", out, fixed = TRUE)
    out <- gsub("0e+00", "0", out, fixed = TRUE)
    return(out)
  }
  formatC(x, format = "f", digits = 2)
}

p_b <- ggplot(cis_mr, aes(estimate, outcome_factor, colour = method_factor, fill = method_factor, shape = method_factor)) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.14,
                linewidth = 0.66, position = position_dodge(width = 0.38)) +
  geom_point(size = 2.5, stroke = 0.42, position = position_dodge(width = 0.38)) +
  facet_grid(. ~ panel_group, scales = "free_x") +
  scale_colour_manual(values = method_colours, name = NULL) +
  scale_fill_manual(values = method_colours, name = NULL) +
  scale_shape_manual(values = c("LD-aware generalized IVW" = 21, "Lead-instrument Wald" = 22), name = NULL) +
  expand_limits(x = 0) +
  scale_x_continuous(n.breaks = 3, labels = format_cis_effect) +
  labs(tag = "B", x = "Effect per unit predicted exon5-6 usage", y = NULL) +
  theme_v9(5.25) +
  theme(
    strip.text = element_text(size = 5.2, face = "bold"),
    axis.text.y = element_text(size = 4.7), panel.spacing = grid::unit(2.0, "mm"),
    legend.position = "bottom", legend.direction = "horizontal",
    legend.text = element_text(size = 4.35), legend.key.width = grid::unit(2.4, "mm"),
    plot.margin = margin(2, 2, 2, 2)
  ) +
  guides(
    colour = "none",
    fill = "none",
    shape = guide_legend(
      nrow = 2,
      byrow = TRUE,
      override.aes = list(
        fill = unname(method_colours),
        colour = unname(method_colours),
        size = 2.35,
        stroke = 0.42
      )
    )
  )

# C. Diagnostic matrix: size encodes criterion strength, colour encodes flags.
diagnostic_order <- c("Minimum F", "Lead-SNP Steiger", "Default PP.H4", "Cochran Q P", "Egger intercept P")
diagnostics[, score := fcase(
  diagnostic == "Minimum instrument F", pmin(value / 10, 1),
  diagnostic == "Steiger direction", value,
  diagnostic == "Exact coloc PP.H4", value,
  diagnostic %chin% c("Cochran Q P", "MR-Egger intercept P"), pmin(value / 0.05, 1),
  default = 0
)]
diagnostics[diagnostic == "Exact coloc PP.H4", diagnostic_label := "Default PP.H4"]
diagnostics[, `:=`(
  diagnostic_factor = factor(diagnostic_label, levels = diagnostic_order),
  outcome_factor = factor(outcome_label, levels = rev(c("AD", "TC", "LDL-C", "non-HDL-C"))),
  status_factor = factor(ifelse(pass, "Not flagged", "Flagged"), levels = c("Not flagged", "Flagged"))
)]

p_c <- ggplot(diagnostics, aes(diagnostic_factor, outcome_factor)) +
  geom_tile(fill = v9_palette[["pale"]], colour = "white", linewidth = 0.62) +
  geom_point(aes(size = score, fill = status_factor), shape = 21, colour = "white", stroke = 0.38) +
  scale_size_continuous(range = c(1.3, 5.0), limits = c(0, 1), guide = "none") +
  scale_fill_manual(values = c("Not flagged" = v9_palette[["teal"]], "Flagged" = v9_palette[["red"]]), name = NULL) +
  labs(tag = "C", x = NULL, y = NULL) +
  theme_v9(5.35) +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_text(size = 4.45, face = "bold", angle = 35, hjust = 1),
    axis.text.y = element_text(size = 4.9), legend.position = "bottom",
    legend.text = element_text(size = 4.6), legend.key.width = grid::unit(2.5, "mm"),
    plot.margin = margin(2, 2, 2, 2)
  ) +
  guides(fill = guide_legend(override.aes = list(size = 3.0)))

# D. Genome-wide lipid-to-AD MR with extended-APOE sensitivity.
global[, exposure_label := factor(
  exposure,
  levels = c("TG", "nonHDL", "LDL", "TC", "HDL"),
  labels = c("TG", "non-HDL-C", "LDL-C", "TC", "HDL-C")
)]
global[, model_factor := factor(model_label, levels = c("Genome-wide", "Extended-APOE excluded"))]
model_colours <- c("Genome-wide" = v9_palette[["graphite"]], "Extended-APOE excluded" = v9_palette[["amber"]])

p_d <- ggplot(global, aes(estimate, exposure_label, colour = model_factor, group = model_factor)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.36, colour = v9_palette[["grey"]]) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.15,
                position = position_dodge(width = 0.38), linewidth = 0.66) +
  geom_point(aes(fill = model_factor), shape = 21, size = 2.7, colour = "white", stroke = 0.40,
             position = position_dodge(width = 0.38)) +
  scale_colour_manual(values = model_colours, guide = "none") +
  scale_fill_manual(values = model_colours, name = NULL) +
  labs(tag = "D", x = "Genome-wide lipid effect on AD", y = NULL) +
  theme_v9(5.65) +
  theme(
    legend.position = "bottom", legend.text = element_text(size = 4.75),
    legend.key.width = grid::unit(2.8, "mm"), plot.margin = margin(2, 2, 2, 2)
  )

# E. Jointly clumped global MVMR estimates and conditional strength.
joint_mvmr[, exposure_label := factor(exposure, levels = c("TG", "LDL", "HDL"), labels = c("TG", "LDL-C", "HDL-C"))]
p_e1 <- ggplot(joint_mvmr, aes(direct_estimate, exposure_label, colour = exposure_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.36, colour = v9_palette[["grey"]]) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.14, linewidth = 0.70) +
  geom_point(aes(fill = exposure_label), shape = 21, size = 2.9, colour = "white", stroke = 0.42) +
  scale_colour_manual(values = lipid_colours, guide = "none") +
  scale_fill_manual(values = lipid_colours, guide = "none") +
  annotate("text", x = -0.121, y = 3.45, label = "Joint MVMR direct effects", hjust = 0,
           size = 1.72, fontface = "bold", colour = v9_palette[["ink"]]) +
  scale_x_continuous(limits = c(-0.125, 0.12), breaks = c(-0.10, -0.05, 0, 0.05, 0.10)) +
  labs(tag = "E", x = "Direct effect on AD", y = NULL) +
  theme_v9(5.25) +
  theme(plot.margin = margin(2, 2, 2, 2))

joint_strength[, exposure_label := factor(exposure, levels = c("TG", "LDL", "HDL"), labels = c("TG", "LDL-C", "HDL-C"))]
p_e2 <- ggplot(joint_strength, aes(assumed_pairwise_sampling_error_correlation, conditional_F,
                                  colour = exposure_label, group = exposure_label)) +
  geom_hline(yintercept = 10, linetype = "dashed", linewidth = 0.38, colour = v9_palette[["grey"]]) +
  geom_line(linewidth = 0.72) +
  geom_point(aes(fill = exposure_label), shape = 21, size = 2.2, colour = "white", stroke = 0.35) +
  scale_colour_manual(values = lipid_colours, name = NULL) +
  scale_fill_manual(values = lipid_colours, guide = "none") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75)) +
  labs(x = "Assumed sampling-error correlation", y = "Conditional F") +
  theme_v9(5.15) +
  theme(
    legend.position = "bottom", legend.direction = "horizontal",
    legend.text = element_text(size = 4.55), legend.key.width = grid::unit(2.5, "mm"),
    plot.margin = margin(2, 2, 2, 2)
  )
p_e <- p_e1 | p_e2

# F. PC-GMM lipid-coefficient attenuation and identification stability.
pc[, lipid_label := factor(lipid, levels = c("TC", "LDL", "nonHDL"), labels = c("TC", "LDL-C", "non-HDL-C"))]
pc[, strength_status := factor(ifelse(all_strength_F_ge_10, "Both F >= 10", "At least one F < 10"),
                               levels = c("Both F >= 10", "At least one F < 10"))]
pc_colours <- c("TC" = v9_palette[["teal"]], "LDL-C" = v9_palette[["blue"]], "non-HDL-C" = v9_palette[["purple"]])

p_f1 <- ggplot(pc, aes(n_pcs, attenuation_fraction, colour = lipid_label, group = lipid_label)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = 1, fill = v9_palette[["pale"]], colour = NA) +
  geom_hline(yintercept = 0, linewidth = 0.36, colour = v9_palette[["graphite"]]) +
  geom_line(linewidth = 0.70, alpha = 0.90) +
  geom_point(aes(shape = strength_status, fill = lipid_label), size = 2.3, stroke = 0.42) +
  scale_colour_manual(values = pc_colours, name = NULL) +
  scale_fill_manual(values = pc_colours, guide = "none") +
  scale_shape_manual(values = c("Both F >= 10" = 16, "At least one F < 10" = 1), name = NULL) +
  scale_x_continuous(breaks = c(3, 5, 10, 15, 20, 25, 30)) +
  scale_y_continuous(limits = c(-1.1, 0.82), breaks = c(-1, -0.5, 0, 0.5)) +
  labs(tag = "F", title = "Lipid-effect attenuation", x = "Retained LD principal components", y = "Lipid-coefficient attenuation\nafter splice adjustment") +
  theme_v9(5.45) +
  theme(
    plot.title = element_text(size = 5.45, face = "bold", hjust = 0),
    legend.position = "bottom", legend.box = "horizontal", legend.text = element_text(size = 4.55),
    legend.key.width = grid::unit(2.6, "mm"), plot.margin = margin(2, 2, 2, 2)
  )

p_f2 <- ggplot(pc, aes(n_pcs, minimum_conditional_F, colour = lipid_label, group = lipid_label)) +
  geom_hline(yintercept = 10, linetype = "dashed", linewidth = 0.45, colour = v9_palette[["red"]]) +
  geom_line(linewidth = 0.70, alpha = 0.90) +
  geom_point(aes(fill = lipid_label), shape = 21, size = 2.15, colour = "white", stroke = 0.36) +
  scale_colour_manual(values = pc_colours, guide = "none") +
  scale_fill_manual(values = pc_colours, guide = "none") +
  scale_x_continuous(breaks = c(3, 5, 10, 15, 20, 25, 30)) +
  labs(title = "Identification strength", x = "Retained LD principal components", y = "Minimum conditional F") +
  theme_v9(5.45) +
  theme(plot.title = element_text(size = 5.45, face = "bold", hjust = 0), plot.margin = margin(2, 2, 2, 2))
p_f <- p_f1 | p_f2

figure_4 <- (p_a | p_b | p_c) / (p_d | p_e) / p_f +
  plot_layout(heights = c(1.00, 0.92, 1.02), widths = c(1.0, 1.15, 0.9)) & panel_tag_theme_v9

save_pub_v9(figure_4, paths$output_dir, "Figure_4_causal_scope_v9", 183, 205)
message("Exported Figure 4 v9 to: ", paths$output_dir)
