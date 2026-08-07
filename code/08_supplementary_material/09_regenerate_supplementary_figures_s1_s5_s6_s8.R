#!/usr/bin/env Rscript
# Figure contract: S1 tests whether the AD-HDL-C correlation persists after
# extended-APOE conditioning; S5 and S6 distinguish exact-event cross-tissue
# consistency from adjacent-junction co-usage; S8 defines the causal scope.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(patchwork)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(ifelse(length(args), args[1], "."), winslash = "/")
out <- file.path(root, "outputs", "supplement_v19", "figures")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

blue <- "#1769AA"; teal <- "#128C7E"; orange <- "#D98200"
purple <- "#7A5195"; red <- "#C44E52"; grey <- "#65737D"
theme_pub <- theme_classic(base_family = "Arial", base_size = 8.2) +
  theme(
    plot.title = element_text(size = 9.2, face = "bold", colour = "#172B3A", margin = margin(b = 5, l = 14)),
    axis.title = element_text(size = 8.2),
    axis.text = element_text(size = 7.2, colour = "#263640"),
    legend.title = element_text(size = 7.5, face = "bold"),
    legend.text = element_text(size = 7),
    legend.position = "bottom",
    panel.grid.major.y = element_line(colour = "#E7ECEF", linewidth = .25),
    plot.margin = margin(8, 10, 8, 12)
  )
tag_theme <- theme(
  plot.tag = element_text(family = "Arial", face = "bold", size = 12),
  plot.tag.position = c(.01, .99)
)

save_fig <- function(plot, number, width = 190, height = 140) {
  stem <- file.path(out, paste0("Supplementary_Figure_S", number))
  ggsave(paste0(stem, ".png"), plot, width = width, height = height, units = "mm", dpi = 300, bg = "white")
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, units = "mm", device = cairo_pdf, bg = "white")
  ggsave(paste0(stem, ".svg"), plot, width = width, height = height, units = "mm", bg = "white")
  ggsave(paste0(stem, ".tiff"), plot, width = width, height = height, units = "mm", dpi = 600, compression = "lzw", bg = "white")
}

# S1: baseline correlation and extended-APOE sensitivity.
ld_all <- read_tsv(file.path(root, "figures", "Figure_01", "data", "Figure_1_ldsc_apoe_conditioning.tsv"), show_col_types = FALSE)
ld <- ld_all %>%
  filter(model == "baseline") %>%
  distinct(trait_label, rg, se, p, lo, hi) %>%
  mutate(trait_label = factor(trait_label, levels = rev(c("HDL-C", "LDL-C", "TG", "TC", "non-HDL-C"))))
p1a <- ggplot(ld, aes(rg, trait_label)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#9AA5AC", linewidth = .4) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = .16, colour = grey, linewidth = .55) +
  geom_point(aes(colour = p < .05), size = 2.6) +
  scale_colour_manual(values = c(`TRUE` = orange, `FALSE` = grey), guide = "none") +
  labs(x = "Genetic correlation (95% CI)", y = NULL, title = "Baseline AD-lipid genetic correlation") +
  theme_pub

apoe <- ld_all %>%
  filter(trait == "HDL", model %in% c("own_lead", "pair_union_leads", "w5Mb")) %>%
  mutate(
    label = recode(model,
      own_lead = "LD-conditioned: trait lead",
      pair_union_leads = "LD-conditioned: union leads",
      w5Mb = "Physical-window sensitivity: 5 Mb"
    ),
    analysis = ifelse(model == "w5Mb", "Physical-window sensitivity", "LD-conditioned"),
    label = factor(label, levels = rev(c(
      "LD-conditioned: trait lead", "LD-conditioned: union leads",
      "Physical-window sensitivity: 5 Mb"
    )))
  )
p1b <- ggplot(apoe, aes(rg, label, colour = analysis)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#9AA5AC", linewidth = .4) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = .14, linewidth = .55) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = c("LD-conditioned" = blue, "Physical-window sensitivity" = teal)) +
  labs(x = "AD-HDL-C genetic correlation (95% CI)", y = NULL, title = "Extended-APOE sensitivity", colour = NULL) +
  theme_pub + theme(legend.justification = "center")
save_fig((p1a | p1b) + plot_annotation(tag_levels = "A") & tag_theme, 1, 190, 105)

# S5: coordinate-identical event across partially overlapping neural tissues.
rep <- read_tsv(file.path(root, "figures", "Figure_03", "data", "Figure_3_exact_replication_matrix.tsv"), show_col_types = FALSE) %>%
  mutate(signed_logp = -log10(p_value) * sign(nes))
p5a <- ggplot(rep, aes(tissue_label, reorder(snp_label, position), fill = signed_logp)) +
  geom_tile(colour = "white", linewidth = .4) +
  scale_fill_gradient2(low = purple, mid = "white", high = orange, midpoint = 0, name = "Risk-aligned\nsigned -log10(P)") +
  labs(x = NULL, y = NULL, title = "Coordinate-identical exon5-6 sQTL") +
  theme_pub + theme(axis.text.x = element_text(angle = 28, hjust = 1))
p5b <- rep %>%
  group_by(tissue_label) %>%
  summarise(median_nes = median(nes), min_p = min(p_value), .groups = "drop") %>%
  ggplot(aes(median_nes, reorder(tissue_label, median_nes))) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#A8B0B5", linewidth = .4) +
  geom_point(aes(size = -log10(min_p)), colour = blue) +
  scale_size(range = c(2.2, 4.8), breaks = c(10, 15, 20)) +
  labs(x = "Median risk-aligned NES", y = NULL, size = "Strongest\n-log10(P)", title = "Cross-tissue direction") +
  theme_pub + theme(legend.key.width = grid::unit(7, "mm"))
save_fig((p5a | p5b) + plot_layout(widths = c(1.28, .72)) + plot_annotation(tag_levels = "A") & tag_theme, 5, 190, 110)

# S6: aggregate adjacent-junction co-usage, not transcript-level PSI.
co <- read_tsv(file.path(root, "figures", "Figure_03", "data", "Figure_3_brain_cousage_summary.tsv"), show_col_types = FALSE) %>%
  mutate(tissue_label = ifelse(is.na(tissue_label), tissue, tissue_label))
p6a <- ggplot(co, aes(spearman_rho, reorder(tissue_label, spearman_rho), size = n_both_nonzero, colour = spearman_rho)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#A8B0B5", linewidth = .4) +
  geom_point(alpha = .9) +
  scale_colour_gradient(low = "#BCD8E8", high = teal) +
  scale_size(range = c(2, 6)) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = "Spearman correlation", y = NULL, size = "Samples with\nboth junctions", colour = "rho", title = "Exon5-6 and exon6-7 co-usage") +
  theme_pub
p6b <- co %>%
  select(tissue_label, n_samples, n_donors, n_both_nonzero) %>%
  pivot_longer(-tissue_label, names_to = "metric", values_to = "count") %>%
  mutate(metric = recode(metric, n_samples = "Samples", n_donors = "Donors", n_both_nonzero = "Both junctions")) %>%
  ggplot(aes(count, reorder(tissue_label, count), colour = metric)) +
  geom_point(size = 2) +
  facet_wrap(~metric, scales = "free_x", nrow = 1) +
  scale_x_continuous(breaks = c(150, 200, 250)) +
  scale_colour_manual(values = c("Samples" = blue, "Donors" = orange, "Both junctions" = teal), guide = "none") +
  labs(x = "Count", y = NULL, title = "Coverage supporting each correlation") +
  theme_pub + theme(
    strip.background = element_blank(), strip.text = element_text(face = "bold"),
    panel.spacing.x = grid::unit(6, "mm")
  )
save_fig((p6a / p6b) + plot_layout(heights = c(1.15, .85)) + plot_annotation(tag_levels = "A") & tag_theme, 6, 190, 145)

# S8: bidirectional systemic MR and local PC-GMM identification diagnostics.
gm <- read_tsv(file.path(root, "figures", "Figure_04", "data", "Figure_4_genomewide_lipid_to_ad.tsv"), show_col_types = FALSE) %>%
  filter(analysis == "genomewide_bidirectional") %>%
  mutate(pair = paste(exposure, "to", outcome), nominal = pvalue < .05)
p8a <- ggplot(gm, aes(estimate, reorder(pair, estimate))) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#A8B0B5", linewidth = .4) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = .15, colour = grey, linewidth = .55) +
  geom_point(aes(colour = nominal), size = 2.3) +
  scale_colour_manual(values = c(`TRUE` = red, `FALSE` = grey), guide = "none") +
  labs(x = "Genome-wide MR estimate (95% CI)", y = NULL, title = "Bidirectional systemic effects") +
  theme_pub

pc <- read_tsv(file.path(root, "figures", "Figure_04", "data", "Figure_4_pc_gmm_dimension_sensitivity.tsv"), show_col_types = FALSE) %>%
  mutate(lipid_label = recode(lipid, TC = "TC", LDL = "LDL-C", nonHDL = "non-HDL-C"))
lipid_colours <- c("TC" = orange, "LDL-C" = blue, "non-HDL-C" = teal)
p8b <- ggplot(pc, aes(n_pcs, indirect_estimate, colour = lipid_label)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "#A8B0B5", linewidth = .4) +
  geom_line(linewidth = .65, na.rm = TRUE) +
  geom_point(aes(shape = all_strength_F_ge_10), size = 1.9, na.rm = TRUE) +
  scale_colour_manual(values = lipid_colours) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "All strength criteria met", `FALSE` = "At least one criterion <10")) +
  labs(x = "Retained LD principal components", y = "Estimated indirect product", colour = NULL, shape = NULL, title = "PC-GMM dimension sensitivity") +
  guides(colour = guide_legend(order = 1, nrow = 1), shape = guide_legend(order = 2, nrow = 1)) +
  theme_pub + theme(legend.box = "vertical")
p8c <- ggplot(pc, aes(n_pcs, minimum_conditional_F, colour = lipid_label)) +
  geom_hline(yintercept = 10, linetype = 2, colour = red, linewidth = .4) +
  geom_line(linewidth = .65, na.rm = TRUE) +
  geom_point(size = 1.8, na.rm = TRUE) +
  scale_colour_manual(values = lipid_colours, guide = "none") +
  labs(x = "Retained LD principal components", y = "Minimum conditional F", title = "Identification strength") +
  theme_pub
save_fig(((p8a | p8b) / p8c) + plot_layout(widths = c(.9, 1.1), heights = c(1.15, .85)) + plot_annotation(tag_levels = "A") & tag_theme, 8, 190, 165)

cat("Updated Supplementary Figures S1, S5, S6 and S8 in", out, "\n")
