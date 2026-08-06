#!/usr/bin/env Rscript
# Figure contract: S4 tests genotype-dependent canonical-cryptic read choice;
# S7 tests whether the local directional estimate survives plausible LD handling.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(patchwork)
  library(scales)
})

root <- normalizePath(commandArgs(trailingOnly = TRUE)[1], winslash = "/")
out <- file.path(root, "outputs", "supplement_v19", "figures")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

blue <- "#1769A6"; teal <- "#148A7B"; orange <- "#D88700"; grey <- "#65737E"; light <- "#DCE7ED"
theme_pub <- theme_classic(base_size = 7.2, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = .35), axis.ticks = element_line(linewidth = .35),
        plot.title = element_text(face = "bold", size = 8.2), strip.text = element_text(face = "bold"),
        legend.title = element_text(size = 6.7), legend.text = element_text(size = 6.5),
        plot.margin = margin(7, 7, 7, 7))
tag_theme <- theme(plot.tag = element_text(family = "Arial", face = "bold", size = 11),
                   plot.tag.position = c(.01, .99))

save_fig <- function(p, n, w = 183, h = 138) {
  stem <- file.path(out, paste0("Supplementary_Figure_S", n))
  ggsave(paste0(stem, ".png"), p, width = w, height = h, units = "mm", dpi = 300, bg = "white")
  ggsave(paste0(stem, ".pdf"), p, width = w, height = h, units = "mm", device = cairo_pdf, bg = "white")
  ggsave(paste0(stem, ".svg"), p, width = w, height = h, units = "mm", bg = "white")
  ggsave(paste0(stem, ".tiff"), p, width = w, height = h, units = "mm", dpi = 600, compression = "lzw", bg = "white")
}

# S4: quantitative-grid validation of the canonical-cryptic acceptor balance.
donor <- read_tsv(file.path(root, "audit", "reviewer_revision", "count_level_acceptor_choice_donors.tsv"), show_col_types = FALSE) %>%
  mutate(genotype = factor(genotype, levels = 0:2, labels = c("0 copies", "1 copy", "2 copies")),
         canonical_percent = 100 * canonical_fraction)
models <- read_tsv(file.path(root, "audit", "reviewer_revision", "count_level_acceptor_choice_models.tsv"), show_col_types = FALSE) %>%
  filter(grepl("Firth|Beta-binomial", analysis)) %>%
  mutate(label = case_when(grepl("Depth-adjusted", analysis) ~ "Depth-adjusted Firth detection",
                           grepl("Beta-binomial", analysis) ~ "Beta-binomial read count",
                           TRUE ~ "Firth donor detection"),
         label = factor(label, levels = c("Beta-binomial read count", "Depth-adjusted Firth detection", "Firth donor detection")))
depth <- read_tsv(file.path(root, "audit", "reviewer_revision", "count_level_acceptor_choice_depth_sensitivity.tsv"), show_col_types = FALSE)

p4a <- ggplot(donor, aes(genotype, canonical_percent, colour = genotype)) +
  geom_jitter(width = .13, height = 0, size = 1.05, alpha = .48) +
  stat_summary(fun = median, geom = "crossbar", width = .48, linewidth = .55, colour = "#243746") +
  scale_colour_manual(values = c("0 copies" = light, "1 copy" = blue, "2 copies" = teal), guide = "none") +
  coord_cartesian(ylim = c(82, 100.5)) +
  labs(x = "rs7080009 AD-risk C-allele dosage", y = "Canonical read fraction (%)", title = "Donor-level local splice choice") + theme_pub

p4b <- ggplot(models, aes(odds_ratio, label)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "#AAB3B9", linewidth = .4) +
  geom_errorbarh(aes(xmin = odds_ratio_95ci_low, xmax = odds_ratio_95ci_high), height = .16, colour = grey, linewidth = .55) +
  geom_point(aes(colour = label), size = 2.5) +
  scale_x_log10(breaks = c(.001, .01, .1, 1), labels = label_number()) +
  scale_colour_manual(values = c("Beta-binomial read count" = orange, "Depth-adjusted Firth detection" = teal, "Firth donor detection" = blue), guide = "none") +
  labs(x = "OR per rs7080009 AD-risk allele (95% CI)", y = NULL, title = "Count and detection models") + theme_pub

p4c <- ggplot(depth, aes(minimum_cluster_depth, odds_ratio_per_risk_allele)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "#AAB3B9", linewidth = .4) +
  geom_line(colour = teal, linewidth = .7) + geom_point(aes(size = donors), colour = teal, alpha = .9) +
  scale_y_log10(breaks = c(.001, .01, .1, 1), labels = label_number()) +
  scale_size(range = c(2, 4.5)) +
  labs(x = "Minimum local cluster depth", y = "Firth OR per risk allele", size = "Donors", title = "Read-depth threshold sensitivity") + theme_pub

save_fig(((p4a | p4b) / p4c) + plot_annotation(tag_levels = "A") & tag_theme, 4, 183, 142)

# S7: asymmetric robustness display for rebuilt ancestry-matched LD.
mr <- read_tsv(file.path(root, "audit", "reviewer_revision", "cis_mr_numerical_sensitivity.tsv"), show_col_types = FALSE)
base <- mr %>% filter(analysis == "ridge", suppressWarnings(as.numeric(parameter)) == 0) %>%
  mutate(lo = estimate - 1.96 * se, hi = estimate + 1.96 * se,
         outcome = factor(outcome, levels = rev(c("AD", "TC", "LDL", "nonHDL"))))
pert <- read_tsv(file.path(root, "audit", "reviewer_revision", "cis_mr_ld_perturbation_sensitivity.tsv"), show_col_types = FALSE) %>%
  mutate(outcome = factor(outcome, levels = rev(c("AD", "TC", "LDL", "nonHDL"))), noise = factor(off_diagonal_noise_sd))
diag <- read_tsv(file.path(root, "audit", "reviewer_revision", "cis_mr_ld_diagnostics.tsv"), show_col_types = FALSE) %>%
  mutate(metric = recode(metric, minimum_eigenvalue = "Minimum eigenvalue", maximum_eigenvalue = "Maximum eigenvalue",
                         condition_number = "Condition number", effective_rank_entropy = "Entropy effective rank",
                         effective_rank_eigenvalue_gt_0.1 = "Eigenvalues > 0.1"),
         metric = factor(metric, levels = rev(c("Minimum eigenvalue", "Maximum eigenvalue", "Condition number", "Entropy effective rank", "Eigenvalues > 0.1"))))

p7a <- ggplot(base, aes(estimate, outcome)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#AAB3B9", linewidth = .4) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = .16, colour = grey, linewidth = .55) +
  geom_point(colour = blue, size = 2.6) + facet_wrap(~outcome, scales = "free_x", nrow = 1) +
  scale_x_continuous(n.breaks = 3, labels = label_number(accuracy = .001), guide = guide_axis(check.overlap = TRUE)) +
  labs(x = "Risk-aligned directional estimate (95% CI)", y = NULL, title = "Rebuilt-LD local estimates") + theme_pub +
  theme(strip.background = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_text(size = 5.7))

p7b <- ggplot(pert, aes(estimate_median, outcome, colour = noise)) +
  geom_errorbarh(aes(xmin = estimate_p025, xmax = estimate_p975), position = position_dodge(width = .46), height = .16, linewidth = .55) +
  geom_point(position = position_dodge(width = .46), size = 2.1) + facet_wrap(~outcome, scales = "free_x", nrow = 1) +
  scale_colour_manual(values = c("0.02" = teal, "0.05" = orange), labels = c("SD 0.02", "SD 0.05")) +
  labs(x = "Estimate under LD perturbation (2.5th-97.5th percentile)", y = NULL, colour = "Off-diagonal noise", title = "LD perturbation sensitivity") + theme_pub +
  theme(strip.background = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "top")

p7c <- ggplot(diag, aes(value, metric)) +
  geom_segment(aes(x = 0, xend = value, yend = metric), colour = "#CBD4D9", linewidth = .7) +
  geom_point(colour = orange, size = 2.5) +
  geom_text(aes(label = formatC(value, digits = 3, format = "fg")), hjust = -.25, size = 2.5, family = "Arial") +
  scale_x_continuous(expand = expansion(mult = c(.02, .2))) +
  labs(x = "Diagnostic value", y = NULL, title = "Local LD matrix diagnostics") + theme_pub

save_fig((p7a / p7b / p7c) + plot_layout(heights = c(1, 1, .8)) + plot_annotation(tag_levels = "A") & tag_theme, 7, 183, 176)

# S9: separated cell-context, disease-state and structural evidence classes.
atlas <- read_tsv(file.path(root, "figures", "Figure_05", "data", "Figure_5_cell_context_atlas.tsv"), show_col_types = FALSE) %>%
  mutate(score = as.numeric(evidence_strength))
dis <- read_tsv(file.path(root, "figures", "Figure_05", "data", "Figure_5_disease_state_rna.tsv"), show_col_types = FALSE) %>%
  mutate(label = paste(source_label, cell_label, sep = " | "))
st <- read_tsv(file.path(root, "figures", "Figure_05", "data", "Figure_5_ec2_structure.tsv"), show_col_types = FALSE)

p9a <- ggplot(atlas, aes(evidence_layer, reorder(context, score), size = score, colour = evidence_class)) +
  geom_point(alpha = .9) + scale_size(range = c(2, 6), breaks = 2:4, labels = c("Contextual", "Moderate", "Direct")) +
  labs(x = NULL, y = NULL, size = "Evidence", colour = NULL, title = "Neural context") + theme_pub +
  theme(axis.text.x = element_text(angle = 28, hjust = 1), legend.position = "bottom", plot.margin = margin(9, 6, 6, 14))
p9b <- ggplot(dis, aes(estimate, reorder(label, estimate), colour = estimate > 0)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#A8B0B5") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = .15, na.rm = TRUE) + geom_point(size = 2.2) +
  scale_colour_manual(values = c(`TRUE` = orange, `FALSE` = blue), guide = "none") +
  labs(x = "Disease-state estimate (95% CI)", y = NULL, title = "Disease-state RNA") + theme_pub +
  theme(plot.margin = margin(9, 6, 6, 14))
p9c <- ggplot(st, aes(residue, pLDDT)) +
  annotate("rect", xmin = 114, xmax = 232, ymin = -Inf, ymax = Inf, fill = "#DDF0EC", alpha = .6) +
  geom_line(linewidth = .65, colour = grey) + geom_vline(xintercept = 150.5, colour = "#7A4FA3", linewidth = .8) +
  annotate("text", x = 153, y = max(st$pLDDT, na.rm = TRUE), label = "AA150/151", hjust = 0, size = 2.7, colour = "#7A4FA3") +
  labs(x = "TSPAN14 residue", y = "AlphaFold pLDDT", title = "Canonical boundary in EC2") + theme_pub +
  theme(plot.margin = margin(9, 6, 6, 14))
save_fig(((p9a | p9b) / p9c) + plot_layout(heights = c(1.15, .9)) + plot_annotation(tag_levels = "A") & tag_theme, 9, 183, 152)

cat("Updated Supplementary Figures S4, S7 and S9 in", out, "\n")
