#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(ifelse(length(args), args[1], "."), winslash = "/")
out <- file.path(root, "outputs/supplement_v19/figures")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
read_current <- function(path) read_tsv(file.path(root, path), show_col_types = FALSE, progress = FALSE)

blue <- "#1769AA"; teal <- "#128C7E"; purple <- "#7A5195"; red <- "#C44E52"; grey <- "#56636D"
theme_pub <- theme_classic(base_family = "Arial", base_size = 8.2) + theme(
  plot.title = element_text(size = 9.2, face = "bold", color = "#172B3A", margin = margin(b = 5, l = 14)),
  axis.title = element_text(size = 8.2), axis.text = element_text(size = 7.2, color = "#263640"),
  legend.title = element_text(size = 7.5, face = "bold"), legend.text = element_text(size = 7),
  legend.position = "bottom", panel.grid.major.y = element_line(color = "#E7ECEF", linewidth = .25)
)
tag_theme <- theme(plot.tag = element_text(family = "Arial", face = "bold", size = 12), plot.tag.position = c(.01, .99))

atlas <- read_current("figures/Figure_05/data/Figure_5_cell_context_atlas.tsv") %>%
  mutate(score = as.numeric(evidence_strength))
p9a <- ggplot(atlas, aes(evidence_layer, reorder(context, score), size = score, color = evidence_class)) +
  geom_point(alpha = .9) +
  scale_size(range = c(2, 6), breaks = 2:4, labels = c("Contextual", "Moderate", "Direct")) +
  labs(x = NULL, y = NULL, size = "Evidence", color = NULL, title = "Neural cell-context evidence") +
  guides(size = guide_legend(order = 1, nrow = 1), color = guide_legend(order = 2, nrow = 2, byrow = TRUE)) +
  theme_pub +
  theme(
    axis.text.x = element_text(angle = 28, hjust = 1), legend.justification = "center",
    legend.box = "vertical", legend.margin = margin(2, 8, 2, 8),
    legend.box.margin = margin(0, 6, 0, 6), plot.margin = margin(8, 10, 8, 16)
  )

dis <- read_current("figures/Figure_05/data/Figure_5_disease_state_rna.tsv") %>%
  mutate(label = paste(source_label, cell_label, sep = " | "))
p9b <- ggplot(dis, aes(estimate, reorder(label, estimate), color = estimate > 0)) +
  geom_vline(xintercept = 0, linetype = 2, color = "#A8B0B5") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = .15, na.rm = TRUE) +
  geom_point(size = 2.2) +
  scale_color_manual(values = c(`TRUE` = red, `FALSE` = blue), guide = "none") +
  labs(x = "Disease-state estimate (95% CI)", y = NULL, title = "Disease-state RNA sensitivity") +
  theme_pub + theme(plot.margin = margin(8, 12, 8, 10))

st <- read_current("figures/Figure_05/data/Figure_5_ec2_structure.tsv")
p9c <- ggplot(st, aes(residue, pLDDT)) +
  annotate("rect", xmin = 114, xmax = 232, ymin = -Inf, ymax = Inf, fill = "#DDF0EC", alpha = .6) +
  geom_line(linewidth = .65, color = grey) +
  geom_vline(xintercept = 150.5, color = purple, linewidth = .8) +
  annotate("text", x = 153, y = max(st$pLDDT, na.rm = TRUE), label = "AA150/151", hjust = 0, size = 2.7, color = purple) +
  labs(x = "TSPAN14 residue", y = "AlphaFold pLDDT", title = "Exact splice boundary within the EC2 region") +
  theme_pub + theme(plot.margin = margin(8, 10, 8, 14))

figure <- ((p9a | p9b) / p9c) +
  plot_layout(widths = c(1.08, 1), heights = c(1.12, 1)) +
  plot_annotation(tag_levels = "A") & tag_theme
stem <- file.path(out, "Supplementary_Figure_S9")
ggsave(paste0(stem, ".png"), figure, width = 195, height = 175, units = "mm", dpi = 300, bg = "white")
ggsave(paste0(stem, ".pdf"), figure, width = 195, height = 175, units = "mm", device = cairo_pdf, bg = "white")
ggsave(paste0(stem, ".svg"), figure, width = 195, height = 175, units = "mm", bg = "white")
ggsave(paste0(stem, ".tiff"), figure, width = 195, height = 175, units = "mm", dpi = 600, compression = "lzw", bg = "white")
cat("Updated", stem, "\n")
