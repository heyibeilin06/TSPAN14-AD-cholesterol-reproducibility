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
reg <- read_tsv(file.path(root, "figures/Figure_01/data/Figure_1_regional_screen.tsv"), show_col_types = FALSE) %>%
  mutate(focus = ifelse(grepl("TSPAN14", locus), "TSPAN14", "Other loci"))

blue <- "#1769AA"; orange <- "#D98200"; red <- "#C44E52"; grey <- "#65737D"
theme_pub <- theme_classic(base_family = "Arial", base_size = 8.2) + theme(
  plot.title = element_text(size = 9.2, face = "bold", color = "#172B3A", margin = margin(b = 5, l = 14)),
  axis.title = element_text(size = 8.2), axis.text = element_text(size = 7.2, color = "#263640"),
  legend.title = element_text(size = 7.5, face = "bold"), legend.text = element_text(size = 7),
  legend.position = "bottom", panel.grid.major.y = element_line(color = "#E7ECEF", linewidth = .25),
  plot.margin = margin(8, 10, 8, 12)
)
tag_theme <- theme(plot.tag = element_text(family = "Arial", face = "bold", size = 12), plot.tag.position = c(.01, .99))

p2a <- ggplot(reg, aes(trait_label, reorder(locus, midpoint), size = PP.H4, color = PP.H4)) +
  geom_point(alpha = .95) +
  scale_size(range = c(1.5, 7), limits = c(0, 1), breaks = c(0, .25, .5, .75, 1)) +
  scale_color_gradient(low = "#DDE6EA", high = blue, limits = c(0, 1), guide = "none") +
  labs(x = NULL, y = NULL, size = "PP.H4", title = "All screened locus-trait pairs") +
  guides(size = guide_legend(nrow = 1, title.position = "left")) +
  theme_pub +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.justification = "center", legend.margin = margin(2, 6, 2, 6))

p2b <- reg %>%
  arrange(PP.H4) %>%
  mutate(pair = paste(locus, trait_label, sep = " | ")) %>%
  ggplot(aes(PP.H4, reorder(pair, PP.H4), color = focus)) +
  geom_segment(aes(x = 0, xend = PP.H4, yend = reorder(pair, PP.H4)), color = "#C9D2D8") +
  geom_point(size = 2.5) +
  geom_vline(xintercept = .8, linetype = 2, color = orange) +
  scale_color_manual(values = c("TSPAN14" = red, "Other loci" = grey)) +
  labs(x = "Posterior probability of a shared signal", y = NULL, title = "Regional evidence ranking", color = NULL) +
  theme_pub + theme(legend.justification = "center")

figure <- (p2a | p2b) + plot_layout(widths = c(1.02, 1)) + plot_annotation(tag_levels = "A") & tag_theme
stem <- file.path(out, "Supplementary_Figure_S2")
ggsave(paste0(stem, ".png"), figure, width = 190, height = 135, units = "mm", dpi = 300, bg = "white")
ggsave(paste0(stem, ".pdf"), figure, width = 190, height = 135, units = "mm", device = cairo_pdf, bg = "white")
ggsave(paste0(stem, ".svg"), figure, width = 190, height = 135, units = "mm", bg = "white")
ggsave(paste0(stem, ".tiff"), figure, width = 190, height = 135, units = "mm", dpi = 600, compression = "lzw", bg = "white")
cat("Updated", stem, "\n")
