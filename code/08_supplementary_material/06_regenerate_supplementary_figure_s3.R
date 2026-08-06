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

blue <- "#1769AA"; teal <- "#128C7E"; orange <- "#D98200"; purple <- "#7A5195"; red <- "#C44E52"
theme_pub <- theme_classic(base_family = "Arial", base_size = 8.2) + theme(
  plot.title = element_text(size = 9.2, face = "bold", color = "#172B3A", margin = margin(b = 5, l = 14)),
  axis.title = element_text(size = 8.2), axis.text = element_text(size = 7.2, color = "#263640"),
  legend.title = element_text(size = 7.5, face = "bold"), legend.text = element_text(size = 7),
  legend.position = "bottom", plot.margin = margin(6, 8, 6, 8),
  panel.grid.major.y = element_line(color = "#E7ECEF", linewidth = .25)
)
tag_theme <- theme(plot.tag = element_text(family = "Arial", face = "bold", size = 12), plot.tag.position = c(.01, .99))

vf <- read_current("figures/Figure_01/data/Figure_1_variant_fingerprint.tsv") %>%
  filter(rank_by_product_pip <= 15)
fine_mapping <- read_current("tables/supplementary/source_data/Table_S05.tsv")
cs <- fine_mapping %>% filter(analysis_block == "credible sets")
ec <- read_current("figures/Figure_02/data/Figure_2_exact_event_coloc.tsv")

p3a <- ggplot(vf, aes(trait_label, reorder(snp, mean_pip), size = pip, color = pip)) +
  geom_point() + scale_size(range = c(1, 6)) +
  scale_color_gradient(low = "#DCE7EC", high = purple) +
  labs(x = NULL, y = NULL, size = "PIP", color = "PIP", title = "Trait-specific posterior inclusion probabilities") +
  theme_pub + theme(axis.text.x = element_text(angle = 30, hjust = 1))

p3b <- ggplot(cs, aes(n_variants, reorder(paste(trait, credible_set, sep = " | "), n_variants), color = top_pip)) +
  geom_segment(aes(x = 0, xend = n_variants, yend = reorder(paste(trait, credible_set, sep = " | "), n_variants)), color = "#CED6DB") +
  geom_point(size = 3) + scale_color_gradient(low = teal, high = red) +
  labs(x = "Variants in credible set", y = NULL, color = "Top PIP", title = "Multiple-signal credible sets") + theme_pub

p3c <- ggplot(ec, aes(pph4, reorder(trait, pph4))) +
  geom_col(width = .62, fill = blue) + geom_vline(xintercept = .8, linetype = 2, color = orange) +
  coord_cartesian(xlim = c(.8, 1)) +
  labs(x = "Default-prior exact-event PP.H4", y = NULL, title = "Exact exon5-6 colocalization") + theme_pub

figure <- ((p3a | p3b) / p3c) + plot_annotation(tag_levels = "A") & tag_theme
stem <- file.path(out, "Supplementary_Figure_S3")
ggsave(paste0(stem, ".png"), figure, width = 190, height = 170, units = "mm", dpi = 300, bg = "white")
ggsave(paste0(stem, ".pdf"), figure, width = 190, height = 170, units = "mm", device = cairo_pdf, bg = "white")
ggsave(paste0(stem, ".svg"), figure, width = 190, height = 170, units = "mm", bg = "white")
ggsave(paste0(stem, ".tiff"), figure, width = 190, height = 170, units = "mm", dpi = 600, compression = "lzw", bg = "white")
cat("Updated", stem, "\n")
