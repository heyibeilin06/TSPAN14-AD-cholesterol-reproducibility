script_arg <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)])
script_dir <- dirname(normalizePath(script_arg, winslash = "/", mustWork = TRUE))
source(file.path(script_dir, "00_figure_style_and_io.R"))

paths <- parse_v9_args()

ldsc <- read_v9_source(
  paths$source_dir, "Figure_1_ldsc_apoe_conditioning.tsv",
  c("display_order", "display_label", "analysis_group", "model", "rg", "lo", "hi", "p")
)
regional <- read_v9_source(
  paths$source_dir, "Figure_1_regional_screen.tsv",
  c("signal_id", "locus", "trait_label", "chromosome", "midpoint", "PP.H3", "PP.H4", "top_snp_pph4", "focus_class")
)
evidence <- read_v9_source(
  paths$source_dir, "Figure_1_evidence_matrix.tsv",
  c("signal_id", "row_label", "locus", "focus_class", "PP.H4", "PP.H3", "conditional_top_snp_share",
    "coloc_susie", "hyprcoloc", "exact_sqtl_pp_h4", "supported_fraction")
)
fingerprint <- read_v9_source(
  paths$source_dir, "Figure_1_variant_fingerprint.tsv",
  c("rank_by_product_pip", "snp", "mean_pip", "n_cs_memberships", "functional_prior_score",
    "functional_anchor", "trait_label", "pip")
)

assert_unique_v9(ldsc, "display_label", "LDSC panel")
assert_unique_v9(regional, "signal_id", "regional skyline")
assert_unique_v9(evidence, "signal_id", "evidence matrix")
assert_unique_v9(fingerprint, c("snp", "trait_label"), "variant fingerprint")

# -----------------------------------------------------------------------------
# A. Genetic-correlation forest with current extended-APOE sensitivity analyses.
# -----------------------------------------------------------------------------
ldsc[, row_factor := factor(display_label, levels = rev(display_label[order(display_order)]))]
ldsc[, display_class := fcase(
  model == "baseline" & display_label == "AD–HDL-C", "Baseline AD–HDL-C",
  model == "baseline", "Other baseline traits",
  model %chin% c("own_lead", "pair_union_leads"), "LD-conditioned",
  default = "Physical-window sensitivity"
)]
ldsc[, p_label := paste0("P = ", format_p_v9(p))]

class_colours <- c(
  "Baseline AD–HDL-C" = v9_palette[["amber"]],
  "Other baseline traits" = v9_palette[["grey"]],
  "LD-conditioned" = v9_palette[["blue"]],
  "Physical-window sensitivity" = v9_palette[["teal"]]
)

p_a <- ggplot(ldsc, aes(rg, row_factor)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.34, colour = v9_palette[["graphite"]]) +
  geom_hline(yintercept = 3.5, linewidth = 0.34, colour = v9_palette[["grid"]]) +
  geom_errorbar(
    aes(xmin = lo, xmax = hi, colour = display_class),
    orientation = "y", width = 0.12, linewidth = 0.68
  ) +
  geom_point(aes(fill = display_class), shape = 21, size = 3.3, colour = "white", stroke = 0.45) +
  geom_text(aes(x = 0.292, label = p_label), hjust = 1, size = 1.75, colour = v9_palette[["graphite"]]) +
  annotate("text", x = -0.167, y = 8.42, label = "Genome-wide LDSC", hjust = 0,
           size = 1.82, fontface = "bold", colour = v9_palette[["graphite"]]) +
  annotate("text", x = -0.167, y = 3.35, label = "Extended-APOE sensitivity", hjust = 0,
           size = 1.82, fontface = "bold", colour = v9_palette[["graphite"]]) +
  scale_colour_manual(values = class_colours, guide = "none") +
  scale_fill_manual(values = class_colours, guide = "none") +
  scale_x_continuous(limits = c(-0.17, 0.30), breaks = c(-0.1, 0, 0.1, 0.2), expand = c(0, 0)) +
  labs(tag = "A", x = "Genetic correlation with AD", y = NULL) +
  theme_v9(6.25) +
  theme(axis.text.y = element_text(size = 5.65, lineheight = 0.9), plot.margin = margin(3, 3, 2, 2))

# -----------------------------------------------------------------------------
# B. Mirrored regional posterior skyline: shared PP.H4 versus distinct PP.H3.
# -----------------------------------------------------------------------------
setorder(regional, chromosome, midpoint, trait_label)
regional[, x_index := .I]
regional[, signal_label := paste0(sub("^chr[0-9]+ ", "", locus), "–", trait_label)]
regional[, label_y := fifelse(PP.H4 >= 0.80, PP.H4, -PP.H3)]

chromosome_bands <- regional[, .(
  xmin = min(x_index) - 0.45,
  xmax = max(x_index) + 0.45,
  xmid = mean(range(x_index))
), by = chromosome]
chromosome_bands[, band := rep(c("white", "shade"), length.out = .N)]

label_data <- regional[focus_class %chin% c("TSPAN14", "MS4A comparator")]
label_data[, `:=`(label_x = as.numeric(x_index), label_y_plot = PP.H4, label_hjust = 0.5)]
label_data[focus_class == "TSPAN14" & trait_label == "LDL-C", `:=`(label_x = x_index - 0.15, label_y_plot = 0.84, label_hjust = 0.5)]
label_data[focus_class == "TSPAN14" & trait_label == "TC", `:=`(label_x = x_index - 0.05, label_y_plot = 1.055, label_hjust = 1.0)]
label_data[focus_class == "TSPAN14" & trait_label == "non-HDL-C", `:=`(label_x = x_index + 0.05, label_y_plot = 1.00, label_hjust = 0.0)]
label_data[focus_class == "MS4A comparator", `:=`(label_x = x_index - 0.05, label_y_plot = 0.72, label_hjust = 0.0)]

p_b <- ggplot() +
  geom_rect(
    data = chromosome_bands,
    aes(xmin = xmin, xmax = xmax, ymin = -1.08, ymax = 1.08, fill = band),
    inherit.aes = FALSE, colour = NA
  ) +
  geom_hline(yintercept = 0, linewidth = 0.38, colour = v9_palette[["ink"]]) +
  geom_hline(yintercept = c(-0.80, 0.80), linewidth = 0.30, linetype = "dashed", colour = v9_palette[["grey"]]) +
  geom_segment(
    data = regional,
    aes(x = x_index, xend = x_index, y = 0, yend = PP.H4, colour = trait_label),
    linewidth = 0.75, lineend = "round"
  ) +
  geom_point(
    data = regional,
    aes(x_index, PP.H4, size = top_snp_pph4, fill = trait_label),
    shape = 21, colour = "white", stroke = 0.42
  ) +
  geom_segment(
    data = regional,
    aes(x = x_index, xend = x_index, y = 0, yend = -PP.H3),
    linewidth = 0.55, colour = v9_palette[["red"]], alpha = 0.72
  ) +
  geom_point(
    data = regional,
    aes(x_index, -PP.H3),
    shape = 21, size = 1.8, fill = v9_palette[["red_soft"]], colour = "white", stroke = 0.32
  ) +
  geom_text(
    data = label_data,
    aes(label_x, label_y_plot, label = signal_label, colour = trait_label, hjust = label_hjust),
    size = 1.72, fontface = "bold"
  ) +
  annotate("text", x = 0.65, y = 1.04, label = "Default-prior PP.H4", hjust = 0,
           size = 1.8, fontface = "bold", colour = v9_palette[["graphite"]]) +
  annotate("text", x = 0.65, y = -1.04, label = "Distinct configurations (PP.H3)", hjust = 0,
           size = 1.8, fontface = "bold", colour = v9_palette[["graphite"]]) +
  scale_fill_manual(values = c(white = "white", shade = v9_palette[["pale"]], v9_trait_colours), guide = "none") +
  scale_colour_manual(values = v9_trait_colours, name = NULL) +
  scale_size_continuous(range = c(1.6, 5.0), limits = c(0, 1), guide = "none") +
  scale_x_continuous(
    breaks = chromosome_bands$xmid, labels = paste0("chr", chromosome_bands$chromosome),
    expand = expansion(mult = c(0.025, 0.035))
  ) +
  scale_y_continuous(
    limits = c(-1.10, 1.10), breaks = c(-1, -0.5, 0, 0.5, 1),
    labels = c("1.0", "0.5", "0", "0.5", "1.0")
  ) +
  labs(tag = "B", x = NULL, y = "Posterior probability") +
  theme_v9(6.2) +
  theme(
    axis.line.x = element_blank(), axis.ticks.x = element_blank(),
    legend.position = "bottom", legend.box = "horizontal",
    legend.margin = margin(0, 0, 0, 0), legend.spacing.x = grid::unit(2, "mm"),
    legend.key.width = grid::unit(2.8, "mm"), legend.text = element_text(size = 5.3),
    plot.margin = margin(3, 3, 2, 2)
  ) +
  guides(
    colour = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(linewidth = 1.2, size = 2.6))
  )

# -----------------------------------------------------------------------------
# C. Locus-by-method evidence matrix.
# -----------------------------------------------------------------------------
method_spec <- data.table(
  source_column = c("PP.H4", "PP.H3", "conditional_top_snp_share", "coloc_susie",
                    "hyprcoloc", "exact_sqtl_pp_h4", "supported_fraction"),
  method = c("Regional\nH4", "Regional\nH3", "Top-SNP\nshare",
             "SuSiE\nH4", "HyPrColoc", "Exact exon5–6\nH4", "Evidence\ncoverage"),
  method_group = c(rep("Regional screen", 3), rep("TSPAN14 resolution", 3), "Summary"),
  family = c("Shared signal", "Distinct signal", "Variant concentration",
             "Fine-mapped signal", "Multi-trait signal", "Exact splicing", "Evidence coverage")
)

for (column in method_spec$source_column) set(evidence, j = column, value = as.numeric(evidence[[column]]))

evidence <- evidence[order(-PP.H4, row_label)]
row_levels <- evidence$row_label
evidence[, row_factor := factor(row_label, levels = rev(row_levels))]

matrix_long <- melt(
  evidence,
  id.vars = c("signal_id", "row_factor", "locus", "focus_class"),
  measure.vars = method_spec$source_column,
  variable.name = "source_column", value.name = "value"
)
matrix_long <- merge(matrix_long, method_spec, by = "source_column", all.x = TRUE, sort = FALSE)
matrix_long[, method := factor(method, levels = method_spec$method)]
matrix_long[, method_group := factor(method_group, levels = c("Regional screen", "TSPAN14 resolution", "Summary"))]
matrix_long[, available := !is.na(value)]

matrix_grid <- unique(matrix_long[, .(method, method_group, row_factor, focus_class)])
matrix_grid[, background := ifelse(focus_class == "TSPAN14", "focus", "base")]

family_colours <- c(
  "Shared signal" = v9_palette[["teal"]],
  "Distinct signal" = v9_palette[["red"]],
  "Variant concentration" = v9_palette[["graphite"]],
  "Fine-mapped signal" = v9_palette[["blue"]],
  "Multi-trait signal" = v9_palette[["teal_dark"]],
  "Exact splicing" = v9_palette[["purple"]],
  "Evidence coverage" = v9_palette[["amber"]]
)

p_c <- ggplot() +
  geom_tile(
    data = matrix_grid,
    aes(method, row_factor, fill = background),
    width = 0.94, height = 0.86, colour = "white", linewidth = 0.38
  ) +
  geom_point(
    data = matrix_long[available == TRUE],
    aes(method, row_factor, size = value, fill = family),
    shape = 21, colour = "white", stroke = 0.34
  ) +
  geom_text(
    data = matrix_long[available == FALSE],
    aes(method, row_factor, label = "×"),
    size = 2.25, colour = "#B5BEC3"
  ) +
  scale_fill_manual(
    values = c(base = v9_palette[["pale"]], focus = "#DFEDF5", family_colours),
    guide = "none"
  ) +
  scale_size_continuous(range = c(0.9, 5.1), limits = c(0, 1), breaks = c(0.25, 0.50, 0.75, 1.00),
                        name = "Posterior / relative support") +
  facet_grid(cols = vars(method_group), scales = "free_x", space = "free_x") +
  scale_x_discrete(position = "top") +
  labs(tag = "C", x = NULL, y = NULL) +
  theme_v9(5.75) +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_text(size = 5.25, face = "bold", lineheight = 0.88, margin = margin(b = 1.5)),
    axis.text.y = element_text(size = 5.25),
    strip.text = element_text(size = 5.3, face = "bold", colour = v9_palette[["graphite"]], margin = margin(b = 1.2)),
    legend.position = "bottom", legend.box = "vertical", legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(2, 2, 2, 2)
  ) +
  guides(
    size = guide_legend(nrow = 1, byrow = TRUE)
  )

# -----------------------------------------------------------------------------
# D. Fine-mapping and functional fingerprint for the top 15 variants.
# -----------------------------------------------------------------------------
variant_meta <- unique(fingerprint[, .(
  rank_by_product_pip, snp, mean_pip, n_cs_memberships,
  functional_prior_score, functional_anchor, priority_class
)])
setorder(variant_meta, rank_by_product_pip)
variant_meta[, snp_factor := factor(snp, levels = rev(snp))]

pip_long <- fingerprint[, .(snp, trait_label, pip)]
pip_long <- merge(pip_long, variant_meta[, .(snp, snp_factor)], by = "snp", all.x = TRUE)
mean_rows <- variant_meta[, .(snp, trait_label = "Mean\nPIP", pip = mean_pip, snp_factor)]
pip_long <- rbindlist(list(pip_long, mean_rows), use.names = TRUE)

column_levels <- c("Mean\nPIP", "AD", "TC", "LDL-C", "non-HDL-C", "CS\ncount", "Functional\nscore", "CRISPRi")
pip_long[, column := factor(trait_label, levels = column_levels)]
variant_meta[, `:=`(
  cs_column = factor("CS\ncount", levels = column_levels),
  score_column = factor("Functional\nscore", levels = column_levels),
  anchor_column = factor("CRISPRi", levels = column_levels)
)]

fingerprint_grid <- CJ(
  snp_factor = factor(levels(variant_meta$snp_factor), levels = levels(variant_meta$snp_factor)),
  column = factor(column_levels, levels = column_levels), unique = TRUE
)

p_d <- ggplot() +
  geom_tile(
    data = fingerprint_grid,
    aes(column, snp_factor), width = 0.94, height = 0.84,
    fill = v9_palette[["pale"]], colour = "white", linewidth = 0.36
  ) +
  geom_point(
    data = pip_long,
    aes(column, snp_factor, size = pip, fill = pip),
    shape = 21, colour = "white", stroke = 0.32
  ) +
  geom_point(
    data = variant_meta[n_cs_memberships > 0],
    aes(cs_column, snp_factor, size = n_cs_memberships / 4),
    shape = 21, fill = v9_palette[["blue_soft"]], colour = "white", stroke = 0.32
  ) +
  geom_point(
    data = variant_meta[!is.na(functional_prior_score)],
    aes(score_column, snp_factor, size = pmin(functional_prior_score / 18, 1)),
    shape = 21, fill = v9_palette[["amber_soft"]], colour = "white", stroke = 0.32
  ) +
  geom_point(
    data = variant_meta[functional_anchor == TRUE],
    aes(anchor_column, snp_factor), shape = 23, size = 3.5,
    fill = v9_palette[["amber"]], colour = "white", stroke = 0.48
  ) +
  scale_fill_gradient(low = v9_palette[["teal_soft"]], high = v9_palette[["teal_dark"]],
                      limits = c(0, 0.40), oob = squish, name = "PIP") +
  scale_size_continuous(range = c(0.8, 4.3), limits = c(0, 1), guide = "none") +
  scale_x_discrete(
    position = "top", drop = FALSE,
    labels = c("Mean\nPIP", "AD", "TC", "LDL", "non-HDL", "CS", "Score", "Anchor")
  ) +
  labs(tag = "D", x = NULL, y = NULL) +
  theme_v9(5.75) +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_text(size = 4.90, face = "bold", lineheight = 0.88),
    axis.text.y = element_text(size = 5.20),
    legend.position = "bottom", legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(2, 2, 2, 2)
  ) +
  guides(fill = guide_colorbar(
    title.position = "left", barwidth = grid::unit(14, "mm"),
    barheight = grid::unit(2.2, "mm"), display = "rectangles"
  ))

top_row <- p_a + p_b + plot_layout(widths = c(0.80, 1.38))
bottom_row <- p_c + p_d + plot_layout(widths = c(1.30, 0.95))
figure_1 <- top_row / bottom_row + plot_layout(heights = c(0.86, 1.25)) & panel_tag_theme_v9

stem <- "Figure_1_APOE_aware_TSPAN14_prioritization_v9"
save_pub_v9(figure_1, paths$output_dir, stem, width_mm = 183, height_mm = 168)
assert_clean_svg_v9(file.path(paths$output_dir, paste0(stem, ".svg")))

citation_ledger <- data.table(
  manuscript_result = c(
    "Genome-wide AD-lipid correlation and extended-APOE sensitivity",
    "Non-APOE regional screen prioritizing TSPAN14",
    "Multi-method evidence convergence across screened loci",
    "TSPAN14 credible-set and functional-anchor interpretation"
  ),
  figure_citation = c("Fig. 1A", "Fig. 1B", "Fig. 1C", "Fig. 1D")
)
fwrite(citation_ledger, file.path(paths$output_dir, "Figure_1_manuscript_citation_ledger.tsv"), sep = "\t")

cat("Figure 1 v9 exported to", paths$output_dir, "\n")
