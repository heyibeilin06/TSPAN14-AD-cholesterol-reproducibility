suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) return(default)
  if (idx == length(args)) stop("Missing value for ", flag, call. = FALSE)
  args[[idx + 1L]]
}

root <- normalizePath(value_after("--project-root", "."), winslash = "/", mustWork = TRUE)
out_dir <- file.path(root, "outputs", "main_figures_v9", "source_data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_required <- function(path, columns) {
  if (!file.exists(path)) stop("Missing required input: ", path, call. = FALSE)
  x <- fread(path, sep = "\t", na.strings = c("NA", ""))
  missing <- setdiff(columns, names(x))
  if (length(missing)) stop("Missing columns in ", basename(path), ": ", paste(missing, collapse = ", "), call. = FALSE)
  x
}

assert_unique <- function(x, keys, label) {
  if (x[, anyDuplicated(.SD), .SDcols = keys] > 0) {
    stop("Duplicate keys in ", label, ": ", paste(keys, collapse = ", "), call. = FALSE)
  }
}

write_clean <- function(x, filename) {
  x <- copy(as.data.table(x))
  char_cols <- names(x)[vapply(x, is.character, logical(1))]
  if (length(char_cols)) {
    for (column in char_cols) {
      set(x, j = column, value = gsub("[\r\n]+", " ", x[[column]], perl = TRUE))
      if (any(grepl("(^|[^A-Za-z])(C:/|C:\\\\Users|D:/|D:\\\\)", x[[column]], perl = TRUE), na.rm = TRUE)) {
        stop("Local absolute path detected in output column ", column, call. = FALSE)
      }
    }
  }
  fwrite(x, file.path(out_dir, filename), sep = "\t", quote = FALSE, na = "NA")
}

v6_source <- file.path(root, "outputs", "main_figures_v6", "source_data")
apoe_source <- file.path(root, "outputs", "mentor_revision", "apoe_conditional")

ldsc_old <- read_required(
  file.path(v6_source, "Figure_1_ldsc_context.tsv"),
  c("window", "comparison", "rg", "se", "p")
)
apoe <- read_required(
  file.path(apoe_source, "apoe_sensitivity_summary.tsv"),
  c("analysis_class", "model", "comparison", "rg", "se", "p", "scope")
)
atlas <- read_required(
  file.path(v6_source, "Figure_1_locus_evidence_atlas.tsv"),
  c("plot_order", "signal_id", "locus", "trait", "chromosome", "midpoint", "PP.H3", "PP.H4",
    "top_snp", "top_snp_pph4", "coloc_susie", "hyprcoloc", "exact_sqtl_pp_h4")
)
fingerprint <- read_required(
  file.path(v6_source, "Figure_1_variant_fingerprint.tsv"),
  c("rank_by_product_pip", "snp", "mean_pip", "n_cs_memberships", "functional_prior_score",
    "functional_anchor", "priority_class", "trait", "pip")
)
prior_table <- read_required(
  file.path(root, "tables", "supplementary", "source_data", "Table_S04.tsv"),
  c("analysis_block", "trait_pair", "pp_h4_p12_div_10")
)

# Panel A: five baseline estimates plus the current extended-APOE sensitivity analyses.
baseline <- ldsc_old[window == "baseline", .(comparison, rg, se, p)]
if (nrow(baseline) != 5L) stop("Expected five baseline LDSC estimates", call. = FALSE)

required_models <- c("baseline", "own_lead", "pair_union_leads", "w5Mb")
observed_models <- apoe[comparison == "AD-HDL", unique(model)]
if (!all(required_models %chin% observed_models)) {
  stop("APOE summary lacks one or more required models: ", paste(setdiff(required_models, observed_models), collapse = ", "), call. = FALSE)
}

sensitivity <- apoe[comparison == "AD-HDL" & model %chin% c("own_lead", "pair_union_leads", "w5Mb"),
  .(comparison, rg, se, p, analysis_class, model, scope)]
baseline[, `:=`(analysis_class = "Genome-wide baseline", model = "baseline", scope = "Unmodified genome-wide statistics")]
ldsc <- rbindlist(list(baseline, sensitivity), use.names = TRUE, fill = TRUE)
ldsc[, trait := sub("^AD-", "", comparison)]
ldsc[, trait_label := fcase(
  trait == "HDL", "HDL-C", trait == "LDL", "LDL-C", trait == "nonHDL", "non-HDL-C", default = trait
)]
ldsc[, display_label := fcase(
  model == "baseline", paste0("AD–", trait_label),
  model == "own_lead", "AD–HDL-C | LD-conditioned, trait leads",
  model == "pair_union_leads", "AD–HDL-C | LD-conditioned, union leads",
  model == "w5Mb", "AD–HDL-C | ±5-Mb exclusion"
)]
ldsc[, analysis_group := ifelse(model == "baseline", "Genome-wide LDSC", "Extended-APOE sensitivity")]
ldsc[, `:=`(lo = rg - 1.96 * se, hi = rg + 1.96 * se)]
ldsc[, display_order := match(display_label, c(
  "AD–HDL-C", "AD–LDL-C", "AD–TG", "AD–TC", "AD–non-HDL-C",
  "AD–HDL-C | LD-conditioned, trait leads",
  "AD–HDL-C | LD-conditioned, union leads",
  "AD–HDL-C | ±5-Mb exclusion"
))]
setorder(ldsc, display_order)
assert_unique(ldsc, c("display_label"), "Figure 1 LDSC source")

# Panels B and C: full non-APOE screen and method-specific evidence atlas.
atlas[, trait_label := fcase(
  trait == "HDL", "HDL-C", trait == "LDL", "LDL-C", trait == "nonHDL", "non-HDL-C", default = trait
)]
atlas[, row_label := paste(locus, trait_label, sep = " | ")]
atlas[, focus_class := fcase(
  grepl("TSPAN14", locus), "TSPAN14",
  grepl("MS4A", locus), "MS4A comparator",
  default = "Other screened locus"
)]
atlas[, posterior_class := cut(
  PP.H4, breaks = c(-Inf, 0.50, 0.80, 0.95, Inf),
  labels = c("<0.50", "0.50–0.79", "0.80–0.94", "≥0.95"), right = FALSE
)]
atlas[, conditional_top_snp_share := fifelse(PP.H4 >= 0.80, top_snp_pph4, NA_real_)]
atlas[, tested_layers := 1L + as.integer(!is.na(coloc_susie)) + as.integer(!is.na(hyprcoloc)) + as.integer(!is.na(exact_sqtl_pp_h4))]
atlas[, supported_layers := as.integer(PP.H4 >= 0.80) +
        as.integer(!is.na(coloc_susie) & coloc_susie >= 0.80) +
        as.integer(!is.na(hyprcoloc) & hyprcoloc >= 0.80) +
        as.integer(!is.na(exact_sqtl_pp_h4) & exact_sqtl_pp_h4 >= 0.80)]
atlas[, supported_fraction := supported_layers / tested_layers]
atlas[, supported_tested_label := paste0(supported_layers, "/", tested_layers)]
assert_unique(atlas, c("signal_id"), "Figure 1 regional screen")

regional <- atlas[, .(
  plot_order, signal_id, locus, trait, trait_label, chromosome, midpoint,
  PP.H3, PP.H4, top_snp, top_snp_pph4, posterior_class, focus_class
)]
evidence <- atlas[, .(
  plot_order, signal_id, locus, trait_label, row_label, focus_class,
  PP.H4, PP.H3, conditional_top_snp_share, coloc_susie, hyprcoloc,
  exact_sqtl_pp_h4, tested_layers, supported_layers, supported_fraction,
  supported_tested_label
)]

regional_prior <- unique(prior_table[
  analysis_block == "coloc prior sensitivity" & trait_pair %chin% c("AD-TC", "AD-LDL-C", "AD-non-HDL-C"),
  .(
    signal_id = fcase(
      trait_pair == "AD-TC", "chr10 TSPAN14 | TC",
      trait_pair == "AD-LDL-C", "chr10 TSPAN14 | LDL",
      trait_pair == "AD-non-HDL-C", "chr10 TSPAN14 | nonHDL"
    ),
    regional_pp_h4_conservative = pp_h4_p12_div_10,
    interpretation = fifelse(
      pp_h4_p12_div_10 >= 0.80,
      "Most prior-robust regional convergence",
      "High default-prior support; prior-sensitive"
    )
  )
])
assert_unique(regional_prior, "signal_id", "Figure 1 regional prior sensitivity")

# Panel D: top 15 variants with full trait-specific PIP records.
variant_rank <- unique(fingerprint[, .(snp, rank_by_product_pip)])
setorder(variant_rank, rank_by_product_pip)
keep_snps <- variant_rank[rank_by_product_pip <= 15, snp]
anchors <- c("rs1870137", "rs1870138", "rs7080009")
if (!all(anchors %chin% keep_snps)) stop("A required functional anchor falls outside the top 15 variants", call. = FALSE)
fingerprint <- fingerprint[snp %chin% keep_snps]
fingerprint[, trait_label := fcase(
  trait %chin% c("LDL", "LDL-C"), "LDL-C",
  trait %chin% c("nonHDL", "Non-HDL-C", "non-HDL-C"), "non-HDL-C",
  default = trait
)]
assert_unique(fingerprint, c("snp", "trait_label"), "Figure 1 variant fingerprint")

write_clean(ldsc, "Figure_1_ldsc_apoe_conditioning.tsv")
write_clean(regional, "Figure_1_regional_screen.tsv")
write_clean(evidence, "Figure_1_evidence_matrix.tsv")
write_clean(regional_prior, "Figure_1_regional_prior_sensitivity.tsv")
write_clean(fingerprint, "Figure_1_variant_fingerprint.tsv")

cat("Figure 1 v9 source preparation complete\n")
cat("LDSC rows:", nrow(ldsc), "\n")
cat("Regional rows:", nrow(regional), "\n")
cat("Evidence rows:", nrow(evidence), "\n")
cat("Fine-mapping variants:", uniqueN(fingerprint$snp), "\n")
