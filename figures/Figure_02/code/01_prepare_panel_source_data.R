#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop("Missing value for ", flag, call. = FALSE)
  args[[i + 1L]]
}

root <- normalizePath(value_after("--project-root", "."), winslash = "/", mustWork = TRUE)
legacy_dir <- file.path(root, "outputs", "main_figures_v8", "source_data")
out_dir <- file.path(root, "outputs", "main_figures_v9", "source_data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

copy_checked <- function(filename, required) {
  source <- file.path(legacy_dir, filename)
  if (!file.exists(source)) stop("Missing audited source: ", source, call. = FALSE)
  x <- fread(source, sep = "\t", na.strings = c("", "NA"))
  missing <- setdiff(required, names(x))
  if (length(missing)) stop("Missing columns in ", filename, ": ", paste(missing, collapse = ", "), call. = FALSE)
  fwrite(x, file.path(out_dir, filename), sep = "\t", na = "NA")
  x
}

tracks <- copy_checked(
  "Figure_2_locus_tracks_grch38.tsv",
  c("snp", "position_grch38", "p_value", "track", "r2_to_rs7080009", "neg_log10_p")
)
markers <- copy_checked(
  "Figure_2_annotation_markers_grch38.tsv",
  c("snp", "snp_role", "functional_prior_score", "position_grch38", "r2_to_rs7080009")
)
scatter <- copy_checked(
  "Figure_2_colocalization_scatter.tsv",
  c("snp", "trait_neg_log10_p", "sqtl_neg_log10_p", "r2_to_rs7080009", "comparison")
)
transcripts <- copy_checked(
  "Figure_2_gencode_v38_transcripts.tsv",
  c("name", "name2", "txStart", "txEnd", "exonStarts", "exonEnds", "display_name", "is_primary", "track_y")
)
exons <- copy_checked(
  "Figure_2_gencode_v38_exons.tsv",
  c("transcript", "exon_rank", "start_grch38", "end_grch38", "start_mb", "end_mb", "track_y", "is_primary")
)
ccre <- copy_checked(
  "Figure_2_regulatory_elements.tsv",
  c("chromStart_0based", "chromEnd_0based", "ccre_id", "ccre_class", "is_enhancer_like", "functional_anchor", "start_mb", "end_mb")
)

# The audited v8 scatter source predates the non-HDL panel. Reconstruct that
# panel by exact SNP matching between the retained regional tracks.
if (!"nonHDL" %chin% scatter$comparison) {
  nonhdl <- tracks[track == "nonHDL", .(
    snp, trait_neg_log10_p = neg_log10_p,
    r2_to_rs7080009
  )]
  sqtl <- tracks[track == "BA24 exact exon5-6 sQTL", .(
    snp, sqtl_neg_log10_p = neg_log10_p
  )]
  nonhdl <- merge(nonhdl, sqtl, by = "snp", all = FALSE)
  nonhdl[, `:=`(
    comparison = "nonHDL",
    ld_bin = cut(
      r2_to_rs7080009,
      breaks = c(-Inf, 0.2, 0.4, 0.6, 0.8, 0.95, Inf),
      labels = c("<0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-0.95", ">=0.95"),
      right = FALSE
    )
  )]
  scatter <- rbindlist(list(scatter, nonhdl), use.names = TRUE, fill = TRUE)
  fwrite(scatter, file.path(out_dir, "Figure_2_colocalization_scatter.tsv"), sep = "\t", na = "NA")
}

triangulation_path <- file.path(legacy_dir, "Figure_5_exact_event_triangulation.tsv")
if (!file.exists(triangulation_path)) stop("Missing exact-event triangulation source", call. = FALSE)
triangulation <- fread(triangulation_path, sep = "\t")
triangulation <- triangulation[trait %chin% c("AD", "TC", "LDL-C", "non-HDL-C")]
stopifnot(nrow(triangulation) == 4L, all(is.finite(triangulation$pph4)))
fwrite(triangulation, file.path(out_dir, "Figure_2_exact_event_coloc.tsv"), sep = "\t")

functional_path <- file.path(
  root, "outputs", "supplementary_material", "supplementary_table_sources",
  "Table_S7_LD_block_functional_annotation.tsv"
)
functional <- fread(functional_path, sep = "\t", na.strings = c("", "NA"))
functional <- functional[snp %chin% markers$snp]
setnames(functional, "pos_hg38", "position_grch38")
functional[, `:=`(
  ld_r2 = pmin(1, fifelse(is.finite(max_bridge_r2), max_bridge_r2, 0)),
  in_exact_interval = relation_to_target_splice_interval == "inside_target_splice_interval",
  near_exact_interval = relation_to_target_splice_interval %chin% c(
    "inside_target_splice_interval", "upstream_of_target_splice_interval", "downstream_of_target_splice_interval"
  ) & nearest_event_boundary_distance_bp <= 20000,
  multi_trait_cs = n_cs_memberships >= 4,
  replicated_regions = fifelse(is.finite(n_niagads_sqtl_regions), n_niagads_sqtl_regions, 0),
  crispri_anchor = grepl("Laub_core_AD_enhancer_SNP", laub_functional_anchor),
  prime_editing_anchor = snp == "rs7922621"
)]

annotation_long <- rbindlist(list(
  functional[, .(snp, metric = "EUR LD r2", value = ld_r2, metric_type = "continuous")],
  functional[, .(snp, metric = "Exact interval", value = as.numeric(in_exact_interval), metric_type = "binary")],
  functional[, .(snp, metric = "Within 20 kb", value = as.numeric(near_exact_interval), metric_type = "binary")],
  functional[, .(snp, metric = "4-trait CS", value = as.numeric(multi_trait_cs), metric_type = "binary")],
  functional[, .(snp, metric = "sQTL tissues", value = pmin(replicated_regions / 3, 1), metric_type = "continuous")],
  functional[, .(snp, metric = "CRISPRi", value = as.numeric(crispri_anchor), metric_type = "binary")],
  functional[, .(snp, metric = "Prime editing", value = as.numeric(prime_editing_anchor), metric_type = "binary")]
), use.names = TRUE)
annotation_long <- merge(
  annotation_long,
  functional[, .(snp, position_grch38, functional_prior_score, functional_priority_grade)],
  by = "snp", all.x = TRUE
)
setorder(annotation_long, position_grch38, metric)
fwrite(annotation_long, file.path(out_dir, "Figure_2_variant_annotation_matrix.tsv"), sep = "\t")

stopifnot(
  all(c("AD", "TC", "LDL", "nonHDL", "BA24 exact exon5-6 sQTL") %chin% unique(tracks$track)),
  nrow(scatter[comparison %chin% c("AD", "TC", "LDL", "nonHDL")]) > 5000,
  any(markers$snp == "rs7080009"),
  any(transcripts$display_name == "ENST00000429989"),
  nrow(ccre) == 47L,
  nrow(annotation_long) == nrow(functional) * 7L
)

manifest <- data.table(
  artifact = c(
    "Figure_2_locus_tracks_grch38.tsv", "Figure_2_colocalization_scatter.tsv",
    "Figure_2_gencode_v38_transcripts.tsv", "Figure_2_gencode_v38_exons.tsv",
    "Figure_2_regulatory_elements.tsv", "Figure_2_annotation_markers_grch38.tsv",
    "Figure_2_exact_event_coloc.tsv", "Figure_2_variant_annotation_matrix.tsv"
  ),
  role = c(
    "Aligned regional association tracks", "Locus-comparison scatter panels",
    "GENCODE v38 transcript models", "GENCODE v38 exon coordinates",
    "ENCODE candidate cis-regulatory elements", "Published and project variant anchors",
    "Exact-event colocalization posterior labels", "Compact functional annotation matrix"
  ),
  coordinate_build = c(rep("GRCh38", 6), "SNP-aligned summary", "GRCh38")
)
fwrite(manifest, file.path(out_dir, "Figure_2_source_manifest.tsv"), sep = "\t")

message("Prepared Figure 2 v9 source tables in: ", out_dir)
