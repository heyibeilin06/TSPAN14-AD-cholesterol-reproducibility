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
out_dir <- file.path(root, "outputs", "main_figures_v9", "source_data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_checked <- function(path, required) {
  if (!file.exists(path)) stop("Missing audited source: ", path, call. = FALSE)
  x <- fread(path, sep = "\t", na.strings = c("", "NA"))
  missing <- setdiff(required, names(x))
  if (length(missing)) stop("Missing columns in ", basename(path), ": ", paste(missing, collapse = ", "), call. = FALSE)
  x
}

replication <- read_checked(
  file.path(root, "outputs", "main_figures_v8", "source_data", "Figure_3_exact_replication_matrix.tsv"),
  c("snp", "risk_allele", "position", "tissue", "p_value", "nes", "phenotype_id", "event_class")
)
donor <- read_checked(
  file.path(root, "outputs", "mentor_revision", "leafcutter_delta_ier", "02_ba24_sample_raw_ier_and_genotype.tsv"),
  c("donor_id", "normalized_phenotype", "genotype", "sample_id", "target_reads", "cluster_reads", "raw_ier")
)
depth <- read_checked(
  file.path(root, "outputs", "mentor_revision", "leafcutter_delta_ier", "04_depth_sensitivity_delta_ier.tsv"),
  c("minimum_cluster_reads", "n", "alt_allele", "per_alt_allele_delta_ier_percentage_points",
    "bootstrap_95ci_low", "bootstrap_95ci_high", "homozygote_contrast_percentage_points")
)
alignment <- read_checked(
  file.path(root, "outputs", "mentor_revision", "leafcutter_delta_ier", "05_alignment_and_metric_audit.tsv"),
  c("variant_id", "alternative_effect_allele", "primary_per_c_allele_delta_ier_percentage_points", "primary_hc3_p_value")
)
counts <- read_checked(
  file.path(root, "outputs", "p1_public_junction_cousage", "gtex_v8_tspan14_brain_junction_counts.tsv"),
  c("sample_id", "donor_id", "tissue", "ex5_6", "ex6_7")
)
cousage <- read_checked(
  file.path(root, "outputs", "p1_public_junction_cousage", "gtex_v8_tspan14_brain_cousage_summary.tsv"),
  c("tissue", "n_samples", "n_donors", "spearman_rho", "spearman_p", "metric_definition")
)
cluster <- read_checked(
  file.path(root, "outputs", "mentor_revision", "leafcutter_delta_ier", "01_leafcutter_target_cluster.tsv"),
  c("leafcutter_start", "leafcutter_end", "is_target_exon5_6", "event_annotation", "cluster_relation")
)
structure <- read_checked(
  file.path(root, "outputs", "targeted_splice_validation", "03_transcript_structure_audit.tsv"),
  c("transcript_id", "event", "junction_coordinates", "structural_status", "protein_transition")
)
exons <- read_checked(
  file.path(root, "outputs", "main_figures_v8", "source_data", "Figure_2_gencode_v38_exons.tsv"),
  c("transcript", "exon_rank", "start_grch38", "end_grch38", "is_primary")
)

replication[, `:=`(
  neg_log10_p = -log10(p_value),
  tissue_label = fcase(
    tissue == "Brain_Anterior_cingulate_cortex_BA24", "BA24",
    tissue == "Brain_Hippocampus", "Hippocampus",
    tissue == "Brain_Putamen_basal_ganglia", "Putamen",
    tissue == "Brain_Spinal_cord_cervical_c-1", "Spinal cord C1",
    default = tissue
  ),
  snp_label = paste0(snp, " [", risk_allele, "]")
)]
setorder(replication, position, tissue_label)

donor[, `:=`(
  raw_ier_percent = raw_ier * 100,
  genotype_label = factor(genotype, levels = 0:2, labels = c("0 copies", "1 copy", "2 copies")),
  risk_aligned_allele = "GTEx C (reverse-complement of risk-aligned G)"
)]

counts[, `:=`(
  log_ex5_6 = log1p(ex5_6),
  log_ex6_7 = log1p(ex6_7)
)]
cousage <- cousage[tissue != "All brain tissues"]
cousage[, tissue_label := sub("^Brain - ", "", tissue)]

primary_exons <- exons[transcript == "ENST00000429989" & exon_rank %in% 5:7,
                       .(exon_rank, start_grch38, end_grch38)]
splice_events <- rbindlist(list(
  data.table(
    event = "Canonical exon5-6", event_class = "exact",
    start_grch38 = structure[event == "project_exon5_exon6", as.numeric(sub(".*:(\\d+)-.*", "\\1", junction_coordinates))],
    end_grch38 = structure[event == "project_exon5_exon6", as.numeric(sub(".*-(\\d+)$", "\\1", junction_coordinates))],
    protein_transition = "AA150/151", evidence = "Exact replicated sQTL"
  ),
  data.table(
    event = "Cryptic acceptor 1", event_class = "competing",
    start_grch38 = cluster[event_annotation == "published_cryptic_exon_1_junction", leafcutter_start],
    end_grch38 = cluster[event_annotation == "published_cryptic_exon_1_junction", leafcutter_end],
    protein_transition = NA_character_, evidence = "Competing acceptor in LeafCutter cluster"
  ),
  data.table(
    event = "Adjacent exon6-7", event_class = "adjacent",
    start_grch38 = structure[event == "adjacent_exon6_exon7", as.numeric(sub(".*:(\\d+)-.*", "\\1", junction_coordinates))],
    end_grch38 = structure[event == "adjacent_exon6_exon7", as.numeric(sub(".*-(\\d+)$", "\\1", junction_coordinates))],
    protein_transition = "AA192/193", evidence = "Co-used structural context"
  )
))

stopifnot(
  nrow(replication) == 24L,
  uniqueN(replication$snp) == 6L,
  uniqueN(replication$tissue_label) == 4L,
  nrow(donor) == 147L,
  all(donor$raw_ier >= 0 & donor$raw_ier <= 1),
  all(sort(unique(donor$genotype)) == 0:2),
  nrow(depth) == 5L,
  nrow(counts) == 2642L,
  nrow(cousage) == 13L,
  nrow(primary_exons) == 3L,
  nrow(splice_events) == 3L,
  alignment$primary_per_c_allele_delta_ier_percentage_points[1] > 1
)

tables <- list(
  Figure_3_exact_replication_matrix.tsv = replication,
  Figure_3_ba24_donor_ier.tsv = donor,
  Figure_3_depth_sensitivity_delta_ier.tsv = depth,
  Figure_3_alignment_metric_audit.tsv = alignment,
  Figure_3_brain_junction_counts.tsv = counts,
  Figure_3_brain_cousage_summary.tsv = cousage,
  Figure_3_primary_exons.tsv = primary_exons,
  Figure_3_splice_events.tsv = splice_events
)
for (filename in names(tables)) fwrite(tables[[filename]], file.path(out_dir, filename), sep = "\t", na = "NA")

manifest <- data.table(
  artifact = names(tables),
  role = c(
    "Six-variant by four-tissue exact-event replication", "BA24 donor-level raw IER and genotype",
    "Read-depth sensitivity of DeltaIER", "Allele-orientation and metric audit",
    "GTEx brain exon5-6 and exon6-7 junction counts", "Brain-region junction co-usage summary",
    "GENCODE v38 exons 5-7 of ENST00000429989", "Coordinate-audited exact, competing and adjacent splice events"
  )
)
fwrite(manifest, file.path(out_dir, "Figure_3_source_manifest.tsv"), sep = "\t")
message("Prepared Figure 3 v9 source tables in: ", out_dir)
