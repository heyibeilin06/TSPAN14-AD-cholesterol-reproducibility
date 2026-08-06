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

v8_dir <- file.path(root, "outputs", "main_figures_v8", "source_data")
mr_dir <- file.path(root, "outputs", "mentor_revision", "complete_mr")
med_dir <- file.path(root, "outputs", "mentor_revision", "mediation_rescue")

instrument <- read_checked(
  file.path(v8_dir, "Figure_3_instrument_effect_atlas.tsv"),
  c("SNP", "BP", "layer", "outcome", "beta", "se", "p_value", "f_statistic", "z_value")
)
cis_mr <- read_checked(
  file.path(v8_dir, "Figure_3_mr_forest.tsv"),
  c("outcome", "method", "estimate", "se", "pvalue", "n_instruments", "lo", "hi", "outcome_label", "method_label")
)
diagnostics <- read_checked(
  file.path(v8_dir, "Figure_3_diagnostic_matrix.tsv"),
  c("outcome", "diagnostic", "value", "pass", "outcome_label", "diagnostic_label", "status")
)
global_mr <- read_checked(
  file.path(mr_dir, "03_genomewide_bidirectional_mr.tsv"),
  c("exposure", "outcome", "n_instruments", "method", "estimate", "se", "pvalue", "egger_intercept_pvalue")
)[outcome == "AD"]
global_apoe <- read_checked(
  file.path(mr_dir, "04_genomewide_bidirectional_apoe_sensitivity.tsv"),
  c("exposure", "outcome", "n_instruments", "method", "estimate", "se", "pvalue", "egger_intercept_pvalue")
)[outcome == "AD"]
joint_mvmr <- read_checked(
  file.path(mr_dir, "17_global_joint_lipid_mvmr_estimates.tsv"),
  c("model", "exposure", "n_jointly_clumped_instruments", "direct_estimate", "se", "pvalue")
)
joint_strength <- read_checked(
  file.path(mr_dir, "18_global_joint_lipid_mvmr_strength_sensitivity.tsv"),
  c("exposure", "assumed_pairwise_sampling_error_correlation", "conditional_F", "conditional_F_pass_10")
)
pc_sensitivity <- read_checked(
  file.path(med_dir, "14_pc_dimension_mediation_sensitivity.tsv"),
  c("lipid", "n_pcs", "cumulative_ld_variance", "splice_direct_estimate", "splice_conditional_F",
    "lipid_conditional_F", "attenuation_fraction", "all_strength_F_ge_10", "confirmatory_mediation", "status")
)

# Harmonize every SNP to a positive exact-sQTL effect before comparing layers.
orientation <- instrument[layer == "Exact exon5-6 sQTL", .(orientation = sign(beta[1])), by = SNP]
instrument <- merge(instrument, orientation, by = "SNP", all.x = TRUE)
instrument[, `:=`(
  risk_aligned_beta = beta * orientation,
  risk_aligned_z = z_value * orientation,
  display_layer = fcase(
    layer == "Exact exon5-6 sQTL", "Exact sQTL",
    outcome == "AD", "AD",
    outcome == "TC", "TC",
    outcome == "LDL", "LDL-C",
    outcome == "nonHDL", "non-HDL-C",
    default = outcome
  )
)]

global_mr[, model_label := "Genome-wide"]
global_apoe[, model_label := "Extended-APOE excluded"]
global <- rbindlist(list(global_mr, global_apoe), use.names = TRUE, fill = TRUE)
global[, `:=`(lo = estimate - 1.96 * se, hi = estimate + 1.96 * se)]

joint_mvmr[, `:=`(lo = direct_estimate - 1.96 * se, hi = direct_estimate + 1.96 * se)]
pc_sensitivity <- pc_sensitivity[
  status == "completed" & is.finite(attenuation_fraction) &
    is.finite(splice_conditional_F) & is.finite(lipid_conditional_F)
]
pc_sensitivity[, minimum_conditional_F := pmin(splice_conditional_F, lipid_conditional_F)]

stopifnot(
  uniqueN(instrument$SNP) == 5L,
  setequal(unique(instrument$display_layer), c("Exact sQTL", "AD", "TC", "LDL-C", "non-HDL-C")),
  nrow(cis_mr) == 8L,
  nrow(diagnostics) == 20L,
  sum(!diagnostics$pass) == 2L,
  nrow(global) == 10L,
  nrow(joint_mvmr) == 3L,
  nrow(joint_strength) == 12L,
  all(joint_strength$conditional_F > 10),
  uniqueN(pc_sensitivity$lipid) == 3L,
  !any(pc_sensitivity$confirmatory_mediation)
)

tables <- list(
  Figure_4_instrument_effect_atlas.tsv = instrument,
  Figure_4_ld_aware_cis_mr.tsv = cis_mr,
  Figure_4_cis_mr_diagnostics.tsv = diagnostics,
  Figure_4_genomewide_lipid_to_ad.tsv = global,
  Figure_4_global_joint_mvmr.tsv = joint_mvmr,
  Figure_4_global_joint_mvmr_strength.tsv = joint_strength,
  Figure_4_pc_gmm_dimension_sensitivity.tsv = pc_sensitivity
)
for (filename in names(tables)) fwrite(tables[[filename]], file.path(out_dir, filename), sep = "\t", na = "NA")

manifest <- data.table(
  artifact = names(tables),
  role = c(
    "Risk-aligned effects of five correlated cis-sQTL instruments", "LD-aware and lead-instrument cis-MR",
    "Colocalization, directionality, strength and pleiotropy diagnostics", "Genome-wide lipid-to-AD MR with extended-APOE sensitivity",
    "Jointly clumped genome-wide lipid MVMR direct effects", "MVMR conditional-strength sensitivity",
    "Local PC-GMM attenuation and identification across PC dimensions"
  )
)
fwrite(manifest, file.path(out_dir, "Figure_4_source_manifest.tsv"), sep = "\t")
message("Prepared Figure 4 v9 source tables in: ", out_dir)
