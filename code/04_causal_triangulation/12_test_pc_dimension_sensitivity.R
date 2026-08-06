suppressPackageStartupMessages(library(MendelianRandomization))

source("scripts/262_run_pc_mvmr.R")

normal_positive_product_probability <- function(a, se_a, b, se_b) {
  pa <- stats::pnorm(a / se_a)
  pb <- stats::pnorm(b / se_b)
  pa * pb + (1 - pa) * (1 - pb)
}

run_dimension_sensitivity <- function(input_path, ld_path, output_dir) {
  frame <- read.delim(input_path, check.names = FALSE)
  ld <- as.matrix(read.delim(ld_path, row.names = 1, check.names = FALSE))
  if (!identical(as.character(frame$SNP), rownames(ld))) stop("LD rows do not match input SNP order.")
  eigenvalues <- pmax(eigen((ld + t(ld)) / 2, symmetric = TRUE, only.values = TRUE)$values, 0)
  cumulative_variance <- cumsum(eigenvalues) / sum(eigenvalues)
  rows <- list()

  for (lipid in c("TC", "LDL", "nonHDL")) {
    n_splice <- 147L
    n_lipid <- safe_sample_size(frame[[paste0("n_", lipid)]])
    n_ad <- safe_sample_size(frame$n_AD)
    splice_input <- mr_input(
      bx = frame$beta_splice, bxse = frame$se_splice,
      by = frame$beta_AD, byse = frame$se_AD, correlation = ld,
      exposure = "Exact exon5-6 splice usage", outcome = "AD", snps = frame$SNP)
    lipid_input <- mr_input(
      bx = frame[[paste0("beta_", lipid)]], bxse = frame[[paste0("se_", lipid)]],
      by = frame$beta_AD, byse = frame$se_AD, correlation = ld,
      exposure = lipid, outcome = "AD", snps = frame$SNP)
    first_input <- mr_input(
      bx = frame[[paste0("beta_", lipid)]], bxse = frame[[paste0("se_", lipid)]],
      by = frame$beta_splice, byse = frame$se_splice, correlation = ld,
      exposure = lipid, outcome = "Exact exon5-6 splice usage", snps = frame$SNP)
    bx <- cbind(frame$beta_splice, frame[[paste0("beta_", lipid)]])
    bxse <- cbind(frame$se_splice, frame[[paste0("se_", lipid)]])
    colnames(bx) <- colnames(bxse) <- c("Exact exon5-6 splice usage", lipid)
    mv_input <- mr_mvinput(
      bx = bx, bxse = bxse, by = frame$beta_AD, byse = frame$se_AD,
      correlation = ld, exposure = colnames(bx), outcome = "AD", snps = frame$SNP)

    for (r in c(3:15, 20, 25, 30)) {
      fits <- tryCatch({
        list(
          splice = mr_pcgmm(splice_input, nx = n_splice, ny = n_ad, r = r, robust = TRUE),
          lipid = mr_pcgmm(lipid_input, nx = n_lipid, ny = n_ad, r = r, robust = TRUE),
          first = mr_pcgmm(first_input, nx = n_lipid, ny = n_splice, r = r, robust = TRUE),
          mv = mr_mvpcgmm(mv_input, nx = c(n_splice, n_lipid), ny = n_ad, r = r, robust = TRUE)
        )
      }, error = function(error) error)
      if (inherits(fits, "error")) {
        rows[[length(rows) + 1]] <- data.frame(
          lipid = lipid, n_pcs = r, cumulative_ld_variance = cumulative_variance[r],
          lipid_total_estimate = NA_real_, lipid_total_p = NA_real_,
          lipid_instrument_F = NA_real_, lipid_to_splice_estimate = NA_real_,
          lipid_to_splice_se = NA_real_, lipid_to_splice_p = NA_real_,
          first_step_instrument_F = NA_real_, splice_direct_estimate = NA_real_,
          splice_direct_p = NA_real_, lipid_direct_estimate = NA_real_,
          lipid_direct_p = NA_real_, splice_conditional_F = NA_real_,
          lipid_conditional_F = NA_real_, attenuation_fraction = NA_real_,
          indirect_estimate = NA_real_, indirect_se = NA_real_, indirect_p = NA_real_,
          posterior_probability_indirect_positive = NA_real_,
          all_strength_F_ge_10 = FALSE, confirmatory_mediation = FALSE,
          status = conditionMessage(fits))
        next
      }
      splice_index <- which(fits$mv@Exposure == "Exact exon5-6 splice usage")
      lipid_index <- which(fits$mv@Exposure == lipid)
      indirect <- delta_product(
        fits$first@Estimate, fits$first@StdError,
        fits$mv@Estimate[splice_index], fits$mv@StdError[splice_index])
      rows[[length(rows) + 1]] <- data.frame(
        lipid = lipid, n_pcs = r, cumulative_ld_variance = cumulative_variance[r],
        lipid_total_estimate = fits$lipid@Estimate,
        lipid_total_p = fits$lipid@Pvalue,
        lipid_instrument_F = fits$lipid@Fstat,
        lipid_to_splice_estimate = fits$first@Estimate,
        lipid_to_splice_se = fits$first@StdError,
        lipid_to_splice_p = fits$first@Pvalue,
        first_step_instrument_F = fits$first@Fstat,
        splice_direct_estimate = fits$mv@Estimate[splice_index],
        splice_direct_p = fits$mv@Pvalue[splice_index],
        lipid_direct_estimate = fits$mv@Estimate[lipid_index],
        lipid_direct_p = fits$mv@Pvalue[lipid_index],
        splice_conditional_F = fits$mv@CondFstat[splice_index],
        lipid_conditional_F = fits$mv@CondFstat[lipid_index],
        attenuation_fraction = 1 - fits$mv@Estimate[lipid_index] / fits$lipid@Estimate,
        indirect_estimate = indirect$estimate,
        indirect_se = indirect$se,
        indirect_p = indirect$pvalue,
        posterior_probability_indirect_positive = normal_positive_product_probability(
          fits$first@Estimate, fits$first@StdError,
          fits$mv@Estimate[splice_index], fits$mv@StdError[splice_index]),
        all_strength_F_ge_10 = fits$lipid@Fstat >= 10 && fits$first@Fstat >= 10 &&
          all(fits$mv@CondFstat >= 10),
        confirmatory_mediation = fits$lipid@Fstat >= 10 && fits$first@Fstat >= 10 &&
          all(fits$mv@CondFstat >= 10) && fits$first@Pvalue < 0.05 && indirect$pvalue < 0.05,
        status = "completed")
    }
  }

  result <- do.call(rbind, rows)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write.table(result, file.path(output_dir, "14_pc_dimension_mediation_sensitivity.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  completed <- result[result$status == "completed", ]
  summary <- do.call(rbind, lapply(split(completed, completed$lipid), function(x) {
    strong <- x[x$all_strength_F_ge_10, ]
    data.frame(
      lipid = x$lipid[1], tested_pc_dimensions = nrow(x),
      dimensions_all_strength_F_ge_10 = nrow(strong),
      dimensions_confirmatory_mediation = sum(x$confirmatory_mediation),
      first_step_positive_dimensions = sum(x$lipid_to_splice_estimate > 0),
      first_step_p_lt_0_05_dimensions = sum(x$lipid_to_splice_p < 0.05),
      indirect_positive_dimensions = sum(x$indirect_estimate > 0),
      indirect_p_lt_0_05_dimensions = sum(x$indirect_p < 0.05),
      attenuation_min = min(x$attenuation_fraction, na.rm = TRUE),
      attenuation_max = max(x$attenuation_fraction, na.rm = TRUE),
      max_directional_probability = max(x$posterior_probability_indirect_positive, na.rm = TRUE)
    )
  }))
  write.table(summary, file.path(output_dir, "15_pc_dimension_mediation_summary.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  print(summary)
}

main_dimension <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  value_after <- function(flag) {
    index <- match(flag, args)
    if (is.na(index) || index == length(args)) stop("Missing ", flag)
    args[[index + 1]]
  }
  run_dimension_sensitivity(value_after("--input"), value_after("--ld"), value_after("--output-dir"))
}

if (sys.nframe() == 0) main_dimension()
