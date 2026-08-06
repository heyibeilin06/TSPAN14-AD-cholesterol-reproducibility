suppressPackageStartupMessages(library(MendelianRandomization))

safe_sample_size <- function(values) {
  as.integer(round(stats::median(as.numeric(values), na.rm = TRUE)))
}

validate_model_input <- function(frame, lipid) {
  required <- c("beta_splice", "se_splice", "beta_AD", "se_AD",
                paste0("beta_", lipid), paste0("se_", lipid))
  missing <- setdiff(required, names(frame))
  if (length(missing) > 0) stop("Missing columns: ", paste(missing, collapse = ", "))
  values <- frame[, required, drop = FALSE]
  if (any(!is.finite(as.matrix(values)))) stop("All model values must be finite.")
  if (any(frame$se_splice <= 0) || any(frame$se_AD <= 0) ||
      any(frame[[paste0("se_", lipid)]] <= 0)) stop("Standard errors must be positive.")
  as.integer(nrow(frame))
}

slot_or_na <- function(fit, slot_name) {
  if (is.null(fit) || !slot_name %in% methods::slotNames(fit)) return(NA)
  methods::slot(fit, slot_name)
}

scalar_or_na <- function(value) {
  if (length(value) == 0 || !is.finite(as.numeric(value[[1]]))) return(NA_real_)
  as.numeric(value[[1]])
}

delta_product <- function(a, se_a, b, se_b) {
  estimate <- as.numeric(a * b)
  se <- as.numeric(sqrt((b^2 * se_a^2) + (a^2 * se_b^2)))
  list(estimate = estimate, se = se,
       pvalue = 2 * stats::pnorm(-abs(estimate / se)))
}

fit_row <- function(fit, analysis, lipid, threshold, robust, exposure_names) {
  if (inherits(fit, "error")) {
    return(data.frame(
      analysis = analysis, lipid = lipid, pc_variance_threshold = threshold,
      robust_overdispersion = robust, exposure = exposure_names,
      estimate = NA, se = NA, ci_lower = NA, ci_upper = NA, pvalue = NA,
      instrument_F = NA, conditional_F = NA, n_pcs = NA,
      overdispersion = NA, heterogeneity = NA, status = conditionMessage(fit)
    ))
  }
  estimates <- slot_or_na(fit, "Estimate")
  standard_errors <- slot_or_na(fit, "StdError")
  lower <- slot_or_na(fit, "CILower")
  upper <- slot_or_na(fit, "CIUpper")
  pvalues <- slot_or_na(fit, "Pvalue")
  conditional_f <- slot_or_na(fit, "CondFstat")
  instrument_f <- slot_or_na(fit, "Fstat")
  data.frame(
    analysis = analysis, lipid = lipid, pc_variance_threshold = threshold,
    robust_overdispersion = robust, exposure = exposure_names,
    estimate = as.numeric(estimates), se = as.numeric(standard_errors),
    ci_lower = as.numeric(lower), ci_upper = as.numeric(upper),
    pvalue = as.numeric(pvalues),
    instrument_F = if (length(instrument_f) == 1) as.numeric(instrument_f) else NA,
    conditional_F = if (length(conditional_f) == length(exposure_names)) as.numeric(conditional_f) else NA,
    n_pcs = scalar_or_na(slot_or_na(fit, "PCs")),
    overdispersion = scalar_or_na(slot_or_na(fit, "Overdispersion")),
    heterogeneity = paste(as.numeric(slot_or_na(fit, "Heter.Stat")), collapse = ";"),
    status = "completed"
  )
}

run_models <- function(input_path, ld_path, output_dir) {
  frame <- read.delim(input_path, check.names = FALSE)
  ld <- as.matrix(read.delim(ld_path, row.names = 1, check.names = FALSE))
  if (!identical(as.character(frame$SNP), rownames(ld))) stop("LD rows do not match the harmonized SNP order.")
  if (!identical(rownames(ld), colnames(ld))) stop("LD row and column labels differ.")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  thresholds <- c(0.95, 0.99, 0.999)
  lipids <- c("TC", "LDL", "nonHDL")
  results <- list()

  for (lipid in lipids) {
    validate_model_input(frame, lipid)
    n_splice <- 147L
    n_lipid <- safe_sample_size(frame[[paste0("n_", lipid)]])
    n_ad <- safe_sample_size(frame$n_AD)
    for (threshold in thresholds) {
      for (robust in c(FALSE, TRUE)) {
        splice_input <- mr_input(
          bx = frame$beta_splice, bxse = frame$se_splice,
          by = frame$beta_AD, byse = frame$se_AD,
          correlation = ld, exposure = "Exact exon5-6 splice usage",
          outcome = "AD", snps = frame$SNP
        )
        splice_fit <- tryCatch(
          mr_pcgmm(splice_input, nx = n_splice, ny = n_ad,
                   thres = threshold, robust = robust),
          error = function(error) error
        )
        results[[length(results) + 1]] <- fit_row(
          splice_fit, "univariable PC-GMM", lipid, threshold, robust,
          "Exact exon5-6 splice usage"
        )

        lipid_input <- mr_input(
          bx = frame[[paste0("beta_", lipid)]], bxse = frame[[paste0("se_", lipid)]],
          by = frame$beta_AD, byse = frame$se_AD,
          correlation = ld, exposure = lipid, outcome = "AD", snps = frame$SNP
        )
        lipid_fit <- tryCatch(
          mr_pcgmm(lipid_input, nx = n_lipid, ny = n_ad,
                   thres = threshold, robust = robust),
          error = function(error) error
        )
        results[[length(results) + 1]] <- fit_row(
          lipid_fit, "univariable PC-GMM", lipid, threshold, robust, lipid
        )

        bx <- cbind(frame$beta_splice, frame[[paste0("beta_", lipid)]])
        bxse <- cbind(frame$se_splice, frame[[paste0("se_", lipid)]])
        colnames(bx) <- colnames(bxse) <- c("Exact exon5-6 splice usage", lipid)
        mv_input <- mr_mvinput(
          bx = bx, bxse = bxse, by = frame$beta_AD, byse = frame$se_AD,
          correlation = ld, exposure = colnames(bx), outcome = "AD", snps = frame$SNP
        )
        mv_fit <- tryCatch(
          mr_mvpcgmm(mv_input, nx = c(n_splice, n_lipid), ny = n_ad,
                     thres = threshold, robust = robust),
          error = function(error) error
        )
        results[[length(results) + 1]] <- fit_row(
          mv_fit, "multivariable PC-GMM", lipid, threshold, robust, colnames(bx)
        )
      }
    }
  }
  result <- do.call(rbind, results)
  write.table(result, file.path(output_dir, "04_pc_gmm_estimates.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  minimum_rows <- list()
  mediation_rows <- list()
  for (lipid in lipids) {
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
    lipid_to_splice_input <- mr_input(
      bx = frame[[paste0("beta_", lipid)]], bxse = frame[[paste0("se_", lipid)]],
      by = frame$beta_splice, byse = frame$se_splice, correlation = ld,
      exposure = lipid, outcome = "Exact exon5-6 splice usage", snps = frame$SNP)
    bx <- cbind(frame$beta_splice, frame[[paste0("beta_", lipid)]])
    bxse <- cbind(frame$se_splice, frame[[paste0("se_", lipid)]])
    colnames(bx) <- colnames(bxse) <- c("Exact exon5-6 splice usage", lipid)
    mv_input <- mr_mvinput(
      bx = bx, bxse = bxse, by = frame$beta_AD, byse = frame$se_AD,
      correlation = ld, exposure = colnames(bx), outcome = "AD", snps = frame$SNP)
    splice_fit <- mr_pcgmm(splice_input, nx = n_splice, ny = n_ad, r = 3, robust = TRUE)
    lipid_fit <- mr_pcgmm(lipid_input, nx = n_lipid, ny = n_ad, r = 3, robust = TRUE)
    first_step_fit <- mr_pcgmm(lipid_to_splice_input, nx = n_lipid, ny = n_splice, r = 3, robust = TRUE)
    mv_fit <- mr_mvpcgmm(mv_input, nx = c(n_splice, n_lipid), ny = n_ad, r = 3, robust = TRUE)
    minimum_rows[[length(minimum_rows) + 1]] <- fit_row(
      splice_fit, "minimum identified PC-GMM", lipid, NA, TRUE,
      "Exact exon5-6 splice usage")
    minimum_rows[[length(minimum_rows) + 1]] <- fit_row(
      lipid_fit, "minimum identified PC-GMM", lipid, NA, TRUE, lipid)
    minimum_rows[[length(minimum_rows) + 1]] <- fit_row(
      first_step_fit, "minimum identified network first step", lipid, NA, TRUE, lipid)
    minimum_rows[[length(minimum_rows) + 1]] <- fit_row(
      mv_fit, "minimum identified multivariable PC-GMM", lipid, NA, TRUE, colnames(bx))

    splice_index <- which(mv_fit@Exposure == "Exact exon5-6 splice usage")
    lipid_index <- which(mv_fit@Exposure == lipid)
    indirect <- delta_product(
      first_step_fit@Estimate, first_step_fit@StdError,
      mv_fit@Estimate[splice_index], mv_fit@StdError[splice_index])
    attenuation <- 1 - (mv_fit@Estimate[lipid_index] / lipid_fit@Estimate)
    strength_pass <- lipid_fit@Fstat >= 10 && first_step_fit@Fstat >= 10 && all(mv_fit@CondFstat >= 10)
    mediation_rows[[length(mediation_rows) + 1]] <- data.frame(
      pathway = paste(lipid, "-> exact exon5-6 splice usage -> AD"),
      n_regional_variants = nrow(frame), n_pcs = 3,
      lipid_total_estimate = lipid_fit@Estimate,
      lipid_total_se = lipid_fit@StdError,
      lipid_total_p = lipid_fit@Pvalue,
      lipid_instrument_F = lipid_fit@Fstat,
      lipid_to_splice_estimate = first_step_fit@Estimate,
      lipid_to_splice_se = first_step_fit@StdError,
      lipid_to_splice_p = first_step_fit@Pvalue,
      first_step_instrument_F = first_step_fit@Fstat,
      splice_direct_estimate = mv_fit@Estimate[splice_index],
      splice_direct_se = mv_fit@StdError[splice_index],
      splice_direct_p = mv_fit@Pvalue[splice_index],
      lipid_direct_estimate = mv_fit@Estimate[lipid_index],
      lipid_direct_se = mv_fit@StdError[lipid_index],
      lipid_direct_p = mv_fit@Pvalue[lipid_index],
      splice_conditional_F = mv_fit@CondFstat[splice_index],
      lipid_conditional_F = mv_fit@CondFstat[lipid_index],
      coefficient_attenuation_fraction = attenuation,
      indirect_product_estimate = indirect$estimate,
      indirect_product_se = indirect$se,
      indirect_product_p = indirect$pvalue,
      all_strength_criteria_pass = strength_pass,
      mediation_supported = strength_pass && first_step_fit@Pvalue < 0.05 && indirect$pvalue < 0.05,
      interpretation = if (strength_pass && first_step_fit@Pvalue < 0.05 && indirect$pvalue < 0.05)
        "identified locus-restricted mediation sensitivity" else
        "mediation not established; retain as locus-restricted decomposition only"
    )
  }
  write.table(do.call(rbind, minimum_rows), file.path(output_dir, "05_minimum_identified_pc_models.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(do.call(rbind, mediation_rows), file.path(output_dir, "06_locus_network_mediation.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  result
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  value_after <- function(flag) {
    index <- match(flag, args)
    if (is.na(index) || index == length(args)) stop("Missing ", flag)
    args[[index + 1]]
  }
  run_models(value_after("--input"), value_after("--ld"), value_after("--output-dir"))
}

if (sys.nframe() == 0) main()
