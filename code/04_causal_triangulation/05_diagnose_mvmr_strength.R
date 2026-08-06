suppressPackageStartupMessages({
  library(data.table)
  library(MVMR)
})

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) stop("Missing ", flag)
  args[[i + 1]]
}

input_path <- value_after("--local-input")
output_dir <- value_after("--output-dir")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

input <- fread(input_path)
strength_rows <- list()
estimate_rows <- list()
threshold_rows <- list()
loo_rows <- list()

for (trait in unique(input$lipid)) {
  d <- input[lipid == trait]
  formatted <- format_mvmr(
    BXGs = as.matrix(d[, .(splice_beta, lipid_beta)]),
    BYG = d$ad_beta,
    seBXGs = as.matrix(d[, .(splice_se, lipid_se)]),
    seBYG = d$ad_se,
    RSID = d$SNP
  )
  strength <- suppressWarnings(strength_mvmr(formatted, gencov = 0))
  estimate <- ivw_mvmr(formatted)
  pleiotropy <- suppressWarnings(pleiotropy_mvmr(formatted, gencov = 0))
  for (j in seq_len(2)) {
    exposure <- c("exact_exon5_6_sQTL", trait)[j]
    strength_rows[[length(strength_rows) + 1]] <- data.table(
      model = paste0("exact_exon5_6_sQTL + ", trait, " -> AD"),
      exposure = exposure,
      n_instruments = nrow(d),
      conditional_F = as.numeric(strength[1, j]),
      conditional_F_pass_10 = as.numeric(strength[1, j]) >= 10,
      genetic_covariance_assumption = 0,
      covariance_basis = "GTEx BA24 sQTL and GLGC lipid exposure samples treated as non-overlapping"
    )
    estimate_rows[[length(estimate_rows) + 1]] <- data.table(
      model = paste0("exact_exon5_6_sQTL + ", trait, " -> AD"),
      exposure = exposure,
      n_instruments = nrow(d),
      estimate = estimate[j, "Estimate"],
      se = estimate[j, "Std. Error"],
      pvalue = estimate[j, "Pr(>|t|)"],
      conditional_F = as.numeric(strength[1, j]),
      q_statistic = pleiotropy$Qstat,
      q_pvalue = pleiotropy$Qpval,
      inference_status = ifelse(
        as.numeric(strength[1, j]) >= 10,
        "conditionally_strong",
        "weak_conditional_instrument_do_not_interpret_direct_effect"
      )
    )
  }

  d[, splice_p := 2 * pnorm(-abs(splice_beta / splice_se))]
  for (threshold in c(5e-8, 1e-5, 1e-4, 1e-3)) {
    ds <- d[splice_p <= threshold]
    status <- if (nrow(ds) <= 2) "insufficient_for_overidentified_two_exposure_MVMR" else "estimable"
    threshold_rows[[length(threshold_rows) + 1]] <- data.table(
      model = paste0("exact_exon5_6_sQTL + ", trait, " -> AD"),
      sqtl_p_threshold = threshold,
      n_instruments = nrow(ds),
      status = status
    )
  }

  if (nrow(d) > 3) {
    for (omitted in d$SNP) {
      dl <- d[SNP != omitted]
      rl <- format_mvmr(
        BXGs = as.matrix(dl[, .(splice_beta, lipid_beta)]),
        BYG = dl$ad_beta,
        seBXGs = as.matrix(dl[, .(splice_se, lipid_se)]),
        seBYG = dl$ad_se,
        RSID = dl$SNP
      )
      sl <- suppressWarnings(strength_mvmr(rl, gencov = 0))
      el <- ivw_mvmr(rl)
      for (j in seq_len(2)) {
        loo_rows[[length(loo_rows) + 1]] <- data.table(
          model = paste0("exact_exon5_6_sQTL + ", trait, " -> AD"),
          omitted_snp = omitted,
          exposure = c("exact_exon5_6_sQTL", trait)[j],
          estimate = el[j, "Estimate"],
          se = el[j, "Std. Error"],
          pvalue = el[j, "Pr(>|t|)"],
          conditional_F = as.numeric(sl[1, j])
        )
      }
    }
  }
}

fwrite(rbindlist(strength_rows), file.path(output_dir, "12_local_mvmr_conditional_strength.tsv"), sep = "\t")
fwrite(rbindlist(estimate_rows), file.path(output_dir, "13_local_mvmr_official_estimates.tsv"), sep = "\t")
fwrite(rbindlist(threshold_rows), file.path(output_dir, "14_local_mvmr_threshold_sensitivity.tsv"), sep = "\t")
fwrite(rbindlist(loo_rows), file.path(output_dir, "15_local_mvmr_leave_one_out.tsv"), sep = "\t")
