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

input <- fread(value_after("--input"))
output_dir <- value_after("--output-dir")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
traits <- c("LDL", "HDL", "TG")

formatted <- format_mvmr(
  BXGs = as.matrix(input[, paste0("beta_", traits), with = FALSE]),
  BYG = input$beta_AD,
  seBXGs = as.matrix(input[, paste0("se_", traits), with = FALSE]),
  seBYG = input$se_AD,
  RSID = input$SNP
)
estimate <- ivw_mvmr(formatted)

estimate_rows <- rbindlist(lapply(seq_along(traits), function(i) {
  data.table(
    model = "LDL + HDL + TG -> AD",
    exposure = traits[i],
    n_jointly_clumped_instruments = nrow(input),
    direct_estimate = estimate[i, "Estimate"],
    se = estimate[i, "Std. Error"],
    pvalue = estimate[i, "Pr(>|t|)"]
  )
}))
fwrite(estimate_rows, file.path(output_dir, "17_global_joint_lipid_mvmr_estimates.tsv"), sep = "\t")

sensitivity_rows <- list()
for (rho in c(0, 0.25, 0.5, 0.75)) {
  gencov <- lapply(seq_len(nrow(input)), function(i) {
    ses <- as.numeric(input[i, paste0("se_", traits), with = FALSE])
    covariance <- rho * outer(ses, ses)
    diag(covariance) <- ses^2
    covariance
  })
  strength <- strength_mvmr(formatted, gencov = gencov)
  pleiotropy <- pleiotropy_mvmr(formatted, gencov = gencov)
  for (i in seq_along(traits)) {
    sensitivity_rows[[length(sensitivity_rows) + 1]] <- data.table(
      model = "LDL + HDL + TG -> AD",
      exposure = traits[i],
      assumed_pairwise_sampling_error_correlation = rho,
      conditional_F = as.numeric(strength[1, i]),
      conditional_F_pass_10 = as.numeric(strength[1, i]) >= 10,
      q_statistic = pleiotropy$Qstat,
      q_pvalue = pleiotropy$Qpval
    )
  }
}
fwrite(rbindlist(sensitivity_rows), file.path(output_dir, "18_global_joint_lipid_mvmr_strength_sensitivity.tsv"), sep = "\t")
