suppressPackageStartupMessages({
  library(data.table)
  library(MRPRESSO)
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
set.seed(20260720)

presso_rows <- list()
loo_rows <- list()
for (key in unique(input[, paste(exposure, outcome, sep = "__")])) {
  parts <- strsplit(key, "__", fixed = TRUE)[[1]]
  exposure_name <- parts[1]
  outcome_name <- parts[2]
  d_full <- input[exposure == exposure_name & outcome == outcome_name][order(p_x)]
  d <- d_full[seq_len(min(.N, 100))]
  if (nrow(d) < 4) next
  fit <- tryCatch(
    mr_presso(
      BetaOutcome = "beta_y",
      BetaExposure = "beta_x",
      SdOutcome = "se_y",
      SdExposure = "se_x",
      OUTLIERtest = TRUE,
      DISTORTIONtest = TRUE,
      data = as.data.frame(d),
      NbDistribution = 1000,
      SignifThreshold = 0.05
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    presso_rows[[length(presso_rows) + 1]] <- data.table(
      exposure = exposure_name, outcome = outcome_name, n_instruments = nrow(d),
      full_primary_instrument_count = nrow(d_full),
      sensitivity_scope = "100 strongest LD-clumped exposure instruments",
      status = "failed", error = conditionMessage(fit)
    )
  } else {
    main <- as.data.table(fit$`Main MR results`, keep.rownames = "estimate_type")
    global_p <- fit$`MR-PRESSO results`$`Global Test`$Pvalue
    distortion <- fit$`MR-PRESSO results`$`Distortion Test`
    outliers <- fit$`MR-PRESSO results`$`Outlier Test`
    for (i in seq_len(nrow(main))) {
      outlier_resolution_ok <- nrow(d) / 1000 <= 0.05
      presso_rows[[length(presso_rows) + 1]] <- data.table(
        exposure = exposure_name,
        outcome = outcome_name,
        n_instruments = nrow(d),
        full_primary_instrument_count = nrow(d_full),
        sensitivity_scope = "100 strongest LD-clumped exposure instruments",
        estimate_type = ifelse(i == 1, "raw", "outlier_corrected"),
        estimate = main$`Causal Estimate`[i],
        se = main$Sd[i],
        pvalue = main$`P-value`[i],
        global_test_pvalue = global_p,
        n_outliers = if (!outlier_resolution_ok || is.null(outliers)) NA_integer_ else sum(as.numeric(outliers$Pvalue) <= 0.05, na.rm = TRUE),
        outlier_test_status = ifelse(outlier_resolution_ok, "resolved", "insufficient_null_resolution_at_1000_simulations"),
        distortion_pvalue = if (is.null(distortion)) NA_real_ else distortion$Pvalue,
        status = "completed",
        error = ""
      )
    }
  }

  for (omitted in d$SNP) {
    dl <- d[SNP != omitted]
    w <- 1 / dl$se_y^2
    estimate <- sum(w * dl$beta_x * dl$beta_y) / sum(w * dl$beta_x^2)
    loo_rows[[length(loo_rows) + 1]] <- data.table(
      exposure = exposure_name,
      outcome = outcome_name,
      omitted_snp = omitted,
      estimate = estimate
    )
  }
}
fwrite(rbindlist(presso_rows, fill = TRUE), file.path(output_dir, "19_bidirectional_mr_presso.tsv"), sep = "\t")
fwrite(rbindlist(loo_rows), file.path(output_dir, "20_bidirectional_leave_one_out.tsv"), sep = "\t")
