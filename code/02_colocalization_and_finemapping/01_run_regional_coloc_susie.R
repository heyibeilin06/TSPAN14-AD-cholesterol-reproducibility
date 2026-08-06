suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})

args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args0, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else file.path("scripts", "18_run_iteration7_ldaware_coloc_susie.R")
ROOT <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
if (!dir.exists(file.path(ROOT, "results"))) ROOT <- normalizePath(getwd(), mustWork = TRUE)

TARGETS <- data.table(
  axis = c("TSPAN14", "TSPAN14", "TSPAN14", "MS4A_boundary"),
  region_id = c("TC_nonAPOE_004", "LDL_nonAPOE_004", "nonHDL_nonAPOE_004", "TG_nonAPOE_004"),
  trait = c("TC", "LDL", "nonHDL", "TG")
)

LOCUS_DIR <- file.path(ROOT, "results", "tables", "coloc_loci_lead250kb")
LD_DIR <- file.path(ROOT, "results", "ld", "iteration7")
OUT_RESULTS <- file.path(ROOT, "results", "tables", "technical_iteration7_ldaware_coloc_susie_results.tsv")
OUT_SIGNALS <- file.path(ROOT, "results", "tables", "technical_iteration7_ldaware_coloc_susie_signal_summary.tsv")
OUT_VARIANTS <- file.path(ROOT, "results", "tables", "technical_iteration7_ldaware_coloc_susie_variant_inputs.tsv")
OUT_REPORT <- file.path(ROOT, "results", "reports", "preliminary_result_28_iteration7_ldaware_coloc_susie.md")

max_variants <- as.integer(Sys.getenv("ITER7_SUSIE_MAX_VARIANTS", "1200"))
min_variants <- as.integer(Sys.getenv("ITER7_SUSIE_MIN_VARIANTS", "100"))

read_ld <- function(region_id) {
  path <- file.path(LD_DIR, paste0(region_id, ".ld.tsv"))
  if (!file.exists(path)) stop("Missing LD matrix: ", path)
  m <- as.matrix(fread(path), rownames = 1)
  storage.mode(m) <- "double"
  m[!is.finite(m)] <- 0
  diag(m) <- 1
  m <- (m + t(m)) / 2
  m
}

make_dataset <- function(dt, beta_col, varbeta_col, maf_col, n_col, type, s_col = NULL, R) {
  d <- list(
    beta = dt[[beta_col]],
    varbeta = dt[[varbeta_col]],
    snp = dt$snp,
    MAF = dt[[maf_col]],
    N = as.integer(round(stats::median(dt[[n_col]], na.rm = TRUE))),
    type = type,
    LD = R
  )
  if (!is.null(s_col)) d$s <- stats::median(dt[[s_col]], na.rm = TRUE)
  d
}

prep_locus <- function(region_id) {
  locus_path <- file.path(LOCUS_DIR, paste0(region_id, ".tsv"))
  if (!file.exists(locus_path)) stop("Missing locus table: ", locus_path)
  dt <- fread(locus_path)
  ld <- read_ld(region_id)
  keep <- intersect(dt$snp, rownames(ld))
  dt <- dt[match(keep, snp)]
  ld <- ld[keep, keep, drop = FALSE]
  dt <- dt[
    is.finite(ad_beta) & is.finite(ad_varbeta) &
      is.finite(trait_beta_aligned_to_ad_a1) & is.finite(trait_varbeta) &
      is.finite(ad_maf) & is.finite(trait_maf) &
      ad_varbeta > 0 & trait_varbeta > 0 &
      ad_maf > 0 & ad_maf < 1 & trait_maf > 0 & trait_maf < 1
  ]
  ld <- ld[dt$snp, dt$snp, drop = FALSE]
  dt[, min_p_any := pmin(ad_p, trait_p, na.rm = TRUE)]
  if (nrow(dt) > max_variants) {
    dt <- dt[order(min_p_any)][seq_len(max_variants)]
    ld <- ld[dt$snp, dt$snp, drop = FALSE]
  }
  list(dt = dt, ld = ld)
}

extract_signal_summary <- function(fit, axis, region_id, trait) {
  if (is.null(fit$summary)) return(data.table())
  x <- as.data.table(fit$summary)
  if (!nrow(x)) return(data.table())
  x[, `:=`(axis = axis, region_id = region_id, trait = trait)]
  setcolorder(x, c("axis", "region_id", "trait", setdiff(names(x), c("axis", "region_id", "trait"))))
  x
}

run_one <- function(axis, region_id, trait) {
  base <- data.table(axis = axis, region_id = region_id, trait = trait)
  prepared <- tryCatch(prep_locus(region_id), error = function(e) e)
  if (inherits(prepared, "error")) {
    return(list(result = cbind(base, data.table(status = "input_failed", nsnps = NA_integer_, best_susie_pph4 = NA_real_, note = conditionMessage(prepared))),
                signals = data.table(), variants = data.table()))
  }
  dt <- prepared$dt
  ld <- prepared$ld
  variants <- dt[, .(
    axis, region_id, trait, snp, chr, bp_ad_build, bp_trait_build,
    ad_a1, ad_a2, trait_a1, trait_a2, allele_flip,
    ad_beta, ad_varbeta, ad_p, ad_maf, ad_n, ad_s,
    trait_beta_aligned_to_ad_a1, trait_varbeta, trait_p, trait_maf, trait_n,
    min_p_any
  )]
  if (nrow(dt) < min_variants) {
    return(list(result = cbind(base, data.table(status = "too_few_variants", nsnps = nrow(dt), best_susie_pph4 = NA_real_, note = "Too few SNPs after LD intersection/QC")),
                signals = data.table(), variants = variants))
  }
  d1 <- make_dataset(dt, "ad_beta", "ad_varbeta", "ad_maf", "ad_n", "cc", "ad_s", ld)
  d2 <- make_dataset(dt, "trait_beta_aligned_to_ad_a1", "trait_varbeta", "trait_maf", "trait_n", "quant", NULL, ld)
  fit <- tryCatch({
    s1 <- coloc::runsusie(d1)
    s2 <- coloc::runsusie(d2)
    coloc::coloc.susie(s1, s2)
  }, error = function(e) e)
  if (inherits(fit, "error")) {
    return(list(result = cbind(base, data.table(status = "susie_failed", nsnps = nrow(dt), best_susie_pph4 = NA_real_, note = conditionMessage(fit))),
                signals = data.table(), variants = variants))
  }
  signals <- extract_signal_summary(fit, axis, region_id, trait)
  pp4_cols <- intersect(c("PP.H4.abf", "PP.H4"), names(signals))
  best_pph4 <- if (length(pp4_cols)) max(signals[[pp4_cols[1]]], na.rm = TRUE) else NA_real_
  if (!is.finite(best_pph4)) best_pph4 <- NA_real_
  list(
    result = cbind(base, data.table(
      status = "completed",
      nsnps = nrow(dt),
      n_signal_pairs = nrow(signals),
      best_susie_pph4 = best_pph4,
      note = "LD-aware coloc.susie completed with local 1000G EUR LD"
    )),
    signals = signals,
    variants = variants
  )
}

outs <- vector("list", nrow(TARGETS))
for (i in seq_len(nrow(TARGETS))) {
  message("Running ", TARGETS$region_id[i])
  outs[[i]] <- run_one(TARGETS$axis[i], TARGETS$region_id[i], TARGETS$trait[i])
}

results <- rbindlist(lapply(outs, `[[`, "result"), fill = TRUE)
signals <- rbindlist(lapply(outs, `[[`, "signals"), fill = TRUE)
variants <- rbindlist(lapply(outs, `[[`, "variants"), fill = TRUE)

fwrite(results, OUT_RESULTS, sep = "\t")
fwrite(signals, OUT_SIGNALS, sep = "\t")
fwrite(variants, OUT_VARIANTS, sep = "\t")

status_counts <- results[, .N, by = status][order(status)]
lines <- c(
  "# Preliminary Result 28: LD-aware coloc.susie",
  "",
  "## Purpose",
  "",
  "This iteration replaces the previous unavailable/identity-LD SuSiE sensitivity step with local signed LD matrices built from the 1000 Genomes Phase 3 EUR PLINK reference.",
  "",
  "## Targets",
  "",
  paste0("- ", results$axis, " / ", results$region_id, " / ", results$trait, ": ", results$status, ", SNPs=", results$nsnps, ", best PP.H4=", signif(results$best_susie_pph4, 4)),
  "",
  "## Status counts",
  "",
  paste(capture.output(print(status_counts)), collapse = "\n"),
  "",
  "## Interpretation rule",
  "",
  "Treat this as a sensitivity layer over coloc.abf and QTL colocalization, not as a replacement. Strong support requires completed model fitting plus high signal-level PP.H4 concordant with the existing ABF and QTL evidence."
)
writeLines(lines, OUT_REPORT)
cat("Saved iteration 7 LD-aware coloc.susie outputs.\n")
