suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})

# P0 sensitivity analysis: use the rebuilt 1000 Genomes EUR LD matrices to
# test whether local TSPAN14 colocalization is stable to SNP-set size and SuSiE
# complexity. It does not substitute for a chromosome-19-excluded LDSC rerun.
args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args0, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else file.path("scripts", "166_p0_tspan14_coloc_susie_sensitivity.R")
ROOT <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
if (!dir.exists(file.path(ROOT, "results"))) ROOT <- normalizePath(getwd(), mustWork = TRUE)

LOCUS_DIR <- file.path(ROOT, "results", "tables", "coloc_loci_lead250kb")
APOE_MANIFEST <- file.path(ROOT, "results", "tables", "coloc_loci_apoe250kb", "manifest.tsv")
LD_DIR <- Sys.getenv("P0_LD_DIR", unset = "D:/SLM_AD_Lipid_data/p0_ld_1000g_eur/matrices")
OUT_DIR <- file.path(ROOT, "outputs", "p0_reanalysis")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_RESULTS <- file.path(OUT_DIR, "p0_coloc_susie_sensitivity.tsv")
OUT_QC <- file.path(OUT_DIR, "p0_coloc_susie_input_qc.tsv")
OUT_APOE <- file.path(OUT_DIR, "p0_apoe_scope_audit.tsv")

# Prespecified low/high settings keep the full 3-trait sensitivity rerun
# tractable while testing both locus truncation and sparse-effect complexity.
MAX_VARIANT_GRID <- as.integer(strsplit(Sys.getenv("P0_MAX_GRID", unset = "500,1200"), ",", fixed = TRUE)[[1]])
SUSIE_L_GRID <- as.integer(strsplit(Sys.getenv("P0_SUSIE_L_GRID", unset = "5,10"), ",", fixed = TRUE)[[1]])
TARGETS <- data.table(
  trait = c("TC", "LDL", "nonHDL"),
  locus = c("TC_nonAPOE_004", "LDL_nonAPOE_004", "nonHDL_nonAPOE_004")
)
target_filter <- Sys.getenv("P0_TARGET_TRAIT", unset = "")
if (nzchar(target_filter)) TARGETS <- TARGETS[trait == target_filter]

read_ld <- function(locus) {
  path <- file.path(LD_DIR, paste0(locus, ".ld.tsv"))
  if (!file.exists(path)) stop("Missing rebuilt EUR LD matrix: ", path)
  ld <- as.matrix(fread(path), rownames = 1)
  storage.mode(ld) <- "double"
  ld <- (ld + t(ld)) / 2
  ld[!is.finite(ld)] <- 0
  diag(ld) <- 1
  ld
}

prepare_locus <- function(locus, n_keep) {
  x <- fread(file.path(LOCUS_DIR, paste0(locus, ".tsv")))
  x <- x[is.finite(ad_beta) & is.finite(ad_varbeta) & ad_varbeta > 0 & is.finite(ad_maf) & ad_maf > 0 & ad_maf < 1 & is.finite(ad_n) & is.finite(ad_s) & is.finite(trait_beta_aligned_to_ad_a1) & is.finite(trait_varbeta) & trait_varbeta > 0 & is.finite(trait_maf) & trait_maf > 0 & trait_maf < 1 & is.finite(trait_n)]
  ld <- read_ld(locus)
  x <- x[snp %in% rownames(ld)]
  x[, min_p := pmin(ad_p, trait_p, na.rm = TRUE)]
  setorder(x, min_p)
  x <- x[seq_len(min(n_keep, .N))]
  ld <- ld[x$snp, x$snp, drop = FALSE]
  list(data = x, ld = ld)
}

make_dataset <- function(x, kind, ld) {
  if (kind == "AD") {
    list(beta = x$ad_beta, varbeta = x$ad_varbeta, snp = x$snp, MAF = x$ad_maf, N = as.integer(round(median(x$ad_n))), type = "cc", s = median(x$ad_s), LD = ld)
  } else {
    list(beta = x$trait_beta_aligned_to_ad_a1, varbeta = x$trait_varbeta, snp = x$snp, MAF = x$trait_maf, N = as.integer(round(median(x$trait_n))), type = "quant", LD = ld)
  }
}

run_one <- function(trait, locus, n_keep, susie_l) {
  base <- data.table(trait = trait, locus = locus, max_variants = n_keep, susie_L = susie_l, ld_orientation = "signed_to_AD_A1")
  prepared <- tryCatch(prepare_locus(locus, n_keep), error = function(e) e)
  if (inherits(prepared, "error")) return(list(result = cbind(base, data.table(status = "input_failed", best_pp_h4 = NA_real_, n_signal_pairs = NA_integer_, note = conditionMessage(prepared))), qc = cbind(base, data.table(n_snps = NA_integer_, note = conditionMessage(prepared)))))
  x <- prepared$data
  qc <- cbind(base, data.table(n_snps = nrow(x), n_missing_ld = sum(!is.finite(prepared$ld)), min_p = min(x$min_p), note = "Rebuilt 1000G EUR LD"))
  fit <- tryCatch({
    s1 <- runsusie(make_dataset(x, "AD", prepared$ld), L = susie_l)
    s2 <- runsusie(make_dataset(x, trait, prepared$ld), L = susie_l)
    coloc.susie(s1, s2)
  }, error = function(e) e)
  if (inherits(fit, "error")) return(list(result = cbind(base, data.table(status = "susie_failed", best_pp_h4 = NA_real_, n_signal_pairs = NA_integer_, note = conditionMessage(fit))), qc = qc))
  signals <- as.data.table(fit$summary)
  pp_column <- intersect(c("PP.H4.abf", "PP.H4"), names(signals))
  best <- if (length(pp_column) && nrow(signals)) max(signals[[pp_column[1]]], na.rm = TRUE) else NA_real_
  list(result = cbind(base, data.table(status = "completed", best_pp_h4 = best, n_signal_pairs = nrow(signals), note = "coloc.susie sensitivity with rebuilt 1000G EUR LD")), qc = qc)
}

grid <- CJ(target_row = seq_len(nrow(TARGETS)), max_variants = MAX_VARIANT_GRID, susie_L = SUSIE_L_GRID)
outputs <- lapply(seq_len(nrow(grid)), function(i) {
  target <- TARGETS[grid$target_row[i]]
  run_one(target$trait, target$locus, grid$max_variants[i], grid$susie_L[i])
})
results <- rbindlist(lapply(outputs, `[[`, "result"), fill = TRUE)
qc <- rbindlist(lapply(outputs, `[[`, "qc"), fill = TRUE)
if (file.exists(OUT_RESULTS)) {
  existing <- fread(OUT_RESULTS)
  if ("ld_orientation" %in% names(existing)) existing <- existing[ld_orientation == "signed_to_AD_A1"] else existing <- existing[0]
  results <- rbind(existing, results, fill = TRUE)
}
if (file.exists(OUT_QC)) {
  existing <- fread(OUT_QC)
  if ("ld_orientation" %in% names(existing)) existing <- existing[ld_orientation == "signed_to_AD_A1"] else existing <- existing[0]
  qc <- rbind(existing, qc, fill = TRUE)
}
setorder(results, trait, max_variants, susie_L)
setorder(qc, trait, max_variants, susie_L)
results <- unique(results, by = c("trait", "locus", "max_variants", "susie_L"), fromLast = TRUE)
qc <- unique(qc, by = c("trait", "locus", "max_variants", "susie_L"), fromLast = TRUE)
fwrite(results, OUT_RESULTS, sep = "\t")
fwrite(qc, OUT_QC, sep = "\t")

tspan <- fread(file.path(LOCUS_DIR, "TC_nonAPOE_004.tsv"))
tspan_start <- min(tspan$bp_ad_build, na.rm = TRUE)
tspan_end <- max(tspan$bp_ad_build, na.rm = TRUE)
apoe <- fread(APOE_MANIFEST)
apoe_scope <- data.table(
  apoe_window_label = c("retained_APOE250kb_manifest", "plus_or_minus_500kb", "plus_or_minus_1Mb", "plus_or_minus_2Mb"),
  APOE_manifest_rows = nrow(apoe),
  TSPAN14_chromosome = unique(tspan$chr)[1],
  TSPAN14_start_hg38 = tspan_start,
  TSPAN14_end_hg38 = tspan_end,
  TSPAN14_excluded = FALSE,
  interpretation = "TSPAN14 is on chromosome 10 and is geometrically unaffected by any chromosome-19 APOE exclusion window. This scope audit does not rerun chromosome-19-excluded genome-wide LDSC."
)
fwrite(apoe_scope, OUT_APOE, sep = "\t")
cat("Wrote P0 coloc.susie sensitivity outputs to ", OUT_DIR, "\n", sep = "")
