suppressPackageStartupMessages({
  library(data.table)
  library(hyprcoloc)
})

# P0 rerun: every posterior reported in the manuscript is generated as a
# separately named model with its own SNP count and source-table provenance.
args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args0, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else file.path("scripts", "163_p0_tspan14_hyprcoloc_reanalysis.R")
ROOT <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
if (!dir.exists(file.path(ROOT, "results"))) ROOT <- normalizePath(getwd(), mustWork = TRUE)

LOCUS_DIR <- file.path(ROOT, "results", "tables", "coloc_loci_lead250kb")
SQTL_FILE <- file.path(ROOT, "results", "tables", "microglia_full_qtl_target_extracts", "SVZ_eur_splicing_peer0_gene.cis_qtl_nominal_tabixed.candidate_extract.tsv")
SQTL_HARMONIZED_FILE <- Sys.getenv("P0_SQTL_HARMONIZED_FILE", unset = file.path(ROOT, "outputs", "p0_reanalysis", "p0_sqtl_full_harmonized_to_ad_a1.tsv"))
LD_DIR <- Sys.getenv("P0_LD_DIR", unset = "D:/SLM_AD_Lipid_data/p0_ld_1000g_eur/matrices")
OUT_DIR <- file.path(ROOT, "outputs", "p0_reanalysis")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_MANIFEST <- file.path(OUT_DIR, "p0_hyprcoloc_model_manifest.tsv")
OUT_RESULTS <- file.path(OUT_DIR, "p0_hyprcoloc_results.tsv")
OUT_QC <- file.path(OUT_DIR, "p0_hyprcoloc_input_qc.tsv")

TARGET_EVENT <- "chr10:80509471:80512144:clu_4260_+:ENSG00000108219.15"
MAX_VARIANTS <- as.integer(Sys.getenv("P0_HYPR_MAX_VARIANTS", "1200"))

MODEL_SPECS <- list(
  list(model_id = "AD_TC", traits = c("AD", "TC"), locus = "TC_nonAPOE_004", include_sqtl = FALSE),
  list(model_id = "AD_LDL", traits = c("AD", "LDL"), locus = "LDL_nonAPOE_004", include_sqtl = FALSE),
  list(model_id = "AD_nonHDL", traits = c("AD", "nonHDL"), locus = "nonHDL_nonAPOE_004", include_sqtl = FALSE),
  list(model_id = "AD_TC_sQTL", traits = c("AD", "TC", "sQTL"), locus = "TC_nonAPOE_004", include_sqtl = TRUE),
  list(model_id = "AD_TC_LDL_nonHDL", traits = c("AD", "TC", "LDL", "nonHDL"), locus = "TC_nonAPOE_004", include_sqtl = FALSE),
  list(model_id = "AD_TC_LDL_nonHDL_sQTL", traits = c("AD", "TC", "LDL", "nonHDL", "sQTL"), locus = "TC_nonAPOE_004", include_sqtl = TRUE)
)

read_ld <- function(locus) {
  path <- file.path(LD_DIR, paste0(locus, ".ld.tsv"))
  if (!file.exists(path)) stop("Missing P0 EUR LD matrix: ", path)
  ld <- as.matrix(fread(path), rownames = 1)
  storage.mode(ld) <- "double"
  ld <- (ld + t(ld)) / 2
  ld[!is.finite(ld)] <- 0
  diag(ld) <- 1
  ld
}

read_gwas <- function(locus, trait) {
  x <- fread(file.path(LOCUS_DIR, paste0(locus, ".tsv")))
  if (trait == "AD") {
    x[, .(snp, trait = "AD", beta = ad_beta, se = sqrt(ad_varbeta), p = ad_p, source = locus)]
  } else {
    x[, .(snp, trait = trait, beta = trait_beta_aligned_to_ad_a1, se = sqrt(trait_varbeta), p = trait_p, source = locus)]
  }
}

read_sqtl <- function() {
  if (file.exists(SQTL_HARMONIZED_FILE)) {
    x <- fread(SQTL_HARMONIZED_FILE)
    x <- x[alignment_status %in% c("ad_a1_is_miga_alt", "ad_a1_is_miga_ref_flipped") & is.finite(beta_aligned_to_ad_a1) & is.finite(slope_se) & slope_se > 0]
    return(x[, .(snp, trait = "sQTL", beta = beta_aligned_to_ad_a1, se = slope_se, p = pval_nominal, source = basename(SQTL_HARMONIZED_FILE))])
  }
  x <- fread(SQTL_FILE)
  x <- x[phenotype_id == TARGET_EVENT & is.finite(slope) & is.finite(slope_se) & slope_se > 0]
  x[, .(snp = variant_id, trait = "sQTL", beta = slope, se = slope_se, p = pval_nominal, source = basename(SQTL_FILE))]
}

prepare_model <- function(spec) {
  base_traits <- setdiff(spec$traits, "sQTL")
  long <- rbindlist(lapply(base_traits, function(trait) {
    locus <- if (trait == "AD") spec$locus else if (trait == "TC") "TC_nonAPOE_004" else if (trait == "LDL") "LDL_nonAPOE_004" else "nonHDL_nonAPOE_004"
    read_gwas(locus, trait)
  }))
  if (spec$include_sqtl) long <- rbind(long, read_sqtl())
  long <- unique(long[is.finite(beta) & is.finite(se) & se > 0], by = c("snp", "trait"))
  make_wide <- function(value_column, prefix) {
    out <- dcast(long, snp ~ trait, value.var = value_column)
    setnames(out, setdiff(names(out), "snp"), paste0(prefix, setdiff(names(out), "snp")))
    out
  }
  beta <- make_wide("beta", "beta__")
  se <- make_wide("se", "se__")
  p <- make_wide("p", "p__")
  wide <- Reduce(function(a, b) merge(a, b, by = "snp", all = FALSE), list(beta, se, p))
  complete_cols <- c(paste0("beta__", spec$traits), paste0("se__", spec$traits), paste0("p__", spec$traits))
  for (column in complete_cols) wide <- wide[is.finite(get(column))]
  ld <- read_ld(spec$locus)
  keep <- intersect(wide$snp, rownames(ld))
  wide <- wide[match(keep, snp)]
  wide[, min_p := do.call(pmin, c(.SD, list(na.rm = TRUE))), .SDcols = paste0("p__", spec$traits)]
  wide <- wide[order(min_p)]
  if (nrow(wide) > MAX_VARIANTS) wide <- wide[seq_len(MAX_VARIANTS)]
  ld <- ld[wide$snp, wide$snp, drop = FALSE]
  beta_m <- as.matrix(wide[, paste0("beta__", spec$traits), with = FALSE])
  se_m <- as.matrix(wide[, paste0("se__", spec$traits), with = FALSE])
  colnames(beta_m) <- spec$traits
  colnames(se_m) <- spec$traits
  rownames(beta_m) <- wide$snp
  rownames(se_m) <- wide$snp
  list(wide = wide, beta = beta_m, se = se_m, ld = ld)
}

run_model <- function(spec) {
  prepared <- tryCatch(prepare_model(spec), error = function(e) e)
  base <- data.table(model_id = spec$model_id, traits_requested = paste(spec$traits, collapse = ";"), locus = spec$locus, includes_project_sqtl = spec$include_sqtl)
  if (inherits(prepared, "error")) {
    return(list(result = cbind(base, data.table(status = "input_failed", traits_clustered = NA_character_, posterior_prob = NA_real_, regional_prob = NA_real_, candidate_snp = NA_character_, posterior_explained_by_snp = NA_real_, dropped_trait = NA_character_, note = conditionMessage(prepared))), qc = cbind(base, data.table(n_variants = NA_integer_, note = conditionMessage(prepared)))))
  }
  result <- tryCatch(hyprcoloc(effect.est = prepared$beta, effect.se = prepared$se, binary.outcomes = ifelse(spec$traits == "AD", 1, 0), trait.names = spec$traits, snp.id = rownames(prepared$beta), ld.matrix = prepared$ld, snpscores = TRUE), error = function(e) e)
  qc <- cbind(base, data.table(
    n_variants = nrow(prepared$wide),
    min_p = min(prepared$wide$min_p, na.rm = TRUE),
    ld_dimension = nrow(prepared$ld),
    n_missing_beta = sum(!is.finite(prepared$beta)),
    n_missing_se = sum(!is.finite(prepared$se)),
    n_missing_ld = sum(!is.finite(prepared$ld)),
    missing_beta_traits = paste(colnames(prepared$beta)[colSums(!is.finite(prepared$beta)) > 0], collapse = ";"),
    missing_se_traits = paste(colnames(prepared$se)[colSums(!is.finite(prepared$se)) > 0], collapse = ";"),
    note = "Input intersected with rebuilt 1000G EUR LD"
  ))
  if (inherits(result, "error")) {
    return(list(result = cbind(base, data.table(status = "hyprcoloc_failed", traits_clustered = NA_character_, posterior_prob = NA_real_, regional_prob = NA_real_, candidate_snp = NA_character_, posterior_explained_by_snp = NA_real_, dropped_trait = NA_character_, note = conditionMessage(result))), qc = qc))
  }
  res <- as.data.table(result$results)
  if (!nrow(res)) res <- data.table(traits = "None", posterior_prob = NA_real_, regional_prob = NA_real_, candidate_snp = NA_character_, posterior_explained_by_snp = NA_real_, dropped_trait = NA_character_)
  setnames(res, "traits", "traits_clustered", skip_absent = TRUE)
  note <- if (spec$include_sqtl && file.exists(SQTL_HARMONIZED_FILE)) "P0 rerun with rebuilt 1000G EUR LD and MiGA sQTL effects aligned to AD A1" else "P0 rerun with rebuilt 1000G EUR LD"
  list(result = cbind(base, data.table(status = "completed"), res, data.table(note = note)), qc = qc)
}

outputs <- lapply(MODEL_SPECS, run_model)
results <- rbindlist(lapply(outputs, `[[`, "result"), fill = TRUE)
qc <- rbindlist(lapply(outputs, `[[`, "qc"), fill = TRUE)
manifest <- rbindlist(lapply(MODEL_SPECS, function(spec) data.table(model_id = spec$model_id, traits_requested = paste(spec$traits, collapse = ";"), locus = spec$locus, ld_file = file.path(LD_DIR, paste0(spec$locus, ".ld.tsv")), sqtl_source = if (spec$include_sqtl) if (file.exists(SQTL_HARMONIZED_FILE)) SQTL_HARMONIZED_FILE else SQTL_FILE else NA_character_, max_variants = MAX_VARIANTS, script = normalizePath(script_path, mustWork = FALSE))), fill = TRUE)
fwrite(manifest, OUT_MANIFEST, sep = "\t")
fwrite(results, OUT_RESULTS, sep = "\t")
fwrite(qc, OUT_QC, sep = "\t")
cat("Wrote P0 HyPrColoc rerun outputs to ", OUT_DIR, "\n", sep = "")
