suppressPackageStartupMessages({
  library(data.table)
  library(hyprcoloc)
})

args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args0, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else file.path("scripts", "19_run_iteration7_ldaware_hyprcoloc.R")
ROOT <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
if (!dir.exists(file.path(ROOT, "results"))) ROOT <- normalizePath(getwd(), mustWork = TRUE)

TABLES <- file.path(ROOT, "results", "tables")
LOCUS_DIR <- file.path(TABLES, "coloc_loci_lead250kb")
LD_DIR <- file.path(ROOT, "results", "ld", "iteration7")
OUT_RESULTS <- file.path(TABLES, "technical_iteration7_ldaware_hyprcoloc_results.tsv")
OUT_SCORES <- file.path(TABLES, "technical_iteration7_ldaware_hyprcoloc_snp_scores.tsv")
OUT_INPUTS <- file.path(TABLES, "technical_iteration7_ldaware_hyprcoloc_inputs.tsv")
OUT_REPORT <- file.path(ROOT, "results", "reports", "preliminary_result_29_iteration7_ldaware_hyprcoloc.md")

max_variants <- as.integer(Sys.getenv("ITER7_HYPR_MAX_VARIANTS", "1200"))

axis_specs <- list(
  list(
    candidate_axis = "chr10/TSPAN14",
    lead_snp = "rs1902660",
    ld_region = "TC_nonAPOE_004",
    traits = c("AD", "TC", "LDL", "nonHDL"),
    files = c(AD = "TC_nonAPOE_004", TC = "TC_nonAPOE_004", LDL = "LDL_nonAPOE_004", nonHDL = "nonHDL_nonAPOE_004")
  ),
  list(
    candidate_axis = "chr11/MS4A_boundary",
    lead_snp = "rs1582763",
    ld_region = "TG_nonAPOE_004",
    traits = c("AD", "TG"),
    files = c(AD = "TG_nonAPOE_004", TG = "TG_nonAPOE_004")
  )
)

read_ld <- function(region_id) {
  path <- file.path(LD_DIR, paste0(region_id, ".ld.tsv"))
  m <- as.matrix(fread(path), rownames = 1)
  storage.mode(m) <- "double"
  m[!is.finite(m)] <- 0
  m <- (m + t(m)) / 2
  diag(m) <- 1
  m
}

read_trait <- function(region_id, trait) {
  path <- file.path(LOCUS_DIR, paste0(region_id, ".tsv"))
  x <- fread(path)
  if (trait == "AD") {
    x[, .(snp, trait = "AD", beta = ad_beta, se = sqrt(ad_varbeta), p = ad_p)]
  } else {
    x[, .(snp, trait = trait, beta = trait_beta_aligned_to_ad_a1, se = sqrt(trait_varbeta), p = trait_p)]
  }
}

prep_axis <- function(spec) {
  long <- rbindlist(Map(function(region_id, trait) read_trait(region_id, trait), spec$files, names(spec$files)), fill = TRUE)
  long <- unique(long[is.finite(beta) & is.finite(se) & se > 0], by = c("snp", "trait"))
  beta_w <- dcast(long, snp ~ trait, value.var = "beta")
  se_w <- dcast(long, snp ~ trait, value.var = "se")
  p_w <- dcast(long, snp ~ trait, value.var = "p")
  setnames(beta_w, spec$traits, paste0("beta__", spec$traits), skip_absent = TRUE)
  setnames(se_w, spec$traits, paste0("se__", spec$traits), skip_absent = TRUE)
  setnames(p_w, spec$traits, paste0("p__", spec$traits), skip_absent = TRUE)
  wide <- Reduce(function(a, b) merge(a, b, by = "snp", all = FALSE), list(beta_w, se_w, p_w))
  needed <- c(paste0("beta__", spec$traits), paste0("se__", spec$traits))
  for (nm in needed) wide <- wide[is.finite(get(nm))]
  ld <- read_ld(spec$ld_region)
  keep <- intersect(wide$snp, rownames(ld))
  wide <- wide[match(keep, snp)]
  pcols <- paste0("p__", spec$traits)
  wide[, min_p_any := do.call(pmin, c(.SD, list(na.rm = TRUE))), .SDcols = pcols]
  wide <- wide[order(min_p_any)]
  if (nrow(wide) > max_variants) wide <- wide[seq_len(max_variants)]
  ld <- ld[wide$snp, wide$snp, drop = FALSE]
  betas <- as.matrix(wide[, paste0("beta__", spec$traits), with = FALSE])
  ses <- as.matrix(wide[, paste0("se__", spec$traits), with = FALSE])
  colnames(betas) <- spec$traits
  colnames(ses) <- spec$traits
  rownames(betas) <- wide$snp
  rownames(ses) <- wide$snp
  list(wide = wide, ld = ld, betas = betas, ses = ses)
}

scores_to_table <- function(scores, spec) {
  if (is.null(scores) || !length(scores)) return(data.table())
  rbindlist(lapply(names(scores), function(cluster_id) {
    v <- scores[[cluster_id]]
    data.table(
      candidate_axis = spec$candidate_axis,
      lead_snp = spec$lead_snp,
      cluster_id = cluster_id,
      snp = names(v),
      snp_score = as.numeric(v)
    )
  }), fill = TRUE)
}

run_axis <- function(spec) {
  prepared <- prep_axis(spec)
  res <- tryCatch(
    hyprcoloc(
      effect.est = prepared$betas,
      effect.se = prepared$ses,
      binary.outcomes = ifelse(spec$traits == "AD", 1, 0),
      trait.names = spec$traits,
      snp.id = rownames(prepared$betas),
      ld.matrix = prepared$ld,
      snpscores = TRUE
    ),
    error = function(e) e
  )
  inputs <- copy(prepared$wide)
  inputs[, `:=`(
    candidate_axis = spec$candidate_axis,
    lead_snp = spec$lead_snp,
    traits_requested = paste(spec$traits, collapse = ";")
  )]
  setcolorder(inputs, c("candidate_axis", "lead_snp", "traits_requested", setdiff(names(inputs), c("candidate_axis", "lead_snp", "traits_requested"))))
  if (inherits(res, "error")) {
    result <- data.table(
      candidate_axis = spec$candidate_axis,
      lead_snp = spec$lead_snp,
      n_snps = nrow(prepared$wide),
      traits_requested = paste(spec$traits, collapse = ";"),
      hyprcoloc_version = as.character(packageVersion("hyprcoloc")),
      iteration = 7,
      status = "hyprcoloc_failed",
      traits = NA_character_,
      posterior_prob = NA_real_,
      regional_prob = NA_real_,
      candidate_snp = NA_character_,
      posterior_explained_by_snp = NA_real_,
      dropped_trait = NA_character_,
      note = conditionMessage(res)
    )
    return(list(results = result, scores = data.table(), inputs = inputs))
  }
  result <- as.data.table(res$results)
  if (!nrow(result)) {
    result <- data.table(
      traits = "None",
      posterior_prob = NA_real_,
      regional_prob = NA_real_,
      candidate_snp = NA_character_,
      posterior_explained_by_snp = NA_real_,
      dropped_trait = NA_character_
    )
  }
  result[, `:=`(
    candidate_axis = spec$candidate_axis,
    lead_snp = spec$lead_snp,
    n_snps = nrow(prepared$wide),
    traits_requested = paste(spec$traits, collapse = ";"),
    hyprcoloc_version = as.character(packageVersion("hyprcoloc")),
    iteration = 7,
    status = "completed",
    note = "LD-aware HyPrColoc completed with local 1000G EUR LD"
  )]
  setcolorder(result, c("candidate_axis", "lead_snp", "n_snps", "traits_requested", "hyprcoloc_version", "iteration", "status", setdiff(names(result), c("candidate_axis", "lead_snp", "n_snps", "traits_requested", "hyprcoloc_version", "iteration", "status"))))
  scores <- scores_to_table(res$snpscores, spec)
  if (nrow(scores)) scores <- scores[order(candidate_axis, cluster_id, -snp_score)]
  list(results = result, scores = scores, inputs = inputs)
}

outs <- lapply(axis_specs, run_axis)
results <- rbindlist(lapply(outs, `[[`, "results"), fill = TRUE)
scores <- rbindlist(lapply(outs, `[[`, "scores"), fill = TRUE)
inputs <- rbindlist(lapply(outs, `[[`, "inputs"), fill = TRUE)

fwrite(results, OUT_RESULTS, sep = "\t")
fwrite(scores, OUT_SCORES, sep = "\t")
fwrite(inputs, OUT_INPUTS, sep = "\t")

lines <- c(
  "# Preliminary Result 29: LD-aware HyPrColoc",
  "",
  "## Purpose",
  "",
  "This sensitivity analysis reruns the formal multitrait colocalization step with local signed 1000 Genomes Phase 3 EUR LD instead of the iteration-6 identity-LD placeholder.",
  "",
  "## Result summary",
  "",
  paste0("- ", results$candidate_axis, ": status=", results$status, ", traits=", results$traits, ", posterior=", signif(results$posterior_prob, 4), ", candidate=", results$candidate_snp, ", dropped=", results$dropped_trait),
  "",
  "## Interpretation",
  "",
  "A supported main axis should retain a coherent multitrait cluster under LD-aware HyPrColoc. Pairwise LD-aware coloc.susie signals that do not survive this multitrait layer should be treated as boundary or sensitivity evidence rather than as the primary disease-lipid mechanism."
)
writeLines(lines, OUT_REPORT)
cat("Saved iteration 7 LD-aware HyPrColoc outputs.\n")
