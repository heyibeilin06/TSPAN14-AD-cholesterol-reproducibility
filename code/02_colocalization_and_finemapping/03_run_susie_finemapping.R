suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})

args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args0, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else file.path("scripts", "22_run_iteration7_ldaware_susie_finemapping.R")
ROOT <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
if (!dir.exists(file.path(ROOT, "results"))) ROOT <- normalizePath(getwd(), mustWork = TRUE)

TARGETS <- data.table(
  axis = c("TSPAN14", "TSPAN14", "TSPAN14", "MS4A_boundary"),
  region_id = c("TC_nonAPOE_004", "LDL_nonAPOE_004", "nonHDL_nonAPOE_004", "TG_nonAPOE_004"),
  trait = c("TC", "LDL", "nonHDL", "TG")
)

LOCUS_DIR <- file.path(ROOT, "results", "tables", "coloc_loci_lead250kb")
LD_DIR <- file.path(ROOT, "results", "ld", "iteration7")
OUT_SIGNALS <- file.path(ROOT, "results", "tables", "technical_iteration7_ldaware_susie_finemap_signals.tsv")
OUT_VARIANTS <- file.path(ROOT, "results", "tables", "technical_iteration7_ldaware_susie_finemap_variants.tsv")
OUT_REPORT <- file.path(ROOT, "results", "reports", "preliminary_result_33_iteration7_ldaware_susie_finemapping.md")

max_variants <- as.integer(Sys.getenv("ITER7_FINEMAP_MAX_VARIANTS", "1200"))

read_ld <- function(region_id) {
  x <- fread(file.path(LD_DIR, paste0(region_id, ".ld.tsv")))
  rn <- x[[1]]
  x[[1]] <- NULL
  m <- as.matrix(x)
  rownames(m) <- rn
  colnames(m) <- rn
  storage.mode(m) <- "double"
  m[!is.finite(m)] <- 0
  m <- (m + t(m)) / 2
  diag(m) <- 1
  m
}

prep_locus <- function(region_id) {
  dt <- fread(file.path(LOCUS_DIR, paste0(region_id, ".tsv")))
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

make_dataset <- function(dt, outcome, R) {
  if (outcome == "AD") {
    list(
      beta = dt$ad_beta,
      varbeta = dt$ad_varbeta,
      snp = dt$snp,
      MAF = dt$ad_maf,
      N = as.integer(round(stats::median(dt$ad_n, na.rm = TRUE))),
      s = stats::median(dt$ad_s, na.rm = TRUE),
      type = "cc",
      LD = R
    )
  } else {
    list(
      beta = dt$trait_beta_aligned_to_ad_a1,
      varbeta = dt$trait_varbeta,
      snp = dt$snp,
      MAF = dt$trait_maf,
      N = as.integer(round(stats::median(dt$trait_n, na.rm = TRUE))),
      type = "quant",
      LD = R
    )
  }
}

cs_membership <- function(fit, snps) {
  out <- data.table(snp = snps, credible_set = NA_character_)
  if (is.null(fit$sets) || is.null(fit$sets$cs) || !length(fit$sets$cs)) return(out)
  for (nm in names(fit$sets$cs)) {
    idx <- fit$sets$cs[[nm]]
    out[snp %in% snps[idx], credible_set := fifelse(is.na(credible_set), nm, paste(credible_set, nm, sep = ";"))]
  }
  out
}

extract_fit <- function(fit, meta, dt) {
  pip <- fit$pip
  if (is.null(pip)) pip <- rep(NA_real_, nrow(dt))
  cm <- cs_membership(fit, dt$snp)
  variants <- data.table(
    axis = meta$axis,
    region_id = meta$region_id,
    trait = meta$trait,
    outcome = meta$outcome,
    snp = dt$snp,
    chr = dt$chr,
    bp_ad_build = dt$bp_ad_build,
    ad_p = dt$ad_p,
    trait_p = dt$trait_p,
    pip = as.numeric(pip)
  )
  variants <- merge(variants, cm, by = "snp", all.x = TRUE, sort = FALSE)
  setorder(variants, -pip)

  n_cs <- if (!is.null(fit$sets) && !is.null(fit$sets$cs)) length(fit$sets$cs) else 0
  signals <- if (n_cs > 0) {
    rbindlist(lapply(names(fit$sets$cs), function(cs_name) {
      members <- dt$snp[fit$sets$cs[[cs_name]]]
      v <- variants[snp %in% members][order(-pip)]
      data.table(
        axis = meta$axis,
        region_id = meta$region_id,
        trait = meta$trait,
        outcome = meta$outcome,
        credible_set = cs_name,
        n_variants = nrow(v),
        top_snp = v$snp[1],
        top_pip = v$pip[1],
        purity_min_abs_corr = if (!is.null(fit$sets$purity)) fit$sets$purity[cs_name, "min.abs.corr"] else NA_real_,
        purity_mean_abs_corr = if (!is.null(fit$sets$purity)) fit$sets$purity[cs_name, "mean.abs.corr"] else NA_real_
      )
    }), fill = TRUE)
  } else {
    data.table(
      axis = meta$axis,
      region_id = meta$region_id,
      trait = meta$trait,
      outcome = meta$outcome,
      credible_set = NA_character_,
      n_variants = 0L,
      top_snp = variants$snp[1],
      top_pip = variants$pip[1],
      purity_min_abs_corr = NA_real_,
      purity_mean_abs_corr = NA_real_
    )
  }
  list(signals = signals, variants = variants)
}

all_signals <- list()
all_variants <- list()
counter <- 1L
for (i in seq_len(nrow(TARGETS))) {
  prepared <- prep_locus(TARGETS$region_id[i])
  for (outcome in c("AD", TARGETS$trait[i])) {
    message("Fine-mapping ", TARGETS$region_id[i], " / ", outcome)
    d <- make_dataset(prepared$dt, ifelse(outcome == "AD", "AD", "trait"), prepared$ld)
    fit <- coloc::runsusie(d)
    extracted <- extract_fit(
      fit,
      list(axis = TARGETS$axis[i], region_id = TARGETS$region_id[i], trait = TARGETS$trait[i], outcome = outcome),
      prepared$dt
    )
    all_signals[[counter]] <- extracted$signals
    all_variants[[counter]] <- extracted$variants
    counter <- counter + 1L
  }
}

signals <- rbindlist(all_signals, fill = TRUE)
variants <- rbindlist(all_variants, fill = TRUE)
fwrite(signals, OUT_SIGNALS, sep = "\t")
fwrite(variants, OUT_VARIANTS, sep = "\t")

tspan <- signals[axis == "TSPAN14"]
lines <- c(
  "# Preliminary Result 33: LD-aware SuSiE fine-mapping",
  "",
  "## Purpose",
  "",
  "This table-level analysis adds an explicit LD-aware SuSiE fine-mapping layer over the same local 1000 Genomes EUR signed LD matrices used by the colocalization sensitivity analysis.",
  "",
  "## TSPAN14 signal summary",
  "",
  paste0("- ", tspan$region_id, " / ", tspan$outcome, " / ", tspan$credible_set, ": top=", tspan$top_snp, ", PIP=", signif(tspan$top_pip, 4), ", n=", tspan$n_variants),
  "",
  "## Outputs",
  "",
  "- `technical_iteration7_ldaware_susie_finemap_signals.tsv`",
  "- `technical_iteration7_ldaware_susie_finemap_variants.tsv`",
  "",
  "## Boundary",
  "",
  "This is locus-level fine-mapping with a 1000G EUR reference panel and should be interpreted as a sensitivity/refinement layer. It is not a replacement for individual-level genotype fine-mapping."
)
writeLines(lines, OUT_REPORT)
cat("Saved iteration 7 LD-aware SuSiE fine-mapping outputs.\n")
