suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag) {
  index <- match(flag, args)
  if (is.na(index) || index == length(args)) stop("Missing ", flag)
  args[[index + 1]]
}

qtl_path <- value_after("--qtl")
ad_path <- value_after("--ad")
ld_path <- value_after("--ld")
output_dir <- value_after("--output-dir")
qtl_sample_size <- as.integer(value_after("--qtl-sample-size"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

qtl <- fread(qtl_path)
ad <- fread(ad_path)
ld_table <- fread(ld_path)
ld <- as.matrix(ld_table[, -1, with = FALSE])
rownames(ld) <- ld_table[[1]]
colnames(ld) <- names(ld_table)[-1]
storage.mode(ld) <- "double"
setnames(qtl, "SNP", "snp")
merged <- merge(qtl, ad, by = "snp")
merged[, allele_factor := fifelse(A1 == ad_a1 & A2 == ad_a2, 1,
                           fifelse(A1 == ad_a2 & A2 == ad_a1, -1, NA_real_))]
merged[, palindromic := (A1 == "A" & A2 == "T") | (A1 == "T" & A2 == "A") |
                       (A1 == "C" & A2 == "G") | (A1 == "G" & A2 == "C")]
input <- merged[
  is.finite(allele_factor) & !palindromic & snp %in% rownames(ld) &
    is.finite(b) & is.finite(SE) & SE > 0 &
    is.finite(ad_beta) & is.finite(ad_varbeta) & ad_varbeta > 0 &
    is.finite(ad_maf) & ad_maf > 0 & ad_maf < 1
]
input[, sqtl_beta_aligned_to_ad_a1 := b * allele_factor]
input[, sqtl_varbeta := SE^2]
setorder(input, snp)
if (nrow(input) < 100) stop("Too few full cis-sQTL SNPs after QC for SuSiE.")
ld <- ld[input$snp, input$snp, drop = FALSE]
ld <- (ld + t(ld)) / 2
diag(ld) <- 1
d_ad <- list(beta = input$ad_beta, varbeta = input$ad_varbeta, snp = input$snp, MAF = input$ad_maf,
             N = as.integer(round(median(input$ad_n))), type = "cc", s = median(input$ad_s), LD = ld)
d_sqtl <- list(beta = input$sqtl_beta_aligned_to_ad_a1, varbeta = input$sqtl_varbeta, snp = input$snp, MAF = input$ad_maf,
               N = qtl_sample_size, type = "quant", LD = ld)
fit <- coloc.susie(runsusie(d_ad), runsusie(d_sqtl))
summary <- as.data.table(fit$summary)
summary[, `:=`(feature = unique(input$Probe), gene = unique(input$Gene), source_summary_scope = "full_cis_sQTL",
               qtl_sample_size = qtl_sample_size, n_input_snps = nrow(input), analysis = "coloc.susie")]
fwrite(summary, file.path(output_dir, "p1_exact_ba24_ad_sqtl_susie.tsv"), sep = "\t")
cat("Saved full_cis_sQTL LD-aware SuSiE colocalization output.\n")
