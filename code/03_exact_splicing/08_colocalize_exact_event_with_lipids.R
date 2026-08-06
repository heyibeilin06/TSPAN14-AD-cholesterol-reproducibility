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
trait_path <- value_after("--trait")
trait_label <- value_after("--trait-label")
output_dir <- value_after("--output-dir")
qtl_sample_size <- as.integer(value_after("--qtl-sample-size"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

qtl <- fread(qtl_path)
trait <- fread(trait_path)
setnames(qtl, "SNP", "snp")
merged <- merge(qtl, trait, by = "snp")
merged[, allele_factor := fifelse(A1 == ad_a1 & A2 == ad_a2, 1,
                           fifelse(A1 == ad_a2 & A2 == ad_a1, -1, NA_real_))]
merged[, palindromic := (A1 == "A" & A2 == "T") | (A1 == "T" & A2 == "A") |
                       (A1 == "C" & A2 == "G") | (A1 == "G" & A2 == "C")]
input <- merged[
  is.finite(allele_factor) & !palindromic &
    is.finite(b) & is.finite(SE) & SE > 0 &
    is.finite(trait_beta_aligned_to_ad_a1) & is.finite(trait_varbeta) & trait_varbeta > 0 &
    is.finite(trait_maf) & trait_maf > 0 & trait_maf < 1
]
input[, sqtl_beta_aligned_to_ad_a1 := b * allele_factor]
input[, sqtl_varbeta := SE^2]
setorder(input, snp)
if (nrow(input) < 100) stop("Too few full cis-sQTL SNPs after allele QC for coloc.")

d_trait <- list(
  beta = input$trait_beta_aligned_to_ad_a1,
  varbeta = input$trait_varbeta,
  snp = input$snp,
  MAF = input$trait_maf,
  N = as.integer(round(median(input$trait_n))),
  type = "quant"
)
d_sqtl <- list(
  beta = input$sqtl_beta_aligned_to_ad_a1,
  varbeta = input$sqtl_varbeta,
  snp = input$snp,
  MAF = input$trait_maf,
  N = qtl_sample_size,
  type = "quant"
)
fit <- coloc.abf(d_trait, d_sqtl)
summary <- as.data.table(as.list(fit$summary))
summary[, `:=`(
  trait = trait_label,
  feature = unique(input$Probe),
  gene = unique(input$Gene),
  source_summary_scope = "full_cis_sQTL",
  qtl_sample_size = qtl_sample_size,
  n_input_snps = nrow(input),
  analysis = "coloc.abf"
)]
fwrite(summary, file.path(output_dir, paste0("p1_exact_ba24_", trait_label, "_sqtl_coloc.tsv")), sep = "\t")
cat("Saved full_cis_sQTL lipid colocalization output.\n")
