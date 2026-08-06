#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 OUTPUT_DIR TOOL_PREFIX" >&2
  exit 2
fi

output_dir="$1"
tool_prefix="$2"
mkdir -p "$output_dir"

panel_url="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
vcf_url="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr10.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
panel="$output_dir/1000G_phase3_sample_panel.tsv"
samples="$output_dir/1000G_phase3_EUR_samples.txt"
vcf="$output_dir/TSPAN14_functional_variants_1000G_EUR.vcf.gz"

curl -L --fail --silent --show-error "$panel_url" -o "$panel"
awk -F '\t' '$3 == "EUR" {print $1}' "$panel" > "$samples"

"$tool_prefix/bcftools" view \
  -r 10:82251000-82270500 \
  -S "$samples" \
  -i 'POS==82251544 || POS==82269461 || POS==82269611 || POS==82269848' \
  -Ou "$vcf_url" | \
  "$tool_prefix/bcftools" annotate --set-id '%CHROM:%POS:%REF:%FIRST_ALT' -Oz -o "$vcf"
"$tool_prefix/bcftools" index -f "$vcf"
"$tool_prefix/bcftools" query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n' "$vcf" > "$output_dir/TSPAN14_functional_variants_phased_genotypes.tsv"
