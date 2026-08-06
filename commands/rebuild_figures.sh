#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
figures=(01 02 03 04 05 06)
scripts=(02_plot_panels_A_B_C_D.R 02_plot_panels_A_B_C_D.R 02_plot_panels_A_B_C_D.R 02_plot_panels_A_B_C_D_E_F.R 02_plot_panels_A_B_C_D_E.R 02_plot_panels_MODEL.R)
for i in "${!figures[@]}"; do
  id="${figures[$i]}"
  out="$root/reproduced/Figure_$id"
  mkdir -p "$out"
  Rscript "$root/figures/Figure_$id/code/${scripts[$i]}" --project-root "$root" --source-dir "$root/figures/Figure_$id/data" --output-dir "$out"
done
