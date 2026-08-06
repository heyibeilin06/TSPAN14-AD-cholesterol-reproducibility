$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $PSScriptRoot
$jobs = @(
  @('01', '02_plot_panels_A_B_C_D.R'),
  @('02', '02_plot_panels_A_B_C_D.R'),
  @('03', '02_plot_panels_A_B_C_D.R'),
  @('04', '02_plot_panels_A_B_C_D_E_F.R'),
  @('05', '02_plot_panels_A_B_C_D_E.R'),
  @('06', '02_plot_panels_MODEL.R')
)
foreach ($job in $jobs) {
  $id = $job[0]
  $script = Join-Path $root "figures/Figure_$id/code/$($job[1])"
  $source = Join-Path $root "figures/Figure_$id/data"
  $output = Join-Path $root "reproduced/Figure_$id"
  New-Item -ItemType Directory -Force -Path $output | Out-Null
  & Rscript $script --project-root $root --source-dir $source --output-dir $output
  if ($LASTEXITCODE -ne 0) { throw "Figure $id failed" }
}
