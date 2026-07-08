param(
    [string]$ProjectRoot = ""
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

if ([string]::IsNullOrWhiteSpace($ProjectRoot)) {
    $ProjectRoot = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot "..\.."))
} else {
    $ProjectRoot = [System.IO.Path]::GetFullPath($ProjectRoot)
}

if (-not (Test-Path -LiteralPath $ProjectRoot)) {
    throw "Project root not found: $ProjectRoot"
}

$wslRoot = (& wsl.exe --exec wslpath -a $ProjectRoot).Trim()
$quotedRoot = $wslRoot.Replace("'", "'\''")

$cmd = "cd '$quotedRoot' && source scripts/bash/env_common.sh && run_env Rscript scripts/R/06_nextstep_annotation_modules.R && run_env Rscript scripts/R/07_qc_sensitivity.R && run_env Rscript scripts/R/08_make_interpretable_shortlist.R && run_env Rscript scripts/R/09_make_analysis_digest.R && run_env Rscript scripts/R/10_make_summary_tables.R && run_env Rscript scripts/R/11_make_analysis_table_workbook.R"

wsl.exe --exec bash -lc $cmd
if ($LASTEXITCODE -ne 0) { throw "Next-step annotation/QC/candidate synthesis failed" }

Write-Host "DONE. Open: $ProjectRoot\results\analysis_tables"
