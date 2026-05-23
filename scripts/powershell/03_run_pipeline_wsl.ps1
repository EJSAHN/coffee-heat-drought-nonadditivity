<#
Run the coffee RNA-seq pipeline inside WSL.
Examples:
  .\scripts\powershell\03_run_pipeline_wsl.ps1 -Mode setup
  .\scripts\powershell\03_run_pipeline_wsl.ps1 -Mode all -Threads 8
#>
[CmdletBinding()]
param(
    [string]$ProjectRoot = "",
    [ValidateSet("setup","index","metadata","quant","deseq2","next","all")]
    [string]$Mode = "all",
    [int]$Threads = 8
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

function Invoke-WslBash {
    param([string]$Command)
    Write-Host "\n[WSL] $Command" -ForegroundColor Cyan
    & wsl bash -lc $Command
    if ($LASTEXITCODE -ne 0) { throw "WSL command failed with exit code $LASTEXITCODE" }
}

try {
    $wslRoot = (& wsl.exe --exec wslpath -a $ProjectRoot).Trim()
} catch {
    throw "WSL is not available. Install Ubuntu with: wsl --install -d Ubuntu"
}

$quotedRoot = $wslRoot.Replace("'", "'\''")
Write-Host "Windows project root: $ProjectRoot" -ForegroundColor Green
Write-Host "WSL project root:     $wslRoot" -ForegroundColor Green

switch ($Mode) {
    "setup"    { Invoke-WslBash "cd '$quotedRoot' && bash scripts/bash/01_setup_env.sh" }
    "index"    { Invoke-WslBash "cd '$quotedRoot' && bash scripts/bash/02_make_indices_and_tx2gene.sh" }
    "metadata" { Invoke-WslBash "cd '$quotedRoot' && bash scripts/bash/03_make_metadata.sh" }
    "quant"    { Invoke-WslBash "cd '$quotedRoot' && THREADS=$Threads bash scripts/bash/04_fastqc_salmon_resume.sh" }
    "deseq2"   { Invoke-WslBash "cd '$quotedRoot' && bash scripts/bash/05_deseq2_nonadditive.sh" }
    "next"     {
        Invoke-WslBash "cd '$quotedRoot' && source scripts/bash/env_common.sh && run_env Rscript scripts/R/06_nextstep_annotation_modules.R && run_env Rscript scripts/R/07_qc_sensitivity.R && run_env Rscript scripts/R/08_make_interpretable_shortlist.R && run_env Rscript scripts/R/09_make_paper_digest.R && run_env Rscript scripts/R/10_make_final_story_tables.R && run_env Rscript scripts/R/11_make_analysis_table_workbook.R"
    }
    "all"      {
        Invoke-WslBash "cd '$quotedRoot' && bash scripts/bash/01_setup_env.sh"
        Invoke-WslBash "cd '$quotedRoot' && bash scripts/bash/02_make_indices_and_tx2gene.sh"
        Invoke-WslBash "cd '$quotedRoot' && bash scripts/bash/03_make_metadata.sh"
        Invoke-WslBash "cd '$quotedRoot' && THREADS=$Threads bash scripts/bash/04_fastqc_salmon_resume.sh"
        Invoke-WslBash "cd '$quotedRoot' && bash scripts/bash/05_deseq2_nonadditive.sh"
        Invoke-WslBash "cd '$quotedRoot' && source scripts/bash/env_common.sh && run_env Rscript scripts/R/06_nextstep_annotation_modules.R && run_env Rscript scripts/R/07_qc_sensitivity.R && run_env Rscript scripts/R/08_make_interpretable_shortlist.R && run_env Rscript scripts/R/09_make_paper_digest.R && run_env Rscript scripts/R/10_make_final_story_tables.R && run_env Rscript scripts/R/11_make_analysis_table_workbook.R"
    }
}

