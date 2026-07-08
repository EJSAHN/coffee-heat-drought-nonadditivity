[CmdletBinding()]
param(
    [string]$ProjectRoot = "",
    [int]$Threads = 8,
    [switch]$SkipStar
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

try {
    $wslRoot = (& wsl.exe --exec wslpath -a $ProjectRoot).Trim()
} catch {
    throw "WSL is not available. Install Ubuntu with: wsl --install -d Ubuntu"
}

$quotedRoot = $wslRoot.Replace("'", "'\''")
Write-Host "== CL153 QC validation analyses" -ForegroundColor Green
Write-Host "Windows project root: $ProjectRoot" -ForegroundColor Green
Write-Host "WSL project root:     $wslRoot" -ForegroundColor Green
Write-Host "Threads:              $Threads" -ForegroundColor Green

function Invoke-WslBash {
    param([string]$Command)
    Write-Host "`n[WSL] $Command" -ForegroundColor Cyan
    & wsl.exe --exec bash -lc $Command
    if ($LASTEXITCODE -ne 0) { throw "WSL command failed with exit code $LASTEXITCODE" }
}

Invoke-WslBash "cd '$quotedRoot' && source scripts/bash/env_common.sh && run_env Rscript scripts/R/12_cl153_low_mapping_sensitivity.R"
Invoke-WslBash "cd '$quotedRoot' && source scripts/bash/env_common.sh && run_env Rscript scripts/R/13_cl153_stage_local_sensitivity.R"

if (-not $SkipStar) {
    Invoke-WslBash "cd '$quotedRoot' && THREADS=$Threads bash scripts/bash/14_cl153_star_genome_alignment_validation.sh"
} else {
    Write-Host "Skipping STAR genome-alignment validation because -SkipStar was supplied." -ForegroundColor Yellow
}

Write-Host "DONE. Open: $ProjectRoot\results\qc_validation" -ForegroundColor Green
