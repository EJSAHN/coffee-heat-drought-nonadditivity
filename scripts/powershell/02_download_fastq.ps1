<#
Download FASTQ files using the bootstrap manifest logic.
Examples:
  .\scripts\powershell\02_download_fastq.ps1 -Mode test -MaxRuns 2
  .\scripts\powershell\02_download_fastq.ps1 -Mode all
#>
[CmdletBinding()]
param(
    [string]$ProjectRoot = "D:\projects\coffee_multistress_memory",
    [ValidateSet("test","all")]
    [string]$Mode = "test",
    [int]$MaxRuns = 2
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

$bootstrap = Join-Path $ProjectRoot "scripts\powershell\00_bootstrap_coffee_multistress.ps1"
if (-not (Test-Path -LiteralPath $bootstrap)) {
    throw "Bootstrap script not found: $bootstrap"
}

Write-Host "ProjectRoot: $ProjectRoot" -ForegroundColor Cyan
Write-Host "Bootstrap:   $bootstrap" -ForegroundColor Cyan

if ($Mode -eq "test") {
    Write-Host "Running FASTQ test download for first $MaxRuns run(s)." -ForegroundColor Yellow
    & $bootstrap -ProjectRoot $ProjectRoot -DownloadFastq -MaxRuns $MaxRuns -SkipReference -SkipSupplements
} else {
    Write-Host "Running FULL FASTQ download. This can be large." -ForegroundColor Yellow
    & $bootstrap -ProjectRoot $ProjectRoot -DownloadFastq -SkipReference -SkipSupplements
}
