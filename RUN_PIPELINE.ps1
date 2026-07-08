[CmdletBinding()]
param(
    [ValidateSet("setup","download-rich-metadata","index","metadata","quant","deseq2","next","qc-validation","all")]
    [string]$Step = "all",

    # Default: repository root containing this RUN_PIPELINE.ps1 file.
    [string]$ProjectRoot = "",

    [int]$Threads = 8
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

if ([string]::IsNullOrWhiteSpace($ProjectRoot)) {
    if ([string]::IsNullOrWhiteSpace($PSScriptRoot)) {
        $ProjectRoot = (Get-Location).Path
    } else {
        $ProjectRoot = $PSScriptRoot
    }
} else {
    $ProjectRoot = [System.IO.Path]::GetFullPath($ProjectRoot)
}

$ProjectRoot = [System.IO.Path]::GetFullPath($ProjectRoot)

if (-not (Test-Path -LiteralPath $ProjectRoot)) {
    throw "Project root not found: $ProjectRoot"
}

if ($Step -eq "download-rich-metadata") {
    $script = Join-Path $ProjectRoot "scripts\powershell\04_download_rich_metadata.ps1"
    & $script -ProjectRoot $ProjectRoot
    exit $LASTEXITCODE
}

$runner = Join-Path $ProjectRoot "scripts\powershell\03_run_pipeline_wsl.ps1"
& $runner -ProjectRoot $ProjectRoot -Mode $Step -Threads $Threads
exit $LASTEXITCODE
