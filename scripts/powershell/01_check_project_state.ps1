<#
Check project state after coffee bootstrap.
Run from anywhere:
  powershell -ExecutionPolicy Bypass -File .\scripts\powershell\01_check_project_state.ps1
#>
[CmdletBinding()]
param(
    [string]$ProjectRoot = "D:\projects\coffee_multistress_memory"
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

function Test-FileReport {
    param([string]$Path, [string]$Label)
    if (Test-Path -LiteralPath $Path) {
        $item = Get-Item -LiteralPath $Path
        Write-Host ("OK   {0,-38} {1} ({2:N2} MB)" -f $Label, $Path, ($item.Length/1MB)) -ForegroundColor Green
        return $true
    } else {
        Write-Host ("MISS {0,-38} {1}" -f $Label, $Path) -ForegroundColor Red
        return $false
    }
}

Write-Host "\n=== Coffee multistress project check ===" -ForegroundColor Cyan
Write-Host "ProjectRoot: $ProjectRoot"

if (-not (Test-Path -LiteralPath $ProjectRoot)) {
    throw "Project root not found: $ProjectRoot"
}

$manifest = Join-Path $ProjectRoot "data\raw\sra_metadata\fastq_download_manifest.tsv"
$ena = Join-Path $ProjectRoot "data\raw\sra_metadata\coffee_all_ena_runs.tsv"
$template = Join-Path $ProjectRoot "config\sample_metadata_template.tsv"
$arabicaRna = Join-Path $ProjectRoot "data\ref\arabica_Cara_1_0_NCBI\GCF_003713225.1_Cara_1.0_rna.fna.gz"
$canephoraCdna = Join-Path $ProjectRoot "data\ref\canephora_AUK_PRJEB4211_v1_EnsemblPlants\Coffea_canephora.AUK_PRJEB4211_v1.cdna.all.fa.gz"
$bootstrap = Join-Path $ProjectRoot "scripts\powershell\00_bootstrap_coffee_multistress.ps1"

$checks = @()
$checks += Test-FileReport $ena "combined ENA metadata"
$checks += Test-FileReport $manifest "FASTQ manifest"
$checks += Test-FileReport $template "sample metadata template"
$checks += Test-FileReport $arabicaRna "arabica transcript FASTA"
$checks += Test-FileReport $canephoraCdna "canephora cDNA FASTA"
$checks += Test-FileReport $bootstrap "bootstrap script"

if (Test-Path -LiteralPath $ena) {
    $enaRows = Import-Csv -Path $ena -Delimiter "`t"
    Write-Host "\nENA runs by study/species/layout:" -ForegroundColor Cyan
    $enaRows | Group-Object study_accession, scientific_name, library_layout | ForEach-Object {
        Write-Host ("  {0,-60} {1,4} runs" -f $_.Name, $_.Count)
    }
}

if (Test-Path -LiteralPath $manifest) {
    $m = Import-Csv -Path $manifest -Delimiter "`t"
    $runs = @($m | Select-Object -ExpandProperty run_accession -Unique)
    $files = @($m)
    $totalBytes = 0.0
    foreach ($row in $m) {
        try { if ($row.fastq_bytes) { $totalBytes += [double]$row.fastq_bytes } } catch { }
    }
    Write-Host "\nFASTQ manifest summary:" -ForegroundColor Cyan
    Write-Host ("  unique runs: {0}" -f $runs.Count)
    Write-Host ("  FASTQ files:  {0}" -f $files.Count)
    if ($totalBytes -gt 0) { Write-Host ("  estimated compressed FASTQ size: {0:N2} GB" -f ($totalBytes/1GB)) }
}

$fastqs = Get-ChildItem -Path (Join-Path $ProjectRoot "data\raw\fastq") -Recurse -File -Include *.fastq.gz,*.fq.gz -ErrorAction SilentlyContinue
Write-Host "\nDownloaded FASTQ files currently present: $(@($fastqs).Count)" -ForegroundColor Cyan
if (@($fastqs).Count -eq 0) {
    Write-Host "  No FASTQ yet. That is expected from the bootstrap log; use 02_download_fastq.ps1." -ForegroundColor Yellow
}

Write-Host "\nWSL status:" -ForegroundColor Cyan
try {
    & wsl --list --verbose
} catch {
    Write-Host "  WSL command not available or WSL not installed." -ForegroundColor Yellow
}

if ($checks -contains $false) {
    Write-Host "\nSome bootstrap pieces are missing. Re-run the bootstrap before full pipeline." -ForegroundColor Yellow
} else {
    Write-Host "\nBootstrap state looks good. Next: FASTQ test download, then WSL pipeline." -ForegroundColor Green
}
