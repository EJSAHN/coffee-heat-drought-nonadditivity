[CmdletBinding()]
param(
    [string]$ProjectRoot = "D:\projects\coffee_multistress_memory"
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

$projects = @("PRJNA787748", "PRJNA1087442", "PRJNA1087679", "PRJNA1088119", "PRJNA1135679")
$fields = "run_accession,study_accession,secondary_study_accession,experiment_accession,experiment_alias,experiment_title,run_alias,sample_accession,secondary_sample_accession,sample_alias,sample_title,sample_description,scientific_name,library_layout,read_count,base_count,fastq_ftp,fastq_bytes"
$enaDir = Join-Path $ProjectRoot "data\raw\sra_metadata\ena"
New-Item -ItemType Directory -Force -Path $enaDir | Out-Null

foreach ($p in $projects) {
    $url = "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$p&result=read_run&fields=$fields&format=tsv&download=true&limit=0"
    $out = Join-Path $enaDir "$p.rich.tsv"
    curl.exe -L -o $out $url
    if ($LASTEXITCODE -ne 0) { throw "curl failed for $p" }
}

$richFiles = Get-ChildItem $enaDir -Filter "*.rich.tsv" | Sort-Object Name
$outCombined = Join-Path $ProjectRoot "data\raw\sra_metadata\coffee_all_ena_runs_RICH.tsv"
$first = $true
foreach ($f in $richFiles) {
    if ($first) {
        Get-Content $f.FullName | Set-Content $outCombined
        $first = $false
    } else {
        Get-Content $f.FullName | Select-Object -Skip 1 | Add-Content $outCombined
    }
}

Import-Csv $outCombined -Delimiter "`t" |
    Select-Object run_accession, study_accession, scientific_name, sample_accession, sample_alias, sample_title, sample_description, experiment_alias, experiment_title, run_alias |
    Export-Csv (Join-Path $ProjectRoot "data\raw\sra_metadata\metadata_triage_rich.csv") -NoTypeInformation

Write-Host "Wrote rich metadata and metadata_triage_rich.csv"
