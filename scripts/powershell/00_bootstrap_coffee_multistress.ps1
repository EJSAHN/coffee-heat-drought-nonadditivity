<#
Bootstrap script for the coffee heat-drought multistress project.
Creates a reproducible project folder, downloads source metadata, paper supplementary files,
reference files, and optionally raw FASTQ files from ENA.

Default project root: repository root containing this script. Override with -ProjectRoot if needed.

Usage examples:
  Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
  .\00_bootstrap_coffee_multistress.ps1
  .\00_bootstrap_coffee_multistress.ps1 -DownloadFastq -MaxRuns 2
  .\00_bootstrap_coffee_multistress.ps1 -DownloadFastq
#>

param(
    [string]$ProjectRoot = "",
    [switch]$DownloadFastq,
    [int]$MaxRuns = 0,
    [switch]$SkipReference,
    [switch]$SkipSupplements
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

if ([string]::IsNullOrWhiteSpace($ProjectRoot)) {
    $ProjectRoot = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot "..\.."))
} else {
    $ProjectRoot = [System.IO.Path]::GetFullPath($ProjectRoot)
}

New-Item -ItemType Directory -Force -Path $ProjectRoot | Out-Null

function Ensure-Dir {
    param([string]$Path)
    if (-not (Test-Path -LiteralPath $Path)) {
        New-Item -ItemType Directory -Path $Path | Out-Null
    }
}

function Write-Log {
    param([string]$Message)
    $timestamp = Get-Date -Format "yyyy-MM-dd HH:mm:ss"
    $line = "[$timestamp] $Message"
    Write-Host $line
    if ($script:LogFile) { Add-Content -Path $script:LogFile -Value $line }
}

function Download-File {
    param(
        [Parameter(Mandatory=$true)][string]$Url,
        [Parameter(Mandatory=$true)][string]$OutFile
    )
    Ensure-Dir (Split-Path -Parent $OutFile)
    if (Test-Path -LiteralPath $OutFile) {
        $size = (Get-Item -LiteralPath $OutFile).Length
        if ($size -gt 0) {
            Write-Log "exists: $OutFile"
            return
        }
    }
    Write-Log "download: $Url"
    $curl = Get-Command curl.exe -ErrorAction SilentlyContinue
    if ($curl) {
        & curl.exe -L --retry 5 --retry-delay 5 --continue-at - --output "$OutFile" "$Url"
        if ($LASTEXITCODE -ne 0) { throw "curl.exe failed for $Url" }
    } else {
        Invoke-WebRequest -Uri $Url -OutFile $OutFile -UseBasicParsing
    }
}

function Write-TextFile {
    param([string]$Path, [string]$Text)
    Ensure-Dir (Split-Path -Parent $Path)
    Set-Content -Path $Path -Value $Text -Encoding UTF8
}

# -----------------------------------------------------------------------------
# 0. Project folder layout
# -----------------------------------------------------------------------------
$dirs = @(
    "data", "data\accessions", "data\raw", "data\raw\sra_metadata", "data\raw\sra_metadata\ncbi_runinfo",
    "data\raw\sra_metadata\ena", "data\raw\fastq", "data\raw\supplementary",
    "data\ref", "data\ref\arabica_Cara_1_0_NCBI", "data\ref\canephora_AUK_PRJEB4211_v1_EnsemblPlants",
    "data\processed", "data\processed\counts", "data\processed\tximport", "data\processed\annotations",
    "results", "results\qc", "results\multiqc", "results\salmon", "results\deseq2", "results\tables", "results\reports",
    "scripts", "scripts\powershell", "scripts\R", "scripts\python",
    "notebooks", "config", "env", "docs", "docs\literature", "docs\protocols", "manuscript", "logs"
)

foreach ($d in $dirs) { Ensure-Dir (Join-Path $ProjectRoot $d) }
$script:LogFile = Join-Path $ProjectRoot "logs\bootstrap_$(Get-Date -Format 'yyyyMMdd_HHmmss').log"
Write-Log "Project root: $ProjectRoot"

# Copy this bootstrap script into the project for provenance, when running from a file.
try {
    if ($PSCommandPath) {
        Copy-Item -Path $PSCommandPath -Destination (Join-Path $ProjectRoot "scripts\powershell\00_bootstrap_coffee_multistress.ps1") -Force
    }
} catch { }

# -----------------------------------------------------------------------------
# 1. Accessions and project notes
# -----------------------------------------------------------------------------
$bioProjects = @("PRJNA787748", "PRJNA1087442", "PRJNA1087679", "PRJNA1088119", "PRJNA1135679")
$accessionTsv = @"
study_accession	note
PRJNA787748	NCBI SRA BioProject listed in Marques et al. 2024 IJMS data availability for Coffea heat/drought transcriptomes
PRJNA1087442	NCBI SRA BioProject listed in Marques et al. 2024 IJMS data availability for Coffea heat/drought transcriptomes
PRJNA1087679	NCBI SRA BioProject listed in Marques et al. 2024 IJMS data availability for Coffea heat/drought transcriptomes
PRJNA1088119	NCBI SRA BioProject listed in Marques et al. 2024 IJMS data availability for Coffea heat/drought transcriptomes
PRJNA1135679	NCBI SRA BioProject listed in Marques et al. 2024 IJMS data availability for Coffea heat/drought transcriptomes
"@
Write-TextFile -Path (Join-Path $ProjectRoot "data\accessions\bioprojects_coffee_heat_drought.tsv") -Text $accessionTsv

$readme = @"
# Coffee heat-drought non-additivity reanalysis

Working title:
Non-additive heat-drought transcriptional logic and recovery hysteresis in Coffea arabica and Coffea canephora

Base folder:
$ProjectRoot

Main source paper:
Marques et al. 2024, International Journal of Molecular Sciences, DOI: 10.3390/ijms25147995

Primary raw-read accessions listed by the paper:
PRJNA787748, PRJNA1087442, PRJNA1087679, PRJNA1088119, PRJNA1135679

Folder logic:
- data/raw/sra_metadata: NCBI/ENA metadata and FASTQ URL manifests
- data/raw/fastq: optional raw FASTQ download target
- data/raw/supplementary: article supplementary material
- data/ref: reference genome / transcriptome / annotation files
- data/processed: counts, tximport, parsed annotations
- results: tables, QC, DESeq2 outputs
- scripts: PowerShell, R, Python scripts
- config: sample metadata and design files
- manuscript: manuscript drafts

First pass:
1. Run bootstrap without -DownloadFastq to get metadata, supplementary files, and reference files.
2. Open data/raw/sra_metadata/coffee_all_ena_runs.tsv and config/sample_metadata_template.tsv.
3. Confirm which runs belong to the exact 2 species/genotypes x water x temperature/recovery design.
4. Then run with -DownloadFastq, optionally -MaxRuns for a small test.

Recommended test:
.\scripts\powershell\00_bootstrap_coffee_multistress.ps1 -DownloadFastq -MaxRuns 2

Full raw download:
.\scripts\powershell\00_bootstrap_coffee_multistress.ps1 -DownloadFastq
"@
Write-TextFile -Path (Join-Path $ProjectRoot "README_PROJECT.md") -Text $readme

# -----------------------------------------------------------------------------
# 2. Download paper supplementary material
# -----------------------------------------------------------------------------
if (-not $SkipSupplements) {
    $suppUrl = "https://www.mdpi.com/article/10.3390/ijms25147995/s1"
    $suppOut = Join-Path $ProjectRoot "data\raw\supplementary\Marques_2024_IJMS_25_7995_supplementary_material.zip"
    Download-File -Url $suppUrl -OutFile $suppOut
}

# -----------------------------------------------------------------------------
# 3. Download NCBI RunInfo and ENA metadata for each BioProject
# -----------------------------------------------------------------------------
$enaFields = "run_accession,study_accession,secondary_study_accession,experiment_accession,sample_accession,secondary_sample_accession,scientific_name,library_strategy,library_source,library_layout,instrument_platform,instrument_model,read_count,base_count,fastq_ftp,fastq_bytes,submitted_ftp,sra_ftp,sample_title,experiment_title"
$allEnaRows = @()

foreach ($acc in $bioProjects) {
    $ncbiUrl = "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=$acc"
    $ncbiOut = Join-Path $ProjectRoot "data\raw\sra_metadata\ncbi_runinfo\$acc.runinfo.csv"
    Download-File -Url $ncbiUrl -OutFile $ncbiOut

    $enaUrl = "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$acc&result=read_run&fields=$enaFields&format=tsv&download=true&limit=0"
    $enaOut = Join-Path $ProjectRoot "data\raw\sra_metadata\ena\$acc.ena_read_run.tsv"
    Download-File -Url $enaUrl -OutFile $enaOut

    try {
        $rows = Import-Csv -Path $enaOut -Delimiter "`t"
        foreach ($r in $rows) { $allEnaRows += $r }
    } catch {
        Write-Log "WARNING: could not parse ENA metadata for ${acc}: $($_.Exception.Message)"
    }
}

$allEnaPath = Join-Path $ProjectRoot "data\raw\sra_metadata\coffee_all_ena_runs.tsv"
if ($allEnaRows.Count -gt 0) {
    $allEnaRows | Sort-Object study_accession, run_accession | Export-Csv -Path $allEnaPath -Delimiter "`t" -NoTypeInformation
    Write-Log "wrote combined ENA run metadata: $allEnaPath"
} else {
    Write-Log "WARNING: no ENA rows parsed. Check individual metadata files in data/raw/sra_metadata/ena."
}

# -----------------------------------------------------------------------------
# 4. Build FASTQ download manifest from ENA fastq_ftp field
# -----------------------------------------------------------------------------
$manifestRows = New-Object System.Collections.Generic.List[object]
foreach ($r in $allEnaRows) {
    if (-not $r.fastq_ftp) { continue }
    $urls = $r.fastq_ftp -split ";"
    $bytes = @()
    if ($r.fastq_bytes) { $bytes = $r.fastq_bytes -split ";" }
    for ($i = 0; $i -lt $urls.Count; $i++) {
        $u = $urls[$i].Trim()
        if ([string]::IsNullOrWhiteSpace($u)) { continue }
        if ($u -notmatch "^https?://" -and $u -notmatch "^ftp://") { $u = "https://$u" }
        if ($u -match "^ftp://") { $u = $u -replace "^ftp://", "https://" }
        $fname = Split-Path $u -Leaf
        $study = if ($r.study_accession) { $r.study_accession } else { "unknown_study" }
        $target = Join-Path $ProjectRoot ("data\raw\fastq\$study\$fname")
        $b = if ($i -lt $bytes.Count) { $bytes[$i] } else { "" }
        $manifestRows.Add([pscustomobject]@{
            run_accession = $r.run_accession
            study_accession = $r.study_accession
            secondary_study_accession = $r.secondary_study_accession
            scientific_name = $r.scientific_name
            library_layout = $r.library_layout
            source_url = $u
            filename = $fname
            target_path = $target
            fastq_bytes = $b
            sample_title = $r.sample_title
            experiment_title = $r.experiment_title
        }) | Out-Null
    }
}

$manifestPath = Join-Path $ProjectRoot "data\raw\sra_metadata\fastq_download_manifest.tsv"
$manifestRows | Export-Csv -Path $manifestPath -Delimiter "`t" -NoTypeInformation
Write-Log "wrote FASTQ manifest: $manifestPath"

# Sample metadata template; it is intentionally not auto-filled beyond run titles.
$templatePath = Join-Path $ProjectRoot "config\sample_metadata_template.tsv"
$sampleRows = @()
foreach ($r in ($allEnaRows | Sort-Object study_accession, run_accession)) {
    $sampleRows += [pscustomobject]@{
        run_accession = $r.run_accession
        species = $r.scientific_name
        genotype = "TODO_Icatu_or_CL153"
        water = "TODO_WW_or_SWD"
        temperature = "TODO_25_37_42_or_REC14"
        timepoint = "TODO"
        replicate = "TODO"
        include = "yes"
        study_accession = $r.study_accession
        sample_title = $r.sample_title
        experiment_title = $r.experiment_title
    }
}
$sampleRows | Export-Csv -Path $templatePath -Delimiter "`t" -NoTypeInformation
Write-Log "wrote sample metadata template: $templatePath"

# -----------------------------------------------------------------------------
# 5. Download reference files
# -----------------------------------------------------------------------------
if (-not $SkipReference) {
    $referenceRows = @(
        [pscustomobject]@{species="Coffea arabica"; source="NCBI RefSeq"; file="genome_fasta"; url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/713/225/GCF_003713225.1_Cara_1.0/GCF_003713225.1_Cara_1.0_genomic.fna.gz"; dest="data\ref\arabica_Cara_1_0_NCBI\GCF_003713225.1_Cara_1.0_genomic.fna.gz"},
        [pscustomobject]@{species="Coffea arabica"; source="NCBI RefSeq"; file="gff"; url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/713/225/GCF_003713225.1_Cara_1.0/GCF_003713225.1_Cara_1.0_genomic.gff.gz"; dest="data\ref\arabica_Cara_1_0_NCBI\GCF_003713225.1_Cara_1.0_genomic.gff.gz"},
        [pscustomobject]@{species="Coffea arabica"; source="NCBI RefSeq"; file="cds"; url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/713/225/GCF_003713225.1_Cara_1.0/GCF_003713225.1_Cara_1.0_cds_from_genomic.fna.gz"; dest="data\ref\arabica_Cara_1_0_NCBI\GCF_003713225.1_Cara_1.0_cds_from_genomic.fna.gz"},
        [pscustomobject]@{species="Coffea arabica"; source="NCBI RefSeq"; file="transcripts"; url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/713/225/GCF_003713225.1_Cara_1.0/GCF_003713225.1_Cara_1.0_rna.fna.gz"; dest="data\ref\arabica_Cara_1_0_NCBI\GCF_003713225.1_Cara_1.0_rna.fna.gz"},
        [pscustomobject]@{species="Coffea arabica"; source="NCBI RefSeq"; file="protein"; url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/713/225/GCF_003713225.1_Cara_1.0/GCF_003713225.1_Cara_1.0_protein.faa.gz"; dest="data\ref\arabica_Cara_1_0_NCBI\GCF_003713225.1_Cara_1.0_protein.faa.gz"},
        [pscustomobject]@{species="Coffea canephora"; source="Ensembl Plants release 62"; file="genome_fasta"; url="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/fasta/coffea_canephora/dna/Coffea_canephora.AUK_PRJEB4211_v1.dna.toplevel.fa.gz"; dest="data\ref\canephora_AUK_PRJEB4211_v1_EnsemblPlants\Coffea_canephora.AUK_PRJEB4211_v1.dna.toplevel.fa.gz"},
        [pscustomobject]@{species="Coffea canephora"; source="Ensembl Plants release 62"; file="gff"; url="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/gff3/coffea_canephora/Coffea_canephora.AUK_PRJEB4211_v1.62.gff3.gz"; dest="data\ref\canephora_AUK_PRJEB4211_v1_EnsemblPlants\Coffea_canephora.AUK_PRJEB4211_v1.62.gff3.gz"},
        [pscustomobject]@{species="Coffea canephora"; source="Ensembl Plants release 62"; file="cdna"; url="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/fasta/coffea_canephora/cdna/Coffea_canephora.AUK_PRJEB4211_v1.cdna.all.fa.gz"; dest="data\ref\canephora_AUK_PRJEB4211_v1_EnsemblPlants\Coffea_canephora.AUK_PRJEB4211_v1.cdna.all.fa.gz"},
        [pscustomobject]@{species="Coffea canephora"; source="Ensembl Plants release 62"; file="cds"; url="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/fasta/coffea_canephora/cds/Coffea_canephora.AUK_PRJEB4211_v1.cds.all.fa.gz"; dest="data\ref\canephora_AUK_PRJEB4211_v1_EnsemblPlants\Coffea_canephora.AUK_PRJEB4211_v1.cds.all.fa.gz"},
        [pscustomobject]@{species="Coffea canephora"; source="Ensembl Plants release 62"; file="protein"; url="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/fasta/coffea_canephora/pep/Coffea_canephora.AUK_PRJEB4211_v1.pep.all.fa.gz"; dest="data\ref\canephora_AUK_PRJEB4211_v1_EnsemblPlants\Coffea_canephora.AUK_PRJEB4211_v1.pep.all.fa.gz"}
    )
    $refManifestPath = Join-Path $ProjectRoot "data\ref\reference_sources.tsv"
    $referenceRows | Export-Csv -Path $refManifestPath -Delimiter "`t" -NoTypeInformation
    foreach ($rr in $referenceRows) {
        Download-File -Url $rr.url -OutFile (Join-Path $ProjectRoot $rr.dest)
    }
}

# -----------------------------------------------------------------------------
# 6. Optional raw FASTQ download
# -----------------------------------------------------------------------------
if ($DownloadFastq) {
    $rowsToDownload = $manifestRows
    if ($MaxRuns -gt 0) {
        $selectedRuns = $manifestRows | Select-Object -ExpandProperty run_accession -Unique | Select-Object -First $MaxRuns
        $runSet = @{}
        foreach ($x in $selectedRuns) { $runSet[$x] = $true }
        $rowsToDownload = $manifestRows | Where-Object { $runSet.ContainsKey($_.run_accession) }
        Write-Log "FASTQ test mode: downloading files from first $MaxRuns unique run(s)."
    } else {
        Write-Log "FASTQ full mode: downloading all files in manifest. This may require substantial disk space."
    }

    foreach ($m in $rowsToDownload) {
        Download-File -Url $m.source_url -OutFile $m.target_path
    }
} else {
    Write-Log "FASTQ download skipped. Re-run with -DownloadFastq when ready."
}

# -----------------------------------------------------------------------------
# 7. Minimal R scripts and env notes
# -----------------------------------------------------------------------------
$rSummary = @'
# Quick metadata summary after bootstrap
library(readr)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) normalizePath(args[1], winslash = "/", mustWork = TRUE) else normalizePath(getwd(), winslash = "/", mustWork = TRUE)
meta_file <- file.path(root, "data/raw/sra_metadata/coffee_all_ena_runs.tsv")
meta <- read_tsv(meta_file, show_col_types = FALSE)

message("Runs by study_accession:")
print(meta %>% count(study_accession, scientific_name, library_layout, sort = TRUE))

message("Sample/experiment title preview:")
print(meta %>% select(run_accession, study_accession, scientific_name, sample_title, experiment_title) %>% head(30))
'@
Write-TextFile -Path (Join-Path $ProjectRoot "scripts\R\00_summarize_ena_metadata.R") -Text $rSummary

$envNotes = @'
# Suggested local software environment

For metadata-only setup, PowerShell is enough.

For a lightweight local RNA-seq route later, use Salmon + tximport/DESeq2 rather than genome-wide STAR alignment:

mamba create -n coffee-multistress -c conda-forge -c bioconda salmon fastqc multiqc seqkit sra-tools
mamba activate coffee-multistress

R packages later:
BiocManager::install(c("tximport", "DESeq2", "edgeR", "limma", "apeglm", "clusterProfiler"))
install.packages(c("tidyverse", "openxlsx"))
'@
Write-TextFile -Path (Join-Path $ProjectRoot "env\README_conda_and_R.md") -Text $envNotes

Write-Log "DONE. Open: $ProjectRoot\README_PROJECT.md"

