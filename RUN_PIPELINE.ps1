[CmdletBinding()]
param(
    [ValidateSet("setup","download-rich-metadata","index","metadata","quant","deseq2","next","all")]
    [string]$Step = "all",
    [string]$ProjectRoot = "D:\projects\coffee_multistress_memory",
    [int]$Threads = 8
)

$ErrorActionPreference = "Stop"

if ($Step -eq "download-rich-metadata") {
    & .\scripts\powershell\04_download_rich_metadata.ps1 -ProjectRoot $ProjectRoot
    exit $LASTEXITCODE
}

& .\scripts\powershell\03_run_pipeline_wsl.ps1 -ProjectRoot $ProjectRoot -Mode $Step -Threads $Threads
