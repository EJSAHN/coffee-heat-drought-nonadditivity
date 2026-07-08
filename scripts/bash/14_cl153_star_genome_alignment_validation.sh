#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(pwd)"
THREADS="${THREADS:-8}"
RELEASE="${ENSEMBL_PLANTS_RELEASE:-62}"

source scripts/bash/env_common.sh

echo "== CL153 STAR genome-alignment validation"
echo "== Project: ${PROJECT_ROOT}"
echo "== Threads: ${THREADS}"
echo "== Ensembl Plants release requested: ${RELEASE}"

OUT_DIR="results/qc_validation/cl153_star_genome_alignment"
REF_DIR="data/ref/canephora_AUK_PRJEB4211_v1_STAR_validation"
mkdir -p "${OUT_DIR}" "${REF_DIR}"

require_tool() {
    if ! run_env "$1" --version >/dev/null 2>&1; then
        echo "ERROR: $1 (STAR/gffread/samtools) not found in the coffee-rnaseq environment."
        echo "Install the environment first, then rerun the QC validation step:"
        echo "  .\\RUN_PIPELINE.ps1 -Step setup -Threads 8"
        echo "  .\\RUN_PIPELINE.ps1 -Step qc-validation -Threads 8"
        exit 1
    fi
}
require_tool STAR
require_tool gffread
require_tool samtools

download_if_missing() {
    local url="$1"
    local out="$2"
    if [[ -s "$out" ]]; then
        echo "Already present: $out"
    else
        echo "Downloading: $url"
        curl -L --retry 5 --retry-delay 10 -o "$out" "$url"
        if [[ ! -s "$out" ]]; then
            echo "ERROR: failed to download $url"
            exit 1
        fi
    fi
}

# This analysis uses STAR genome alignment with Ensembl Plants release 62 AUK_PRJEB4211_v1.
# It provides a genome-count comparison to the primary transcript-level Salmon workflow.
DNA_GZ="${REF_DIR}/Coffea_canephora.AUK_PRJEB4211_v1.dna.toplevel.fa.gz"
GFF_GZ="${REF_DIR}/Coffea_canephora.AUK_PRJEB4211_v1.${RELEASE}.gff3.gz"

DNA_URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-${RELEASE}/fasta/coffea_canephora/dna/Coffea_canephora.AUK_PRJEB4211_v1.dna.toplevel.fa.gz"
GFF_URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-${RELEASE}/gff3/coffea_canephora/Coffea_canephora.AUK_PRJEB4211_v1.${RELEASE}.gff3.gz"

download_if_missing "$DNA_URL" "$DNA_GZ"
download_if_missing "$GFF_URL" "$GFF_GZ"

DNA_FA="${REF_DIR}/Coffea_canephora.AUK_PRJEB4211_v1.dna.toplevel.fa"
GFF="${REF_DIR}/Coffea_canephora.AUK_PRJEB4211_v1.${RELEASE}.gff3"
GTF="${REF_DIR}/Coffea_canephora.AUK_PRJEB4211_v1.${RELEASE}.gtf"

if [[ ! -s "$DNA_FA" ]]; then
    echo "Decompressing genome FASTA"
    gunzip -c "$DNA_GZ" > "$DNA_FA"
fi
if [[ ! -s "$GFF" ]]; then
    echo "Decompressing GFF3"
    gunzip -c "$GFF_GZ" > "$GFF"
fi
if [[ ! -s "$GTF" ]]; then
    echo "Converting GFF3 to GTF with gffread"
    run_env gffread "$GFF" -T -o "$GTF"
fi

MANIFEST="${OUT_DIR}/cl153_star_manifest.tsv"
run_env python - <<'PY'
from pathlib import Path
import csv, sys, re

root = Path(".")
meta_path = root / "config" / "sample_metadata_curated.tsv"
if not meta_path.exists():
    raise SystemExit(f"Missing {meta_path}")

def norm(s):
    return (s or "").strip()

with open(meta_path, newline='', encoding='utf-8-sig') as f:
    rows = list(csv.DictReader(f, delimiter='\t'))

if not rows:
    raise SystemExit("sample metadata is empty")

cols = set(rows[0].keys())

def get(row, candidates):
    lower_map = {c.lower(): c for c in row.keys()}
    for cand in candidates:
        if cand in row:
            return row[cand]
        if cand.lower() in lower_map:
            return row[lower_map[cand.lower()]]
    return ""

sample_candidates = ["sample_id", "sample", "Sample", "sample_name"]
run_candidates = ["run_accession", "Run", "run", "run_id", "Run accession", "run_accession"]
include_candidates = ["include", "Include"]
species_candidates = ["index_species", "species", "Species"]
genotype_candidates = ["genotype", "Genotype"]
water_candidates = ["water", "Water"]
stage_candidates = ["stage", "Stage"]
rep_candidates = ["replicate", "rep", "Replicate"]

out = []
for row in rows:
    include = norm(get(row, include_candidates)).lower()
    if include and include not in {"yes", "y", "true", "1", "include"}:
        continue
    species = norm(get(row, species_candidates)).lower()
    genotype = norm(get(row, genotype_candidates))
    if not (species == "canephora" or genotype.upper() == "CL153"):
        continue
    sample_id = norm(get(row, sample_candidates))
    run = norm(get(row, run_candidates))
    water = norm(get(row, water_candidates))
    stage = norm(get(row, stage_candidates))
    rep = norm(get(row, rep_candidates))
    if not sample_id:
        sample_id = "_".join(x for x in [genotype, water, stage, f"R{rep}" if rep else ""] if x)
    if not run:
        # Try deriving run from prior salmon quant directory names not possible; fail explicitly.
        raise SystemExit(f"Could not find run accession for sample {sample_id}. Metadata columns: {list(row.keys())}")
    matches1 = sorted(root.glob(f"data/raw/fastq/**/{run}_1.fastq.gz"))
    matches2 = sorted(root.glob(f"data/raw/fastq/**/{run}_2.fastq.gz"))
    if not matches1 or not matches2:
        raise SystemExit(f"Missing FASTQ pair for run {run}, sample {sample_id}")
    out.append({
        "sample_id": sample_id,
        "run_accession": run,
        "water": water,
        "stage": stage,
        "replicate": rep,
        "fastq_1": str(matches1[0]),
        "fastq_2": str(matches2[0]),
    })

if len(out) != 24:
    print(f"WARNING: expected 24 CL153 rows, found {len(out)}", file=sys.stderr)

dest = root / "results/qc_validation/cl153_star_genome_alignment/cl153_star_manifest.tsv"
dest.parent.mkdir(parents=True, exist_ok=True)
with open(dest, "w", newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, delimiter='\t', fieldnames=["sample_id","run_accession","water","stage","replicate","fastq_1","fastq_2"])
    w.writeheader()
    w.writerows(out)
print(dest)
print("n =", len(out))
PY

STAR_INDEX="${REF_DIR}/STAR_index_sjdbOverhang100_sa13"
mkdir -p "$STAR_INDEX"

if [[ ! -s "${STAR_INDEX}/Genome" ]]; then
    echo "== Building STAR index"
    run_env STAR \
        --runThreadN "$THREADS" \
        --runMode genomeGenerate \
        --genomeDir "$STAR_INDEX" \
        --genomeFastaFiles "$DNA_FA" \
        --sjdbGTFfile "$GTF" \
        --sjdbOverhang 100 \
        --genomeSAindexNbases 13 \
        --limitGenomeGenerateRAM 32000000000
else
    echo "STAR index already present: ${STAR_INDEX}"
fi

STAR_OUT="${OUT_DIR}/star_alignments"
mkdir -p "$STAR_OUT"

# STAR uses FIFO files during mapping. Windows-mounted drives (/mnt/c, /mnt/d)
# do not support FIFO files, so keep STAR temporary directories on the native
# Linux filesystem and write only final outputs to the project directory.
STAR_TMP_BASE="${HOME}/star_tmp/coffee_cl153_star_validation"
mkdir -p "$STAR_TMP_BASE"

echo "== STAR-aligning CL153 libraries"
tail -n +2 "$MANIFEST" | while IFS=$'\t' read -r sample run water stage rep fq1 fq2; do
    sample_dir="${STAR_OUT}/${sample}"
    prefix="${sample_dir}/"
    mkdir -p "$sample_dir"
    if [[ -s "${prefix}ReadsPerGene.out.tab" && -s "${prefix}Log.final.out" ]]; then
        echo "Skip existing STAR output: $sample"
        continue
    fi
    echo "STAR run: $sample ($run)"
    tmpdir="${STAR_TMP_BASE}/${sample}_STARtmp"
    rm -rf "$tmpdir"
    run_env STAR \
        --runThreadN "$THREADS" \
        --genomeDir "$STAR_INDEX" \
        --readFilesIn "$fq1" "$fq2" \
        --readFilesCommand zcat \
        --outFileNamePrefix "$prefix" \
        --outTmpDir "$tmpdir" \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts \
        --outSAMattrRGline ID:"$sample" SM:"$sample" \
        --limitBAMsortRAM 12000000000
    if [[ -s "${prefix}Aligned.sortedByCoord.out.bam" ]]; then
        run_env samtools index -@ "$THREADS" "${prefix}Aligned.sortedByCoord.out.bam" || true
    fi
    rm -rf "$tmpdir" || true
done

echo "== Running R summary and DESeq2 sensitivity"
run_env Rscript scripts/R/14_cl153_star_genome_alignment_validation.R

echo "== DONE CL153 STAR genome-alignment validation"
echo "Outputs written to: ${OUT_DIR}"
