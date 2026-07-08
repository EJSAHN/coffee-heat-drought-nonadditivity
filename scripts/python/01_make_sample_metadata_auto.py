#!/usr/bin/env python3
import argparse
import os
import re
import shutil
from pathlib import Path
import pandas as pd


def win_to_wsl_path(p: str) -> str:
    if not isinstance(p, str) or not p:
        return ''
    # Already POSIX
    if p.startswith('/'):
        return p
    # Windows drive path, for example C:\path\to\repo -> /mnt/c/path/to/repo
    m = re.match(r'^([A-Za-z]):\\(.*)$', p)
    if m:
        drive = m.group(1).lower()
        rest = m.group(2).replace('\\', '/')
        return f'/mnt/{drive}/{rest}'
    return p.replace('\\', '/')


def sanitize(x: str) -> str:
    x = re.sub(r'[^A-Za-z0-9_.-]+', '_', str(x).strip())
    x = re.sub(r'_+', '_', x).strip('_')
    return x or 'sample'


def infer_species(scientific_name: str, text: str):
    blob = f'{scientific_name} {text}'.lower()
    if 'arabica' in blob:
        return 'arabica', 'Coffea arabica'
    if 'canephora' in blob or 'conilon' in blob or 'cl153' in blob or 'clone 153' in blob:
        return 'canephora', 'Coffea canephora'
    return 'UNPARSED', scientific_name or 'UNPARSED'


def infer_genotype(index_species: str, text: str):
    blob = text.lower()
    if re.search(r'icatu', blob):
        return 'Icatu'
    if re.search(r'cl\s*153|clone\s*153|conilon', blob):
        return 'CL153'
    # Safe defaults from the study design, but keep these obvious.
    if index_species == 'arabica':
        return 'Icatu'
    if index_species == 'canephora':
        return 'CL153'
    return 'UNPARSED'


def infer_water(text: str):
    blob = text.lower().replace('-', ' ')
    if re.search(r'\bswd\b|severe\s+water\s+deficit|water\s+deficit|drought|dry', blob):
        return 'SWD'
    if re.search(r'\bww\b|well\s+watered|wellwatered|control', blob):
        return 'WW'
    return 'UNPARSED'


def infer_stage(text: str):
    blob = text.lower().replace('-', ' ')
    if re.search(r'rec\s*14|recovery|recover|14\s*d', blob):
        return 'REC14'
    if re.search(r'\b42\b|42c|42\s*°?c', blob):
        return 'T42'
    if re.search(r'\b37\b|37c|37\s*°?c', blob):
        return 'T37'
    if re.search(r'\b25\b|25c|25\s*°?c|control\s+temperature', blob):
        return 'T25'
    return 'UNPARSED'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--manifest', required=True)
    ap.add_argument('--ena', required=True)
    ap.add_argument('--out-auto', required=True)
    ap.add_argument('--out-curated', required=True)
    args = ap.parse_args()

    manifest = Path(args.manifest)
    if not manifest.exists():
        raise SystemExit(f'Missing manifest: {manifest}')

    m = pd.read_csv(manifest, sep='\t', dtype=str).fillna('')
    if m.empty:
        raise SystemExit(f'Manifest is empty: {manifest}')

    rows = []
    for run, g in m.groupby('run_accession', sort=True):
        g = g.sort_values('filename')
        first = g.iloc[0]
        title = ' '.join([str(first.get('sample_title','')), str(first.get('experiment_title',''))])
        sci = str(first.get('scientific_name',''))
        index_species, species = infer_species(sci, title)
        genotype = infer_genotype(index_species, title)
        water = infer_water(title)
        stage = infer_stage(title)
        layout = str(first.get('library_layout','')).upper() or 'UNPARSED'
        study = str(first.get('study_accession',''))
        fastqs = [win_to_wsl_path(x) for x in g['target_path'].tolist() if str(x).strip()]
        fastqs = sorted(fastqs)
        fastq_1 = ''
        fastq_2 = ''
        if layout == 'PAIRED' or len(fastqs) >= 2:
            # Prefer files with _1/_2 or .1/.2; otherwise sorted order.
            r1 = [x for x in fastqs if re.search(r'(_1|_R1|\.1)\.f(ast)?q\.gz$', x, re.I)]
            r2 = [x for x in fastqs if re.search(r'(_2|_R2|\.2)\.f(ast)?q\.gz$', x, re.I)]
            fastq_1 = r1[0] if r1 else (fastqs[0] if fastqs else '')
            fastq_2 = r2[0] if r2 else (fastqs[1] if len(fastqs) > 1 else '')
        else:
            fastq_1 = fastqs[0] if fastqs else ''

        base = '_'.join([genotype, water, stage, run])
        sample_id = sanitize(base)
        rows.append({
            'sample_id': sample_id,
            'run_accession': run,
            'study_accession': study,
            'species': species,
            'index_species': index_species,
            'genotype': genotype,
            'water': water,
            'stage': stage,
            'replicate': 'AUTO',
            'library_layout': layout,
            'fastq_1': fastq_1,
            'fastq_2': fastq_2,
            'fastq_1_exists': str(os.path.exists(fastq_1)).lower() if fastq_1 else 'false',
            'fastq_2_exists': str(os.path.exists(fastq_2)).lower() if fastq_2 else 'false',
            'salmon_index': 'data/ref/salmon_index/arabica_Cara_1_0' if index_species == 'arabica' else ('data/ref/salmon_index/canephora_AUK_PRJEB4211_v1' if index_species == 'canephora' else 'UNPARSED'),
            'include': 'yes',
            'sample_title': str(first.get('sample_title','')),
            'experiment_title': str(first.get('experiment_title','')),
        })

    df = pd.DataFrame(rows)
    # Replicate numbering within parsed condition group.
    df['replicate'] = df.groupby(['genotype','water','stage']).cumcount().add(1).astype(str)

    out_auto = Path(args.out_auto)
    out_auto.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_auto, sep='\t', index=False)
    print(f'Wrote auto metadata: {out_auto} ({len(df)} runs)')

    out_curated = Path(args.out_curated)
    if not out_curated.exists():
        shutil.copyfile(out_auto, out_curated)
        print(f'Created curated metadata template: {out_curated}')
    else:
        print(f'Curated metadata already exists, not overwriting: {out_curated}')

    summary = df.groupby(['index_species','genotype','water','stage'], dropna=False).size().reset_index(name='n')
    print('\nParsed condition summary:')
    print(summary.to_string(index=False))

    bad = df[(df[['index_species','genotype','water','stage']] == 'UNPARSED').any(axis=1)]
    if len(bad):
        print('\nWARNING: Some rows have UNPARSED fields. Edit config/sample_metadata_curated.tsv before quant/DESeq2.')
        print(bad[['run_accession','index_species','genotype','water','stage','sample_title','experiment_title']].head(20).to_string(index=False))


if __name__ == '__main__':
    main()
