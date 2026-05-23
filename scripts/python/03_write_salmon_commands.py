#!/usr/bin/env python3
import argparse
import os
from pathlib import Path
import pandas as pd


def shell_quote(s: str) -> str:
    return "'" + str(s).replace("'", "'\\''") + "'"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--metadata', required=True)
    ap.add_argument('--threads', type=int, default=8)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.metadata, sep='\t', dtype=str).fillna('')
    if 'include' in df.columns:
        df = df[df['include'].str.lower().isin(['yes','y','true','1'])]

    cmds = []
    skipped = []
    for _, r in df.iterrows():
        sid = r.get('sample_id') or r.get('run_accession')
        sid = ''.join(c if c.isalnum() or c in '._-' else '_' for c in sid)
        layout = r.get('library_layout','').upper()
        idx = r.get('salmon_index','')
        fq1 = r.get('fastq_1','')
        fq2 = r.get('fastq_2','')
        outdir = f'results/salmon/{sid}'
        if not idx or idx == 'UNPARSED' or not os.path.isdir(idx):
            skipped.append((sid, 'missing salmon index', idx))
            continue
        if not fq1 or not os.path.exists(fq1):
            skipped.append((sid, 'missing fastq_1', fq1))
            continue
        if (layout == 'PAIRED' or fq2) and (not fq2 or not os.path.exists(fq2)):
            skipped.append((sid, 'missing fastq_2', fq2))
            continue
        cmds.append(f"if [[ -s {shell_quote(outdir + '/quant.sf')} ]]; then echo 'SKIP existing {sid}'; else mkdir -p {shell_quote(outdir)}; salmon quant -i {shell_quote(idx)} -l A " +
                    (f"-1 {shell_quote(fq1)} -2 {shell_quote(fq2)} " if (layout == 'PAIRED' or fq2) else f"-r {shell_quote(fq1)} ") +
                    f"--validateMappings --seqBias --gcBias -p {args.threads} -o {shell_quote(outdir)}; fi")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open('w', encoding='utf-8') as w:
        w.write('#!/usr/bin/env bash\nset -euo pipefail\n')
        for c in cmds:
            w.write(c + '\n')
    os.chmod(out, 0o755)

    print(f'Wrote {len(cmds)} Salmon commands: {out}')
    if skipped:
        skip_path = Path('logs/salmon_skipped.tsv')
        with skip_path.open('w', encoding='utf-8') as w:
            w.write('sample_id\treason\tvalue\n')
            for row in skipped:
                w.write('\t'.join(map(str,row)) + '\n')
        print(f'Skipped {len(skipped)} samples. See {skip_path}')
        for row in skipped[:20]:
            print('SKIP', row)


if __name__ == '__main__':
    main()
