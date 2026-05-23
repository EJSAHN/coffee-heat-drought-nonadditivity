#!/usr/bin/env python3
import argparse
import gzip
import re
from pathlib import Path


def open_text(path):
    path = str(path)
    if path.endswith('.gz'):
        return gzip.open(path, 'rt', encoding='utf-8', errors='replace')
    return open(path, 'rt', encoding='utf-8', errors='replace')


def parse_header(header: str):
    # Header without leading '>'
    first = header.split()[0]
    tx = first.split('|')[0]

    patterns = [
        r'gene:([^\s]+)',
        r'gene=([^\]\s;]+)',
        r'\[gene=([^\]]+)\]',
        r'locus_tag=([^\]\s;]+)',
        r'\[locus_tag=([^\]]+)\]',
        r'GeneID:([0-9]+)',
        r'gene_id[ =]"?([^";\s]+)',
        r'Parent=gene:([^;\s]+)',
        r'Parent=gene-([^;\s]+)',
    ]
    gene = None
    for pat in patterns:
        m = re.search(pat, header)
        if m:
            gene = m.group(1).strip()
            break
    if not gene:
        # NCBI transcript accessions often have gene symbols in [gene=...]. If absent,
        # use the transcript id as a safe fallback; tximport then works at transcript level.
        gene = tx
    return tx, gene, header


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--species', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    with open_text(args.fasta) as fh:
        for line in fh:
            if line.startswith('>'):
                tx, gene, header = parse_header(line[1:].strip())
                rows.append((tx, gene))

    seen = set()
    n = 0
    with out.open('w', encoding='utf-8') as w:
        w.write('TXNAME\tGENEID\n')
        for tx, gene in rows:
            key = (tx, gene)
            if key in seen:
                continue
            seen.add(key)
            w.write(f'{tx}\t{gene}\n')
            n += 1

    print(f'Wrote {n} tx2gene rows for {args.species}: {out}')


if __name__ == '__main__':
    main()
