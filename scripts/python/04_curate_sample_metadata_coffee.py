#!/usr/bin/env python3
import argparse
import re
from pathlib import Path
import pandas as pd


def code_prefix(value: str) -> str:
    m = re.match(r"^\d+", str(value or ""))
    return m.group(0) if m else ""


def set_row(row, code):
    study = row["study_accession"]
    genotype = row["genotype"]
    row["include"] = "no"
    row["water"] = "UNMAPPED"
    row["stage"] = "UNMAPPED"

    if study == "PRJNA787748":
        if genotype == "CL153" and code == "1":
            row["water"], row["stage"], row["include"] = "WW", "T25", "yes"
        elif genotype == "CL153" and code == "3":
            row["water"], row["stage"], row["include"] = "SWD", "T25", "yes"
        elif genotype == "Icatu" and code == "7":
            row["water"], row["stage"], row["include"] = "WW", "T25", "yes"
        elif genotype == "Icatu" and code == "9":
            row["water"], row["stage"], row["include"] = "SWD", "T25", "yes"

    elif study == "PRJNA1087442":
        row["water"] = "WW"
        if code == "31": row["stage"], row["include"] = "T37", "yes"
        elif code == "43": row["stage"], row["include"] = "T42", "yes"
        elif code == "55": row["stage"], row["include"] = "REC14", "yes"

    elif study == "PRJNA1087679":
        row["water"] = "SWD"
        if code in {"32", "33"}: row["stage"], row["include"] = "T37", "yes"
        elif code == "45": row["stage"], row["include"] = "T42", "yes"
        elif code == "57": row["stage"], row["include"] = "REC14", "yes"

    elif study == "PRJNA1088119":
        row["water"] = "WW"
        if code == "25": row["stage"], row["include"] = "T37", "yes"
        elif code == "37": row["stage"], row["include"] = "T42", "yes"
        elif code == "49": row["stage"], row["include"] = "REC14", "yes"

    elif study == "PRJNA1135679":
        row["water"] = "SWD"
        if code in {"26", "27"}: row["stage"], row["include"] = "T37", "yes"
        elif code == "39": row["stage"], row["include"] = "T42", "yes"
        elif code == "51": row["stage"], row["include"] = "REC14", "yes"
    return row


def main():
    ap = argparse.ArgumentParser(description="Curate coffee heat-drought sample metadata from ENA rich metadata aliases.")
    ap.add_argument("--metadata", default="config/sample_metadata_curated.tsv")
    ap.add_argument("--rich", default="data/raw/sra_metadata/metadata_triage_rich.csv")
    ap.add_argument("--out", default="config/sample_metadata_curated.tsv")
    args = ap.parse_args()

    meta = pd.read_csv(args.metadata, sep="\t", dtype=str).fillna("")
    rich = pd.read_csv(args.rich, dtype=str).fillna("")
    rich_by_run = rich.set_index("run_accession").to_dict("index")

    rows = []
    for _, r in meta.iterrows():
        row = r.to_dict()
        rr = rich_by_run.get(row["run_accession"], {})
        code = code_prefix(rr.get("sample_alias", ""))
        rows.append(set_row(row, code))
    out = pd.DataFrame(rows)

    for _, idxs in out[out["include"].str.lower() == "yes"].groupby(["genotype", "water", "stage"]).groups.items():
        sorted_idx = sorted(list(idxs), key=lambda i: out.loc[i, "run_accession"])
        for rep, i in enumerate(sorted_idx, start=1):
            out.loc[i, "replicate"] = str(rep)
            out.loc[i, "sample_id"] = f"{out.loc[i, 'genotype']}_{out.loc[i, 'water']}_{out.loc[i, 'stage']}_R{rep}"
    for i in out[out["include"].str.lower() != "yes"].index:
        out.loc[i, "sample_id"] = f"{out.loc[i, 'genotype']}_EXCLUDE_{out.loc[i, 'run_accession']}"

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)
    included = out[out["include"].str.lower() == "yes"]
    print(included.groupby(["index_species", "genotype", "water", "stage"]).size().reset_index(name="n").to_string(index=False))
    print(f"Included samples: {len(included)}")


if __name__ == "__main__":
    main()
