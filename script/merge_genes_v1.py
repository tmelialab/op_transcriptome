#!/usr/bin/env python3
"""
Merge all quant.genes.sf files (gene-level Salmon output) into 2 matrices:
  1. tpm_genes_matrix.tsv:      Name, Length, EffectiveLength, TPM per sample
  2. numreads_genes_matrix.tsv: Name, Length, EffectiveLength, NumReads per sample

Adapted from merge_quant_v4.py (transcript-level). Same safeguards:
  - Deduplicate gene IDs (second occurrence -> _dup2).
  - Fast C-level pandas to_csv.
  - Write to .tmp + verify counts + atomic rename.
"""

import os
import sys
import time
import subprocess
import pandas as pd
from pathlib import Path

BASE = Path(os.path.expanduser("~/september"))
LIST_FILE = BASE / "list_sample.txt"
QUANT_DIR = BASE / "quant"
OUT_DIR = BASE / "results"
QUANT_FILE = "quant.genes.sf"


def make_unique(ids):
    seen = {}
    out = []
    dups = []
    for x in ids:
        if x not in seen:
            seen[x] = 1
            out.append(x)
        else:
            seen[x] += 1
            new = f"{x}_dup{seen[x]}"
            dups.append((x, new))
            out.append(new)
    return out, dups


def verify(path, nrows_exp, ncols_exp):
    with open(path) as f:
        nlines = sum(1 for _ in f)
    h = subprocess.run(["head", "-1", str(path)], capture_output=True, text=True)
    ncols = h.stdout.count("\t") + 1
    print(f"  {path.name}: {nlines} lines, {ncols} cols, "
          f"{os.path.getsize(path)/1024/1024:.1f} MB", flush=True)
    assert nlines == nrows_exp, f"row mismatch {nlines} != {nrows_exp}"
    assert ncols == ncols_exp, f"col mismatch {ncols} != {ncols_exp}"


def main():
    start = time.time()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("[1/4] Reading sample list...", flush=True)
    samples = []
    with open(LIST_FILE) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            p = line.split("\t")
            if len(p) >= 2:
                samples.append((p[0], p[1]))
    print(f"  Found {len(samples)} samples", flush=True)

    print("[2/4] Reading quant.genes.sf files...", flush=True)
    tpm_cols = {}
    nr_cols = {}
    ref_ids = ref_len = ref_eff = None
    missing = []

    for i, (proj, srr) in enumerate(samples):
        qf = QUANT_DIR / proj / srr / QUANT_FILE
        if not qf.exists():
            missing.append(f"{proj}/{srr}")
            continue
        try:
            df = pd.read_csv(qf, sep="\t")
            if df.empty:
                missing.append(f"{proj}/{srr} (empty)")
                continue
            if ref_ids is None:
                ref_ids, dups = make_unique(df["Name"].tolist())
                if dups:
                    print(f"  NOTE: duplicate ID(s): {dups}", flush=True)
                ref_len = df["Length"].tolist()
                ref_eff = df["EffectiveLength"].tolist()
            if len(df) == len(ref_ids):
                tpm_cols[srr] = df["TPM"].to_numpy()
                nr_cols[srr] = df["NumReads"].to_numpy()
            else:
                tmp = df.drop_duplicates(subset="Name", keep="first").set_index("Name")
                base = [u.split("_dup")[0] if "_dup" in u else u for u in ref_ids]
                tpm_cols[srr] = tmp["TPM"].reindex(base).fillna(0.0).to_numpy()
                nr_cols[srr] = tmp["NumReads"].reindex(base).fillna(0.0).to_numpy()
        except Exception as e:
            missing.append(f"{proj}/{srr} (error: {e})")
        if (i + 1) % 50 == 0:
            print(f"  Read {i+1}/{len(samples)}...", flush=True)

    print(f"  OK: {len(tpm_cols)}, missing: {len(missing)}", flush=True)
    for m in missing[:10]:
        print(f"    {m}", flush=True)
    if not tpm_cols:
        sys.exit(1)

    print("[3/4] Building dataframes...", flush=True)
    cols = list(tpm_cols.keys())
    idx = pd.Index(ref_ids, name="transcript_id")
    meta = pd.DataFrame({"Length": ref_len, "EffectiveLength": ref_eff}, index=idx)
    tpm_out = pd.concat([meta, pd.DataFrame(tpm_cols, index=idx)], axis=1)
    nr_out = pd.concat([meta, pd.DataFrame(nr_cols, index=idx)], axis=1)
    if tpm_out.index.duplicated().any():
        print("  WARN: duplicated index remains, dropping dup rows (keep first)", flush=True)
        tpm_out = tpm_out[~tpm_out.index.duplicated(keep="first")]
        nr_out = nr_out[~nr_out.index.duplicated(keep="first")]
    print(f"  TPM {tpm_out.shape}, NR {nr_out.shape}", flush=True)

    print("[4/4] Writing...", flush=True)
    nrows_exp = len(tpm_out) + 1
    ncols_exp = len(tpm_out.columns) + 1
    tpm_path = OUT_DIR / "tpm_genes_matrix.tsv"
    nr_path = OUT_DIR / "numreads_genes_matrix.tsv"

    tmp = tpm_path.with_suffix(".tsv.tmp")
    tpm_out.to_csv(tmp, sep="\t", float_format="%.4f")
    verify(tmp, nrows_exp, ncols_exp)
    os.replace(tmp, tpm_path)
    print(f"  -> {tpm_path.name} OK", flush=True)
    del tpm_out

    tmp = nr_path.with_suffix(".tsv.tmp")
    nr_out.to_csv(tmp, sep="\t", float_format="%.2f")
    verify(tmp, nrows_exp, ncols_exp)
    os.replace(tmp, nr_path)
    print(f"  -> {nr_path.name} OK", flush=True)

    print(f"\nDone in {time.time()-start:.1f}s: {len(tpm_cols)} samples", flush=True)


if __name__ == "__main__":
    main()
