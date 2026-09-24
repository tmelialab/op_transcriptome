#!/usr/bin/env python3
"""
Merge all quant.sf files into 2 matrices (v3, fixed):
  1. tpm_matrix.tsv:      transcript_id, Length, EffectiveLength, TPM per sample
  2. numreads_matrix.tsv: transcript_id, Length, EffectiveLength, NumReads per sample

Fixes vs v2:
  - Deduplicate transcript IDs (reference contains unassigned_transcript_568 twice
    with different Length). Second occurrence gets suffix _dup2, logged to stdout.
  - Write to .tmp files first, verify row/col counts, then atomic rename.
  - Single-pass matrix build (lower memory).
"""

import os
import sys
import time
import pandas as pd
import numpy as np
from pathlib import Path

BASE = Path(os.path.expanduser("~/september"))
LIST_FILE = BASE / "list_sample.txt"
QUANT_DIR = BASE / "quant"
OUT_DIR = BASE / "results"


def make_unique(ids):
    """Make a list of IDs unique by appending _dup2, _dup3... to repeats."""
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
            parts = line.split("\t")
            if len(parts) >= 2:
                samples.append((parts[0], parts[1]))
    print(f"  Found {len(samples)} samples", flush=True)

    print("[2/4] Reading quant.sf files...", flush=True)
    tpm_dict = {}
    numreads_dict = {}
    ref_ids = None
    ref_len = None
    ref_efflen = None
    dup_reported = False
    missing = []

    for i, (proj, srr) in enumerate(samples):
        qf = QUANT_DIR / proj / srr / "quant.sf"
        if not qf.exists():
            missing.append(f"{proj}/{srr}")
            continue
        try:
            df = pd.read_csv(qf, sep="\t")
            if df.empty:
                missing.append(f"{proj}/{srr} (empty)")
                continue
            if ref_ids is None:
                # First file defines canonical row order + metadata
                raw_ids = df["Name"].tolist()
                uniq_ids, dups = make_unique(raw_ids)
                if dups:
                    print(f"  NOTE: {len(dups)} duplicate transcript ID(s) in reference:", flush=True)
                    for orig, new in dups:
                        print(f"    {orig} -> {new}", flush=True)
                ref_ids = uniq_ids
                ref_len = df["Length"].tolist()
                ref_efflen = df["EffectiveLength"].tolist()
                dup_reported = True
            # Align this sample's values to canonical order.
            # Use groupby-friendly approach: drop duplicate labels keeping first,
            # then reindex. But to preserve dup rows distinctly, map positionally
            # when row count matches.
            if len(df) == len(ref_ids):
                tpm_dict[srr] = df["TPM"].tolist()
                numreads_dict[srr] = df["NumReads"].tolist()
            else:
                tmp = df.drop_duplicates(subset="Name", keep="first").set_index("Name")
                # reindex to de-duplicated canonical ids
                canon_nodup = []
                seen = set()
                for u in ref_ids:
                    base = u.split("_dup")[0] if "_dup" in u else u
                    canon_nodup.append(base)
                tpm_dict[srr] = tmp["TPM"].reindex(canon_nodup).fillna(0.0).tolist()
                numreads_dict[srr] = tmp["NumReads"].reindex(canon_nodup).fillna(0.0).tolist()
        except Exception as e:
            missing.append(f"{proj}/{srr} (error: {e})")

        if (i + 1) % 50 == 0:
            print(f"  Read {i+1}/{len(samples)} files...", flush=True)

    print(f"  Successfully read: {len(tpm_dict)}", flush=True)
    print(f"  Missing/skipped: {len(missing)}", flush=True)
    for m in missing[:10]:
        print(f"    {m}", flush=True)

    if len(tpm_dict) == 0:
        print("ERROR: No quant.sf files found!")
        sys.exit(1)

    print("[3/4] Building matrices...", flush=True)
    cols = list(tpm_dict.keys())
    tpm_mat = np.column_stack([tpm_dict[c] for c in cols])
    nr_mat = np.column_stack([numreads_dict[c] for c in cols])
    print(f"  TPM shape: {tpm_mat.shape}, NumReads shape: {nr_mat.shape}", flush=True)

    header = ["transcript_id", "Length", "EffectiveLength"] + cols

    print("[4/4] Writing output files (via .tmp + verify)...", flush=True)
    tpm_path = OUT_DIR / "tpm_matrix.tsv"
    nr_path = OUT_DIR / "numreads_matrix.tsv"

    def write_matrix(path, mat, fmt):
        tmp = path.with_suffix(path.suffix + ".tmp")
        with open(tmp, "w") as out:
            out.write("\t".join(header) + "\n")
            for i in range(len(ref_ids)):
                row = [ref_ids[i], str(ref_len[i]), f"{ref_efflen[i]:.4f}" if fmt == "%.4f" else f"{ref_efflen[i]:.2f}"]
                if fmt == "%.4f":
                    row += [f"{v:.4f}" for v in mat[i]]
                else:
                    row += [f"{v:.2f}" for v in mat[i]]
                out.write("\t".join(row) + "\n")
        # verify
        with open(tmp) as f:
            nlines = sum(1 for _ in f)
        import subprocess
        ncol = subprocess.run(["head", "-1", str(tmp)], capture_output=True, text=True)
        ncols = ncol.stdout.count("\t") + 1
        print(f"  {tmp.name}: {nlines} lines, {ncols} cols, {os.path.getsize(tmp)/1024/1024:.1f} MB", flush=True)
        assert nlines == len(ref_ids) + 1, f"row count mismatch: {nlines} != {len(ref_ids)+1}"
        assert ncols == len(header), f"col count mismatch: {ncols} != {len(header)}"
        os.replace(tmp, path)
        print(f"  -> {path.name} OK", flush=True)

    write_matrix(tpm_path, tpm_mat, "%.4f")
    write_matrix(nr_path, nr_mat, "%.2f")

    elapsed = time.time() - start
    print(f"\nDone in {elapsed:.1f} seconds", flush=True)
    print(f"Total: {len(ref_ids)} transcripts x {len(cols)} samples", flush=True)


if __name__ == "__main__":
    main()
