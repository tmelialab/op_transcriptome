#!/usr/bin/env python3
"""
Merge all quant.sf files from Salmon quantification into matrix files.
Output:
  1. tpm_matrix.tsv       - TPM values (transcript × sample)
  2. numreads_matrix.tsv  - NumReads values (transcript × sample)
  3. combined_matrix.tsv  - TPM + NumReads side by side (transcript × sample×2)
"""

import os
import sys
import time
import pandas as pd
import numpy as np
from pathlib import Path

# Paths
BASE = Path(os.path.expanduser("~/september"))
LIST_FILE = BASE / "list_sample.txt"
QUANT_DIR = BASE / "quant"
OUT_DIR = BASE / "results"

def main():
    start = time.time()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Read sample list
    print("[1/5] Reading sample list...")
    samples = []
    with open(LIST_FILE) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                proj, srr = parts[0], parts[1]
                samples.append((proj, srr))
    print(f"  Found {len(samples)} samples")

    # Read all quant.sf files
    print("[2/5] Reading quant.sf files...")
    tpm_dict = {}
    numreads_dict = {}
    missing = []
    empty = []

    for i, (proj, srr) in enumerate(samples):
        qf = QUANT_DIR / proj / srr / "quant.sf"
        if not qf.exists():
            missing.append(f"{proj}/{srr}")
            continue
        try:
            df = pd.read_csv(qf, sep="\t")
            if df.empty:
                empty.append(f"{proj}/{srr}")
                continue
            tpm_dict[srr] = df.set_index("Name")["TPM"]
            numreads_dict[srr] = df.set_index("Name")["NumReads"]
        except Exception as e:
            missing.append(f"{proj}/{srr} (error: {e})")

        if (i + 1) % 50 == 0:
            print(f"  Read {i+1}/{len(samples)} files...")

    print(f"  Successfully read: {len(tpm_dict)}")
    print(f"  Missing: {len(missing)}")
    if missing:
        print(f"    Examples: {missing[:5]}")
    print(f"  Empty: {len(empty)}")

    if len(tpm_dict) == 0:
        print("ERROR: No quant.sf files found!")
        sys.exit(1)

    # Build TPM matrix
    print("[3/5] Building TPM matrix...")
    tpm_df = pd.DataFrame(tpm_dict)
    tpm_df.index.name = "transcript_id"
    tpm_df = tpm_df.fillna(0.0)
    print(f"  Shape: {tpm_df.shape[0]} transcripts × {tpm_df.shape[1]} samples")

    # Build NumReads matrix
    print("[4/5] Building NumReads matrix...")
    numreads_df = pd.DataFrame(numreads_dict)
    numreads_df.index.name = "transcript_id"
    numreads_df = numreads_df.fillna(0.0)
    print(f"  Shape: {numreads_df.shape[0]} transcripts × {numreads_df.shape[1]} samples")

    # Build Combined matrix (TPM and NumReads side by side)
    print("[5/5] Building combined matrix...")
    # Create paired columns: SRR_TPM, SRR_NumReads
    combined_cols = []
    for srr in tpm_df.columns:
        combined_cols.append(f"{srr}_TPM")
        combined_cols.append(f"{srr}_NumReads")

    combined_data = np.zeros((tpm_df.shape[0], len(combined_cols)))
    for i, srr in enumerate(tpm_df.columns):
        combined_data[:, i*2] = tpm_df[srr].values
        combined_data[:, i*2+1] = numreads_df[srr].values

    combined_df = pd.DataFrame(combined_data, index=tpm_df.index, columns=combined_cols)
    combined_df.index.name = "transcript_id"

    # Write output files
    print("Writing output files...")
    tpm_path = OUT_DIR / "tpm_matrix.tsv"
    numreads_path = OUT_DIR / "numreads_matrix.tsv"
    combined_path = OUT_DIR / "combined_matrix.tsv"

    tpm_df.to_csv(tpm_path, sep="\t", float_format="%.4f")
    print(f"  TPM matrix: {tpm_path} ({os.path.getsize(tpm_path) / 1024 / 1024:.1f} MB)")

    numreads_df.to_csv(numreads_path, sep="\t", float_format="%.2f")
    print(f"  NumReads matrix: {numreads_path} ({os.path.getsize(numreads_path) / 1024 / 1024:.1f} MB)")

    combined_df.to_csv(combined_path, sep="\t", float_format="%.4f")
    print(f"  Combined matrix: {combined_path} ({os.path.getsize(combined_path) / 1024 / 1024:.1f} MB)")

    elapsed = time.time() - start
    print(f"\nDone in {elapsed:.1f} seconds")
    print(f"Total: {tpm_df.shape[0]} transcripts × {tpm_df.shape[1]} samples")

if __name__ == "__main__":
    main()
