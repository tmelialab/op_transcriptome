#!/usr/bin/env python3
"""
Merge all quant.sf files into 2 matrices as requested by dosen:
  1. tpm_matrix.tsv:     transcript_id, Length, EffectiveLength, TPM per sample
  2. numreads_matrix.tsv: transcript_id, Length, EffectiveLength, NumReads per sample
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

def main():
    start = time.time()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Read sample list
    print("[1/4] Reading sample list...")
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
    print("[2/4] Reading quant.sf files...")
    tpm_dict = {}
    numreads_dict = {}
    length_col = None
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
            # Set index to Name (transcript_id)
            df = df.set_index("Name")
            tpm_dict[srr] = df["TPM"]
            numreads_dict[srr] = df["NumReads"]
            # Save Length and EffectiveLength from first file (same for all)
            if length_col is None:
                length_col = df["Length"]
                efflength_col = df["EffectiveLength"]
        except Exception as e:
            missing.append(f"{proj}/{srr} (error: {e})")

        if (i + 1) % 50 == 0:
            print(f"  Read {i+1}/{len(samples)} files...")

    print(f"  Successfully read: {len(tpm_dict)}")
    print(f"  Missing/skipped: {len(missing)}")
    if missing:
        print(f"    Examples: {missing[:5]}")

    if len(tpm_dict) == 0:
        print("ERROR: No quant.sf files found!")
        sys.exit(1)

    # Build matrices
    print("[3/4] Building matrices...")
    tpm_df = pd.DataFrame(tpm_dict)
    numreads_df = pd.DataFrame(numreads_dict)

    # Add Length and EffectiveLength columns (from first sample, they're the same across all)
    first_sample = list(tpm_dict.keys())[0]
    qf_first = None
    for proj, srr in samples:
        if srr == first_sample:
            qf_first = QUANT_DIR / proj / srr / "quant.sf"
            break
    meta_df = pd.read_csv(qf_first, sep="\t").set_index("Name")
    length_series = meta_df["Length"]
    efflength_series = meta_df["EffectiveLength"]

    # Build TPM matrix: transcript_id, Length, EffectiveLength, sample1_TPM, sample2_TPM, ...
    print("  TPM matrix...")
    tpm_out = pd.DataFrame({
        "Length": length_col if length_col is not None else length_col,
        "EffectiveLength": df["EffectiveLength"] if "EffectiveLength" in df.columns else 0,
    })
    # Actually, let me rebuild properly
    # Use the first quant.sf to get Length and EffectiveLength
    first_srr = list(tpm_dict.keys())[0]
    first_proj = [p for p, s in samples if s == first_srr][0]
    first_df = pd.read_csv(QUANT_DIR / first_proj / first_srr / "quant.sf", sep="\t").set_index("Name")

    length_series = first_df["Length"]
    efflength_series = first_df["EffectiveLength"]

    # Build TPM matrix
    tpm_matrix = pd.DataFrame(tpm_dict)
    tpm_matrix.index.name = "transcript_id"

    # Add Length and EffectiveLength columns at the front
    tpm_out = pd.concat([length_series, efflength_series, tpm_matrix], axis=1)
    tpm_out.columns = ["Length", "EffectiveLength"] + list(tpm_matrix.columns)
    tpm_out = tpm_out.fillna(0.0)

    # Build NumReads matrix
    numreads_matrix = pd.DataFrame(numreads_dict)
    numreads_matrix.index.name = "transcript_id"

    nr_out = pd.concat([length_series, efflength_series, numreads_matrix], axis=1)
    nr_out.columns = ["Length", "EffectiveLength"] + list(numreads_matrix.columns)
    nr_out = nr_out.fillna(0.0)

    print(f"  TPM matrix: {tpm_out.shape[0]} transcripts × {tpm_out.shape[1]} columns")
    print(f"  NumReads matrix: {nr_out.shape[0]} transcripts × {nr_out.shape[1]} columns")

    # Write output files
    print("[4/4] Writing output files...")
    tpm_path = OUT_DIR / "tpm_matrix.tsv"
    nr_path = OUT_DIR / "numreads_matrix.tsv"

    tpm_out.to_csv(tpm_path, sep="\t", float_format="%.4f")
    print(f"  TPM: {tpm_path} ({os.path.getsize(tpm_path) / 1024 / 1024:.1f} MB)")

    nr_out.to_csv(nr_path, sep="\t", float_format="%.2f")
    print(f"  NumReads: {nr_path} ({os.path.getsize(nr_path) / 1024 / 1024:.1f} MB)")

    elapsed = time.time() - start
    print(f"\nDone in {elapsed:.1f} seconds")
    print(f"Total: {tpm_out.shape[0]} transcripts × {len(samples)} samples")
    print(f"\nFiles: {tpm_path}, {nr_path}")

if __name__ == "__main__":
    main()
