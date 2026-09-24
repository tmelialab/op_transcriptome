#!/usr/bin/env python3
"""
Analyze merged quantification results.
"""

import os
import pandas as pd
import numpy as np
from pathlib import Path

BASE = Path(os.path.expanduser("~/september"))
RESULTS = BASE / "results"
LIST_FILE = BASE / "list_sample.txt"

def main():
    print("=" * 70)
    print("  RNA-seq QUANTIFICATION RESULTS ANALYSIS")
    print("=" * 70)

    # Load matrices
    print("\n[1] Loading matrices...")
    tpm = pd.read_csv(RESULTS / "tpm_matrix.tsv", sep="\t", index_col=0)
    numreads = pd.read_csv(RESULTS / "numreads_matrix.tsv", sep="\t", index_col=0)
    print(f"  TPM matrix: {tpm.shape[0]} transcripts × {tpm.shape[1]} samples")
    print(f"  NumReads matrix: {numreads.shape[0]} transcripts × {numreads.shape[1]} samples")

    # Read sample metadata
    meta = {}
    with open(LIST_FILE) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                meta[parts[1]] = {"proj": parts[0], "type": parts[2] if len(parts) > 2 else "unknown", "strand": parts[3] if len(parts) > 3 else "unknown"}

    # =========================================================================
    print("\n" + "=" * 70)
    print("  GLOBAL BASIC STATISTICS")
    print("=" * 70)

    # TPM stats
    all_tpm = tpm.values.flatten()
    all_tpm_nonzero = all_tpm[all_tpm > 0]
    print(f"\n--- TPM (all values) ---")
    print(f"  Total values:     {len(all_tpm):,}")
    print(f"  Non-zero:        {len(all_tpm_nonzero):,} ({100*len(all_tpm_nonzero)/len(all_tpm):.1f}%)")
    print(f"  Mean:            {np.mean(all_tpm):.4f}")
    print(f"  Median:          {np.median(all_tpm):.4f}")
    print(f"  Max:             {np.max(all_tpm):.4f}")
    print(f"  Std:             {np.std(all_tpm):.4f}")

    print(f"\n--- TPM (non-zero only) ---")
    print(f"  Mean:            {np.mean(all_tpm_nonzero):.4f}")
    print(f"  Median:          {np.median(all_tpm_nonzero):.4f}")
    print(f"  25th percentile: {np.percentile(all_tpm_nonzero, 25):.4f}")
    print(f"  75th percentile: {np.percentile(all_tpm_nonzero, 75):.4f}")

    # NumReads stats
    all_nr = numreads.values.flatten()
    all_nr_nonzero = all_nr[all_nr > 0]
    print(f"\n--- NumReads (all values) ---")
    print(f"  Total values:     {len(all_nr):,}")
    print(f"  Non-zero:        {len(all_nr_nonzero):,} ({100*len(all_nr_nonzero)/len(all_nr):.1f}%)")
    print(f"  Mean:            {np.mean(all_nr):.2f}")
    print(f"  Median:          {np.median(all_nr):.2f}")
    print(f"  Max:             {np.max(all_nr):.2f}")

    # =========================================================================
    print("\n" + "=" * 70)
    print("  PER-SAMPLE STATISTICS")
    print("=" * 70)

    # Per sample: how many genes expressed (TPM > 0)
    genes_per_sample = (tpm > 0).sum(axis=0)
    tpm_sum_per_sample = tpm.sum(axis=0)

    print(f"\n--- Expressed genes (TPM > 0) per sample ---")
    print(f"  Mean:   {genes_per_sample.mean():.0f} genes")
    print(f"  Median: {genes_per_sample.median():.0f} genes")
    print(f"  Min:    {genes_per_sample.min()} genes (sample: {genes_per_sample.idxmin()})")
    print(f"  Max:    {genes_per_sample.max()} genes (sample: {genes_per_sample.idxmax()})")
    print(f"  Std:    {genes_per_sample.std():.0f}")

    print(f"\n--- Total TPM per sample ---")
    print(f"  Mean:   {tpm_sum_per_sample.mean():.1f}")
    print(f"  Median: {tpm_sum_per_sample.median():.1f}")
    print(f"  Min:    {tpm_sum_per_sample.min():.1f} (sample: {tpm_sum_per_sample.idxmin()})")
    print(f"  Max:    {tpm_sum_per_sample.max():.1f} (sample: {tpm_sum_per_sample.idxmax()})")

    # =========================================================================
    print("\n" + "=" * 70)
    print("  TOP 20 HIGHEST-EXPRESSED TRANSCRIPTS (mean TPM)")
    print("=" * 70)

    mean_tpm = tpm.mean(axis=1).sort_values(ascending=False)
    print(f"\n{'Rank':<6}{'Transcript':<25}{'Mean TPM':>12}{'Median TPM':>12}{'Max TPM':>12}{'% Samples >0':>15}")
    print("-" * 82)
    for i, (tid, val) in enumerate(mean_tpm.head(20).items()):
        row_tpm = tpm.loc[tid]
        pct_expressed = 100 * (row_tpm > 0).sum() / len(row_tpm)
        print(f"{i+1:<6}{tid:<25}{val:>12.2f}{row_tpm.median():>12.2f}{row_tpm.max():>12.2f}{pct_expressed:>14.1f}%")

    # =========================================================================
    print("\n" + "=" * 70)
    print("  TOP 20 TRANSCRIPTS BY READ COUNT (mean NumReads)")
    print("=" * 70)

    mean_nr = numreads.mean(axis=1).sort_values(ascending=False)
    print(f"\n{'Rank':<6}{'Transcript':<25}{'Mean Reads':>12}{'Median Reads':>13}{'Max Reads':>12}")
    print("-" * 68)
    for i, (tid, val) in enumerate(mean_nr.head(20).items()):
        row_nr = numreads.loc[tid]
        print(f"{i+1:<6}{tid:<25}{val:>12.1f}{row_nr.median():>13.1f}{row_nr.max():>12.1f}")

    # =========================================================================
    print("\n" + "=" * 70)
    print("  TPM vs NumReads CORRELATION (per sample)")
    print("=" * 70)

    tpm_sample_sum = tpm.sum(axis=0)
    nr_sample_sum = numreads.sum(axis=0)
    corr = tpm_sample_sum.corr(nr_sample_sum)
    print(f"\n  Pearson correlation (total TPM vs total NumReads): {corr:.4f}")
    print(f"  (Should be close to 1.0 since TPM is normalized from NumReads)")

    # =========================================================================
    print("\n" + "=" * 70)
    print("  SUMMARY PER PROJECT")
    print("=" * 70)

    # Group by project
    proj_stats = []
    for srr in tpm.columns:
        if srr in meta:
            proj = meta[srr]["proj"]
            proj_stats.append({
                "project": proj,
                "sample": srr,
                "genes_expressed": int((tpm[srr] > 0).sum()),
                "total_tpm": tpm[srr].sum(),
                "mean_tpm": tpm[srr].mean(),
            })

    proj_df = pd.DataFrame(proj_stats)
    proj_summary = proj_df.groupby("project").agg(
        n_samples=("sample", "count"),
        avg_genes_expressed=("genes_expressed", "mean"),
        avg_total_tpm=("total_tpm", "mean"),
        avg_mean_tpm=("mean_tpm", "mean"),
    ).sort_values("n_samples", ascending=False)

    print(f"\n{'Project':<20}{'#Samples':>10}{'Avg Genes >0':>14}{'Avg Total TPM':>15}{'Avg Mean TPM':>14}")
    print("-" * 73)
    for proj, row in proj_summary.iterrows():
        print(f"{proj:<20}{int(row['n_samples']):>10}{row['avg_genes_expressed']:>14.0f}{row['avg_total_tpm']:>15.1f}{row['avg_mean_tpm']:>14.4f}")

    # =========================================================================
    print("\n" + "=" * 70)
    print("  EXPRESSION DISTRIBUTION (TPM bins)")
    print("=" * 70)

    bins = [-0.001, 0, 0.01, 0.1, 1, 5, 10, 50, 100, 500, 1000, float("inf")]
    labels = ["=0", "0-0.01", "0.01-0.1", "0.1-1", "1-5", "5-10", "10-50", "50-100", "100-500", "500-1000", ">1000"]

    # Per sample distribution (average)
    avg_dist = []
    for srr in tpm.columns:
        counts = pd.cut(tpm[srr], bins=bins, labels=labels, right=False).value_counts()
        avg_dist.append(counts)

    avg_dist_df = pd.DataFrame(avg_dist).mean(axis=0)
    total_transcripts = len(tpm)

    print(f"\n{'TPM Range':<15}{'Avg Count':>12}{'% of Transcripts':>18}")
    print("-" * 45)
    for label in labels:
        count = avg_dist_df[label]
        pct = 100 * count / total_transcripts
        print(f"{label:<15}{count:>12.0f}{pct:>17.1f}%")

    # =========================================================================
    print("\n" + "=" * 70)
    print("  OUTLIER SAMPLES (active gene count)")
    print("=" * 70)

    mean_genes = genes_per_sample.mean()
    std_genes = genes_per_sample.std()
    threshold_low = mean_genes - 2 * std_genes
    threshold_high = mean_genes + 2 * std_genes

    outliers_low = genes_per_sample[genes_per_sample < threshold_low]
    outliers_high = genes_per_sample[genes_per_sample > threshold_high]

    print(f"\n  Mean genes expressed: {mean_genes:.0f}")
    print(f"  Threshold low (< mean-2*std):  {threshold_low:.0f}")
    print(f"  Threshold high (> mean+2*std): {threshold_high:.0f}")

    if len(outliers_low) > 0:
        print(f"\n  ⚠️  LOW outliers ({len(outliers_low)} samples):")
        for srr, val in outliers_low.sort_values().items():
            proj = meta.get(srr, {}).get("proj", "?")
            print(f"    {srr} ({proj}): {val:.0f} genes")
    else:
        print(f"\n  ✅ No low outliers")

    if len(outliers_high) > 0:
        print(f"\n  ⚠️  HIGH outliers ({len(outliers_high)} samples):")
        for srr, val in outliers_high.sort_values(ascending=False).items():
            proj = meta.get(srr, {}).get("proj", "?")
            print(f"    {srr} ({proj}): {val:.0f} genes")
    else:
        print(f"\n  ✅ No high outliers")

    # =========================================================================
    print("\n" + "=" * 70)
    print("  DONE!")
    print("=" * 70)
    print(f"\nOutput files in: {RESULTS}/")
    print(f"  tpm_matrix.tsv      ({os.path.getsize(RESULTS / 'tpm_matrix.tsv') / 1024 / 1024:.1f} MB)")
    print(f"  numreads_matrix.tsv ({os.path.getsize(RESULTS / 'numreads_matrix.tsv') / 1024 / 1024:.1f} MB)")
    print(f"  combined_matrix.tsv ({os.path.getsize(RESULTS / 'combined_matrix.tsv') / 1024 / 1024:.1f} MB)")

if __name__ == "__main__":
    main()
