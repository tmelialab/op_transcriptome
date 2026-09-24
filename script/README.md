# script/

Pipeline + analysis scripts. This folder is the source of truth — edit here,
then copy to the cluster project dir (`cp script/* ~/september/scripts/`,
re-runnable). All compute scripts are idempotent: existing outputs are
skipped, so re-submits are safe.

## Index builds (run once each)

| File | What it does | Idempotency note |
|---|---|---|
| `build_star_index.sbatch` | `STAR --runMode genomeGenerate` for EG11 (`--sjdbOverhang 100`, `--genomeSAindexNbases 14`), 48 CPU / 240G | Re-runnable; output dir created with `mkdir -p`. Re-running rebuilds the index. |
| `build_salmon_index.sbatch` | Fixes the NCBI GTF (strips empty `transcript_id ""` that aborts `gffread ≥ 0.12`, drops the one `"?"`-strand trans-spliced `rps12` line), builds transcripts with `gffread`, concatenates genome as decoys (`EG11_gentrome.fa`), then `salmon index` | Fixed GTF is only written if missing (`-s` guard). Rest re-runs deterministically. |

## Mapping + quant

| File | What it does | Idempotency note |
|---|---|---|
| `map_chunk.sbatch START END` | STAR 2-pass (`--twopassMode Basic`, `--alignIntronMax 500000`) over `list_sample.txt` lines START–END, 8 CPU / 48G. Writes `$OUT/<SRR>Aligned.sortedByCoord.out.bam` + `samtools index` + `mapping_statistics.txt` + `num_lines.txt` | Skips samples whose BAM already exists (`-s` check). Safe to resubmit after timeouts. |
| `quant_chunk.sbatch START END` | `salmon quant -l A` over lines START–END, 8 CPU / 16G, GTF passed via `-g` | Skips samples whose `quant.sf` already exists. Safe to resubmit. |
| `mapx86.sbatch` | Legacy per-sample STAR job array (`1-407`), same flags as chunks | Same skip-if-BAM-exists guard. Kept for reference; chunks are preferred under tight submit limits. |
| `quant.sbatch` | Legacy per-sample salmon job array (`1-407`) | Same skip-if-`quant.sf`-exists guard. Kept for reference. |
| `submit_driver.sh` | Submits map chunks (25 samples) then quant chunks (50 samples) in passes (max 4 passes), never more than `MAXJOBS=4` of your jobs queued; tracks submissions in `logs/driver_submitted.txt` | Fully re-runnable: only chunks with missing outputs and not yet submitted are queued. Kill/restart anytime (run in `screen`). Clears the tracking file after each fully-completed pass. Uses `squeue -u "$USER"` — no hardcoded username. |

Parser contract (all mapping/quant scripts): each `list_sample.txt` line is
split into exactly 4 tab-separated fields `PROJ SAMPLEID TYPE STRAND`.
`TYPE=single` reads `<SAMPLE>.fastq.gz`; anything else reads
`<SAMPLE>_1/_2.fastq.gz`.

## Merge + analysis (login node, stdlib + pandas/numpy)

| File | What it does | Status |
|---|---|---|
| `merge_quant_v4.py` | Merges all `quant/<PROJ>/<SRR>/quant.sf` into `results/tpm_matrix.tsv` + `results/numreads_matrix.tsv` (transcript_id, Length, EffectiveLength + one column per sample). Dedups transcript IDs (`_dup2` suffix), pure-C `to_csv`, writes via `.tmp` + row/col verification + atomic rename; tolerates missing samples (listed, not fatal) | **Canonical transcript-level merge — use this.** |
| `merge_genes_v1.py` | Same safeguards for `quant.genes.sf` → `tpm_genes_matrix.tsv` + `numreads_genes_matrix.tsv` | Use for gene-level (needs gene-level salmon outputs). |
| `merge_quant.py`, `merge_quant_v2.py`, `merge_quant_v3.py` | Earlier merge iterations | History only; superseded by v4 (v3's pure-Python writer was too slow). |
| `analyze_quant.py` | Read-only report over the merged matrices: global TPM/NumReads stats, per-sample expressed-gene counts, top-20 transcripts, TPM↔NumReads correlation, per-PROJNA summary, TPM-bin distribution, ±2σ outlier samples | No inputs modified. Expects `results/tpm_matrix.tsv`, `numreads_matrix.tsv` (+ `combined_matrix.tsv` size line if present). |

All Python scripts locate data via `~/september` (`os.path.expanduser`),
so they work for any username with the documented layout.

## Sample + fetch helpers

| File | What it does |
|---|---|
| `list_sample.txt` | 407 samples, 4 tab-separated columns, no header. `#`/blank lines skipped. Tabs required; trailing newline on every line required. |
| `PRJNA778743.urls.txt` / `PRJNA778743.md5` | Upstream URLs + checksums to (re)fetch the PRJNA778743 fastq set. |

## `orig_dosen/`

Unmodified originals (`mapx86.slurm`, `quant.sh`, `genome_index.pjsub`)
from a different scheduler/site (different queues, `--mem` units, absolute
`/home/…` paths, 3-column sample parsing). Historical reference only —
do not submit as-is. The ported equivalents at this level already contain
the fixes: 4-column parsing, corrected output prefix/path, `zcat`
piping, `samtools index` on the real BAM path.

## Site adaptation checklist (per fresh cluster)

1. The two conda setup lines in each `.sbatch`
   (`source …/conda.sh`, `conda activate …`) → your env.
2. `#SBATCH -p`, `--mem`, `--time` → your partitions/limits
   (current values: `short`, map 8 CPU/48G, quant 8 CPU/16G, indexes
   48 CPU/240G and 16 CPU/64G, 23 h walltime).
3. Submit from the project root (`~/september`) so relative `logs/`
   directives resolve; `mkdir -p logs` first.
