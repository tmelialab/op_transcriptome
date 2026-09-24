# op_transcriptome

RNA-seq mapping + quantification pipeline (STAR 2-pass mapping, salmon
quantification) with merge + summary analysis scripts.

Reference: *Elaeis guineensis* EG11 (`GCF_000442705.2`).
Sample list: 407 samples (`script/list_sample.txt`).

This repo is the source of truth for these scripts. Everything needed to
reproduce a run lives here; large/site-local data (fastq, indexes, BAMs,
`quant.sf` outputs) is intentionally not committed.

## Layout

```
script/
  build_star_index.sbatch     STAR genome index build (EG11)
  build_salmon_index.sbatch   salmon gentrome + decoy index build (EG11)
  map_chunk.sbatch            STAR mapping, serial per-chunk (START END)
  quant_chunk.sbatch          salmon quant, serial per-chunk (START END)
  mapx86.sbatch               STAR mapping, legacy per-sample job array
  quant.sbatch                salmon quant, legacy per-sample job array
  submit_driver.sh            driver: submits map/quant chunks, keeps jobs ≤ limit
  list_sample.txt             407 samples, 4 tab-separated columns
  analyze_quant.py            read-only summary of merged matrices
  merge_quant.py              merge quant.sf (v1, history)
  merge_quant_v2.py           merge quant.sf (v2, history)
  merge_quant_v3.py           merge quant.sf (v3, history)
  merge_quant_v4.py           merge quant.sf (v4, canonical transcript-level)
  merge_genes_v1.py           merge quant.genes.sf (gene-level)
  PRJNA778743.md5             checksums for PRJNA778743 refetch
  PRJNA778743.urls.txt        URLs for PRJNA778743 refetch
  orig_dosen/                 unmodified originals (historical reference only)
```

## Requirements

- Slurm cluster with STAR, salmon, samtools, gffread available
  (here via one conda env; adjust the two conda setup lines in each
  `.sbatch` to your cluster).
- Python 3 with pandas + numpy (merge/analysis scripts).
- Reference files (not in git, fetch once):
  `GCF_000442705.2_EG11_genomic.fna` and
  `GCF_000442705.2_EG11_genomic.gtf` (NCBI, *E. guineensis* EG11).

## Setup (idempotent — safe to re-run)

All commands below are safe to run more than once: directories use
`mkdir -p`, copies overwrite with identical content, and every pipeline
script skips outputs that already exist.

```bash
# 1. Project layout (run from your home on the cluster)
mkdir -p ~/september/{scripts,logs,reference,data/rna,star_index,salmon,map,quant,results}

# 2. Install scripts + sample list from this repo
cp script/*.sbatch script/*.sh script/*.py ~/september/scripts/
cp script/list_sample.txt ~/september/list_sample.txt

# 3. Put the EG11 reference files in place (fetch once, keep)
#    ~/september/reference/GCF_000442705.2_EG11_genomic.fna
#    ~/september/reference/GCF_000442705.2_EG11_genomic.gtf

# 4. Stage fastq as ~/september/data/rna/<PROJ>/<SAMPLE>[_1,_2].fastq.gz
#    (single-end: <SAMPLE>.fastq.gz; paired-end: <SAMPLE>_1/_2.fastq.gz)
#    Verify with md5sum where checksums are available.

# 5. Adjust the two conda setup lines in each .sbatch to your cluster, then
#    submit every job from ~/september so relative logs/ paths resolve:
cd ~/september
```

## Run order

```bash
cd ~/september
sbatch scripts/build_star_index.sbatch        # STAR index (once)
sbatch scripts/build_salmon_index.sbatch      # salmon index (once)
```

Mapping + quant run in chunks (chunk scripts are the supported path; the
per-sample arrays `mapx86.sbatch` / `quant.sbatch` are kept for reference):

```bash
# manual chunk, e.g. samples 1-25 (re-runnable: finished samples are skipped)
sbatch scripts/map_chunk.sbatch 1 25
sbatch scripts/quant_chunk.sbatch 1 50
```

or via the driver (persistent `screen` session; survives Ctrl-C and
re-runs — already-submitted chunks are tracked, finished samples skipped):

```bash
screen -S rnaseq
bash ~/september/scripts/submit_driver.sh
# detach with Ctrl-A D; reattach with: screen -r rnaseq
```

Merge + analysis (run on the login node; re-runnable, outputs are
atomically replaced only after verification):

```bash
python3 ~/september/scripts/merge_quant_v4.py   # transcript-level matrices
python3 ~/september/scripts/analyze_quant.py    # read-only summary report
# gene-level (needs salmon --geneMap / quant.genes.sf outputs):
python3 ~/september/scripts/merge_genes_v1.py
```

Results land in `~/september/results/`:
`tpm_matrix.tsv`, `numreads_matrix.tsv`
(`tpm_genes_matrix.tsv`, `numreads_genes_matrix.tsv` for gene-level).

## Sample list format

`list_sample.txt`: one sample per line, 4 tab-separated columns, no header:

```
<PROJ>  <SAMPLEID>  <TYPE>  <STRAND>
```

- `TYPE`: `single` | `paired`
- `STRAND`: `unstranded` | `fr-firststrand` | `fr-secondstrand`
  (documentary for STAR/salmon here; salmon runs with `-l A` auto-detect)
- Lines starting with `#` and blank lines are skipped.

Map chunks cover 25 samples, quant chunks 50. The driver keeps at most 4
of your jobs in the system (site submit limit is ~5, so pending
`(AssocMaxJobsLimit)` is normal). Jobs that hit the walltime are simply
re-queued by the next driver pass — check `logs/driver_submitted.txt`
before submitting anything by hand to avoid duplicates.
