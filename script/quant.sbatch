#!/bin/bash
# NOTE: cluster-portable paths ($HOME/september layout, logs/ relative to submit dir).
# Before first submit: adjust the conda setup lines below to your cluster and
# submit from the project root (mkdir -p logs first) so relative log paths resolve.
#SBATCH -p short
#SBATCH --job-name=salmon_quant_EG11
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=23:00:00
#SBATCH --array=1-407%60
#SBATCH --output=logs/quant_%A_%a.out
#SBATCH --error=logs/quant_%A_%a.err

# Salmon quantification (adapted from originals/quant.sh into a per-sample job array)
# - FIX vs original: parse 4 columns (PROJ SAMPLEID TYPE STRAND)
# - libtype "-l A" (auto-detect, as in the original); STRAND column is documentary
# - idempotent: existing quant.sf files are skipped (safe to resubmit)
set -e

: "${SLURM_ARRAY_TASK_ID:?must run as a job array}"

source /mgpfs/apps/bioinformatics/apps/miniforge3/24.3.0-0/etc/profile.d/conda.sh
conda activate training_qc

HOME_DIR="$HOME/september/data"
SALMON_INDEX="$HOME/september/salmon/salmon_index"
QUANT_DIR="$HOME/september/quant"
GTF_FILE="$HOME/september/reference/GCF_000442705.2_EG11_genomic.gtf"
LIST="$HOME/september/list_sample.txt"

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")
if [[ -z "$LINE" || "$LINE" =~ ^# ]]; then
  echo "task ${SLURM_ARRAY_TASK_ID}: empty/comment line, skip"
  exit 0
fi

IFS=$'\t' read -r PROJ SAMPLEID TYPE STRAND <<< "$LINE"
echo "task ${SLURM_ARRAY_TASK_ID}: $PROJ $SAMPLEID $TYPE $STRAND (salmon $(salmon --version 2>&1 | head -1))"

OUT="${QUANT_DIR}/${PROJ}/${SAMPLEID}"
if [[ -s "${OUT}/quant.sf" ]]; then
  echo "quant.sf already exists, skip: $OUT"
  exit 0
fi

if [[ "$TYPE" == "single" ]]; then
  R1="${HOME_DIR}/rna/${PROJ}/${SAMPLEID}.fastq.gz"
  salmon quant -i "$SALMON_INDEX" -l A -r "$R1" -o "$OUT" -p 8 -g "$GTF_FILE"
else
  R1="${HOME_DIR}/rna/${PROJ}/${SAMPLEID}_1.fastq.gz"
  R2="${HOME_DIR}/rna/${PROJ}/${SAMPLEID}_2.fastq.gz"
  salmon quant -i "$SALMON_INDEX" -l A -1 "$R1" -2 "$R2" -o "$OUT" -p 8 -g "$GTF_FILE"
fi

echo "done task ${SLURM_ARRAY_TASK_ID}; stop $(date)"