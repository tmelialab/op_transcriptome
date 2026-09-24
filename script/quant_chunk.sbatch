#!/bin/bash
# NOTE: cluster-portable paths ($HOME/september layout, logs/ relative to submit dir).
# Before first submit: adjust the conda setup lines below to your cluster and
# submit from the project root (mkdir -p logs first) so relative log paths resolve.
#SBATCH -p short
#SBATCH --job-name=salmon_quant_chunk
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=23:00:00
#SBATCH --output=logs/quant_chunk_%j.out
#SBATCH --error=logs/quant_chunk_%j.err

# Salmon quant per-CHUNK serial (fallback limit AssocMaxSubmitJobLimit).
# Usage: sbatch quant_chunk.sbatch START END
# Idempoten: quant.sf yang sudah ada di-skip (aman resubmit).
set -e

START=${1:?isi START}; END=${2:?isi END}

source /mgpfs/apps/bioinformatics/apps/miniforge3/24.3.0-0/etc/profile.d/conda.sh
conda activate training_qc

HOME_DIR="$HOME/september/data"
SALMON_INDEX="$HOME/september/salmon/salmon_index"
QUANT_DIR="$HOME/september/quant"
GTF_FILE="$HOME/september/reference/GCF_000442705.2_EG11_genomic.gtf"
LIST="$HOME/september/list_sample.txt"

echo "chunk $START-$END | $(salmon --version 2>&1 | head -1) | start $(date)"

for N in $(seq "$START" "$END"); do
  LINE=$(sed -n "${N}p" "$LIST")
  if [[ -z "$LINE" || "$LINE" =~ ^# ]]; then echo "[$N] skip (kosong/komentar)"; continue; fi
  IFS=$'\t' read -r PROJ SAMPLEID TYPE STRAND <<< "$LINE"
  echo "[$N] $PROJ $SAMPLEID $TYPE $STRAND"

  OUT="${QUANT_DIR}/${PROJ}/${SAMPLEID}"
  if [[ -s "${OUT}/quant.sf" ]]; then echo "  quant.sf ada, skip"; continue; fi

  if [[ "$TYPE" == "single" ]]; then
    R1="${HOME_DIR}/rna/${PROJ}/${SAMPLEID}.fastq.gz"
    salmon quant -i "$SALMON_INDEX" -l A -r "$R1" -o "$OUT" -p 8 -g "$GTF_FILE"
  else
    R1="${HOME_DIR}/rna/${PROJ}/${SAMPLEID}_1.fastq.gz"
    R2="${HOME_DIR}/rna/${PROJ}/${SAMPLEID}_2.fastq.gz"
    salmon quant -i "$SALMON_INDEX" -l A -1 "$R1" -2 "$R2" -o "$OUT" -p 8 -g "$GTF_FILE"
  fi
done

echo "done chunk $START-$END; stop $(date)"