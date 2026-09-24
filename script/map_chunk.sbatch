#!/bin/bash
# NOTE: cluster-portable paths ($HOME/september layout, logs/ relative to submit dir).
# Before first submit: adjust the conda setup lines below to your cluster and
# submit from the project root (mkdir -p logs first) so relative log paths resolve.
#SBATCH -p short
#SBATCH --job-name=star_map_chunk
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=23:00:00
#SBATCH --output=logs/map_chunk_%j.out
#SBATCH --error=logs/map_chunk_%j.err

# Serial per-chunk STAR mapping (fallback for AssocMaxSubmitJobLimit:
# large arrays are rejected, 1 chunk job = 1 slot). Usage: sbatch map_chunk.sbatch START END
# Idempotent: existing BAMs are skipped (safe to resubmit).
set -e

START=${1:?provide START}; END=${2:?provide END}

source /mgpfs/apps/bioinformatics/apps/miniforge3/24.3.0-0/etc/profile.d/conda.sh
conda activate training_qc

export PARALLEL=8
export OMP_NUM_THREADS=${PARALLEL}

GTF_FILE="$HOME/september/reference/GCF_000442705.2_EG11_genomic.gtf"
HOME_DIR="$HOME/september/data"
GENOME="$HOME/september/star_index/EG11_STAR"
LIST="$HOME/september/list_sample.txt"

echo "chunk $START-$END | STAR $(STAR --version) | start $(date)"

for N in $(seq "$START" "$END"); do
  LINE=$(sed -n "${N}p" "$LIST")
  if [[ -z "$LINE" || "$LINE" =~ ^# ]]; then echo "[$N] skip (empty/comment line)"; continue; fi
  IFS=$'\t' read -r PROJ SAMPLEID TYPE STRAND <<< "$LINE"
  echo "[$N] $PROJ $SAMPLEID $TYPE $STRAND"

  OUTDIR="$HOME/september/map/${PROJ}/${SAMPLEID}"
  PREFIX="${OUTDIR}/${SAMPLEID}"
  BAM="${PREFIX}Aligned.sortedByCoord.out.bam"
  mkdir -p "$OUTDIR"
  cd "$OUTDIR"

  if [[ -s "$BAM" ]]; then echo "  BAM exists, skip"; continue; fi

  if [[ "$TYPE" == "single" ]]; then
    R1="${HOME_DIR}/rna/${PROJ}/${SAMPLEID}.fastq.gz"
    STAR --runMode alignReads --runThreadN "$PARALLEL" --genomeDir "$GENOME" \
         --readFilesIn "$R1" --outFileNamePrefix "$PREFIX" --sjdbGTFfile "$GTF_FILE" \
         --twopassMode Basic --alignIntronMax 500000 --outSAMtype BAM SortedByCoordinate \
         --outSAMattributes Standard --readFilesCommand zcat
  else
    R1="${HOME_DIR}/rna/${PROJ}/${SAMPLEID}_1.fastq.gz"
    R2="${HOME_DIR}/rna/${PROJ}/${SAMPLEID}_2.fastq.gz"
    STAR --runMode alignReads --runThreadN "$PARALLEL" --genomeDir "$GENOME" \
         --readFilesIn "$R1" "$R2" --outFileNamePrefix "$PREFIX" --sjdbGTFfile "$GTF_FILE" \
         --twopassMode Basic --alignIntronMax 500000 --outSAMtype BAM SortedByCoordinate \
         --outSAMattributes Standard --readFilesCommand zcat
  fi

  samtools index "$BAM"
  samtools view -F 0x904 -c "$BAM" > mapping_statistics.txt
  zcat "$R1" | wc -l > num_lines.txt
done

echo "done chunk $START-$END; stop $(date)"