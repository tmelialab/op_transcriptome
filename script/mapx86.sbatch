#!/bin/bash
# NOTE: cluster-portable paths ($HOME/september layout, logs/ relative to submit dir).
# Before first submit: adjust the conda setup lines below to your cluster and
# submit from the project root (mkdir -p logs first) so relative log paths resolve.
#SBATCH -p short
#SBATCH --job-name=star_map_EG11
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=23:00:00
#SBATCH --array=1-407%40
#SBATCH --output=logs/map_%A_%a.out
#SBATCH --error=logs/map_%A_%a.err

# STAR mapping -> BRIN (adaptasi mapx86.slurm dari BU SCC, jadi job array per-sampel)
# - serial 407 sampel tidak muat 24 jam -> 1 task = 1 sampel, maks 40 jalan bareng
# - FIX vs original: baca 4 kolom (PROJ SAMPLEID TYPE STRAND); versi original baca 3 kolom
#   sehingga TYPE="single\tunstranded" dan SEMUA sampel single-end salah masuk cabang paired
# - FIX vs original: prefix & path BAM dibetulkan (original: "${OUT}Aligned..." tanpa slash
#   + "-o /${OUT}..." jadi path absolut dari root)
# - idempoten: BAM yang sudah ada di-skip (aman resubmit)
set -e

: "${SLURM_ARRAY_TASK_ID:?harus jalan sebagai job array}"

source /mgpfs/apps/bioinformatics/apps/miniforge3/24.3.0-0/etc/profile.d/conda.sh
conda activate training_qc

export PARALLEL=8
export OMP_NUM_THREADS=${PARALLEL}

GTF_FILE="$HOME/september/reference/GCF_000442705.2_EG11_genomic.gtf"
HOME_DIR="$HOME/september/data"
GENOME="$HOME/september/star_index/EG11_STAR"
LIST="$HOME/september/list_sample.txt"

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")
if [[ -z "$LINE" || "$LINE" =~ ^# ]]; then
  echo "task ${SLURM_ARRAY_TASK_ID}: baris kosong/komentar, skip"
  exit 0
fi

IFS=$'\t' read -r PROJ SAMPLEID TYPE STRAND <<< "$LINE"
echo "task ${SLURM_ARRAY_TASK_ID}: $PROJ $SAMPLEID $TYPE $STRAND"

OUTDIR="$HOME/september/map/${PROJ}/${SAMPLEID}"
PREFIX="${OUTDIR}/${SAMPLEID}"
BAM="${PREFIX}Aligned.sortedByCoord.out.bam"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

if [[ -s "$BAM" ]]; then
  echo "BAM sudah ada, skip mapping: $BAM"
  exit 0
fi

if [[ "$TYPE" == "single" ]]; then
  R1="${HOME_DIR}/rna/${PROJ}/${SAMPLEID}.fastq.gz"
  # shellcheck disable=SC2086
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

# index + stats
samtools index "$BAM"
samtools view -F 0x904 -c "$BAM" > mapping_statistics.txt
zcat "$R1" | wc -l > num_lines.txt

echo "done task ${SLURM_ARRAY_TASK_ID}; stop $(date)"