#!/bin/bash
# NOTE: cluster-portable paths ($HOME/september layout, logs/ relative to submit dir).
# Before first submit: adjust the conda setup lines below to your cluster and
# submit from the project root (mkdir -p logs first) so relative log paths resolve.
#SBATCH -p short
#SBATCH --job-name=star_index_EG11
#SBATCH -N 1
#SBATCH --cpus-per-task=48
#SBATCH --mem=240G
#SBATCH --time=23:00:00
#SBATCH --output=logs/star_index_%j.out
#SBATCH --error=logs/star_index_%j.err

# STAR genome index build for Elaeis guineensis EG11 (GCF_000442705.2)
# Port of originals/genome_index.pjsub (Fugaku #PJM + spack) to Slurm + conda
set -e

source /mgpfs/apps/bioinformatics/apps/miniforge3/24.3.0-0/etc/profile.d/conda.sh
conda activate training_qc

GENOME=$HOME/september/reference/GCF_000442705.2_EG11_genomic.fna
GTF=$HOME/september/reference/GCF_000442705.2_EG11_genomic.gtf
INDEX=$HOME/september/star_index/EG11_STAR

echo "=========================================================="
echo "STAR version: $(STAR --version)"
echo "Starting   : $(date)"
echo "Node       : $(hostname)"
echo "=========================================================="

mkdir -p "$INDEX"
STAR --runThreadN 48 \
     --runMode genomeGenerate \
     --genomeDir "$INDEX" \
     --genomeFastaFiles "$GENOME" \
     --sjdbGTFfile "$GTF" \
     --sjdbOverhang 100 \
     --genomeSAindexNbases 14 \
     --limitGenomeGenerateRAM 100000000000

echo "=== index contents ==="
ls -la "$INDEX"
echo "done. Stop time is $(date)"