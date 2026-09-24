#!/bin/bash
# NOTE: cluster-portable paths ($HOME/september layout, logs/ relative to submit dir).
# Before first submit: adjust the conda setup lines below to your cluster and
# submit from the project root (mkdir -p logs first) so relative log paths resolve.
#SBATCH -p short
#SBATCH --job-name=salmon_index_EG11
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=23:00:00
#SBATCH --output=logs/salmon_index_%j.out
#SBATCH --error=logs/salmon_index_%j.err

# Salmon index for Elaeis guineensis EG11 (GCF_000442705.2)
set -e

source /mgpfs/apps/bioinformatics/apps/miniforge3/24.3.0-0/etc/profile.d/conda.sh
conda activate training_qc

GENOME=$HOME/september/reference/GCF_000442705.2_EG11_genomic.fna
GTF=$HOME/september/reference/GCF_000442705.2_EG11_genomic.gtf
SALDIR=$HOME/september/salmon

echo "gffread: $(which gffread)"
echo "salmon : $(salmon --version 2>&1 | head -1)"
echo "Starting: $(date)"

cd "$SALDIR"

# FIX: GTF NCBI menaruh transcript_id "" (kosong) di baris gene -> gffread >=0.12
# abort "no valid ID". Strip atribut kosong tsb (teruji OK); exon/CDS/transcript utuh.
FIXED_GTF="$SALDIR/EG11.fixed.gtf"
if [[ ! -s "$FIXED_GTF" ]]; then
  # 1) strip transcript_id "" (gffread>=0.12 abort "no valid ID")
  # 2) buang 1 baris strand "?" (rps12 kloroplas trans-splicing) -> gffread abort
  #    "Error parsing strand"; STAR index asli tidak terpengaruh (sudah COMPLETED)
  sed 's/transcript_id ""; //g' "$GTF" | awk -F'\t' '/^#/ || $7 != "?"' > "$FIXED_GTF"
  echo "fixed GTF: $(du -h "$FIXED_GTF" | cut -f1)"
fi

gffread "$FIXED_GTF" -g "$GENOME" -w EG11.transcript.fasta
grep "^>" "$GENOME" | cut -d " " -f 1 | sed 's/^>//' > decoys.txt
cat EG11.transcript.fasta "$GENOME" > EG11_gentrome.fa

echo "transcripts: $(grep -c '^>' EG11.transcript.fasta)"
echo "decoys     : $(wc -l < decoys.txt)"
echo "gentrome   : $(du -h EG11_gentrome.fa | cut -f1)"

salmon index -t EG11_gentrome.fa -d decoys.txt -p 48 -i salmon_index
echo "=== index contents ==="
ls -la "$SALDIR/salmon_index"
echo "done. Stop time is $(date)"