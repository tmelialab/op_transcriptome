#!/bin/bash
#SBATCH -p mem2      # Specifying queue
#SBATCH -n 12              # Specifying the number of CPUs (1-56)
#SBATCH --mem 27G  # Specifying memory usage (5-1500G)
#SBATCH -t 24:00:00   # Specifying elapsed time limits (1-24hours)


############################################################################
# SETTING VARS
############################################################################

# Now let's keep track of some information just in case anything goes wrong
start=$(date +%s%N)
echo "=========================================================="
echo "Starting on : $(date)"
echo "Running on node : $(hostname)"
echo "Current directory : $(pwd)"
echo "=========================================================="



#export PARALLEL=12                 # Specify the number of thread
#export OMP_NUM_THREADS=${PARALLEL} # Environment variable setting


GTF_FILE="/home/ra022310/l00039/genome/EG5.1_Genes.V3.gtf"
HOME_DIR="/home/ra022310/tisha"
GENOME="/home/ra022310/tisha/genome/EG5.1_STAR"

# execute job
IFS=$'\t'
while read PROJ SAMPLEID TYPE; do

        if [[ "$PROJ" =~ \#.* ]];then
                echo "skipped comment line"
        else
        		echo $PROJ
        		echo $SAMPLEID
        		echo $TYPE

                OUT=/home/ra022310/tisha/map/dna/${PROJ}/${SAMPLEID}/${SAMPLEID}   

                mkdir -p $OUT             
                cd $OUT

                if [[ "$TYPE" = "single" ]];then
                	    R1=${HOME_DIR}/data/rna/${PROJ}/${SAMPLEID}.fastq.gz

                        /home/ra022310/l00039/tools/STAR-2.7.11b/bin/Linux_x86_64_static/STAR \
                        --runMode alignReads \
                        --runThreadN $PARALLEL \
                        --genomeDir $GENOME \
                        --readFilesIn $R1 \
                        --outFileNamePrefix $OUT \
                        --sjdbGTFfile $GTF_FILE \
                        --twopassMode Basic \
                        --alignIntronMax 500000 \
                        --outSAMtype BAM SortedByCoordinate \
                        --outSAMattributes Standard \
                        --readFilesCommand zcat
                else
                	R1=${HOME_DIR}/data/rna/${PROJ}/${SAMPLEID}_1.fastq.gz
                	R2=${HOME_DIR}/data/rna/${PROJ}/${SAMPLEID}_2.fastq.gz

                        /home/ra022310/l00039/tools/STAR-2.7.11b/bin/Linux_x86_64_static/STAR \
                        --runMode alignReads \
                        --runThreadN $PARALLEL \
                        --genomeDir $GENOME \
                        --readFilesIn $R1 $R2 \
                        --outFileNamePrefix $OUT \
                        --sjdbGTFfile $GTF_FILE \
                        --twopassMode Basic \
                        --alignIntronMax 500000 \
                        --outSAMtype BAM SortedByCoordinate \
                        --outSAMattributes Standard \
                        --readFilesCommand zcat
                fi

                #Buat index file
                samtools index ${OUT}Aligned.sortedByCoord.out.bam -o /${OUT}Aligned.sortedByCoord.out.bam.bai
                #Simpan statistics mapping
                samtools view -F 0x904 -c ${OUT}Aligned.sortedByCoord.out.bam > mapping_statistics.txt
                zcat $R1 | wc -l > num_lines.txt

                #Buag bigwig file

       fi       
done < /home/ra022310/l00039/script/list_sample_mapx86_slurm.txt

echo "done running"
echo "=========================================================="


# print out some diagnostic stuff
end=$(date +%s%N)
duration=$(((end - start) / 1000000000 /60)) #in mins

echo "Stop time is $(date)"
echo "Duration: ${duration}"



# Usage: STAR  [options]... --genomeDir /path/to/genome/index/   --readFilesIn R1.fq R2.fq
# STAR version=2.7.10a
# ### Run Parameters
# runMode                         alignReads
#     string: type of the run.
#                                 alignReads             ... map reads
#                                 genomeGenerate         ... generate genome files
#                                 inputAlignmentsFromBAM ... input alignments from BAM. Presently only works with --outWigType and --bamRemoveDuplicates options.
#                                 liftOver               ... lift-over of GTF files (--sjdbGTFfile) between genome assemblies using chain file(s) from --genomeChainFiles.
#                                 soloCellFiltering  </path/to/raw/count/dir/>   </path/to/output/prefix>    ... STARsolo cell filtering ("calling") without remapping, followed by the path to raw count directory and output (filtered) prefix

# runThreadN                      1
#     int: number of threads to run STAR

# ### Genome Parameters
# genomeDir                   ./GenomeDir/
#     string: path to the directory where genome files are stored (for --runMode alignReads) or will be generated (for --runMode generateGenome)

# genomeLoad                NoSharedMemory
#     string: mode of shared memory usage for the genome files. Only used with --runMode alignReads.
#                           LoadAndKeep     ... load genome into shared and keep it in memory after run
#                           LoadAndRemove   ... load genome into shared but remove it after run
#                           LoadAndExit     ... load genome into shared memory and exit, keeping the genome in memory for future runs
#                           Remove          ... do not map anything, just remove loaded genome from memory
#                           NoSharedMemory  ... do not use shared memory, each job will have its own private copy of the genome

# genomeFastaFiles            -
#     string(s): path(s) to the fasta files with the genome sequences, separated by spaces. These files should be plain text FASTA files, they *cannot* be zipped.
#                             Required for the genome generation (--runMode genomeGenerate). Can also be used in the mapping (--runMode alignReads) to add extra (new) sequences to the genome (e.g. spike-ins).


    

# ### Splice Junctions Database
# sjdbGTFfile                             -
#     string: path to the GTF file with annotations

# ### Read Parameters
# readFilesIn                 Read1 Read2
#     string(s): paths to files that contain input read1 (and, if needed,  read2)
# 
#readFilesCommand             -
#    string(s): command line to execute for each of the input file. This command should generate FASTA or FASTQ text and send it to stdout
#               For example: zcat - to uncompress .gz files, bzcat - to uncompress .bz2 files, etc.

# ### Output: general
# outFileNamePrefix               ./
#     string: output files name prefix (including full or relative path). Can only be defined on the command line.

# ### Output: SAM and BAM
# outSAMtype                      SAM
#     strings: type of SAM/BAM output
#                                 1st word:
#                                 BAM  ... output BAM without sorting
#                                 SAM  ... output SAM without sorting
#                                 None ... no SAM/BAM output
#                                 2nd, 3rd:
#                                 Unsorted           ... standard unsorted
#                                 SortedByCoordinate ... sorted by coordinate. This option will allocate extra memory for sorting which can be specified by --limitBAMsortRAM.

# outSAMattributes                Standard
#     string: a string of desired SAM attributes, in the order desired for the output SAM. Tags can be listed in any combination/order.
#                                 ***Presets:
#                                 None        ... no attributes
#                                 Standard    ... NH HI AS nM
#                                 All         ... NH HI AS nM NM MD jM jI MC ch                                                                    
#                                 ***Alignment:
#                                 NH          ... number of loci the reads maps to: =1 for unique mappers, >1 for multimappers. Standard SAM tag.
#                                 HI          ... multiple alignment index, starts with --outSAMattrIHstart (=1 by default). Standard SAM tag.
#                                 AS          ... local alignment score, +1/-1 for matches/mismateches, score* penalties for indels and gaps. For PE reads, total score for two mates. Stadnard SAM tag.
#                                 nM          ... number of mismatches. For PE reads, sum over two mates.
#                                 NM          ... edit distance to the reference (number of mismatched + inserted + deleted bases) for each mate. Standard SAM tag.
#                                 MD          ... string encoding mismatched and deleted reference bases (see standard SAM specifications). Standard SAM tag.
#                                 jM          ... intron motifs for all junctions (i.e. N in CIGAR): 0: non-canonical; 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5: AT/AC, 6: GT/AT. If splice junctions database is used, and a junction is annotated, 20 is added to its motif value.
#                                 jI          ... start and end of introns for all junctions (1-based).
#                                 XS          ... alignment strand according to --outSAMstrandField.
#                                 MC          ... mate's CIGAR string. Standard SAM tag.
#                                 ch          ... marks all segment of all chimeric alingments for --chimOutType WithinBAM output.
#                                 cN          ... number of bases clipped from the read ends: 5' and 3'
#                                 ***Variation:
#                                 vA          ... variant allele
#                                 vG          ... genomic coordinate of the variant overlapped by the read.
#                                 vW          ... 1 - alignment passes WASP filtering; 2,3,4,5,6,7 - alignment does not pass WASP filtering. Requires --waspOutputMode SAMtag.
#                                 ***STARsolo:
#                                 CR CY UR UY ... sequences and quality scores of cell barcodes and UMIs for the solo* demultiplexing.
#                                 GX GN       ... gene ID and gene name for unique-gene reads.
#                                 gx gn       ... gene IDs and gene names for unique- and multi-gene reads.
#                                 CB UB       ... error-corrected cell barcodes and UMIs for solo* demultiplexing. Requires --outSAMtype BAM SortedByCoordinate.
#                                 sM          ... assessment of CB and UMI.
#                                 sS          ... sequence of the entire barcode (CB,UMI,adapter).
#                                 sQ          ... quality of the entire barcode.
#                                 ***Unsupported/undocumented:
#                                 ha          ... haplotype (1/2) when mapping to the diploid genome. Requires genome generated with --genomeTransformType Diploid .                               
#                                 rB          ... alignment block read/genomic coordinates.
#                                 vR          ... read coordinate of the variant.

# outSAMmapqUnique        255
#     int: 0 to 255: the MAPQ value for unique mappers


# ### Output Wiggle
# outWigType          None
#     string(s): type of signal output, e.g. "bedGraph" OR "bedGraph read1_5p". Requires sorted BAM: --outSAMtype BAM SortedByCoordinate .
#                     1st word:
#                     None       ... no signal output
#                     bedGraph   ... bedGraph format
#                     wiggle     ... wiggle format
#                     2nd word:
#                     read1_5p   ... signal from only 5' of the 1st read, useful for CAGE/RAMPAGE etc
#                     read2      ... signal from only 2nd read

# outWigNorm              RPM
#     string: type of normalization for the signal
#                         RPM    ... reads per million of mapped reads
#                         None   ... no normalization, "raw" counts


# ### 2-pass Mapping
# twopassMode                 None
#     string: 2-pass mapping mode.
#                             None        ... 1-pass mapping
#                             Basic       ... basic 2-pass mapping, with all 1st pass junctions inserted into the genome indices on the fly

