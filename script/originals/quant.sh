#!/bin/bash
#quantifying using salmon

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



#export PARALLEL=12             # Specify the number of thread
#export OMP_NUM_THREADS=${PARALLEL} # Environment variable setting



HOME_DIR="/home/ra022310/tisha"
SALMON_INDEX="/home/ra022310/l00039/genome/salmon_index"
#QUANT_DIR="/home/ra022310/l00039/quant"
#GTF_FILE="/home/ra022310/l00039/genome/EG5.1_Genes.V3.gtf"
QUANT_DIR="/home/ra022310/tisha/quant_stringtie"
GTF_FILE="/home/ra022310/tisha/stringtie/op_trans-merge.gtf"

# execute job
IFS=$'\t'
while read PROJ SAMPLEID TYPE STRAND; do

        if [[ "$PROJ" =~ \#.* ]];then
                echo "skipped comment line"
        else
        	echo $PROJ
        	echo $SAMPLEID
        	echo $TYPE

                OUT=${QUANT_DIR}/${PROJ}/${SAMPLEID}                


                if [[ "$TYPE" = "single" ]];then

                	R1=${HOME_DIR}/data/rna/${PROJ}/${SAMPLEID}.fastq.gz

                        salmon quant \
                        -i $SALMON_INDEX \
                        -l A \
                        -r $R1 \
                        -o $OUT \
                        -p 8 \
                        -g $GTF_FILE


                else
                	R1=${HOME_DIR}/data/rna/${PROJ}/${SAMPLEID}_1.fastq.gz
                	R2=${HOME_DIR}/data/rna/${PROJ}/${SAMPLEID}_2.fastq.gz

                        salmon quant \
                        -i $SALMON_INDEX \
                        -l A \
                        -1 $R1 \
                        -2 $R2 \
                        -o $OUT \
                        -p 8 \
                        -g $GTF_FILE
                fi
       fi       
done < /home/ra022310/l00039/script/list_sample.txt

echo "done running"
echo "=========================================================="


# print out some diagnostic stuff
end=$(date +%s%N)
duration=$(((end - start) / 1000000000 /60)) #in mins

echo "Stop time is $(date)"
echo "Duration: ${duration}"

