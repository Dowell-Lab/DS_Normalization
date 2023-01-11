#!/bin/bash 
#SBATCH --job-name=hisat2 # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=allenma@colorado.edu # Where to send mail
#SBATCH --nodes=1 # Run on a single node
#SBATCH --ntasks=32     # Number of CPU (processer cores i.e. tasks) In this example I use 1. I only need one, since none of the commands I run are parallelized.
#SBATCH --mem=15gb # Memory limit
#SBATCH --time=23:00:00 # Time limit hrs:min:sec
#SBATCH --output=/scratch/Shares/dowell/down/temp_genomehg38/e_and_o/hisat2.%j.out # Standard output
#SBATCH --error=/scratch/Shares/dowell/down/temp_genomehg38/e_and_o/hisat2.%j.err # Standard error log


WKDIR=/scratch/Shares/dowell/down/temp_genomehg38/

mkdir -p ${WKDIR}/sams/

module load hisat2/2.1.0

###******NEED SPECIAL CASE FOR NON-PAIREDEND READS***********

if [[ ${TYPE} == "PAIRED" ]]
then

    hisat2 \
        -p 32 \
        -t \
        --rg-id ${RGID} \
        --rg PL:ILLUMINA \
        --rg PU:${PU} \
        --rg SM:${SM} \
        -x /scratch/Shares/dowell/genomes/hg38/HISAT2/genome \
        -1 ${WKDIR}/fastqs/${FASTQPRE}.end1.fq \
        -2 ${WKDIR}/fastqs/${FASTQPRE}.end2.fq \
        -S ${WKDIR}/sams/${SAM}

elif [[ ${TYPE} == "UNPAIRED" ]]
then
    hisat2 \
        -t \
        -p 32 \
        --rg-id ${RGID} \
        --rg PL:ILLUMINA \
        --rg PU:${PU} \
        --rg SM:${SM} \
        -x /scratch/Shares/dowell/genomes/hg38/HISAT2/genome \
        -U ${WKDIR}/fastqs/${FASTQPRE}.fastq \
        -S ${WKDIR}/sams/${SAM}

fi

