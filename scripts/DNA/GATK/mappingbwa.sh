#!/bin/bash 
#SBATCH --job-name=bwa # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=allenma@colorado.edu # Where to send mail
#SBATCH --nodes=1 # Run on a single node
#SBATCH --ntasks=32     # Number of CPU (processer cores i.e. tasks) In this example I use 1. I only need one, since none of the commands I run are parallelized.
#SBATCH --mem=15gb # Memory limit
#SBATCH --time=23:00:00 # Time limit hrs:min:sec
#SBATCH --output=/scratch/Shares/dowell/down/temp_genomehg38/e_and_o/bwa.%j.out # Standard output
#SBATCH --error=/scratch/Shares/dowell/down/temp_genomehg38/e_and_o/bwa.%j.err # Standard error log

INDIR=/Shares/down/genome/read_mapping_hg38/
WKDIR=/scratch/Shares/dowell/down/temp_genomehg38/

mkdir -p ${WKDIR}/sams/

module load bwa/0.7.15

###******NEED SPECIAL CASE FOR NON-PAIREDEND READS***********

fasta=/scratch/Shares/dowell/genomes/hg38/hg38.fa
#bwa index $fasta
echo $RGLINE

if [[ ${TYPE} == "PAIRED" ]]
then

    bwa mem \
        -t 32 \
	-R $RGLINE \
	$fasta \
        ${INDIR}/fastqs/${FASTQPRE}.end1.fq \
        ${INDIR}/fastqs/${FASTQPRE}.end2.fq \
        > ${WKDIR}/sams/${SAM}

elif [[ ${TYPE} == "UNPAIRED" ]]
then
      bwa mem \
        -t 32 \
        -R $RGLINE \
	$fasta \
        ${INDIR}/fastqs/${FASTQPRE}.fastq \
        > ${WKDIR}/sams/${SAM}    

fi

