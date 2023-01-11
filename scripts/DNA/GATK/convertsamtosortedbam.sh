#!/bin/bash 
#SBATCH --job-name=convert # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=allenma@colorado.edu # Where to send mail
#SBATCH --nodes=1 # Run on a single node
#SBATCH --ntasks=1     # Number of CPU (processer cores i.e. tasks) In this example I use 1. I only need one, since none of the commands I run are parallelized.
#SBATCH --mem=30gb # Memory limit
#SBATCH --time=23:00:00 # Time limit hrs:min:sec
#SBATCH --output=/scratch/Shares/dowell/down/temp_genomehg38/e_and_o/convert.%j.out # Standard output
#SBATCH --error=/scratch/Shares/dowell/down/temp_genomehg38/e_and_o/convert.%j.err # Standard error log

WKDIR=/scratch/Shares/dowell/down/temp_genomehg38/
module load samtools/1.8

mkdir -p ${WKDIR}/unsorted_bams/
mkdir -p ${WKDIR}/sorted_bams/
mkdir -p ${WKDIR}/flagstats/

echo "samtools view"
    samtools view \
    -bS \
    -o ${WKDIR}/unsorted_bams/${rootname}.bam \
    ${WKDIR}/sams/${SAM}

echo "samtools sort"
    samtools sort -o ${WKDIR}/sorted_bams/${rootname}.bam ${WKDIR}/unsorted_bams/${rootname}.bam

echo "samtools index"
    samtools index ${WKDIR}/sorted_bams/${rootname}.bam

echo "samtools flagstat"
    samtools flagstat ${WKDIR}/sorted_bams/${rootname}.bam > ${WKDIR}/flagstats/${rootname}.bam.flagstat

