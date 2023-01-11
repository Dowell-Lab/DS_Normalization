#!/bin/bash 
#SBATCH --job-name=picardmerge # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=allenma@colorado.edu # Where to send mail
#SBATCH --nodes=1 # Run on a single node
#SBATCH --ntasks=16     # Number of CPU (processer cores i.e. tasks) In this example I use 1. I only need one, since none of the commands I run are parallelized.
#SBATCH --mem=20gb # Memory limit
#SBATCH --time=20:00:00 # Time limit hrs:min:sec
#SBATCH --output=/scratch/Shares/dowell/down/temp_genomehg38/e_and_o/picardmerge.%j.out # Standard output
#SBATCH --error=/scratch/Shares/dowell/down/temp_genomehg38/e_and_o/picardmerge.%j.err # Standard error log

WKDIR=/scratch/Shares/dowell/down/temp_genomehg38/

mkdir -p ${WKDIR}/merged_files/
mkdir -p ${WKDIR}/merged_files_validation/
mkdir -p ${WKDIR}/merged_files_markdups/

echo ${SLURM_JOBID}
echo $filestomerge
echo $rootname

#java -Djava.io.tmpdir=/scratch/Users/allenma/tmp/${SLURM_JOBID} -XX:ParallelGCThreads=16 -jar /opt/picard/2.6.0/picard-2.6.0.jar MergeSamFiles ${filestomerge} O=${WKDIR}/merged_files/${rootname}.bam SORT_ORDER=coordinate
#java -Djava.io.tmpdir=/scratch/Users/allenma/tmp/${SLURM_JOBID} -XX:ParallelGCThreads=16 -jar /opt/picard/2.6.0/picard-2.6.0.jar ValidateSamFile I=${WKDIR}/merged_files/${rootname}.bam O=${WKDIR}/merged_files_validation/${rootname}.bam.val  MODE=SUMMARY
java -Djava.io.tmpdir=/scratch/Users/allenma/tmp/${SLURM_JOBID} -XX:ParallelGCThreads=16 -jar /opt/picard/2.6.0/picard-2.6.0.jar MarkDuplicates I=${WKDIR}/merged_files/${rootname}.bam O=${WKDIR}/merged_files_markdups/${rootname}.bam M=${WKDIR}/merged_files_markdups/${rootname}.dupmetrics.txt 

