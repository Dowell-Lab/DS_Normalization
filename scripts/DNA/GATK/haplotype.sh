#!/bin/bash 
#SBATCH --job-name=haplo # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=allenma@colorado.edu # Where to send mail
#SBATCH --nodes=1 # Run on a single node
#SBATCH --ntasks=16     # Number of CPU (processer cores i.e. tasks) In this example I use 1. I only need one, since none of the commands I run are parallelized.
#SBATCH -p long
#SBATCH --mem=110gb # Memory limit
#SBATCH --time=200:00:00 # Time limit hrs:min:sec
#SBATCH --output=/scratch/Shares/dowell/down/temp_genomehg38/e_and_o/haplo.%j.out # Standard output
#SBATCH --error=/scratch/Shares/dowell/down/temp_genomehg38/e_and_o/haplo.%j.err # Standard error log

WKDIR=/scratch/Shares/dowell/down/temp_genomehg38/

mkdir -p ${WKDIR}/gvcfs/diploid/
mkdir -p ${WKDIR}/gvcfs/triploid/
mkdir -p /scratch/Users/allenma/tmp/${SLURM_JOBID}/

echo ${SLURM_JOBID}


java -Djava.io.tmpdir=/scratch/Users/allenma/tmp/${SLURM_JOBID} -XX:ParallelGCThreads=16 -XX:MaxHeapSize=100g\
	-jar /opt/gatk/4.1.0/gatk-package-4.1.0.0-local.jar HaplotypeCaller \
	-R $humanfasta \
	-I ${WKDIR}/base_recal_bams/${rootname}.bam \
	-O ${WKDIR}/gvcfs/diploid/${rootname}.g.vcf \
	-ERC GVCF --dbsnp $dbsnpnumvcf

