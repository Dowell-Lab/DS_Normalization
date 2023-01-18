#!/bin/bash 
#SBATCH --job-name=baserecal # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=allenma@colorado.edu # Where to send mail
#SBATCH --nodes=1 # Run on a single node
#SBATCH --ntasks=16     # Number of CPU (processer cores i.e. tasks) In this example I use 1. I only need one, since none of the commands I run are parallelized.
#SBATCH -p long
#SBATCH --mem=20gb # Memory limit
#SBATCH --time=80:00:00 # Time limit hrs:min:sec
#SBATCH --output=/scratch/Shares/dowell/down/temp_genomehg38/e_and_o/baserecal.%j.out # Standard output
#SBATCH --error=/scratch/Shares/dowell/down/temp_genomehg38/e_and_o/baserecal.%j.err # Standard error log

WKDIR=/scratch/Shares/dowell/down/temp_genomehg38/

mkdir -p ${WKDIR}/baserecalltable/
mkdir -p ${WKDIR}/base_recal_bams/
mkdir -p ${WKDIR}/base_recal_files/
mkdir -p /scratch/Users/allenma/tmp/${SLURM_JOBID}/

echo ${SLURM_JOBID}
module load R/3.6.1

java -Djava.io.tmpdir:/scratch/Users/allenma/tmp/${SLURM_JOBID} -XX:ParallelGCThreads=16 \
        -jar /opt/gatk/4.1.0/gatk-package-4.1.0.0-local.jar CreateSequenceDictionary \
        -R $humanfasta

java -Djava.io.tmpdir:/scratch/Users/allenma/tmp/${SLURM_JOBID} -XX:ParallelGCThreads=16 \
	-jar /opt/gatk/4.1.0/gatk-package-4.1.0.0-local.jar BaseRecalibrator \
	-R $humanfasta \
	-I ${WKDIR}/merged_files_markdups/${rootname}.bam \
	--known-sites ${dbsnpnumvcf} \
	--known-sites $goldstandardindelsvcf \
	-O ${WKDIR}/baserecalltable/${rootname}.recal.table

java -Djava.io.tmpdir:/scratch/Users/allenma/tmp/${SLURM_JOBID} -XX:ParallelGCThreads=16 \
	-jar /opt/gatk/4.1.0/gatk-package-4.1.0.0-local.jar ApplyBQSR \
	-R $humanfasta \
	-I ${WKDIR}/merged_files_markdups/${rootname}.bam \
	-bqsr ${WKDIR}/baserecalltable/${rootname}.recal.table \
	-O ${WKDIR}/base_recal_bams/${rootname}.bam

java -Djava.io.tmpdir=/scratch/Users/allenma/tmp/${SLURM_JOBID} -XX:ParallelGCThreads=16 \
	-jar /opt/gatk/4.1.0/gatk-package-4.1.0.0-local.jar AnalyzeCovariates \
	-bqsr ${WKDIR}/baserecalltable/${rootname}.recal.table \
	-plots ${WKDIR}/base_recal_files/${rootname}.AnalyzeCovariates.pdf