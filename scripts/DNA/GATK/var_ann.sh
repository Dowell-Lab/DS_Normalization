#!/bin/bash 
#SBATCH --job-name=ann # Job name
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=allenma@colorado.edu # Where to send mail
#SBATCH --nodes=1 # Run on a single node
#SBATCH --ntasks=4     # Number of CPU (processer cores i.e. tasks) In this example I use 1. I only need one, since none of the commands I run are parallelized.
#SBATCH --mem=110gb # Memory limit
#SBATCH --time=20:00:00 # Time limit hrs:min:sec
#SBATCH --output=/scratch/Shares/dowell/down/temp_genomehg38/e_and_o/VariantAnnotator.%j.out # Standard output
#SBATCH --error=/scratch/Shares/dowell/down/temp_genomehg38/e_and_o/VariantAnnotator.%j.err # Standard error log

WKDIR=/scratch/Shares/dowell/down/temp_genomehg38/

mkdir -p ${WKDIR}/ann_Fam_vcfs/
mkdir -p /scratch/Users/allenma/tmp/${SLURM_JOBID}/

echo ${SLURM_JOBID}



java -Djava.io.tmpdir=/scratch/Users/allenma/tmp/${SLURM_JOBID} -XX:ParallelGCThreads=4 -XX:MaxHeapSize=100g \
        -jar /opt/gatk/4.1.0/gatk-package-4.1.0.0-local.jar IndexFeatureFile \
	-F ${WKDIR}/full_Fam_vcf/${rootname}

java -Djava.io.tmpdir=/scratch/Users/allenma/tmp/${SLURM_JOBID} -XX:ParallelGCThreads=4 -XX:MaxHeapSize=100g \
        -jar /opt/gatk/4.1.0/gatk-package-4.1.0.0-local.jar VariantAnnotator \
	-R $humanfasta \
	-variant ${WKDIR}/full_Fam_vcf/${rootname} \
	-O ${WKDIR}/ann_Fam_vcfs/${rootname} \
	--dbsnp $dbsnpnumvcf \
#	--resource:ann1000G $ann1000G \
#	-E ann1000G.InbreedingCoeff \
#	-E ann1000G.AF \
#	--resource-allele-concordance 

