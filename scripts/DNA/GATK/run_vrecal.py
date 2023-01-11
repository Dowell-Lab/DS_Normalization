import os
dbdir="/scratch/Shares/dowell/down/temp_genomehg38/tempdb/"
humanfasta="/scratch/Shares/dowell/genomes/hg38/hg38.fa"

if __name__=="__main__":
    for path_and_filename in os.listdir(dbdir):
	print (path_and_filename)
	Fam=path_and_filename.split("/")[-1]
	print (Fam)
	for path_and_filename in os.listdir(dbdir+Fam+"/"):
		chromosome=path_and_filename.split("/")[-1]
		print (chromosome) 
        	#Fam humanfasta
		#os.system('sbatch --export=Fam="' + Fam + '",humanfasta="' + humanfasta + '",chromosome="' + chromosome + '" Famvcf.sh')


#not done
#last time i did this I uesd the below two script... but I had to vary the 
#/Shares/down/genome/GATKhg19/step7_VQSR/diploidvcfsSNPscripts/FamA.varrecal.sh
#/Shares/down/genome/GATKhg19/step7_VQSR/diploidvcfsSNP100.0_vcfsINDEL99.9applytovcfscripts/FamNA.applyvarrecall.sh 
#the first script takes in a vcf and a list of tranches and a mode (INDEL or SNP)
#and creates a tranche file and a rcal file
#I started with 5 family vcf files
#this makes 10 recal and tranch files (one for each mode)
#the second script takes in the vsf file, rcal file, tranche file and a  fileter tranche level and a mode
# and outputs a vcf
#so I made a crap tone of shell scripts, a crap ton of vcf files!


#!/bin/sh
#SBATCH --job-name=FamA  # Job name
#SBATCH --nodes=1
#SBATCH --ntasks=4# Number of CPU (processer cores i.e. tasks) In this example I use 32 for bowtie2. I only need one, since none of the commands I run are parallelized.
#SBATCH --mem=48Gb
#SBATCH --time=4:00:00
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Mary.A.Allen@colorado.edu
#SBATCH -p short



#SBATCH --error=/scratch/Shares/dowell/Down/genome/GATKhg19/step7_VQSR/diploidvcfsSNPscripts/FamA.err
#SBATCH --output=/scratch/Shares/dowell/Down/genome/GATKhg19/step7_VQSR/diploidvcfsSNPscripts/FamA.out



#java -Djava.io.tmpdir=/tmp/${PBS_JOBID}/ -XX:ParallelGCThreads=4 -jar /opt/gatk/3.3-0/GenomeAnalysisTK.jar -T VariantRecalibrator -R /scratch/Shares/dowell/Down/genome/GATKhg19/ftp.broadinstitute.org/bundle/2.8/hg19/ucsc.hg19.fasta -input /scratch/Shares/dowell/Down/genome/GATKhg19/step5_Gentotypegvcf/diploidvcfs/FamA.vcf -recalFile /scratch/Shares/dowell/Down/genome/GATKhg19/step7_VQSR/diploidvcfsSNPrecal/FamA.recal -tranchesFile /scratch/Shares/dowell/Down/genome/GATKhg19/step7_VQSR/diploidvcfsSNPtranche/FamA.tranches -nt 4 -resource:hapmap,known=false,training=true,truth=true,prior=15.0 /scratch/Shares/dowell/Down/genome/GATKhg19/ftp.broadinstitute.org/bundle/2.8/hg19/hapmap_3.3.hg19.vcf -resource:omni,known=false,training=true,truth=true,prior=12.0 /scratch/Shares/dowell/Down/genome/GATKhg19/ftp.broadinstitute.org/bundle/2.8/hg19/1000G_omni2.5.hg19.vcf -resource:1000G,known=false,training=true,truth=false,prior=10.0 /scratch/Shares/dowell/Down/genome/GATKhg19/ftp.broadinstitute.org/bundle/2.8/hg19/1000G_phase1.snps.high_confidence.hg19.vcf -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /scratch/Shares/dowell/Down/genome/GATKhg19/ftp.broadinstitute.org/bundle/2.8/hg19/dbsnp_138.hg19.vcf -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0

#!/bin/sh
#SBATCH --job-name=FamA  # Job name
#SBATCH --nodes=1
#SBATCH --ntasks=4# Number of CPU (processer cores i.e. tasks) In this example I use 32 for bowtie2. I only need one, since none of the commands I run are parallelized.
#SBATCH --mem=48Gb
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Mary.A.Allen@colorado.edu
#SBATCH -p long



#SBATCH --error=/scratch/Shares/dowell/Down/genome/GATKhg19/step7_VQSR/diploidvcfsSNP99.0_vcfsINDEL99.9applytovcfscripts/FamA.err
#SBATCH --output=/scratch/Shares/dowell/Down/genome/GATKhg19/step7_VQSR/diploidvcfsSNP99.0_vcfsINDEL99.9applytovcfscripts/FamA.out



#java -Djava.io.tmpdir=/tmp/${SLURM__JOBID}/ -XX:ParallelGCThreads=4 -jar /opt/gatk/3.3-0/GenomeAnalysisTK.jar -T ApplyRecalibration -R /scratch/Shares/dowell/Down/genome/GATKhg19/ftp.broadinstitute.org/bundle/2.8/hg19/ucsc.hg19.fasta -input /scratch/Shares/dowell/Down/genome/GATKhg19/step7_VQSR/diploidvcfsSNP99.0_vcfs/FamA.vcf -recalFile /scratch/Shares/dowell/Down/genome/GATKhg19/step7_VQSR/diploidvcfsSNP99.0_vcfsINDELrecal/FamA.recal -tranchesFile /scratch/Shares/dowell/Down/genome/GATKhg19/step7_VQSR/diploidvcfsSNP99.0_vcfsINDELtranche/FamA.tranches -nt 4 -o /scratch/Shares/dowell/Down/genome/GATKhg19/step7_VQSR/diploidvcfsSNP99.0_vcfsINDEL99.9_vcfs/FamA.vcf --ts_filter_level 99.9 -mode INDEL

