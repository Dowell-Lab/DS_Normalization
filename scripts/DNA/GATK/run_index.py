import os
humanfasta="/scratch/Shares/dowell/genomes/hg38/hg38.fa"
ftpdir="/scratch/Shares/dowell/down/temp_genomehg38/ftp.broadinstitute.org/bundle/hg38/"
dbsnpnumvcf=ftpdir+"dbsnp_146.hg38.vcf.gz"
indir="/scratch/Shares/dowell/down/temp_genomehg38/ann_Fam_vcfs/"
ann1000G="/scratch/Shares/dowell/down/temp_genomehg38/ftp.broadinstitute.org/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz"

if __name__=="__main__":
    for path_and_filename in os.listdir(indir):
	Fam=path_and_filename.split("/")[-1]
	infiles = ""
	print (Fam)
	if Fam.endswith("vcf"):
		os.system('sbatch --export=rootname="' + Fam + '",dbsnpnumvcf="' + dbsnpnumvcf + '",humanfasta="' + humanfasta + '",ann1000G="' + ann1000G + '" index_vcf.sh')


