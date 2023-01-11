import os
dbdir="/scratch/Shares/dowell/down/temp_genomehg38/tempdb/"
humanfasta="/scratch/Shares/dowell/genomes/hg38/hg38.fa"
allchrs = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', "chrM"]

if __name__=="__main__":
    for path_and_filename in os.listdir(dbdir):
	Fam=path_and_filename.split("/")[-1]
	infiles = ""
	for chromosome in allchrs:
		#/scratch/Shares/dowell/down/temp_genomehg38/Fam_vcf/FamA_chr1/FamA_chr1.vcf.gz
		infilename="/scratch/Shares/dowell/down/temp_genomehg38/Fam_vcf/"+Fam+"_"+chromosome+"/"+Fam+"_"+chromosome+".vcf.gz "
        	infiles = infiles+"-I "+infilename+" "
	print (Fam)
	#print (infiles)
	os.system('sbatch --export=Fam="' + Fam + '",infiles="' + infiles + '" Famfullvcf.sh')

