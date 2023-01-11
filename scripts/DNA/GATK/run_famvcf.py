import os
dbdir="/scratch/Shares/dowell/down/temp_genomehg38/tempdb/"
humanfasta="/scratch/Shares/dowell/genomes/hg38/hg38.fa"

if __name__=="__main__":
    for path_and_filename in os.listdir(dbdir):
	print (path_and_filename)
	Fam=path_and_filename.split("/")[-1]
	for subpath_and_filename in os.listdir(dbdir+Fam):
		chromosome=subpath_and_filename.split("/")[-1]
        	print (chromosome)
		os.system('sbatch --export=Fam="' + Fam + '",humanfasta="' + humanfasta +'",chromosome="' + chromosome + '" Famvcf.sh')


