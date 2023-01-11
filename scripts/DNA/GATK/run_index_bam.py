import os
indir="/scratch/Shares/dowell/down/temp_genomehg38/base_recal_bams/"

if __name__=="__main__":
    for path_and_filename in os.listdir(indir):
	pathandfilename=path_and_filename.split("/")[-1]
	if pathandfilename.endswith("bam"):
		os.system('sbatch --export=pathandfilename="' + indir+path_and_filename + '" index_bam.sh')


