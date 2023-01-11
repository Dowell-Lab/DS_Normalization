import os
indir="/scratch/Shares/dowell/down/temp_genomehg38/scripts/cohorts/"

if __name__=="__main__":
    for path_and_filename in os.listdir(indir):
        print (path_and_filename)
	rootname= path_and_filename.split("/")[-1]
	chromosome, Fam, typeofgvcf, justmap = rootname.split("_")
	print (chromosome, Fam) 
        #chromosome sample_map Fam
	os.system('sbatch --export=Fam="' + Fam + '",sample_map="' + indir+path_and_filename + '",chromosome="' + chromosome + '" joinhap.sh')


