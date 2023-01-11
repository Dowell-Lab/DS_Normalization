import run_hisat2
from run_hisat2 import indv_params
import os

indir="/scratch/Shares/dowell/down/temp_genomehg38/sorted_bams/"
outdir1="/scratch/Shares/dowell/down/temp_genomehg38/merged_bams/"
outdir2="/scratch/Shares/dowell/down/temp_genomehg38/validate_bams/"
outdir3="/scratch/Shares/dowell/down/temp_genomehg38/merge_markdup_bams/"


def main():
    peoplesam = {}
    for filepre, pu, sm, rg_id, sam, read_type in indv_params:
        bam = sam.replace(".s", ".b")
        if sm in peoplesam.keys():
            oldlist = peoplesam[sm]
            oldlist = oldlist +[bam]
            peoplesam[sm] = oldlist
        else:
            peoplesam[sm] = [bam]
    for person in peoplesam.keys():
        filestomerge= [indir+bamfile for bamfile in peoplesam[person]]
	filestomerge = " I=".join(filestomerge)
	filestomerge = " I="+filestomerge
	rootname=person
	os.system('sbatch --export=filestomerge="' + filestomerge + '",rootname="' + rootname + '" picardmergeandprep.sh')


if __name__=="__main__":
	main()



