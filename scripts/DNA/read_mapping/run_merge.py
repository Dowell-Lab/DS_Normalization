import run_hisat2
from run_hisat2 import indv_params
import os

if __name__=="__main__":
    peoplesam = {}
    for filepre, pu, sm, rg_id, sam, read_type in indv_params:
	if sm in peoplesam.keys():
            oldlist = peoplesam[sm]
            oldlist = oldlist +[sam]
            peoplesam[sm] = oldlist
        else:
            peoplesam[sm] = [sam]
    print (peoplesam)
    filestomerge= 
    os.system('sbatch --export=filestomerge="' + filestomerge + '",rootname="' + rootname + '" picardmergeandprep.sh')



