import samples
from samples import indv_params
import os

if __name__=="__main__":
    peoplesam = {}
    for filepre, pu, sm, rg_id, sam, read_type in indv_params:
    	rootname= sam[0:-4]
	print (rootname)
	os.system('sbatch --export=SAM="' + sam + '",rootname="' + rootname + '" convertsamtosortedbam.sh')



