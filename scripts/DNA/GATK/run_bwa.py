import os
import samples
from samples import *

if __name__=="__main__":
    for fastq_pre, pu, sm, rg_id, sam, read_type in indv_params:
	rgline=r'@RG\t'+r'ID:'+rg_id+r'\tSM:'+sm+r'\tPL:ILLUMINA'+r'\tPU:'+pu
	#@RG\tID:Seq01p\tSM:Seq01\tPL:ILLUMINA\tPI:330
        os.system('sbatch --export=FASTQPRE="' + fastq_pre + '",SAM="' + sam + '",RGLINE="' + rgline + 
              '",TYPE="' + read_type + '" mappingbwa.sh')

