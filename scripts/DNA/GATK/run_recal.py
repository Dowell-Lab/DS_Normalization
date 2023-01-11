import run_hisat2
from run_hisat2 import indv_params
import os

humanfasta="/scratch/Shares/dowell/genomes/hg38/hg38.fa"
ftpdir="/scratch/Shares/dowell/down/temp_genomehg38/ftp.broadinstitute.org/bundle/hg38/"
dbsnpnumvcf=ftpdir+"dbsnp_146.hg38.vcf.gz"
goldstandardindelsvcf=ftpdir+"Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

rerun=["Eli", "Elizabeth", "Douglas"]

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
#    for person in peoplesam.keys():
    for person in rerun:
	rootname=person
	os.system('sbatch --export=rootname="' + rootname + '",humanfasta="' + humanfasta +'",dbsnpnumvcf="' + dbsnpnumvcf +'",goldstandardindelsvcf="' + goldstandardindelsvcf +'" baserecal.sh')

if __name__=="__main__":
	main()



