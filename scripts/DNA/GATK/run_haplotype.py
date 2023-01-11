import run_hisat2
from run_hisat2 import indv_params
import os

humanfasta="/scratch/Shares/dowell/genomes/hg38/hg38.fa"
#/scratch/Shares/dowell/down/temp_genomehg38/ftp.broadinstitute.org/bundle/hg38/
ftpdir="/scratch/Shares/dowell/down/temp_genomehg38/ftp.broadinstitute.org/bundle/hg38/"
dbsnpnumvcf=ftpdir+"dbsnp_146.hg38.vcf.gz"
goldstandardindelsvcf=ftpdir+"Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
indir="/scratch/Shares/dowell/down/temp_genomehg38/base_recal_bams/"

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
	os.system('sbatch --export=filestomerge="' + filestomerge + '",rootname="' + rootname + '",humanfasta="' + humanfasta +'",dbsnpnumvcf="' + dbsnpnumvcf +'",goldstandardindelsvcf="' + goldstandardindelsvcf +'" haplotype.sh')
	os.system('sbatch --export=filestomerge="' + filestomerge + '",rootname="' + rootname + '",humanfasta="' + humanfasta +'",dbsnpnumvcf="' + dbsnpnumvcf +'",goldstandardindelsvcf="' + goldstandardindelsvcf +'" haplotype_tri.sh')

if __name__=="__main__":
	main()



