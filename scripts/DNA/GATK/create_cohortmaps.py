import pandas

cohortsdir="/scratch/Shares/dowell/down/temp_genomehg38/scripts/cohorts/"
digvcfdir="/scratch/Shares/dowell/down/temp_genomehg38/gvcfs/diploid/"
trigvcfdir="/scratch/Shares/dowell/down/temp_genomehg38/gvcfs/triploid/"
pedfile="/scratch/Shares/dowell/down/temp_genomehg38/scripts/ped.txt"
allchrs = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', "chrM"]

def create_cohorts():
	df = pandas.read_table(pedfile, names=["fam", "person", "father", "mother", "gender", "affected"], sep="\t")	
	print(df["affected"].unique())
	fams = df["fam"].unique().tolist()
	for fam in fams:
		minidf = df[df["fam"]==fam]
		for chromosome in allchrs:
			if chromosome=="chr21":
				minidf1 = minidf[minidf["affected"]==2]
				minidf1["filename"] = trigvcfdir+minidf1["person"]+".g.vcf"
				minidf1 = minidf1[["person", "filename"]]
				minidf2 = minidf[minidf["affected"]!=2]
				minidf2["filename"] = digvcfdir+minidf2["person"]+".g.vcf"
				minidf2 = minidf2[["person", "filename"]]
				totalminidf = minidf1.append(minidf2, ignore_index=True)
				totalminidf.to_csv(cohortsdir+chromosome+"_"+fam+"_all.sample_map", sep="\t", header=False, index=False)
			else:
				minidf3 = minidf
				minidf3["filename"] = digvcfdir+minidf3["person"]+".g.vcf"
				minidf3 = minidf3[["person", "filename"]]
				minidf3.to_csv(cohortsdir+chromosome+"_"+fam+"_all.sample_map", sep="\t", header=False, index=False)
#	for chromsome in chromosomes:
#		for fam in fams.keys:
#			if chromosome!="chr21":
#				wf = open(mapdir+chromosome+"_"+fam+".sample_map", "w")
#				peopledata = 
				
			
create_cohorts()		
