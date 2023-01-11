# Line plot of MFC for Sample Depth/Composition
library("reshape2")
library("plyr")
library("ggplot2")
library("miscTools")
library("stringr")
library("DESeq2")
library("tidyverse")
theme_set(
  theme_classic(base_size = 18)
)

#pointers to original files Mary made
RNAindir<-"/Shares/down/mixed/RNAGRO01132021/data/RNAseqfiles/"
GROseq_indir<-"/Shares/down/mixed/RNAGRO01132021/data/GROfiles/"
RNAcoveragedat<-"res_featureCounts_gene_idfull_143138.coverage.csv"

#Below is the bed that got used for RNA-seq. 
#This is not the bed file that got used for GRO-seq exactly becuase to be acuurate in the GROS-seq I had to remove any gene bodys/tss that were in the data more than once. 
ori_worldbed <- "/scratch/Shares/dowell/genomes/hg38/hg38_refseq.bed"
#Below are the bed file I used for GRO-seq anaysis
#genes.bed and tss.bed came from /Users/allenma/humangtf.ipynb
#in breif i keep the longest isoform with uniq body start and stop, body(+1000, -500) start and stop, tss(-500, +500) start and stop
#I also removed gene <3000bp and genes whose body or tss <0 in coordinates
#below I will only keep genes that were anayized by both GR0-seq and RNA-seq
beddir<-"/Shares/down/mixed/RNAGRO01132021/data/bedfiles/"
GROann <- paste(beddir,"masterannotation.bedlike",sep="")
RNAann <-paste(RNAindir, "res_featureCounts_gene_idfull_143138.annotation.csv", sep="")

#read the annotation tables
annotationori <-read.table(ori_worldbed, sep="\t", col.names=c("chr", "start", "stop", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts"))
annotationGRO <-read.table(GROann, header=TRUE, sep="\t")
row.names(annotationGRO)<-annotationGRO$name
annotationRNA <-read.csv(RNAann)
row.names(annotationRNA)<-annotationRNA$GeneID

#keep only the genes in both GRO and RNA annotation files
annotationmerge <-merge(annotationRNA, annotationGRO, all = TRUE,by="row.names")
row.names(annotationmerge)<-annotationmerge$GeneID

#keep only the genes on the main chromosomes
minichrs <-c("chr1", "chr2", "chr3","chr4", "chr5", "chr6", "chr7", "chr8", "chr9","chr10","chr11", "chr12", "chr13","chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")
detach("package:plyr")
library(dplyr)
annotationmerge <-annotationmerge %>%
  arrange(GeneID) %>%
  filter(chr %in% minichrs)

#add a column that has TRUE if the gene is not on chr21 for use as a control gene
chr21genes <- annotationmerge[annotationmerge$chr=="chr21",]
annotationmerge$controlgenes <- TRUE
annotationmerge$controlgenes18 <- annotationmerge$chr!="chr18"
annotationmerge$controlgenes19 <- annotationmerge$chr!="chr19"
annotationmerge$controlgenes20 <- annotationmerge$chr!="chr20"
annotationmerge$controlgenes21 <- annotationmerge$chr!="chr21"
annotationmerge$controlgenes21shuffle <- sample(annotationmerge$chr!="chr21")


annotationmergeno21 <- annotationmerge %>% filter(controlgenes21)
annotationmergeno21shuffle <- annotationmerge %>% filter(controlgenes21shuffle)
lesschrs <- c("chr21","chr22")

###Initial DESeq2 Run with NO adjustments####
####UNCORRECTED SIM COMPARISONS####
# I'm purposely picking a high dispersion/noise file because we can see just how well our methods work to correct it
# Note that since this is simulated data, my gut feeling is that the corrections won't be quite as good as with real data
RNAcountdat <- read.csv("/Users/sahu0957/trisomy_normalization/RNAGRO01132021/data/RNAseqfiles/BIOREPEATreplicatefinalsims_a0.03_p4_r1_s1_n3.csv", 
                        sep=",", row.names=1)

#RNAcountdat$Gene_id <- rownames(RNAcountdat)
library(dplyr)
RNAcountdat<- RNAcountdat %>%
  filter(Gene_id %in% annotationmerge$name) %>%
  arrange(Gene_id) %>%
  column_to_rownames('Gene_id')

# Make "fake" metadata in the same structure as the real data. Just so I can recycle code...
names <- c("filename","tech_rep","biological_rep","Family","Person","ploidy","samplegroup")
RNAmetadata <- read.table(text = "",
                          col.names = names,
                          stringsAsFactors = FALSE)
RNAmetadata <- as.data.frame(RNAmetadata,stringsAsFactors=FALSE)
colnames(RNAmetadata) <- names
for (i in 1:(ncol(RNAcountdat)/4)){
  newmetadata1 <- c(paste0("Elvis_rep",i,"_counts"),"0",as.character(i),"E","Elvis","ploidy_typical",paste0("Elvis_",i))
  RNAmetadata[nrow(RNAmetadata)+1,] <- newmetadata1
}

for (i in 1:(ncol(RNAcountdat)/4)){
  newmetadata1 <- c(paste0("EMom_rep",i,"_counts"),"0",as.character(i),"E","EMom","ploidy_typical",paste0("EMom_",i))
  RNAmetadata[nrow(RNAmetadata)+1,] <- newmetadata1
}

for (i in 1:(ncol(RNAcountdat)/4)){
  newmetadata1 <- c(paste0("EDad_rep",i,"_counts"),"0",as.character(i),"E","EDad","ploidy_typical",paste0("EDad_",i))
  RNAmetadata[nrow(RNAmetadata)+1,] <- newmetadata1
}

for (i in 1:(ncol(RNAcountdat)/4)){
  newmetadata2 <- c(paste0("Eddie_rep",i,"_counts"),"0",as.character(i),"E","Eddie","ploidy_trisomy",paste0("Eddie_",i))
  RNAmetadata[nrow(RNAmetadata)+1,] <- newmetadata2
}

# Don't do an interaction term with sims, since they are generated independent of batch info... 
# We could always make sims from each replicate instead of sampling a normal distribution around their means
# This would mirror biological replication. 
RNAdds <- DESeqDataSetFromMatrix(countData = RNAcountdat, colData = RNAmetadata, design = ~Person)
RNAdds <-DESeq(RNAdds, betaPrior = FALSE) # Earlier DESeq2 versions include a betaPrior for the intercept term... wild!

person1 = "Eddie"
person2 = "Elvis"
RNAddsres<-results(RNAdds,  contrast=c("Person", person1, person2))
resdata <- as.data.frame(RNAddsres)
fullresdata <- merge(resdata,annotationmerge,by.x=0,by.y=3)
fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]
medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]

detach("package:dplyr")
library(plyr)
medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]
library(plyr)
medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians
common_fullresdata_rna <- merge(medresdata,common_ids,by.x=1,by.y=1)

signif_medresdata <- medresdata[medresdata$padj<.01,]
notsignif_medresdata <- medresdata[medresdata$padj>=.01,]
ggplot() + 
  geom_violin(data=fullresdata,trim=TRUE,aes(x=chr, y=log2FoldChange)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  #  geom_hline(yintercept=log2(0.66667),linetype="dashed",color="red") +
  geom_hline(yintercept=0) +
  ylab(paste0("Log2FC")) +
  geom_hline(yintercept = log2(1.5), linetype="dotted", color="blue" ) +
  xlab("Chromosome") +
  theme_classic(base_size = 24) +
  geom_jitter(data = signif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=1.5, alpha=0.85,color="red",fill="red") +
  geom_jitter(data = notsignif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=1.5, alpha=0.85,color="black",fill="black") +
  stat_summary(data=fullresdata,aes(x=chr, y=log2FoldChange),fun.y=median, geom="crossbar", size=0.25, color="orange") +
  theme(legend.title = element_blank())

ggsave("rnaseq_deseq2_realsims_new_violins_asympt03_extrapois4_med140_med100.svg",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "svg")
ggsave("rnaseq_deseq2_realsims_new_violins_asympt03_extrapois4_med140_med100.png",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "png")

###Flat run done #####
###looping runs with NO correction:RNA-seq#### 
# We can see whether there's a pattern associated with more/less dispersion by looping DESeq2 runs through all of our simulated files
# each count file is in the same format, where the sim info is in the name of the file and the counts in the file itself
# so we can parse out the parameters of the simulation first and then run DESeq2
# a = asymptotic dispersion
# p = extra-poisson noise
# n = number of replicates in a single experiment
# r = number of experiments
# s = scalar for depth
countfiles <- list.files(path = "/Users/sahu0957/trisomy_normalization/RNAGRO01132021/data/RNAseqfiles/",pattern = "AdjustAll")
countfiles
median_colnames <- c("chr","med","asympt_Disp","extraPoisson","rep","scale","n","sizefactor")
median_scaleddata <- read.table(text = "",
                                col.names = median_colnames,
                                stringsAsFactors = FALSE)
library(data.table)
library(plyr)
library(dplyr)

for (k in countfiles){
  detach("package:plyr")
  library(dplyr)
  disp_info <- str_split(k, "_")
  asympt_Dispersion = substring(disp_info[[1]][3], 2)
  extraPoisson_noise= substring(disp_info[[1]][4], 2)
  rep_number= substring(disp_info[[1]][5], 2)
  scale_number=substring(disp_info[[1]][6], 2)
  n_number <- substring(disp_info[[1]][7], 2, nchar(disp_info[[1]][7])-4)
  sizefactor_type <- "DEFAULT"
  asympt_Dispersion
  extraPoisson_noise
  rep_number
  scale_number
  n_number
  RNAcountdat <- read.csv(paste0("/Users/sahu0957/trisomy_normalization/RNAGRO01132021/data/RNAseqfiles/",
                                 k), 
                          sep=",", row.names=1)
  #Remove genes that are not in annotation (which was filtered to remove genes not in both and not on main chromosome above) 
  #put the data frames in the same order
  
  
  RNAcountdat<- RNAcountdat %>%
    filter(Gene_id %in% annotationmerge$name) %>%
    arrange(Gene_id) %>%
    column_to_rownames('Gene_id')
  #read in metadata about the files
  names <- c("filename","tech_rep","biological_rep","Family","Person","ploidy","samplegroup")
  RNAmetadata <- read.table(text = "",
                            col.names = names,
                            stringsAsFactors = FALSE)
  RNAmetadata <- as.data.frame(RNAmetadata,stringsAsFactors=FALSE)
  colnames(RNAmetadata) <- names
  for (i in 1:(ncol(RNAcountdat)/2)){
    newmetadata1 <- c(paste0("Elvis_rep",i,"_counts"),"0",as.character(i),"E","Elvis","ploidy_typical",paste0("Elvis_",i))
    RNAmetadata[nrow(RNAmetadata)+1,] <- newmetadata1
  }
  
  for (i in 1:(ncol(RNAcountdat)/2)){
    newmetadata2 <- c(paste0("Eddie_rep",i,"_counts"),"0",as.character(i),"E","Eddie","ploidy_trisomy",paste0("Eddie_",i))
    RNAmetadata[nrow(RNAmetadata)+1,] <- newmetadata2
  }
  
  RNAdds <- DESeqDataSetFromMatrix(countData = RNAcountdat, colData = RNAmetadata, design = ~Person)
  RNAdds <-DESeq(RNAdds,fitType = c("parametric"))
  
  detach("package:dplyr")
  library(plyr)
  
  person1 = "Eddie"
  person2 = "Elvis"
  RNAddsres<-results(RNAdds,  contrast=c("Person", person1, person2))
  resdata <- as.data.frame(RNAddsres)
  fullresdata <- merge(resdata,annotationmerge,by.x=0,by.y=3)
  fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]
  medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]
  medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
  medians$asympt_Disp <- asympt_Dispersion
  medians$extraPoisson <- extraPoisson_noise
  medians$rep <- rep_number
  medians$scale <- scale_number
  medians$n <- n_number
  medians$sizefactor <- sizefactor_type
  median_scaleddata <- rbind(median_scaleddata,medians)
  
  # With size factor adjustment
  sizefactor_type <- "ADJUSTED"
  controlgenes <-ifelse(rownames(dds) %in% anuplodygenes, FALSE, TRUE)
  
  RNAdds <- DESeqDataSetFromMatrix(countData = RNAcountdat, 
                                   colData = RNAmetadata, design = ~Person)
  RNAdds <-estimateSizeFactors(RNAdds, controlGenes=controlgenes) #adds a column to colData(ddsCollapsed) called sizeFactor if you use normmatrix then normalizationFactors(ddsCollapsed) instead.
  
  ### Correcting for ploidy should yield a distribution tightly around 1.0
  RNAdds <- estimateDispersionsGeneEst(RNAdds) #adds mu to assays(ddsCollapsed) mu-has one column per sample and fills in 
  #elementMetadata(ddsCollapsed) with baseMean, baseVar, allZero, dispGeneEst, dispGeneIter
  RNAdds<-estimateDispersionsFit(RNAdds,fitType="parametric") # and fills in dispersionFunction(ddsCollapsed) and adds a column in  elementMetadata(ddsCollapsed) called dispFit
  RNAdds <- estimateDispersionsMAP(RNAdds) #adds a attr(,"dispPriorVar") to dispersionFunction(ddsCollapsed) 
  #and adds 4 columns to elementMetadata(ddsCollapsed): dispersion, dispIter, dispOutlier, dispMAP
  RNAdds <- nbinomWaldTest(RNAdds,betaPrior = FALSE) #adds both H and cooks to assays(ddsCollapsed)  both are for each sample you have.
  #Adds a ton of columns to elementMetadata(ddsCollapsed) 
  #including Intercept and all values they expect you to compare 
  
  person1 = "Eddie"
  person2 = "Elvis"
  RNAddsres<-results(RNAdds,  contrast=c("Person", person1, person2))
  resdata <- as.data.frame(RNAddsres)
  fullresdata <- merge(resdata,annotationmerge,by.x=0,by.y=3)
  fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]
  medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]
  medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
  medians$asympt_Disp <- asympt_Dispersion
  medians$extraPoisson <- extraPoisson_noise
  medians$rep <- rep_number
  medians$scale <- scale_number
  medians$n <- n_number
  medians$sizefactor <- sizefactor_type
  median_scaleddata <- rbind(median_scaleddata,medians)
  
  medians
  
}
medians
#saveit_sscaled_rna <- median_scaleddata
#median_scaleddata <- saveit_sscaled_rna
median_scaleddata
saveit_sscaled
write.table(median_scaleddata,file ="/Users/sahu0957/tmp.csv",sep="\t",row.names = FALSE)
median_scaleddata <- read.table(file="/Users/sahu0957/tmp.csv",sep="\t",header = TRUE)
file.remove("/Users/sahu0957/tmp.csv")

# To visualize, we first "marginalize" on one of the parameters so we're comparing like with like on all samples
#test_mediandata <- as.data.frame(median_scaleddata[median_scaleddata$chr=="chr21",])
test_mediandata <- as.data.frame(median_scaleddata[median_scaleddata$sizefactor=="DEFAULT",])
test_mediandata <- test_mediandata[which(test_mediandata$asympt_Disp==.01), ]
test_mediandata <- test_mediandata[which(test_mediandata$n==3), ]
test_mediandata <- test_mediandata[which(test_mediandata$extraPoisson==1), ]
#test_mediandata <- test_mediandata[which(test_mediandata$scale>.01), ]
test_mediandata$scale <- as.numeric(test_mediandata$scale)
View(median_scaleddata)
test_mediandata
ggplot(test_mediandata, aes(x=extraPoisson,y=med,color=chr)) +
  geom_point() +
  geom_smooth(method="gam",span=1,level=.95) +
  labs(colour = "Chromosome",title = "Depth, a=.01, b=1") +
  geom_hline(yintercept = 1.5,linetype="dashed") +
  geom_hline(yintercept = 1.0, linetype="dashed") +
  ylab("Median Fold Change") +
  xlab("Depth (Scale)") +
  #  annotate("segment", x = 0.001, xend = 0.001, y = 1.2, yend = 1.28, colour = "red", size=0.5, alpha=1, arrow=arrow()) +
  #  annotate("segment", x = 0.16, xend = 0.16, y = 1.12, yend = 1.2, colour = "red", size=0.5, alpha=1, arrow=arrow()) +
  #xlab(expression(paste(alpha))) +
  theme_classic(base_size=18)

ggsave("rnaseq_deseq2_bothchrs_varyscale_chr22.svg",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "svg")
ggsave("rnaseq_deseq2_bothchrs_varyscale_chr22.png",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "png")

#### START CORRECTING NEWER SIM DATA####

# Again, deliberately picking a real noisy file to try this out on.
# These sims have NOT had any repeats removed. I'll need to remake these sims to reflect that later
RNAcountdat <- read.csv("/Users/sahu0957/trisomy_normalization/RNAGRO01132021/data/RNAseqfiles/LAGRANGEreplicatefinalsims_a0.08_p8_r1_s1_n3.csv",
                        sep=",", row.names=1)
#Remove genes that are not in annotation (which was filtered to remove genes not in both and not on main chromosome above) 
#put the data frames in the same order
RNAcountdat<- RNAcountdat %>%
  filter(Gene_id %in% annotationmerge$name) %>%
  arrange(Gene_id) %>%
  column_to_rownames('Gene_id')
#read in metadata about the files
names <- c("filename","tech_rep","biological_rep","Family","Person","ploidy","samplegroup")

RNAmetadata <- read.table(text = "",
                          col.names = names,
                          stringsAsFactors = FALSE)
RNAmetadata <- as.data.frame(RNAmetadata,stringsAsFactors=FALSE)
colnames(RNAmetadata) <- names
head(RNAcountdat)
# We have to make fake metadata here for the design matrix. Note that sims incorporate NO batch effects!
for (i in 1:(ncol(RNAcountdat)/4)){
  newmetadata1 <- c(paste0("Elvis_rep",i,"_counts"),"0",as.character(i),"E","Elvis","ploidy_typical",paste0("Elvis_",i))
  RNAmetadata[nrow(RNAmetadata)+1,] <- newmetadata1
}

for (i in 1:(ncol(RNAcountdat)/4)){
  newmetadata1 <- c(paste0("EDad_rep",i,"_counts"),"0",as.character(i),"E","EDad","ploidy_typical",paste0("EDad_",i))
  RNAmetadata[nrow(RNAmetadata)+1,] <- newmetadata1
}

for (i in 1:(ncol(RNAcountdat)/4)){
  newmetadata1 <- c(paste0("EMom_rep",i,"_counts"),"0",as.character(i),"E","EMom","ploidy_typical",paste0("EMom_",i))
  RNAmetadata[nrow(RNAmetadata)+1,] <- newmetadata1
}

for (i in 1:(ncol(RNAcountdat)/4)){
  newmetadata2 <- c(paste0("Eddie_rep",i,"_counts"),"0",as.character(i),"E","Eddie","ploidy_trisomy21",paste0("Eddie_",i))
  RNAmetadata[nrow(RNAmetadata)+1,] <- newmetadata2
}

RNAbed <-ori_worldbed
filetable <- RNAmetadata
minichrs <-c("chr1", "chr2", "chr3","chr4", "chr5", "chr6", "chr7", "chr8", "chr9","chr10","chr11", "chr12", "chr13","chr14", "chr15", "chr16", "chr17", "chr18", "chr19","chr20","chr21","chr22")

masterannotationdf <- annotationmerge
masterannotationdf_only21 <- masterannotationdf %>% filter(chr=="chr21")
head(masterannotationdf_only21)
anuplodygenes <-masterannotationdf_only21[["name"]] #Only chr21 genes
baseploidy <- 2
alt_ploidy <-3

#run Deseq on RNA-seq with 21
ddsFull <- DESeqDataSetFromMatrix(countData = RNAcountdat, colData = RNAmetadata, design = ~ Person) #use this for all samples
#ddsFull$type <- relevel(ddsFull$type, "WT")
#dds <- collapseReplicates( ddsFull,groupby = ddsFull$samplegroup,run = ddsFull$samplegroup )
dds <- ddsFull

ddsnocorr <- DESeq(ddsFull)
lesschrs
person1 = "Eddie" #T21
person2 = "Elvis" #D21
RNAddsres<-results(ddsnocorr,  contrast=c("Person", person1, person2))
resdata <- as.data.frame(RNAddsres)
fullresdata <- merge(resdata,annotationmerge,by.x=0,by.y=3)
common_fullresdata_rna <- merge(fullresdata,common_ids,by.x=1,by.y=1)

fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]
medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]


library(plyr)

medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians

signif_medresdata <- medresdata[medresdata$padj<.01,]
notsignif_medresdata <- medresdata[medresdata$padj>=.01,]
ggplot() + 
  geom_violin(data=fullresdata,trim=TRUE,aes(x=chr, y=log2FoldChange)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  #  geom_hline(yintercept=log2(0.66667),linetype="dashed",color="red") +
  geom_hline(yintercept=0) +
  geom_hline(yintercept=log2(1.5),color="blue",linetype="dashed") +
  ylab(paste0("Log2FC")) +
  xlab("Chromosome") +
  theme_classic(base_size = 24) +
  geom_jitter(data = signif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=1.5, alpha=0.85,color="red",fill="red") +
  geom_jitter(data = notsignif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=1.5, alpha=0.85,color="black",fill="black") +
  stat_summary(data=fullresdata,aes(x=chr, y=log2FoldChange),fun.y=median, geom="crossbar", size=0.25, color="orange") +
  theme(legend.title = element_blank())

ggsave("simseq_lowreps_08disp_8poiss_uncorrected.svg",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "svg")
ggsave("simseq_lowreps_08disp_8poiss_uncorrected.png",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "png")

