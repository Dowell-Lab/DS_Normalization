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
RNAcountdat <- read.csv("/Users/sahu0957/trisomy_normalization/RNAGRO01132021/data/RNAseqfiles/NEWSCALEreplicatefinalsims_a0.0", 
                        sep=",", row.names=1)
head(RNAcountdat)
nrow(RNAcountdat)
#RNAcountdat$Gene_id <- rownames(RNAcountdat)
library(dplyr)
nrow(annotationmerge)
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
for (i in 1:(ncol(RNAcountdat)/2)){
  newmetadata1 <- c(paste0("Elvis_rep",i,"_counts"),"0",as.character(i),"E","Elvis","ploidy_typical",paste0("Elvis_",i))
  RNAmetadata[nrow(RNAmetadata)+1,] <- newmetadata1
}

for (i in 1:(ncol(RNAcountdat)/2)){
  newmetadata2 <- c(paste0("Eddie_rep",i,"_counts"),"0",as.character(i),"E","Eddie","ploidy_trisomy",paste0("Eddie_",i))
  RNAmetadata[nrow(RNAmetadata)+1,] <- newmetadata2
}

# Don't do an interaction term with sims, since they are generated independent of batch info... 
# We could always make sims from each replicate instead of sampling a normal distribution around their means
# This would mirror biological replication. 
RNAmetadata
RNAdds <- DESeqDataSetFromMatrix(countData = RNAcountdat, colData = RNAmetadata, design = ~Person)
RNAdds <-DESeq(RNAdds, betaPrior = FALSE) # Earlier DESeq2 versions include a betaPrior for the intercept term... wild!

person1 = "Eddie"
person2 = "Elvis"
RNAddsres<-results(RNAdds,  contrast=c("Person", person1, person2))
resdata <- as.data.frame(RNAddsres)
fullresdata <- merge(resdata,annotationnew,by.x=0,by.y=4)
fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]
medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]

detach("package:dplyr")
library(plyr)
medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians

medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]
View(medresdata)
library(plyr)
medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians

common_fullresdata_rna <- merge(medresdata,common_ids,by.x=1,by.y=1)

nrow(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$padj<.01,])
nrow(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$chr=="chr21",])
(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$padj<.01 & common_fullresdata_rna$chr=="chr21",])
(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$padj<.01 & common_fullresdata_rna$chr=="chr21" & common_fullresdata_rna$log2FoldChange<0,])

nrow(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$log2FoldChange<.48 & common_fullresdata_rna$chr=="chr21",])


nrow(common_fullresdata_rna[common_fullresdata_rna$padj<.01,])

signif_medresdata <- medresdata[medresdata$padj<.01,]
notsignif_medresdata <- medresdata[medresdata$padj>=.01,]
View(signif_medresdata)
signif_medresdata
nrow(notsignif_medresdata)
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

ggsave("rnaseq_deseq2_realsims_new_violins_med101_med154.svg",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "svg")
ggsave("rnaseq_deseq2_realsims_new_violins_med101_med154.png",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "png")




ggplot(fullresdata, aes(x=chr, y=log2FoldChange)) + 
  geom_violin(trim=TRUE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_hline(yintercept=log2(1.5),linetype="dashed",color="blue") +
  #  geom_hline(yintercept=log2(0.66667),linetype="dashed",color="red") +
  geom_hline(yintercept=0) +
  ylab(paste0("log2FC(Sim-T21/Sim-D21)")) +
  xlab("Chromosome") +
  theme_classic(base_size = 24) +
  geom_jitter(shape=23, position=position_jitter(0.1,seed = 1),size=0.9, alpha=0.5,color="black",fill="black") +
  scale_colour_manual(values=c("black","black"), guide=FALSE) +
  scale_fill_manual(values=c("black","black"), guide = FALSE) +
  stat_summary(fun.y=median, geom="crossbar", size=0.25, color="orange") +
  #  stat_summary(fun.y=mean, geom="crossbar", size=0.25, color="green") +
  theme(legend.title = element_blank())

ggsave("rnaseq_deseq2_realsims_violins_withfcshrinkage_med101_med117.svg",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "svg")
ggsave("rnaseq_deseq2_realsims_violins_withfcshrinkage_med101_med117.png",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "png")

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
  fullresdata <- merge(resdata,annotationnew,by.x=0,by.y=4)
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
  fullresdata <- merge(resdata,annotationnew,by.x=0,by.y=4)
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
test_mediandata <- as.data.frame(median_scaleddata[median_scaleddata$chr=="chr21",])
test_mediandata <- as.data.frame(median_scaleddata[median_scaleddata$sizefactor=="DEFAULT",])
View(test_mediandata)
test_mediandata <- test_mediandata[which(test_mediandata$asympt_Disp==.01), ]
test_mediandata <- test_mediandata[which(test_mediandata$n==3), ]
test_mediandata <- test_mediandata[which(test_mediandata$extraPoisson==1), ]
test_mediandata <- test_mediandata[which(test_mediandata$scale>.01), ]
test_mediandata$scale <- as.numeric(test_mediandata$scale)

test_mediandata
ggplot(test_mediandata, aes(x=scale,y=med,color=chr)) +
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

ggsave("rnaseq_deseq2_bothchrs_varyscale.svg",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "svg")
ggsave("rnaseq_deseq2_bothchrs_varyscale.png",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "png")


###looping runs with NO correction: GRO-seq#### 
# We can see whether there's a pattern associated with more/less dispersion by looping DESeq2 runs through all of our simulated files
# each count file is in the same format, where the sim info is in the name of the file and the counts in the file itself
# so we can parse out the parameters of the simulation first and then run DESeq2
# a = asymptotic dispersion
# p = extra-poisson noise
# n = number of replicates in a single experiment
# r = number of experiments
# s = scalar for depth
countfiles <- list.files(path = "/Users/sahu0957/trisomy_normalization/RNAGRO01132021/data/RNAseqfiles/",pattern = "GROSIMS")
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
  asympt_Dispersion = substring(disp_info[[1]][2], 10)
  extraPoisson_noise= substring(disp_info[[1]][3], 2)
  rep_number= substring(disp_info[[1]][4], 2)
  scale_number=substring(disp_info[[1]][5], 2)
  n_number <- substring(disp_info[[1]][6], 2, nchar(disp_info[[1]][6])-4)
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
  fullresdata <- merge(resdata,annotationnew,by.x=0,by.y=4)
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
  fullresdata <- merge(resdata,annotationnew,by.x=0,by.y=4)
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
#saveit_sscaled <- median_scaleddata
median_scaleddata <- saveit_sscaled
median_scaleddata
saveit_sscaled
write.table(median_scaleddata,file ="/Users/sahu0957/tmp.csv",sep="\t",row.names = FALSE)
median_scaleddata <- read.table(file="/Users/sahu0957/tmp.csv",sep="\t",header = TRUE)
file.remove("/Users/sahu0957/tmp.csv")

# To visualize, we first "marginalize" on one of the parameters so we're comparing like with like on all samples
test_mediandata <- as.data.frame(median_scaleddata[median_scaleddata$chr=="chr21",])
test_mediandata <- as.data.frame(median_scaleddata[median_scaleddata$sizefactor=="DEFAULT",])
View(test_mediandata)
test_mediandata <- test_mediandata[which(test_mediandata$asympt_Disp==.01), ]
test_mediandata <- test_mediandata[which(test_mediandata$n==3), ]
test_mediandata <- test_mediandata[which(test_mediandata$extraPoisson==1), ]
test_mediandata <- test_mediandata[which(test_mediandata$scale>.01), ]
test_mediandata$scale <- as.numeric(test_mediandata$scale)

test_mediandata
ggplot(test_mediandata, aes(x=scale,y=med,color=chr)) +
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

ggsave("groseq_deseq2_bothchrs_varyscale.svg",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "svg")
ggsave("groseq_deseq2_bothchrs_varyscale.png",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "png")

ggplot() + 
  geom_violin(data=test_mediandata,trim=TRUE,aes(x=chr, y=log2FoldChange)) +
  #  geom_hline(yintercept=log2(0.66667),linetype="dashed",color="red") +
  geom_hline(yintercept=0) +
  ylab(paste0("Log2FC")) +
  geom_hline(yintercept = log2(1.5), linetype="dotted", color="blue" ) +
  theme_classic(base_size = 30) +
  #  geom_jitter(data = signif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=2.5, alpha=0.85,color="red",fill="red") +
  #  geom_jitter(data = notsignif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=2.5, alpha=0.85,color="black",fill="black") +
  stat_summary(data=fullresdata,aes(x=chr, y=log2FoldChange),fun.y=median, geom="crossbar", size=0.25, color="orange") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.title = element_blank()) +
  theme(axis.title.x = element_blank())

ggsave("groseq_deseq2_realsims_varyscale.svg",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "svg")
ggsave("groseq_deseq2_realsims_varyscale.png",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "png")


test_mediandata <- median_scaleddata[median_scaleddata$extraPoisson==1,]

ggplot(test_mediandata, aes(x=extraPoisson,y=med,color=chr)) +
  geom_smooth(method="gam") +
  geom_point() +
  geom_hline(yintercept = 1.5,linetype="dashed") +
  ylab("Median Fold Change") +
  xlab("extra-Poisson Noise") +
  #xlab(expression(paste(alpha))) +
  theme_classic(base_size=18)

ggsave("rnaseq_deseq2_realsims_3reps_varydisp.svg",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "svg")
ggsave("rnaseq_deseq2_realsims_3reps_varydisp.png",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "png")

#####LOOPING DONE!#############

####START REAL DATA####
### Now we can pick a couple values that look really poor, and use our methods to correct them
# Dispersion/noise is higher when chr21 is included. Scaling factors are incorrect. Low expression genes
# are very noisy
# 1. Use NormFactors, and SizeFactors with controlgenes
# 2. Remove low expression genes from MFC calculations

########RNA-SEQ#######
####Uncorrected Real RNA-seq data####
RNAbed <-ori_worldbed
filetable <- RNAmetadata
minichrs <-c("chr1", "chr2", "chr3","chr4", "chr5", "chr6", "chr7", "chr8", "chr9","chr10","chr11", "chr12", "chr13","chr14", "chr15", "chr16", "chr17", "chr18", "chr19","chr20","chr21","chr22")
annotationmerge
masterannotationdf <- annotationmerge

library(dplyr)
masterannotationdf_only21 <- masterannotationdf %>% filter(chr=="chr21")
anuplodygenes <-masterannotationdf_only21[["name"]] #Only chr21 genes
baseploidy <- 2
alt_ploidy <-3

#this is loading the metadata
RNAmetadata=read.table("/Shares/down/mixed/RNAGRO01132021/scripts/Deseq2_analysis/RNAinfo.txt", sep="\t", header=TRUE)
RNAcoveragedat

RNAcountdat <- read.csv("/Users/sahu0957/trisomy_normalization/RNAGRO01132021/data/RNAseqfiles/res_featureCounts_gene_idfull_143138.coverage.csv", 
                        sep=",", row.names=1)
head(RNAcountdat)
RNAcountdat$Gene_id <- rownames(RNAcountdat)
RNAcountdat<- RNAcountdat %>%
#  filter(Gene_id %in% annotationmerge$name) %>%
  arrange(Gene_id)

RNAcountdat <- RNAcountdat[,-c(ncol(RNAcountdat))]

RNAmetadata=read.table("/Shares/down/mixed/RNAGRO01132021/scripts/Deseq2_analysis/RNAinfo.txt", sep="\t", header=TRUE)
RNAmetadata$samplegroup <-paste(RNAmetadata$Person, "_", RNAmetadata$biological_rep, sep="")
head(RNAcountdat)


#this is where you would remove samples if they are bad
#RNAcountdat <- read.csv(RNAcoveragedat,sep="\t",skip=1)
RNAcountdat <- read.csv("/scratch/Users/sahu0957/ds_normalization/RNA/nextflowhg38_20220214/counts/featurecounts_rna_withmultis.sense.txt",
                        sep="\t", skip=1)

RNAcountdat <- RNAcountdat[!(duplicated(RNAcountdat$Geneid)),]
rownames(RNAcountdat) <- RNAcountdat$Geneid
RNAcountdat<- RNAcountdat %>%
#  filter(Geneid %in% masterannotationdf$name) %>%
  arrange(Geneid)


#run Deseq on RNA-seq with 21
RNAcountdat <- RNAcountdat[,-c(1:6)]

RNAcountdat <- RNAcountdat[,!(colnames(RNAcountdat) %like% "Ethan")]
RNAmetadata <- RNAmetadata[!(RNAmetadata$Person %like% "Ethan"),]
ddsFull <- DESeqDataSetFromMatrix(countData = RNAcountdat, colData = RNAmetadata, design = ~biological_rep + Person)

# Can run DESeq2 unaltered here to compare results
ddsFull <- DESeq(ddsFull)
RNAddsres<-results(ddsFull,contrast=c("Person","Ethan","Eric"))
withmultisresdata <- as.data.frame(RNAddsres)

fullresdata <- merge(resdata,annotationnew,by.x=0,by.y=4)
fullresdata <- merge(fullresdata,common_ids,by.x=1,by.y=1)
fullresdata <- na.omit(fullresdata)

RNA_uncorrected_fullresdata <- fullresdata
nrow(RNA_uncorrected_fullresdata)

merged_RNAGRO_fullresdata <- merge(RNA_uncorrected_fullresdata,GRO_uncorrected_fullresdata,by.x=1,by.y=1)
merged_RNAGRO_fullresdata <- merged_RNAGRO_fullresdata[!(duplicated(merged_RNAGRO_fullresdata$V2.y)),]
qPCR_candidates <- merged_RNAGRO_fullresdata[merged_RNAGRO_fullresdata$log2FoldChange.x<log2(1.3) & 
                                               merged_RNAGRO_fullresdata$log2FoldChange.x>log2(1.1) & 
                                               merged_RNAGRO_fullresdata$log2FoldChange.y<log2(1.3) & 
                                               merged_RNAGRO_fullresdata$log2FoldChange.y>log2(1.1) &
                                               merged_RNAGRO_fullresdata$baseMean.x > 50 & 
                                               merged_RNAGRO_fullresdata$baseMean.y > 50 &
                                               merged_RNAGRO_fullresdata$chr.x == "chr21", ]

write.csv(qPCR_candidates,file="/Users/sahu0957/trisomy_normalization/qpcr_candidates.csv",sep = "\t")
write.csv(merged_RNAGRO_fullresdata,file="/Users/sahu0957/trisomy_normalization/rna_gro_uncorrected_withmultimaps_results.csv",sep = "\t")

common_fullresdata_test <- merge(fullresdata,common_ids,by.x=1,by.y=1)

#Total chr21 genes:
nrow(common_fullresdata_test[!(duplicated(common_fullresdata_test$V2)) & common_fullresdata_test$chr=="chr21" & is.na(common_fullresdata_test$log2FoldChange),])

#Total chr21 genes Up:
nrow(common_fullresdata_test[!(duplicated(common_fullresdata_test$V2)) & common_fullresdata_test$chr=="chr21" & common_fullresdata_test$log2FoldChange>=.5849 & !(is.na(common_fullresdata_test$log2FoldChange)),])

#Total chr21 genes Down:
nrow(common_fullresdata_test[!(duplicated(common_fullresdata_test$V2)) & common_fullresdata_test$chr=="chr21" & common_fullresdata_test$log2FoldChange<.5849 & !(is.na(common_fullresdata_test$log2FoldChange)),])


nrow(common_fullresdata_test[!(duplicated(common_fullresdata_test$V2)) & common_fullresdata_test$log2FoldChange<=.5849 & common_fullresdata_test$chr=="chr21",])

# DEGs
nrow(common_fullresdata_test[!(duplicated(common_fullresdata_test$V2)) & common_fullresdata_test$padj<.01 & !(is.na(common_fullresdata_test$padj)),])

# chr21 DEGs
nrow(common_fullresdata_test[!(duplicated(common_fullresdata_test$V2)) & common_fullresdata_test$padj<.01 & common_fullresdata_test$chr=="chr21" & !(is.na(common_fullresdata_test$padj)),])

nrow(common_fullresdata_test[!(duplicated(common_fullresdata_test$V2)) & common_fullresdata_test$log2FoldChange<.5849 & common_fullresdata_test$chr=="chr21",])
nrow(common_fullresdata_test[!(duplicated(common_fullresdata_test$V2)) & common_fullresdata_test$log2FoldChange>=.5849 & common_fullresdata_test$chr=="chr21",])

nrow(fullresdata)

lesschrs
fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]

medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]

library(plyr)
medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians

signif_medresdata <- medresdata[medresdata$padj<.01,]
notsignif_medresdata <- medresdata[medresdata$padj>=.01,]
nrow(signif_medresdata)
nrow(notsignif_medresdata)
medians
fullresdata$chr <- factor(fullresdata$chr, levels=minichrs)

CDF <- ecdf(fullresdata[fullresdata$chr=="chr21",]$log2FoldChange)
View(CDF)
CDF
# draw CDF Plot
sd(na.omit(fullresdata$log2FoldChange))
sd(fullresdata[fullresdata$chr=="chr21",]$log2FoldChange)
ggplot(fullresdata[fullresdata$chr=="chr21",], aes(x = log2FoldChange)) +
  stat_ecdf() +
  geom_vline(xintercept = log2(1.5),color="red",linetype="dotted") +
  geom_vline(xintercept=log2(1.0),color="blue",linetype="dotted") +
  geom_vline(xintercept=log2(1.5)-(2*0.933),color="purple",linetype="dotted") +
  theme_classic()

ggsave("rnaseq_ethaneric_cdf2.svg",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "svg")
ggsave("rnaseq_ethaneric_cdf2.png",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "png")


ggplot() + 
  geom_violin(data=fullresdata,trim=TRUE,aes(x=chr, y=log2FoldChange)) +
  #  geom_hline(yintercept=log2(0.66667),linetype="dashed",color="red") +
  geom_hline(yintercept=0) +
  ylab(paste0("Log2FC")) +
  geom_hline(yintercept = log2(1.5), linetype="dotted", color="blue" ) +
  theme_classic(base_size = 30) +
  geom_jitter(data = signif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=2.5, alpha=0.85,color="red",fill="red") +
  geom_jitter(data = notsignif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=2.5, alpha=0.85,color="black",fill="black") +
  stat_summary(data=fullresdata,aes(x=chr, y=log2FoldChange),fun.y=median, geom="crossbar", size=0.25, color="orange") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.title = element_blank()) +
  theme(axis.title.x = element_blank())

ggsave("rnaseq_elieric_noethan_med1082747_med1139804_28signifonchr20and21.svg",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "svg")
ggsave("rnaseq_elieric_noethan_med1082747_med1139804_28signifonchr20and21.png",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "png")


#Volcano Plot
fullresdata <- merge(resdata,annotationnew,by.x=0,by.y=4)
signif_fullresdata <- fullresdata[fullresdata$padj<.01,]
notsignif_fullresdata <- fullresdata[fullresdata$padj>=.01,]
fullresdata_chr21 <- fullresdata[fullresdata$chr=="chr21",]

fullresdata <- na.omit(fullresdata)
head(fullresdata)
ggplot() + 
  geom_point(data=fullresdata,aes(x=log2FoldChange, y=-log10(pvalue)),size=1,alpha=0.5) +
  geom_point(data=signif_fullresdata,aes(x=log2FoldChange, y=-log10(pvalue)),size=1,alpha=0.5,color="red") +
  geom_point(data=fullresdata_chr21,aes(x=log2FoldChange, y=-log10(pvalue)),size=1,alpha=0.5,color="green") +
  ylab(paste0("-log10(pvalue)")) +
  geom_vline(xintercept = log2(1.5), linetype="dotted", color="blue" ) +
  geom_vline(xintercept = log2(1), linetype="dotted", color="black" ) +
  xlab("log2FC") +
  theme_classic(base_size = 24) +
  theme(legend.title = element_blank())

ggsave("rnaseq_ethaneric_volcano.svg",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "svg")
ggsave("rnaseq_ethaneric_volcano.png",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "png")

#  Quintile of FC violin plots
fullresdata <- merge(resdata,annotationnew,by.x=0,by.y=4)
fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]
fullresdata[is.na(fullresdata)] = 0

quantile(fullresdata$baseMean,probs=seq(0,1,.2))

quantile1 <- fullresdata[fullresdata$baseMean<=2.787199e-01,]
quantile1$quantile <- 1
medians1 <- ddply(quantile1, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians1$quantile <- 1
medians1
quantile1
quantile2 <- fullresdata[fullresdata$baseMean <= 6.657754e+00 & fullresdata$baseMean> 2.787199e-01,]
quantile2$quantile <- 2
medians2 <- ddply(quantile2, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians2$quantile <- 2

quantile3 <- fullresdata[fullresdata$baseMean<= 1.617177e+02 & fullresdata$baseMean> 6.657754e+00,]
quantile3$quantile <- 3
medians3 <- ddply(quantile3, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians3$quantile <- 3

quantile4 <- fullresdata[fullresdata$baseMean<= 6.485463e+02 & fullresdata$baseMean> 1.617177e+02,]
quantile4$quantile <- 4
medians4 <- ddply(quantile4, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians4$quantile <- 4

quantile5 <- fullresdata[fullresdata$baseMean<= 2.512809e+07 & fullresdata$baseMean> 6.485463e+02,]
quantile5$quantile <- 5
medians5 <- ddply(quantile5, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians5$quantile <- 5

quantileall <- fullresdata
quantileall$quantile <- 6
medians6 <- ddply(quantileall, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians6$quantile <- 6

newfullres <- rbind(quantile1,quantile2,quantile3,quantile4,quantile5,quantileall)
newmedians <- rbind(medians1,medians2,medians3,medians4,medians5,medians6)
newmedians

# Facet Grid of Violins
colnames(newfullres)
newfullres <- newfullres[newfullres$chr %in% lesschrs,]

ggplot() + 
  geom_violin(data=newfullres,aes(x=chr, y=log2FoldChange)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_hline(yintercept=0) +
  ylab(paste0("Log2FC")) +
  geom_hline(yintercept = log2(1.5), linetype="dotted", color="blue" ) +
  xlab("Chromosome") +
  theme_classic(base_size = 24) +
  #geom_jitter(aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=2.5, alpha=0.85,color="black",fill="black") +
  #stat_summary(data=newfullres,aes(x=chr, y=log2FoldChange),fun.y=median, geom="crossbar", size=0.25, color="orange") +
  theme(legend.title = element_blank()) + 
  facet_wrap(~quantile, ncol = 5) +
  geom_hline(data = newmedians, aes(yintercept = log2(med)),
             colour = "orange")

ggsave("blehrnaseq_ethaneric_facets.svg",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "svg")
ggsave("blehrnaseq_ethaneric_facets.png",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "png")

# Fetch all genes dispersion unchanged
#this is loading the metadata
RNAmetadata=read.table("/Shares/down/mixed/RNAGRO01132021/scripts/Deseq2_analysis/RNAinfo.txt", sep="\t", header=TRUE)
RNAcountdat <- read.csv("/scratch/Users/sahu0957/ds_normalization/RNA/nextflowhg38_20220214/counts/featurecounts_rna_withmultis.sense.txt",
                        sep="\t", skip=1)

RNAcountdat <- RNAcountdat[!(duplicated(RNAcountdat$Geneid)),]
RNAcountdat<- RNAcountdat %>%
  filter(Geneid %in% masterannotationdf$name) %>%
  arrange(Geneid) %>%
  column_to_rownames('Geneid')

head(RNAcountdat)
#run Deseq on RNA-seq with 21
RNAcountdat <- RNAcountdat[,-c(1:5)]
RNAcountdat <- RNAcountdat[,!(colnames(RNAcountdat) %like% "Eric")]
RNAmetadata <- RNAmetadata[!(RNAmetadata$Person %like% "Eric"),]
ddsFull <- DESeqDataSetFromMatrix(countData = RNAcountdat, colData = RNAmetadata, design = ~biological_rep + Person)
ddsFull <- DESeq(ddsFull)
annotationmergeno21
metrics_info_unchanged <- as.data.frame(mcols(ddsFull))
metrics_info_unchanged
metrics_info_unchanged <- metrics_info_unchanged[!(rownames(metrics_info_unchanged) %in% rownames(annotationmergeno21)),]
nrow(metrics_info_unchanged)
# Fetch all gene dispersion leaving chr21

#RNAcountdat <- RNAcountdat[(rownames(RNAcountdat) %in% annotationmergeno21$GeneID),]
ddsFull <- DESeqDataSetFromMatrix(countData = RNAcountdat, colData = RNAmetadata, design = ~biological_rep + Person)
ddsFull <- DESeq(ddsFull)

metrics_info_no21 <- mcols(ddsFull)
nrow(mcols(ddsFull))
# Fetch all gene dispersions leaving out Ethan
RNAmetadata=read.table("/Shares/down/mixed/RNAGRO01132021/scripts/Deseq2_analysis/RNAinfo.txt", sep="\t", header=TRUE)
RNAcountdat <- read.csv("/scratch/Users/sahu0957/ds_normalization/RNA/nextflowhg38_20220214/counts/featurecounts_rna_withmultis.sense.txt",
                        sep="\t", skip=1)

RNAcountdat <- RNAcountdat[!(duplicated(RNAcountdat$Geneid)),]
RNAcountdat<- RNAcountdat %>%
  filter(Geneid %in% masterannotationdf$name) %>%
  arrange(Geneid) %>%
  column_to_rownames('Geneid')

#run Deseq on RNA-seq with Ethan removed
RNAcountdat <- RNAcountdat[,-c(1:5)]
RNAcountdat <- RNAcountdat[,!(colnames(RNAcountdat) %like% "Ethan")]
RNAmetadata <- RNAmetadata[!(RNAmetadata$Person %like% "Ethan"),]
ddsFullnew <- DESeqDataSetFromMatrix(countData = RNAcountdat, colData = RNAmetadata, design = ~biological_rep + Person)
ddsFullnew <- DESeq(ddsFullnew)
metrics_info_noethan <- mcols(ddsFullnew)
metrics_info_noethan <- metrics_info_noethan[!(rownames(metrics_info_noethan) %in% rownames(annotationmergeno21)),]

nrow(metrics_info_no21)
nrow(metrics_info_noethan)
nrow(metrics_info_unchanged)
as.data.frame(metrics_info_no21)
tmpmerge <- merge(as.data.frame(metrics_info_unchanged),as.data.frame(metrics_info_no21),by.x=0,by.y=0)
dispersion_nochr21_noethan <- merge(as.data.frame(metrics_info_unchanged),as.data.frame(metrics_info_noethan),by.x=0,by.y=0)
View(dispersion_nochr21_noethan)
colnames(dispersion_nochr21_noethan)
View(dispersion_nochr21_noethan)
dispersion_nochr21_noethan[dispersion_nochr21_noethan$dispMAP.x]
as.data.frame(dispersion_nochr21_noethan$dispMAP.x)

eq <- function(x,y) {
  m <- lm(y ~ x)
  as.character(
    as.expression(
      substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                 list(a = format(coef(m)[1], digits = 4),
                      b = format(coef(m)[2], digits = 4),
                      r2 = format(summary(m)$r.squared, digits = 3)))
    )
  )
}

eq(dispersion_nochr21_noethan$dispMAP.x,dispersion_nochr21_noethan$dispMAP.y)

ggplot(df,aes(x = wt, y = hp)) + 
  geom_point() + 
  geom_smooth(method = "lm", se=FALSE) 
View(dispersion_nochr21_noethan)
ggplot() + 
  geom_point(data=dispersion_nochr21_noethan,aes(x=as.numeric(as.character(dispFit.x)), y=as.numeric(as.character(dispFit.y)))) +
  geom_smooth(method="lm",color="blue") +
  geom_abline(color="red") +
  ylab(paste0("T21 Removed")) +
  xlab("D21 Brother Removed") +
  scale_x_log10() +
  scale_y_log10() +
  geom_text(x = 0, y = 0, label = eq(dispersion_nochr21_noethan$dispMAP.y,dispersion_nochr21_noethan$dispMAP.y), parse = TRUE) +
  ggtitle("Gene-wise Dispersion Fitted Estimates") +
  theme_classic(base_size = 24)
  #geom_jitter(aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=2.5, alpha=0.85,color="black",fill="black") +
  #stat_summary(data=newfullres,aes(x=chr, y=log2FoldChange),fun.y=median, geom="crossbar", size=0.25, color="orange") +
 # theme(legend.title = element_blank()) 

ggsave("rnaseq_genewise_dispersion_brothersremoved.svg",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "svg")
ggsave("rnaseq_genewise_dispersion_brothersremoved.png",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "png")


####Correcting Real Data ####


controlgenes <-ifelse(rownames(dds) %in% anuplodygenes, FALSE, TRUE) # Only normalize on non-aneuploidy genes

samplepinfo<-as.data.frame(colData(dds))

ploidy_trisomy21 = ifelse(rownames(dds) %in% anuplodygenes, alt_ploidy, baseploidy)
ploidy_typical <- rep(baseploidy, nrow(dds))
ploidy_trisomy21

normFactors <- matrix(do.call(cbind, mget(paste0(dds$ploidy))),ncol=ncol(dds),nrow=nrow(dds),dimnames=list(1:nrow(dds),1:ncol(dds)))
normFactors <- normFactors/baseploidy

unique(normFactors)


ddsCollapsed<-DESeq(dds)
ddsCollapsed_normfactor <-estimateSizeFactors(dds,normMatrix=normFactors, controlGenes=controlgenes) #adds a column to colData(ddsCollapsed) called sizeFactor if you use normmatrix then normalizationFactors(ddsCollapsed) instead.
ddsCollapsed_normfactor <- estimateDispersionsGeneEst(ddsCollapsed_normfactor) #adds mu to assays(ddsCollapsed) mu-has one column per sample and fills in 
#elementMetadata(ddsCollapsed) with baseMean, baseVar, allZero, dispGeneEst, dispGeneIter
ddsCollapsed_normfactor<-estimateDispersionsFit(ddsCollapsed_normfactor) # and fills in dispersionFunction(ddsCollapsed) and adds a column in  elementMetadata(ddsCollapsed) called dispFit
ddsCollapsed_normfactor <- estimateDispersionsMAP(ddsCollapsed_normfactor) #adds a attr(,"dispPriorVar") to dispersionFunction(ddsCollapsed) 
#and adds 4 columns to elementMetadata(ddsCollapsed): dispersion, dispIter, dispOutlier, dispMAP
ddsCollapsed_normfactor <- nbinomWaldTest(ddsCollapsed_normfactor) #adds both H and cooks to assays(ddsCollapsed)  both are for each sample you have.
#Adds a ton of columns to elementMetadata(ddsCollapsed) 
#including Intercept and all values they expect you to compare 
#SE_Intercept (SE stands for standard error)  and SE for every thing they expect you to compare. 
#Wald static for what they expect you to compare. 
#betaConv  betaIter         deviance  maxCooks

person1 = "Ethan"
person2 = "Eric"
RNAddsres<-results(ddsCollapsed_normfactor,  contrast=c("Person", person1, person2))
resdata <- as.data.frame(RNAddsres)
fullresdata <- merge(resdata,annotationnew,by.x=0,by.y=4)
fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]
quantile(fullresdata$baseMean, probs = seq(0,1,0.2))
# Select a value which removes the bottom 2 quintiles of expression
#fullresdata <- fullresdata[fullresdata$baseMean > 12,] # Remove genes with very low expression
medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]
nrow(fullresdata)
nrow(fullresdata[fullresdata$padj<0.05 & fullresdata$chr == "chr21",])
medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians



####REMOVING REPEATS####
# We noticed that a lot of samples had specific genes that always appeared near 1.0, especially rRNA. Could it be
# that repeats are causing multi-mapping issues for these genes? To test, we removed reads over RNA repeats in genome

RNAbed <-ori_worldbed
filetable <- RNAmetadata
minichrs <-c("chr1", "chr2", "chr3","chr4", "chr5", "chr6", "chr7", "chr8", "chr9","chr10","chr11", "chr12", "chr13","chr14", "chr15", "chr16", "chr17", "chr18", "chr19","chr20","chr21","chr22")

masterannotationdf <- annotationmerge
masterannotationdf_only21 <- masterannotationdf %>% filter(chr=="chr21")
head(masterannotationdf_only21)
anuplodygenes <-masterannotationdf_only21[["name"]] #Only chr21 genes
baseploidy <- 2
alt_ploidy <-3

#this is loading the metadata
RNAmetadata=read.table("/Shares/down/mixed/RNAGRO01132021/scripts/Deseq2_analysis/RNAinfo.txt", sep="\t", header=TRUE)

#this is where you would remove samples if they are bad
RNAcountdat <- read.csv("/scratch/Users/sahu0957/ds_normalization/RNA/nextflowhg38_20220214/counts/featurecounts_rna_rmsked.sense.txt",
                        sep="\t", skip=1)

nrow(RNAcountdat)

# filter to only canonical chrs
head(RNAcountdat)
#RNAcountdat <- RNAcountdat[RNAcountdat$Chr %in% minichrs,]

colnames(RNAcountdat)
# colnames are different due to my featurecounts script. fixing that here by removing the path. Kinda hacky, go remove in script!
names <- colnames(RNAcountdat)
names <- gsub(".*removed..","",names)
colnames(RNAcountdat) <- names
colnames(RNAcountdat)

# We've got 1300 duplicate entries... How was this filtered out before?
RNAcountdat <- RNAcountdat[!(duplicated(RNAcountdat$Geneid)),]
# Filter to genes common in both RNA and GRO seq

library(tibble)
RNAcountdat<- RNAcountdat %>%
#  filter(Geneid %in% annotationmerge$GeneID) %>%
  arrange(Geneid) %>%
  column_to_rownames('Geneid')
dim(RNAcountdat)
RNAcountdat <- RNAcountdat[,-c(1:5)]
#run Deseq on RNA-seq with 21
colnames(RNAcountdat)
ddsFull <- DESeqDataSetFromMatrix(countData = RNAcountdat, colData = RNAmetadata, design = ~ biological_rep + Person) #use this for all samples
dds <- ddsFull

# You can see the results of unfiltered DESeq2 by uncommenting this, but we have more correcting to do first!
ddsFull <- DESeq(ddsFull)

person1 = "Ethan"
person2 = "Eric"
RNAddsres<-results(ddsFull,  contrast=c("Person", person1, person2))
nomultisresdata <- as.data.frame(RNAddsres)
withmultisresdata
nomultisresdata
withandnoresdata <- merge(withmultisresdata,nomultisresdata,by.x=0,by.y=0)
withandnoresdata_annot <- merge(withandnoresdata,annotationnew,by.x=1,by.y=4)
View(withandnoresdata_annot)
colnames(withandnoresdata_annot)

ggplot() + 
  geom_point(data=withandnoresdata_annot[withandnoresdata_annot$chr=="chr21",],aes(x=log2FoldChange.x, y=log2FoldChange.y),size=1) +
  geom_hline(yintercept=log2(1.5),linetype="dashed",color="red") +
  geom_vline(xintercept=log2(1.5),linetype="dashed",color="red") +
  geom_hline(yintercept=log2(1),linetype="dashed",color="blue") +
  geom_vline(xintercept=log2(1),linetype="dashed",color="blue") +
  ylab(paste0("Log2FC No Repeats")) +
  xlab("Log2FC With Repeats") +
  theme_classic(base_size = 24) +
  theme(legend.title = element_blank())

ggsave("rnaseq_withwithoutmultis.svg",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "svg")
ggsave("rnaseq_withwithoutmultis.png",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "png")


annotationnew
fullresdata <- merge(resdata,annotationnew,by.x=0,by.y=4)
lesschrs <- c("chr21","chr22")
#fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]
#fullresdata <- na.omit(fullresdata)

#fullresdata <- merge(resdata,annotationnew,by.x=0,by.y=4)

common_fullresdata_test <- merge(fullresdata,common_ids,by.x=1,by.y=1)

#Total chr21 genes:
nrow(common_fullresdata_test[!(duplicated(common_fullresdata_test$V2)) & common_fullresdata_test$chr=="chr21",])

#Total chr21 genes Up:
nrow(common_fullresdata_test[!(duplicated(common_fullresdata_test$V2)) & common_fullresdata_test$chr=="chr21" & common_fullresdata_test$log2FoldChange>=.5849 & !(is.na(common_fullresdata_test$log2FoldChange)),])

#Total chr21 genes Down:
nrow(common_fullresdata_test[!(duplicated(common_fullresdata_test$V2)) & common_fullresdata_test$chr=="chr21" & common_fullresdata_test$log2FoldChange<.5849 & !(is.na(common_fullresdata_test$log2FoldChange)),])


nrow(common_fullresdata_test[!(duplicated(common_fullresdata_test$V2)) & common_fullresdata_test$log2FoldChange<=.5849 & common_fullresdata_test$chr=="chr21",])

# DEGs
nrow(common_fullresdata_test[!(duplicated(common_fullresdata_test$V2)) & common_fullresdata_test$padj<.01 & !(is.na(common_fullresdata_test$padj)),])

# chr21 DEGs
nrow(common_fullresdata_test[!(duplicated(common_fullresdata_test$V2)) & common_fullresdata_test$padj<.01 & common_fullresdata_test$chr=="chr21" & !(is.na(common_fullresdata_test$padj)),])



nrow(fullresdata[!(duplicated(fullresdata$Row.names)) & fullresdata$log2FoldChange < 0.5849 & fullresdata$chr=="chr21",])


medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]

medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians

signif_medresdata <- medresdata[medresdata$padj<.01,]
notsignif_medresdata <- medresdata[medresdata$padj>=.01,]
nrow(signif_medresdata)
nrow(notsignif_medresdata)
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

ggsave("rnaseq_unchanged.svg",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "svg")
ggsave("rnaseq_unchanged.png",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "png")


masterannotationdf_only21 <- masterannotationdf %>% filter(chr=="chr21")
head(masterannotationdf_only21)
anuplodygenes <-masterannotationdf_only21[["name"]] #Only chr21 genes
anuplodygenes
baseploidy <- 2
alt_ploidy <-3


controlgenes <-ifelse(rownames(dds) %in% anuplodygenes, FALSE, TRUE) # Only normalize on non-aneuploidy genes
controlgenes
samplepinfo<-as.data.frame(colData(dds))

ploidy_trisomy21 = ifelse(rownames(dds) %in% anuplodygenes, alt_ploidy, baseploidy)
ploidy_trisomy21
ploidy_typical <- rep(baseploidy, nrow(dds))

normFactors <- matrix(do.call(cbind, mget(paste0(dds$ploidy))),ncol=ncol(dds),nrow=nrow(dds),dimnames=list(1:nrow(dds),1:ncol(dds)))
normFactors <- normFactors/baseploidy
unique(normFactors)

ddsCollapsed<-DESeq(dds,betaPrior = FALSE)
ddsCollapsed_normfactor <-estimateSizeFactors(dds,normMatrix=normFactors, controlGenes=controlgenes) #adds a column to colData(ddsCollapsed) called sizeFactor if you use normmatrix then normalizationFactors(ddsCollapsed) instead.
confused <- as.data.frame(normalizationFactors(ddsCollapsed_normfactor))
head(confused)
head(masterannotationdf)
conffused <- merge(confused,masterannotationdf,by.x=0,by.y=0)


### Correcting for ploidy should yield a distribution tightly around 1.0
ddsCollapsed_normfactor <- estimateDispersionsGeneEst(ddsCollapsed_normfactor) #adds mu to assays(ddsCollapsed) mu-has one column per sample and fills in 
#elementMetadata(ddsCollapsed) with baseMean, baseVar, allZero, dispGeneEst, dispGeneIter
ddsCollapsed_normfactor<-estimateDispersionsFit(ddsCollapsed_normfactor) # and fills in dispersionFunction(ddsCollapsed) and adds a column in  elementMetadata(ddsCollapsed) called dispFit
ddsCollapsed_normfactor <- estimateDispersionsMAP(ddsCollapsed_normfactor) #adds a attr(,"dispPriorVar") to dispersionFunction(ddsCollapsed) 
#and adds 4 columns to elementMetadata(ddsCollapsed): dispersion, dispIter, dispOutlier, dispMAP
ddsCollapsed_normfactor <- nbinomWaldTest(ddsCollapsed_normfactor,betaPrior = FALSE) #adds both H and cooks to assays(ddsCollapsed)  both are for each sample you have.
#Adds a ton of columns to elementMetadata(ddsCollapsed) 
#including Intercept and all values they expect you to compare 
#SE_Intercept (SE stands for standard error)  and SE for every thing they expect you to compare. 
#Wald static for what they expect you to compare. 
#betaConv  betaIter         deviance  maxCooks

person1 = "Ethan"
person2 = "Eric"
RNAddsres<-results(ddsCollapsed_normfactor,  contrast=c("Person", person1, person2))
resdata <- as.data.frame(RNAddsres)
fullresdata <- merge(resdata,annotationnew,by.x=0,by.y=4)
fullresdata <- fullresdata[fullresdata$chr %in% minichrs,]
quantile(fullresdata$baseMean, probs = seq(0,1,0.2))
#View(fullresdata)
#fullresdata <- fullresdata[fullresdata$baseMean > 12,] # Remove genes with very low expression
#fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]
#medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]

dosagelowgenes <- medresdata[medresdata$log2FoldChange < 0.13 & medresdata$chr=="chr21",]
#write.csv(dosagelowgenes,'/Users/sahu0957/ds_deseq_normalization/dosagecompensated_genes_norepeats_ethaneli.csv',sep = '\t')

common_fullresdata_rna <- merge(fullresdata,common_ids,by.x=1,by.y=1)
#View(common_fullresdata_rna)
# Total chr21 genes 
nrow(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$chr=="chr21" & common_fullresdata_rna$log2FoldChange<0 & !(is.na(common_fullresdata_rna$log2FoldChange)),])

nrow(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$chr=="chr21" & !(is.na(common_fullresdata_rna$padj)),])
nrow(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$log2FoldChange<0 & common_fullresdata_rna$chr=="chr21" & !(is.na(common_fullresdata_rna$padj)) & common_fullresdata_rna$baseMean>60,])
nrow(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$chr=="chr21" & common_fullresdata_rna$baseMean<60 & !(is.na(common_fullresdata_rna$padj)) & common_fullresdata_rna$log2FoldChange<0,])
nrow(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$chr=="chr21" & common_fullresdata_rna$log2FoldChange>0 & !(is.na(common_fullresdata_rna$log2FoldChange)),]) 


nrow(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$chr=="chr21" & common_fullresdata_rna$log2FoldChange<0 & common_fullresdata_rna$baseMean>=60 & common_fullresdata_rna$padj<.01 & !(is.na(common_fullresdata_rna$padj)),])

nrow(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$log2FoldChange<.48 & common_fullresdata_rna$chr=="chr21",])


nrow(common_fullresdata_rna[common_fullresdata_rna$padj<.01,])

medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians
colnames(fullresdata)

signif_medresdata <- (common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$padj<.01,])

notsignif_medresdata <- (common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$padj>=.01,])
nrow(signif_medresdata)
nrow(notsignif_medresdata)
ggplot() + 
  geom_violin(data=fullresdata,trim=TRUE,aes(x=chr, y=log2FoldChange)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  #  geom_hline(yintercept=log2(0.66667),linetype="dashed",color="red") +
  geom_hline(yintercept=0) +
  ylab(paste0("Log2FC")) +
  xlab("Chromosome") +
  theme_classic(base_size = 24) +
  geom_jitter(data = signif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=1.5, alpha=0.85,color="red",fill="red") +
  geom_jitter(data = notsignif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=1.5, alpha=0.85,color="black",fill="black") +
  stat_summary(data=fullresdata,aes(x=chr, y=log2FoldChange),fun.y=median, geom="crossbar", size=0.25, color="orange") +
  theme(legend.title = element_blank())

ggsave("rnaseq_changed.svg",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "svg")
ggsave("rnaseq_changed.png",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "png")



ddsFull <- DESeq(ddsFull)

person1 = "Ethan"
person2 = "Eric"
RNAddsres<-results(ddsFull,  contrast=c("Person", person1, person2))
resdata <- as.data.frame(RNAddsres)
fullresdata <- merge(resdata,annotationnew,by.x=0,by.y=4)
fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]
medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]

medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians

masterannotationdf_only21 <- masterannotationdf %>% filter(chr=="chr21")
head(masterannotationdf_only21)
anuplodygenes <-masterannotationdf_only21[["name"]] #Only chr21 genes
anuplodygenes
baseploidy <- 2
alt_ploidy <-3


controlgenes <-ifelse(rownames(dds) %in% anuplodygenes, FALSE, TRUE) # Only normalize on non-aneuploidy genes
controlgenes
samplepinfo<-as.data.frame(colData(dds))

ploidy_fortrisomy = ifelse(rownames(dds) %in% anuplodygenes, alt_ploidy, baseploidy)
ploidy_fortrisomy
ploidy_typical <- rep(baseploidy, nrow(dds))

dds$ploidy[11]
ploidy_typical
View(mget(paste0(dds$ploidy[1])))

normFactors <- matrix(do.call(cbind, mget(paste0(dds$ploidy))),ncol=ncol(dds),nrow=nrow(dds),dimnames=list(1:nrow(dds),1:ncol(dds)))
normFactors <- normFactors/baseploidy
unique(normFactors)

ddsCollapsed<-DESeq(dds,betaPrior = FALSE)
ddsCollapsed_normfactor <-estimateSizeFactors(dds,normMatrix=normFactors, controlGenes=controlgenes) #adds a column to colData(ddsCollapsed) called sizeFactor if you use normmatrix then normalizationFactors(ddsCollapsed) instead.

### Correcting for ploidy should yield a distribution tightly around 1.0
ddsCollapsed_normfactor <- estimateDispersionsGeneEst(ddsCollapsed_normfactor) #adds mu to assays(ddsCollapsed) mu-has one column per sample and fills in 
#elementMetadata(ddsCollapsed) with baseMean, baseVar, allZero, dispGeneEst, dispGeneIter
ddsCollapsed_normfactor<-estimateDispersionsFit(ddsCollapsed_normfactor) # and fills in dispersionFunction(ddsCollapsed) and adds a column in  elementMetadata(ddsCollapsed) called dispFit
ddsCollapsed_normfactor <- estimateDispersionsMAP(ddsCollapsed_normfactor) #adds a attr(,"dispPriorVar") to dispersionFunction(ddsCollapsed) 
#and adds 4 columns to elementMetadata(ddsCollapsed): dispersion, dispIter, dispOutlier, dispMAP
ddsCollapsed_normfactor <- nbinomWaldTest(ddsCollapsed_normfactor,betaPrior = FALSE) #adds both H and cooks to assays(ddsCollapsed)  both are for each sample you have.
#Adds a ton of columns to elementMetadata(ddsCollapsed) 
#including Intercept and all values they expect you to compare 
#SE_Intercept (SE stands for standard error)  and SE for every thing they expect you to compare. 
#Wald static for what they expect you to compare. 
#betaConv  betaIter         deviance  maxCooks

person1 = "Ethan"
person2 = "Eric"
RNAddsres<-results(ddsCollapsed_normfactor,  contrast=c("Person", person1, person2))
resdata <- as.data.frame(RNAddsres)
fullresdata <- merge(resdata,annotationnew,by.x=0,by.y=4)
fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]
quantile(fullresdata$baseMean, probs = seq(0,1,0.2))
fullresdata <- fullresdata[fullresdata$baseMean > 60,] # Remove genes with very low expression

medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]

dosagelowgenes <- medresdata[medresdata$log2FoldChange < 0.13 & medresdata$chr=="chr21",]
common_fullresdata_rna <- merge(medresdata,common_ids,by.x=1,by.y=1)

nrow(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$padj<.01,])
nrow(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$chr=="chr21",])
nrow(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$padj<.01 & common_fullresdata_rna$chr=="chr21",])
nrow(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$log2FoldChange<.48 & common_fullresdata_rna$chr=="chr21",])


nrow(common_fullresdata_rna[common_fullresdata_rna$padj<.01,])

medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians
colnames(fullresdata)

signif_medresdata <- medresdata[medresdata$padj<.01,]
notsignif_medresdata <- medresdata[medresdata$padj>=.01,]
nrow(signif_medresdata)
nrow(notsignif_medresdata)
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

ggsave("rnaseq_unchanged.svg",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "svg")
ggsave("rnaseq_unchanged.png",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "png")

#### Change the null hypothesis- Real RNA data ####
minichrs <-c("chr1", "chr2", "chr3","chr4", "chr5", "chr6", "chr7", "chr8", "chr9","chr10","chr11", "chr12", "chr13","chr14", "chr15", "chr16", "chr17", "chr18", "chr19","chr20","chr21","chr22")

#this is loading the metadata
RNAmetadata=read.table("/Shares/down/mixed/RNAGRO01132021/scripts/Deseq2_analysis/RNAinfo.txt", sep="\t", header=TRUE)
filetable <- RNAmetadata

#this is where you would remove samples if they are bad
RNAcountdat <- read.csv("/scratch/Users/sahu0957/ds_normalization/RNA/nextflowhg38_20220214/counts/featurecounts_rna_rmsked.sense.txt",
                        sep="\t", skip=1)

nrow(RNAcountdat)

# filter to only canonical chrs
head(RNAcountdat)
#RNAcountdat <- RNAcountdat[RNAcountdat$Chr %in% minichrs,]

colnames(RNAcountdat)
# colnames are different due to my featurecounts script. fixing that here by removing the path. Kinda hacky, go remove in script!
names <- colnames(RNAcountdat)
names <- gsub(".*removed..","",names)
colnames(RNAcountdat) <- names
colnames(RNAcountdat)

# We've got 1300 duplicate entries... How was this filtered out before?
RNAcountdat <- RNAcountdat[!(duplicated(RNAcountdat$Geneid)),]
# Filter to genes common in both RNA and GRO seq

library(tibble)
annotationmerge
head(RNAcountdat)
library(dplyr)
RNAcountdat<- RNAcountdat %>%
  filter(Geneid %in% annotationmerge$GeneID) %>%
  arrange(Geneid) %>%
  column_to_rownames('Geneid')
dim(RNAcountdat)
RNAcountdat <- RNAcountdat[,-c(1:5)]
#run Deseq on RNA-seq with 21
colnames(RNAcountdat)
ddsFull <- DESeqDataSetFromMatrix(countData = RNAcountdat, colData = RNAmetadata, design = ~ biological_rep + Person) #use this for all samples
dds <- ddsFull

# You can see the results of unfiltered DESeq2 by uncommenting this, but we have more correcting to do first!
ddsFull <- DESeq(ddsFull,)

person1 = "Ethan"
person2 = "Eric"
?results
RNAddsres<-results(ddsFull,  contrast=c("Person", person1, person2),altHypothesis="lessAbs",alpha=.1,lfcThreshold=log2(1.5))
resdata <- as.data.frame(RNAddsres)
fullresdata <- merge(resdata,annotationnew,by.x=0,by.y=4)

# Make a violin
lesschrs

fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]

medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]

common_fullresdata_rna <- merge(medresdata,common_ids,by.x=1,by.y=1)
View(common_fullresdata_rna)
nrow(common_fullresdata_rna[common_fullresdata_rna$padj<.01,])

medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians
colnames(fullresdata)

signif_medresdata <- medresdata[medresdata$padj<.1,]
notsignif_medresdata <- medresdata[medresdata$padj>=.1,]
View(medresdata)
nrow(signif_medresdata)
nrow(notsignif_medresdata)
ggplot() + 
  geom_violin(data=fullresdata,trim=TRUE,aes(x=chr, y=log2FoldChange)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  #  geom_hline(yintercept=log2(0.66667),linetype="dashed",color="red") +
  geom_hline(yintercept=0) +
  ylab(paste0("Log2FC")) +
  geom_hline(yintercept = log2(1.5), linetype="dotted", color="blue" ) +
  xlab("Chromosome") +
  theme_classic(base_size = 24) +
  geom_jitter(data = notsignif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=1.5, alpha=0.85,color="black",fill="black") +
  geom_jitter(data = signif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=1.5, alpha=0.85,color="red",fill="red") +
  stat_summary(data=fullresdata,aes(x=chr, y=log2FoldChange),fun.y=median, geom="crossbar", size=0.25, color="orange") +
  theme(legend.title = element_blank())

ggsave("rnaseq_newalthypothesis_21signif.svg",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "svg")
ggsave("rnaseq_newalthypothesis_21signif.png",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "png")


###REPEAT-REMOVED DATA DONE!####
###COUNT ONLY EXONS RNA-seq REAL DATA####
# Okay, last little thing- for RNA-seq, we've been counting over whole genes to keep it consistent with
# GRO-seq. We should go ahead and just do exons to be consistent. Keep the other corrections too

RNAbed <-ori_worldbed
filetable <- RNAmetadata
minichrs <-c("chr1", "chr2", "chr3","chr4", "chr5", "chr6", "chr7", "chr8", "chr9","chr10","chr11", "chr12", "chr13","chr14", "chr15", "chr16", "chr17", "chr18", "chr19","chr20","chr21","chr22")

#this is loading the metadata
RNAmetadata=read.table("/Shares/down/mixed/RNAGRO01132021/scripts/Deseq2_analysis/RNAinfo.txt", sep="\t", header=TRUE)
RNAmetadata
#this is where you would remove samples if they are bad
RNAcountdat <- read.csv("/scratch/Users/sahu0957/ds_normalization/RNA/nextflowhg38_20220214/counts/featurecounts_rna_exonsonly_nomultis.sense.txt",
                        sep="\t", skip=1)

# filter to only canonical chrs
nrow(RNAcountdat)

RNAcountdat <- RNAcountdat[RNAcountdat$Chr %in% minichrs,]

colnames(RNAcountdat)
# colnames are different due to my featurecounts script. fixing that here by removing the path. Kinda hacky, go remove in script!
names <- colnames(RNAcountdat)
names <- gsub(".*removed..","",names)
colnames(RNAcountdat) <- names
colnames(RNAcountdat)
# Aggregate exons sums per Geneid
RNAcountdat <- aggregate(. ~ Geneid + Chr,data=RNAcountdat,FUN="sum")
# Remove duplicated entries (some NRs appear on multiple chromosomes due to bad annotation)
RNAcountdat <- RNAcountdat[!(duplicated(RNAcountdat$Geneid)),]

# Assign ploidy for genes
masterannotationdf <- RNAcountdat
masterannotationdf_only21 <- masterannotationdf %>% filter(Chr=="chr21")
head(masterannotationdf_only21)
anuplodygenes <-masterannotationdf_only21[["Geneid"]] #Only chr21 genes
baseploidy <- 2
alt_ploidy <-3
anuplodygenes
length(anuplodygenes)
# Assign and sort rownames
RNAcountdat<- RNAcountdat %>%
  arrange(Geneid) %>%
  column_to_rownames('Geneid')
dim(RNAcountdat)

RNAcountdat <- RNAcountdat[,-c(1:5)]
#run Deseq on RNA-seq with 21
colnames(RNAcountdat)
RNAmetadata
ddsFull <- DESeqDataSetFromMatrix(countData = RNAcountdat, colData = RNAmetadata, design = ~ biological_rep + Person) #use this for all samples
dds <- ddsFull

# You can see the results of unfiltered DESeq2 by uncommenting this, but we have more correcting to do first!
ddsFull <- DESeq(ddsFull)

person1 = "Ethan"
person2 = "Eric"
RNAddsres<-results(ddsFull,  contrast=c("Person", person1, person2))
resdata <- as.data.frame(RNAddsres)
fullresdata <- merge(resdata,annotationnew,by.x=0,by.y=4)
fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]
medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]

medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians


controlgenes <-ifelse(rownames(dds) %in% anuplodygenes, FALSE, TRUE) # Only normalize on non-aneuploidy genes

samplepinfo<-as.data.frame(colData(dds))
anuplodygenes
ploidy_fortrisomy
ploidy_fortrisomy = ifelse(rownames(dds) %in% anuplodygenes, alt_ploidy, baseploidy)
ploidy_typical <- rep(baseploidy, nrow(dds))

normFactors <- matrix(do.call(cbind, mget(paste0(dds$ploidy))),ncol=ncol(dds),nrow=nrow(dds),dimnames=list(1:nrow(dds),1:ncol(dds)))
normFactors <- normFactors/baseploidy
unique(normFactors)

ddsCollapsed<-DESeq(dds,betaPrior = FALSE)
ddsCollapsed_normfactor <-estimateSizeFactors(dds,normMatrix=normFactors, controlGenes=controlgenes) #adds a column to colData(ddsCollapsed) called sizeFactor if you use normmatrix then normalizationFactors(ddsCollapsed) instead.
ddsCollapsed_normfactor <- estimateDispersionsGeneEst(ddsCollapsed_normfactor) #adds mu to assays(ddsCollapsed) mu-has one column per sample and fills in 
#elementMetadata(ddsCollapsed) with baseMean, baseVar, allZero, dispGeneEst, dispGeneIter
ddsCollapsed_normfactor<-estimateDispersionsFit(ddsCollapsed_normfactor) # and fills in dispersionFunction(ddsCollapsed) and adds a column in  elementMetadata(ddsCollapsed) called dispFit
ddsCollapsed_normfactor <- estimateDispersionsMAP(ddsCollapsed_normfactor) #adds a attr(,"dispPriorVar") to dispersionFunction(ddsCollapsed) 
#and adds 4 columns to elementMetadata(ddsCollapsed): dispersion, dispIter, dispOutlier, dispMAP
ddsCollapsed_normfactor <- nbinomWaldTest(ddsCollapsed_normfactor,betaPrior = FALSE) #adds both H and cooks to assays(ddsCollapsed)  both are for each sample you have.
#Adds a ton of columns to elementMetadata(ddsCollapsed) 
#including Intercept and all values they expect you to compare 
#SE_Intercept (SE stands for standard error)  and SE for every thing they expect you to compare. 
#Wald static for what they expect you to compare. 
#betaConv  betaIter         deviance  maxCooks

person1 = "Ethan"
person2 = "Eric"
RNAddsres<-results(ddsCollapsed_normfactor,  contrast=c("Person", person1, person2))
resdata <- as.data.frame(RNAddsres)
fullresdata <- merge(resdata,annotationnew,by.x=0,by.y=4)
fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]
quantile(fullresdata$baseMean, probs = seq(0,1,0.2))
fullresdata <- fullresdata[fullresdata$baseMean > 60,] # Remove genes with very low expression

medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]
dosagelowgenes <- medresdata[medresdata$log2FoldChange < 0.13 & medresdata$chr=="chr21",]
dosagelowgenes <- dosagelowgenes[!(duplicated(dosagelowgenes$Row.names)),]
write.csv(dosagelowgenes,'/Users/sahu0957/ds_deseq_normalization/dosagecompensated_genes_norepeats_exonsonly_ethaneric.csv',sep = '\t')
medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians
colnames(fullresdata)
ggplot(fullresdata, aes(x=chr, y=log2FoldChange)) + 
  geom_violin(trim=TRUE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  #  geom_hline(yintercept=log2(0.66667),linetype="dashed",color="red") +
  geom_hline(yintercept=0) +
  ylab(paste0("Log2FC")) +
  geom_hline(yintercept = log2(1.5), linetype="dotted", color="blue" ) +
  xlab("Chromosome") +
  theme_classic(base_size = 24) +
  geom_jitter(shape=23, position=position_jitter(0.1,seed = 1),size=0.9, alpha=0.5) +
  theme(legend.title = element_blank())

### TRIAL WITH PUBLIC DATA ####
# Okay, time to bring it all together
RNAbed <-ori_worldbed
filetable <- RNAmetadata
minichrs <-c("chr1", "chr2", "chr3","chr4", "chr5", "chr6", "chr7", "chr8", "chr9","chr10","chr11", "chr12", "chr13","chr14", "chr15", "chr16", "chr17", "chr18", "chr19","chr20","chr21","chr22")

#Metadata for public run
RNAmetadata = data.frame(filename=c('SRR5874663.sorted.bam','SRR5874664.sorted.bam','SRR5874665.sorted.bam',
                                    'SRR5874666.sorted.bam','SRR5874667.sorted.bam','SRR5874668.sorted.bam'),
                         biological_rep=c("A","B","C","A","B","C"),
                         ploidy=c("ploidy_trisomy21","ploidy_trisomy21","ploidy_trisomy21","ploidy_typical","ploidy_typical","ploidy_typical"),
                         sample=c("DS","DS","DS","control","control","control"))
RNAmetadata
ploidy_trisomy21
# Okay, we'll have to make a new metadata table
#this is where you would remove samples if they are bad
RNAcountdat <- read.csv("/scratch/Users/sahu0957/ds_normalization/RNA/public_data/ipscs/nextflowhg38_nomultimaps_20220308/counts/featurecounts_exonsonly_nomultimaps.sense.txt",
                        sep="\t", skip=1)
head(RNAcountdat)
# filter to only canonical chrs
nrow(RNAcountdat)

RNAcountdat <- RNAcountdat[RNAcountdat$Chr %in% minichrs,]

colnames(RNAcountdat)
# colnames are different due to my featurecounts script. fixing that here by removing the path. Kinda hacky, go remove in script!
names <- colnames(RNAcountdat)
names <- gsub(".*SRR","SRR",names)
names
colnames(RNAcountdat) <- names
colnames(RNAcountdat)
# Aggregate exons sums per Geneid
RNAcountdat <- aggregate(. ~ Geneid + Chr,data=RNAcountdat,FUN="sum")
# Remove duplicated entries (some NRs appear on multiple chromosomes due to bad annotation)
RNAcountdat <- RNAcountdat[!(duplicated(RNAcountdat$Geneid)),]

# Assign ploidy for genes
masterannotationdf <- RNAcountdat
masterannotationdf_only21 <- masterannotationdf %>% filter(Chr=="chr21")
head(masterannotationdf_only21)
anuplodygenes <-masterannotationdf_only21[["Geneid"]] #Only chr21 genes
baseploidy <- 2
alt_ploidy <-3
anuplodygenes
length(anuplodygenes)
# Assign and sort rownames
RNAcountdat<- RNAcountdat %>%
  arrange(Geneid) %>%
  column_to_rownames('Geneid')
dim(RNAcountdat)

RNAcountdat <- RNAcountdat[,-c(1:5)]
#run Deseq on RNA-seq with 21
colnames(RNAcountdat)
RNAmetadata
ddsFull <- DESeqDataSetFromMatrix(countData = RNAcountdat, colData = RNAmetadata, design = ~ biological_rep + sample) #use this for all samples
dds <- ddsFull

# You can see the results of unfiltered DESeq2 by uncommenting this, but we have more correcting to do first!
ddsFull <- DESeq(ddsFull)

person1 = "DS"
person2 = "control"
RNAddsres<-results(ddsFull,  contrast=c("sample", person1, person2))
resdata <- as.data.frame(RNAddsres)
fullresdata <- merge(resdata,annotationnew,by.x=0,by.y=4)
fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]
medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]

medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians


controlgenes <-ifelse(rownames(dds) %in% anuplodygenes, FALSE, TRUE) # Only normalize on non-aneuploidy genes

samplepinfo<-as.data.frame(colData(dds))
anuplodygenes
ploidy_fortrisomy
ploidy_trisomy21
ploidy_fortrisomy = ifelse(rownames(dds) %in% anuplodygenes, alt_ploidy, baseploidy)
ploidy_typical <- rep(baseploidy, nrow(dds))

normFactors <- matrix(do.call(cbind, mget(paste0(dds$ploidy))),ncol=ncol(dds),nrow=nrow(dds),dimnames=list(1:nrow(dds),1:ncol(dds)))
normFactors <- normFactors/baseploidy
unique(normFactors)

ddsCollapsed<-DESeq(dds,betaPrior = FALSE)
ddsCollapsed_normfactor <-estimateSizeFactors(dds,normMatrix=normFactors, controlGenes=controlgenes) #adds a column to colData(ddsCollapsed) called sizeFactor if you use normmatrix then normalizationFactors(ddsCollapsed) instead.
ddsCollapsed_normfactor <- estimateDispersionsGeneEst(ddsCollapsed_normfactor) #adds mu to assays(ddsCollapsed) mu-has one column per sample and fills in 
#elementMetadata(ddsCollapsed) with baseMean, baseVar, allZero, dispGeneEst, dispGeneIter
ddsCollapsed_normfactor<-estimateDispersionsFit(ddsCollapsed_normfactor) # and fills in dispersionFunction(ddsCollapsed) and adds a column in  elementMetadata(ddsCollapsed) called dispFit
ddsCollapsed_normfactor <- estimateDispersionsMAP(ddsCollapsed_normfactor) #adds a attr(,"dispPriorVar") to dispersionFunction(ddsCollapsed) 
#and adds 4 columns to elementMetadata(ddsCollapsed): dispersion, dispIter, dispOutlier, dispMAP
ddsCollapsed_normfactor <- nbinomWaldTest(ddsCollapsed_normfactor,betaPrior = FALSE) #adds both H and cooks to assays(ddsCollapsed)  both are for each sample you have.
#Adds a ton of columns to elementMetadata(ddsCollapsed) 
#including Intercept and all values they expect you to compare 
#SE_Intercept (SE stands for standard error)  and SE for every thing they expect you to compare. 
#Wald static for what they expect you to compare. 
#betaConv  betaIter         deviance  maxCooks

person1 = "DS"
person2 = "control"
RNAddsres<-results(ddsCollapsed_normfactor,  contrast=c("sample", person1, person2))
resdata <- as.data.frame(RNAddsres)
fullresdata <- merge(resdata,annotationnew,by.x=0,by.y=4)
fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]
quantile(fullresdata$baseMean, probs = seq(0,1,0.2))
fullresdata <- fullresdata[fullresdata$baseMean > 200,] # Remove genes with very low expression

medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]
dosagelowgenes <- medresdata[medresdata$log2FoldChange < 0.13 & medresdata$chr=="chr21",]
nrow(dosagelowgenes)
head(dosagelowgenes)
dosagelowgenes <- dosagelowgenes[!(duplicated(dosagelowgenes$Row.names)),]
head(RNAcountdat)
#write.csv(dosagelowgenes,'/Users/sahu0957/ds_deseq_normalization/dosagecompensated_genes_norepeats_exonsonly_ethaneric.csv',sep = '\t')
#View(dosagelowgenes)
medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians
colnames(fullresdata)
ggplot(fullresdata, aes(x=chr, y=log2FoldChange)) + 
  geom_violin(trim=TRUE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  #  geom_hline(yintercept=log2(0.66667),linetype="dashed",color="red") +
  geom_hline(yintercept=0) +
  ylab(paste0("Log2FC")) +
  geom_hline(yintercept = log2(1.5), linetype="dotted", color="blue" ) +
  xlab("Chromosome") +
  theme_classic(base_size = 24) +
  geom_jitter(shape=23, position=position_jitter(0.1,seed = 1),size=0.9, alpha=0.5) +
  theme(legend.title = element_blank())

# Didn't really improve MFC values, but definitely removed false positives!

#### START CORRECTING SIM DATA####

# Again, deliberately picking a real noisy file to try this out on.
RNAcountdat <- read.csv("/Users/sahu0957/trisomy_normalization/RNAGRO01132021/data/RNAseqfiles/NEWSCALEreplicatefinalsims_a0.075_p1_r1_s1_n3.csv",
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

# We have to make fake metadata here for the design matrix. Note that sims incorporate NO batch effects!
for (i in 1:(ncol(RNAcountdat)/2)){
  newmetadata1 <- c(paste0("Elvis_rep",i,"_counts"),"0",as.character(i),"E","Elvis","ploidy_typical",paste0("Elvis_",i))
  RNAmetadata[nrow(RNAmetadata)+1,] <- newmetadata1
}

for (i in 1:(ncol(RNAcountdat)/2)){
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

person1 = "Eddie" #T21
person2 = "Elvis" #D21
RNAddsres<-results(ddsnocorr,  contrast=c("Person", person1, person2))
resdata <- as.data.frame(RNAddsres)
fullresdata <- merge(resdata,annotationnew,by.x=0,by.y=4)
#fullresdata <- fullresdata[fullresdata$chr %in% minichrs,]
#medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]

common_fullresdata_rna <- merge(fullresdata,common_ids,by.x=1,by.y=1)
log2(1.5)
log2(1.3)
log2(1.7)
nrow(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$chr=="chr21",])
nrow(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$chr=="chr21" & common_fullresdata_rna$log2FoldChange<0.58,])

nrow(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$chr=="chr21" & common_fullresdata_rna$log2FoldChange>0.3785 & common_fullresdata_rna$log2FoldChange<0.7655,])
nrow(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$chr=="chr21" & is.na(common_fullresdata_rna$log2FoldChange),])
#View(common_fullresdata_rna)
nrow(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$chr=="chr21" & !(is.na(common_fullresdata_rna$padj)),])
(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$padj<.01 & common_fullresdata_rna$chr=="chr21" & !(is.na(common_fullresdata_rna$padj)) & common_fullresdata_rna$baseMean>60,])
nrow(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$chr=="chr21" & common_fullresdata_rna$baseMean>60,])
nrow(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$chr=="chr21" & common_fullresdata_rna$log2FoldChange>0 & !(is.na(common_fullresdata_rna$log2FoldChange)),]) 


View(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$chr=="chr21" & common_fullresdata_rna$log2FoldChange<0 & common_fullresdata_rna$baseMean>60 & common_fullresdata_rna$padj>=.01 & !(is.na(common_fullresdata_rna$padj)),])

nrow(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$log2FoldChange<.48 & common_fullresdata_rna$chr=="chr21",])


nrow(common_fullresdata_rna[common_fullresdata_rna$padj<.01,])






medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians


controlgenes <-ifelse(rownames(dds) %in% anuplodygenes, FALSE, TRUE) # Only normalize on non-aneuploidy genes

samplepinfo<-as.data.frame(colData(dds))
ploidy_trisomy21 = ifelse(rownames(dds) %in% anuplodygenes, alt_ploidy, baseploidy)
ploidy_typical <- rep(baseploidy, nrow(dds))
normFactors <- matrix(do.call(cbind, mget(paste0(dds$ploidy))),ncol=ncol(dds),nrow=nrow(dds),dimnames=list(1:nrow(dds),1:ncol(dds)))
normFactors <- normFactors/baseploidy
unique(normFactors)


ddsCollapsed<-DESeq(dds)
ddsCollapsed_normfactor <-estimateSizeFactors(dds,normMatrix=normFactors, controlGenes=controlgenes) #adds a column to colData(ddsCollapsed) called sizeFactor if you use normmatrix then normalizationFactors(ddsCollapsed) instead.
ddsCollapsed_normfactor <- estimateDispersionsGeneEst(ddsCollapsed_normfactor) #adds mu to assays(ddsCollapsed) mu-has one column per sample and fills in 
#elementMetadata(ddsCollapsed) with baseMean, baseVar, allZero, dispGeneEst, dispGeneIter
ddsCollapsed_normfactor<-estimateDispersionsFit(ddsCollapsed_normfactor) # and fills in dispersionFunction(ddsCollapsed) and adds a column in  elementMetadata(ddsCollapsed) called dispFit
ddsCollapsed_normfactor <- estimateDispersionsMAP(ddsCollapsed_normfactor) #adds a attr(,"dispPriorVar") to dispersionFunction(ddsCollapsed) 
#and adds 4 columns to elementMetadata(ddsCollapsed): dispersion, dispIter, dispOutlier, dispMAP
ddsCollapsed_normfactor <- nbinomWaldTest(ddsCollapsed_normfactor) #adds both H and cooks to assays(ddsCollapsed)  both are for each sample you have.

person1 = "Eddie"
person2 = "Elvis"
RNAddsres<-results(ddsCollapsed_normfactor,  contrast=c("Person", person1, person2))
resdata <- as.data.frame(RNAddsres)
fullresdata <- merge(resdata,annotationnew,by.x=0,by.y=4)
fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]
quantile(fullresdata$baseMean, probs = seq(0,1,0.1))
fullresdata <- fullresdata[fullresdata$baseMean > 20,] # Remove genes with very low expression

medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]

medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians

nrow(fullresdata)
nrow(fullresdata[fullresdata$padj<.01,])
nrow(fullresdata[fullresdata$padj<.01 & fullresdata$chr=="chr21",])
nrow(fullresdata[fullresdata$chr=="chr21",])



medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]
medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians
nrow(medresdata[medresdata$chr=="chr21" & medresdata$log2FoldChange<0,])
signif_medresdata <- medresdata[medresdata$padj<.01,]
notsignif_medresdata <- medresdata[medresdata$padj>=.01,]
nrow(medresdata[medresdata$chr=="chr21",])
nrow(signif_medresdata[signif_medresdata$chr=="chr21",])
nrow(notsignif_medresdata)
ggplot() + 
  geom_violin(data=fullresdata,trim=TRUE,aes(x=chr, y=log2FoldChange)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  #  geom_hline(yintercept=log2(0.66667),linetype="dashed",color="red") +
  geom_hline(yintercept=0) +
  ylab(paste0("Log2FC")) +
  xlab("Chromosome") +
  theme_classic(base_size = 24) +
  geom_jitter(data = signif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=1.5, alpha=0.85,color="red",fill="red") +
  geom_jitter(data = notsignif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=1.5, alpha=0.85,color="black",fill="black") +
  stat_summary(data=fullresdata,aes(x=chr, y=log2FoldChange),fun.y=median, geom="crossbar", size=0.25, color="orange") +
  theme(legend.title = element_blank())

ggsave("simseq_changed.svg",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "svg")
ggsave("simseq_changed.png",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "png")

common_fullresdata_rna <- merge(fullresdata,common_ids,by.x=1,by.y=1)
log2(1.5)
nrow(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$chr=="chr21" & common_fullresdata_rna$log2FoldChange<0.58,])
nrow(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$chr=="chr21" & is.na(common_fullresdata_rna$log2FoldChange),])

nrow(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$chr=="chr21" & !(is.na(common_fullresdata_rna$padj)),])
(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$padj<.01 & common_fullresdata_rna$chr=="chr21" & !(is.na(common_fullresdata_rna$padj)) & common_fullresdata_rna$baseMean>60,])
nrow(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$chr=="chr21" & common_fullresdata_rna$baseMean>60,])
nrow(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$chr=="chr21" & common_fullresdata_rna$log2FoldChange>0 & !(is.na(common_fullresdata_rna$log2FoldChange)),]) 


View(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$chr=="chr21" & common_fullresdata_rna$log2FoldChange<0 & common_fullresdata_rna$baseMean>60 & common_fullresdata_rna$padj>=.01 & !(is.na(common_fullresdata_rna$padj)),])

nrow(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$log2FoldChange<.48 & common_fullresdata_rna$chr=="chr21",])


nrow(common_fullresdata_rna[common_fullresdata_rna$padj<.01,])

####RNA-seq alternative hypothesis correction####
RNAcountdat <- read.csv("/Users/sahu0957/trisomy_normalization/RNAGRO01132021/data/RNAseqfiles/AdjustAll_replicatefinalsims_a0.03_p8.0_r1_s1_n3.csv", 
                        sep=",", row.names=1)
head(RNAcountdat)  
#Remove genes that are not in annotation (which was filtered to remove genes not in both and not on main chromosome above) 
#put the data frames in the same order
annotationmerge
library(dplyr)
RNAcountdat<- RNAcountdat %>%
  filter(Gene_id %in% annotationmerge$name) %>%
  arrange(Gene_id)  %>%
  column_to_rownames('Gene_id')
RNAcountdat
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
RNAdds <-DESeq(RNAdds)

person1 = "Eddie"
person2 = "Elvis"
RNAddsres<-results(RNAdds,  contrast=c("Person", person1, person2),lfcThreshold = log2(1.5),altHypothesis = "lessAbs")
resdata <- as.data.frame(RNAddsres)
fullresdata <- merge(resdata,annotationnew,by.x=0,by.y=4)
fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]
medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]
medresdata
medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians

signif_medresdata <- medresdata[medresdata$padj<.1,]
notsignif_medresdata <- medresdata[medresdata$padj>=.1,]
View(medresdata)
nrow(signif_medresdata)
nrow(notsignif_medresdata)
ggplot() + 
  geom_violin(data=fullresdata,trim=TRUE,aes(x=chr, y=log2FoldChange)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  #  geom_hline(yintercept=log2(0.66667),linetype="dashed",color="red") +
  geom_hline(yintercept=0) +
  ylab(paste0("Log2FC")) +
#  ylim() +
  geom_hline(yintercept = log2(1.5), linetype="dotted", color="blue" ) +
  xlab("Chromosome") +
  theme_classic(base_size = 24) +
  geom_jitter(data = notsignif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=1.5, alpha=0.85,color="black",fill="black") +
  geom_jitter(data = signif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=1.5, alpha=0.85,color="red",fill="red") +
  stat_summary(data=fullresdata,aes(x=chr, y=log2FoldChange),fun.y=median, geom="crossbar", size=0.25, color="orange") +
  theme(legend.title = element_blank())

ggsave("rnasims_adjustalthypothesis_0signif.svg",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "svg")
ggsave("rnasims_adjustalthypothesis_0signif.png",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "png")

####SIM CORRECTION DONE!####

#There is an increase in varience when chr21 included...
#two things seem different in chr21 removal vs. shuffles asymptDisp and dispPriorVar
#dispersion = asymptDisp + extraPois / mean
#So as your samples have more dispersion, the FC estimate gets pushed towards 0. Our sims are using a parameterized fit for dispersion estimate,
#so all chr21 genes have higher dispersion, but I think it's conceivable that some oddballs are just driving down the median FC estimate too
#Worth exploring!



### Various trial things here... looking at hyperparameters of the data

RNAcountdatdropchr21 <-RNAcountdat %>%
  rownames_to_column('gene') %>%
  filter(gene %in% annotationmergeno21$name) %>%
  arrange(gene) %>%
  column_to_rownames('gene')

RNAcountdatdropchr21shuffle <-RNAcountdat %>%
  rownames_to_column('gene') %>%
  filter(gene %in% annotationmergeno21shuffle$name) %>%
  arrange(gene) %>%
  column_to_rownames('gene')

RNAcountdat <- read.csv("/Users/sahu0957/trisomy_normalization/RNAGRO01132021/data/RNAseqfiles/res_featureCounts_gene_idfull_143138.coverage.csv", 
                        sep=",", row.names=1)
head(RNAcountdat)
RNAcountdat$Gene_id <- rownames(RNAcountdat)
RNAcountdat<- RNAcountdat %>%
  filter(Gene_id %in% annotationmerge$name) %>%
  arrange(Gene_id) %>%
  column_to_rownames('Gene_id')
nrow(RNAcountdat)
nrow(chr21genes)
#View(chr21genes)
RNAmetadata=read.table("/Shares/down/mixed/RNAGRO01132021/scripts/Deseq2_analysis/RNAinfo.txt", sep="\t", header=TRUE)
RNAmetadata$samplegroup <-paste(RNAmetadata$Person, "_", RNAmetadata$biological_rep, sep="")
head(RNAcountdat)




RNAdds <- DESeqDataSetFromMatrix(countData = RNAcountdat, colData = RNAmetadata, design = ~ biological_rep+ Person)
RNAdds <-DESeq(RNAdds)

RNAdds_no21dds <- DESeqDataSetFromMatrix(countData = RNAcountdatdropchr21, colData = RNAmetadata, design = ~ biological_rep+ Person)
RNAdds_no21dds <-DESeq(RNAdds_no21dds)

RNAdds_shuffleremove21info <- function(whichshuffle) {
  RNAcountdatdrop_sample <- RNAcountdat[sample(nrow(RNAcountdat), nrow(RNAcountdatdropchr21)), ]
  RNAdds_shuffleno21 <- DESeqDataSetFromMatrix(countData = RNAcountdatdrop_sample, colData = RNAmetadata, design = ~ biological_rep+ Person)
  RNAdds_shuffleno21dds <-estimateSizeFactors(RNAdds_shuffleno21) 
  RNAdds_shuffleno21dds <- estimateDispersionsGeneEst(RNAdds_shuffleno21dds) 
  RNAdds_shuffleno21dds<-estimateDispersionsFit(RNAdds_shuffleno21dds) 
  RNAdds_shuffleno21dds <- estimateDispersionsMAP(RNAdds_shuffleno21dds) 
  RNAdds_shuffleno21dds <- nbinomWaldTest(RNAdds_shuffleno21dds)
  RNAdds_shuffleno21dds@dispersionFunction
  colnames(disp_mat)=c("asymptDisp", "extraPois", "varLogDispEsts", "dispPriorVar")
  info <-c(attributes(RNAdds_shuffleno21dds@dispersionFunction)$coefficients[1], attributes(RNAdds_shuffleno21dds@dispersionFunction)$coefficients[2], attributes(RNAdds_no21dds@dispersionFunction)$varLogDispEsts, attributes(RNAdds_shuffleno21dds@dispersionFunction)$dispPriorVar)
  return(info)}

disp_mat <- as.data.frame(matrix(, nrow = 0, ncol = 4))
disp_mat
colnames(disp_mat)=c("asymptDisp", "extraPois", "varLogDispEsts", "dispPriorVar")
disp_mat
disp_mat
for (i in 1:25)
{rdf=as.data.frame(rbind(RNAdds_shuffleremove21info(i)))
colnames(rdf)=c("asymptDisp", "extraPois", "varLogDispEsts", "dispPriorVar")
disp_mat = rbind(disp_mat, rdf)}

disp_mat
disp_mat$sample <-"shuffled removal chr21 genes"
disp_mat
realdispinfo = c(attributes(RNAdds@dispersionFunction)$coefficients[1], attributes(RNAdds@dispersionFunction)$coefficients[2], attributes(RNAdds@dispersionFunction)$varLogDispEsts, attributes(RNAdds@dispersionFunction)$dispPriorVar, "real")
disp_nochr21info = c(attributes(RNAdds_no21dds@dispersionFunction)$coefficients[1], attributes(RNAdds_no21dds@dispersionFunction)$coefficients[2], attributes(RNAdds_no21dds@dispersionFunction)$varLogDispEsts, attributes(RNAdds_no21dds@dispersionFunction)$dispPriorVar, "no chr21")
rdf = as.data.frame(rbind(realdispinfo, disp_nochr21info))
rdf
colnames(rdf)=c("asymptDisp", "extraPois", "varLogDispEsts", "dispPriorVar", "sample")
disp_mat <-rbind(disp_mat, rdf)
disp_mat
disp_mat$asymptDisp <- as.numeric(disp_mat$asymptDisp)
disp_mat$extraPois <- as.numeric(disp_mat$extraPois)
disp_mat$varLogDispEsts <- as.numeric(disp_mat$varLogDispEsts)
disp_mat$dispPriorVar <- as.numeric(disp_mat$dispPriorVar)

minidisp_mat <- disp_mat %>% filter(sample=="shuffled removal chr21 genes")
realdisp_mat <-disp_mat %>% filter(sample=="real")
nochr21disp_mat <-disp_mat %>% filter(sample=="no chr21")

disp_mat

ggplot(disp_mat, aes(x=sample, y=asymptDisp)) + 
  geom_boxplot() +
  theme_classic(base_size=18)

ggsave("rnaseq_asymptdisp.svg",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "svg")
ggsave("rnaseq_asymptdisp.png",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "png")


ggplot(disp_mat, aes(x=sample, y=extraPois)) + 
  geom_boxplot() +
  theme_classic(base_size=18)

ggsave("rnaseq_extrapois.svg",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "svg")
ggsave("rnaseq_extrapois.png",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "png")

ggplot(disp_mat, aes(x=sample, y=varLogDispEsts)) + 
  geom_boxplot() +
  theme_classic(base_size=18)

ggplot(disp_mat, aes(x=sample, y=dispPriorVar)) + 
  geom_boxplot() +
  theme_classic(base_size=18)

#There is an increase in varience when chr21 included...
#the results of this change are that all of the dispersion estimates are affected by the changes on chr21
emd_RNAdds
emd_RNAdds <-elementMetadata(RNAdds)
emd_RNAdds
rownames(emd_RNAdds) <- rownames(RNAdds)
emd_RNAdds_no21dds <-elementMetadata(RNAdds_no21dds)
rownames(emd_RNAdds_no21dds) <- rownames(RNAdds_no21dds)
emd_RNAdds$names <-rownames(emd_RNAdds)
emd <- merge(emd_RNAdds, emd_RNAdds_no21dds, all = FALSE,by="row.names", suffixes = c(".alldata",".nochr21"))
emd
rownames(emd)<-emd$names
emd
head(annotationmerge)
emd <- merge(annotationmerge, emd, all = FALSE, by.x=1,by.y=0)
emd
rownames(emd) <- emd$GeneID
emd <-as.data.frame(emd)
View(emd)
colnames(emd)
ggplot(emd, aes(x=dispGeneEst.nochr21, y=dispGeneEst.alldata)) + 
  geom_point() +
  geom_abline() +
  theme_classic()
emd

emd


person1 = "Ethan"
person2 = "Eric"

RNAseq_alteredgenes <- function(person1, person2){
  RNAdds <- DESeqDataSetFromMatrix(countData = RNAcountdat, colData = RNAmetadata, design = ~ biological_rep+ Person)
  RNAdds <-DESeq(RNAdds)
  RNAddsres<-results(RNAdds,  contrast=c("Person", person1, person2))
  RNAdds_no21dds <- DESeqDataSetFromMatrix(countData = RNAcountdatdropchr21, colData = RNAmetadata, design = ~ biological_rep+ Person)
  RNAdds_no21dds <-DESeq(RNAdds_no21dds)
  RNAdds_no21ddsres<-results(RNAdds_no21dds,contrast=c("Person", person1, person2))
  RNAdds_no21ddsshuffle <- DESeqDataSetFromMatrix(countData = RNAcountdatdropchr21shuffle, colData = RNAmetadata, design = ~ biological_rep+ Person)
  RNAdds_no21ddsshuffle <-DESeq(RNAdds_no21dds)
  RNAdds_no21ddsshuffleres<-results(RNAdds_no21ddsshuffle,contrast=c("Person", person1, person2))
  #filter all the results for genes that are in all three versions
  # world<-
  RNAddsresSig <- subset(RNAddsres[order(RNAddsres$pvalue),], padj < 0.1)
  RNAdds_no21ddsresSig <- subset(RNAdds_no21ddsres[order(RNAdds_no21ddsres$pvalue),], padj < 0.1)
  
  RNAdds_no21ddsresSig <- as.data.frame(RNAdds_no21ddsresSig)
  RNAdds_no21ddsresSig_up <- RNAdds_no21ddsresSig %>% filter(log2FoldChange>0)
  RNAdds_no21ddsresSig_down <- RNAdds_no21ddsresSig %>% filter(log2FoldChange<=0)
  RNAddsresSig_no21_up <- RNAddsresSig_no21 %>% filter(log2FoldChange>0)
  RNAddsresSig_no21_down <- RNAddsresSig_no21 %>% filter(log2FoldChange<=0)
  info1 = c(nrow(RNAdds_no21ddsresSig_up),nrow(RNAddsresSig_no21_up))
  info2 = c(nrow(RNAdds_no21ddsresSig_down), nrow(RNAddsresSig_no21_down))
  rdfup=as.data.frame(rbind(info1))
  colnames(rdfup)=c("RNAseq_no_21", "RNAseq")
  rownames(rdfup)=c(paste(person1, "_vs_", person2,"__up" ,sep=""))
  rdfdown=as.data.frame(rbind(info2))
  colnames(rdfdown)=c("RNAseq_no_21", "RNAseq")
  rdfdown$RNAseq_no_21 = -1*rdfdown$RNAseq_no_21
  rdfdown$RNAseq = -1*rdfdown$RNAseq
  rownames(rdfdown)=c(paste(person1, "_vs_", person2,"__down" ,sep=""))
  rdf=rbind(rdfup, rdfdown)
  return(rdf)}

diffexpdf <- rbind(RNAseq_alteredgenes("Eli", "Eric"), RNAseq_alteredgenes("Elizabeth", "Eric"), RNAseq_alteredgenes("Eli", "Ethan"),RNAseq_alteredgenes("Elizabeth", "Ethan"), RNAseq_alteredgenes("Elizabeth", "Eli"), RNAseq_alteredgenes("Eric", "Ethan"),  
                   RNAseq_alteredgenes("Eddie", "Ethan")) 
diffexpdf$compairison <-rownames(diffexpdf)
diffexpdf$diff <- diffexpdf$RNAseq - diffexpdf$RNAseq_no_21
diffexpdf_long =  gather(diffexpdf, sample, ngenes, -compairison, -diff) %>%separate(compairison, into=c("comparison", "direction"), sep="__")
diffexpdf_long$comparison <- factor(diffexpdf_long$comparison, levels=c("Eli_vs_Eric", "Elizabeth_vs_Eric", "Elizabeth_vs_Eli", "Eli_vs_Ethan", "Elizabeth_vs_Ethan", "Eric_vs_Ethan","Eddie_vs_Ethan"))

ggplot(diffexpdf_long, aes(x=comparison, y=ngenes, color=sample)) + geom_point()+ theme(axis.text.x = element_text(angle = 90))


ggplot(diffexpdf_long, aes(x=comparison, y=diff, color=direction)) + geom_point()+ theme(axis.text.x = element_text(angle = 90))


person1 = "Eli"
person2 = "Eric"

RNAmetadata

RNAcountdat <- read.csv("/Users/sahu0957/trisomy_normalization/RNAGRO01132021/data/RNAseqfiles/res_featureCounts_gene_idfull_143138.coverage.csv", 
                        sep=",", row.names=1)
head(RNAcountdat)
RNAcountdat$Gene_id <- rownames(RNAcountdat)
RNAcountdat<- RNAcountdat %>%
  filter(Gene_id %in% annotationmerge$name) %>%
  arrange(Gene_id) %>%
  column_to_rownames('Gene_id')
nrow(RNAcountdat)
nrow(chr21genes)
#View(chr21genes)
RNAmetadata=read.table("/Shares/down/mixed/RNAGRO01132021/scripts/Deseq2_analysis/RNAinfo.txt", sep="\t", header=TRUE)
RNAmetadata$samplegroup <-paste(RNAmetadata$Person, "_", RNAmetadata$biological_rep, sep="")
head(RNAcountdat)

RNAtestcountdat <- RNAcountdat
# Inform DESeq2 of copy number

newRNAtestcount <- as.data.frame(cbind(rownames(RNAtestcountdat),0))
head(newRNAtestcount)

chr21genes <- na.omit(chr21genes)
chr21genes
for (i in colnames(RNAtestcountdat)){
  #i <- "Eddie_repA_tech1.sorted.bam"
  copynum <- ifelse(i %like% "Eddie" | i %like% "Ethan",1.5,1)
  copynum
  i
  tmpcol <- as.data.frame(cbind(rownames(RNAtestcountdat),RNAtestcountdat[,i]))
  tmpchr21 <- tmpcol[tmpcol$V1 %in% chr21genes$Row.names,]
  tmpnochr21 <- tmpcol[!(tmpcol$V1 %in% chr21genes$Row.names),]
  tmpchr21$V2 <- round(as.numeric(tmpchr21$V2)*copynum)
  tmpnochr21$V2 <- round(as.numeric(tmpnochr21$V2)*1)
  tmpdf <- rbind(tmpchr21,tmpnochr21)
  rownames(tmpdf) <- tmpdf$V1
  colnames(tmpdf) <- c("genes",i)
  newRNAtestcount <- merge(newRNAtestcount,tmpdf, by.x=1,by.y=1)
}
rownames(newRNAtestcount) <- newRNAtestcount$V1
head(newRNAtestcount)
RNAtestcountdat <- newRNAtestcount[,-c(1,2)]
head(RNAtestcountdat)
RNAtestcountdat <- RNAtestcountdat
#[,colnames(RNAtestcountdat) %like% person1 | colnames(RNAtestcountdat) %like% person2] 
RNAtestmetadata <- RNAmetadata
#[RNAmetadata$filename %like% person1 | RNAmetadata$filename %like% person2,]
RNAtestcountdatdropchr21 <- RNAcountdatdropchr21
#[,colnames(RNAcountdatdropchr21) %like% person1 | colnames(RNAcountdatdropchr21) %like% person2]
RNAtestcountdatdropchr21shuffle <- RNAcountdatdropchr21shuffle
#[,colnames(RNAcountdatdropchr21shuffle) %like% person1 | colnames(RNAcountdatdropchr21shuffle) %like% person2]
RNAddstest_no21dds <- DESeqDataSetFromMatrix(countData = RNAtestcountdatdropchr21, colData = RNAtestmetadata, design = ~ biological_rep+ Person)
RNAddstest_no21dds <-DESeq(RNAddstest_no21dds)
RNAddstest_no21ddsres<-results(RNAddstest_no21dds,contrast=c("Person", person1, person2))

#plotDispEsts(RNAddstest)
RNAddstest <- DESeqDataSetFromMatrix(countData = RNAtestcountdat, colData = RNAtestmetadata, design = ~ biological_rep+ Person)
RNAddstest <- estimateSizeFactors(RNAddstest)
RNAddstest <- estimateDispersionsGeneEst(RNAddstest)
mcols(RNAddstest)

estima
#RNAddstest <- DESeq(RNAddstest)
RNAddstest
count_for_dispersions_df <- as.data.frame(counts(RNAddstest, normalized=TRUE))
count_for_dispersions_df <- count_for_dispersions_df[,colnames(count_for_dispersions_df) %like% "Ethan" | colnames(count_for_dispersions_df) %like% "Eric"]
count_for_dispersions_df$mean <- rowMeans(count_for_dispersions_df)

dispersion_dataframe <- as.data.frame(cbind(rownames(RNAddstest),mcols(RNAddstest)))
head(count_for_dispersions_df)
head(dispersion_dataframe)
head(annotationmerge)
dispersion_final <- merge(dispersion_dataframe,annotationmerge,by.x=1,by.y=1)
dispersion_final2 <- merge(dispersion_final,count_for_dispersions_df,by.x=1,by.y=0)
dispersion_final2$dispersion <- as.numeric(as.character(dispersion_final2$dispGeneEst))
dispersion_final2 <- dispersion_final2[dispersion_final2$chr %in% lesschrs,]
View(dispersion_final2)
?estimateDispersions
lesschrs <- c("chr21","chr22")
ggplot(dispersion_final2, aes(x=chr, y=log(dispersion))) + 
  geom_violin(trim=TRUE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  #  geom_hline(yintercept=log2(0.66667),linetype="dashed",color="red") +
  geom_hline(yintercept=0) +
  ylab(paste0("MLE Gene-wise Dispersion")) +
  xlab("Chromosome") +
  theme_classic(base_size = 24) +
  geom_jitter(shape=23, position=position_jitter(0.1,seed = 1),size=0.9, alpha=0.5) +
  theme(legend.title = element_blank())
#View(fullresdata)

medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]
medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians




dispersion_final2$color <- ifelse(dispersion_final2$chr=="chr21","red",ifelse(dispersion_final2$chr=="chr20","green","black"))
dispersion_final2$V2 <- as.numeric(as.character(dispersion_final2$V2)) 

chr22dispersion <- dispersion_final2[dispersion_final2$chr=="chr22",]
chr21dispersion <- dispersion_final2[dispersion_final2$chr=="chr21",]





ggplot(dispersion_final2, aes(x=mean, y=V2)) + 
  scale_x_log10() +
  geom_point(color="gray",alpha=0.5) + ##change opacity of points +
  scale_y_log10() +
  geom_point(data=chr22dispersion, aes(x=mean,y=V2),color="blue",alpha=0.3) +
  geom_point(data=chr21dispersion, aes(x=mean,y=V2),color="#249c24",alpha=0.6) +
  geom_hline(yintercept=0.0235,linetype="dotted") +
  theme_classic(base_size = 18) + ##use the classic theme template 
  xlab("Mean of Normalized Counts") + ##label the x-axis
  ylab("Dispersion Estimate") + ##label the y-axis
  ggtitle("Dispersions") ##to add title

ggsave("dispersion_bychr.svg",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "svg")
ggsave("dispersion_bychr.png",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "png")

dev.copy(svg,'/Users/sahu0957/ds_deseq_normalization/dispersion_estimates_noadjustments.svg')
dev.off()

RNAcountdat <- read.csv("/Users/sahu0957/trisomy_normalization/RNAGRO01132021/data/RNAseqfiles/res_featureCounts_gene_idfull_143138.coverage.csv", 
                        sep=",", row.names=1)
head(RNAcountdat)
RNAcountdat$Gene_id <- rownames(RNAcountdat)
RNAcountdat<- RNAcountdat %>%
  filter(Gene_id %in% annotationmerge$name) %>%
  arrange(Gene_id) %>%
  column_to_rownames('Gene_id')
nrow(RNAcountdat)
nrow(chr21genes)
#View(chr21genes)
RNAmetadata=read.table("/Shares/down/mixed/RNAGRO01132021/scripts/Deseq2_analysis/RNAinfo.txt", sep="\t", header=TRUE)
RNAmetadata$samplegroup <-paste(RNAmetadata$Person, "_", RNAmetadata$biological_rep, sep="")
head(RNAcountdat)



colnames(RNAcountdat)
nrow(RNAcountdat)
#RNAcountdat <- RNAcountdat[!(rownames(RNAcountdat) %in% chr21genes$Row.names),]
rownames(RNAcountdat)
chr21genes
RNAcountdatnoT21 <- RNAcountdat[,!(colnames(RNAcountdat) %like% "Ethan")]

RNAmetadata <- RNAmetadata[RNAmetadata$Person!="Ethan",]
RNAcountdatnoT21
RNAddstestnoT21 <- DESeqDataSetFromMatrix(countData = RNAcountdatnoT21, colData = RNAmetadata, design = ~ biological_rep+ Person)




#sizeFactors(RNAddstest) <- sizeFactors(RNAddstest_no21dds)
RNAddstestnoT21 <-DESeq(RNAddstestnoT21)
noEthanDispersions <- dispersions(RNAddstestnoT21)


RNAcountdat <- read.csv("/Users/sahu0957/trisomy_normalization/RNAGRO01132021/data/RNAseqfiles/res_featureCounts_gene_idfull_143138.coverage.csv", 
                        sep=",", row.names=1)
head(RNAcountdat)
RNAcountdat$Gene_id <- rownames(RNAcountdat)
RNAcountdat<- RNAcountdat %>%
  filter(Gene_id %in% annotationmerge$name) %>%
  arrange(Gene_id) %>%
  column_to_rownames('Gene_id')
nrow(RNAcountdat)
nrow(chr21genes)
#View(chr21genes)
RNAmetadata=read.table("/Shares/down/mixed/RNAGRO01132021/scripts/Deseq2_analysis/RNAinfo.txt", sep="\t", header=TRUE)
RNAmetadata$samplegroup <-paste(RNAmetadata$Person, "_", RNAmetadata$biological_rep, sep="")
head(RNAcountdat)



RNAddstest <- DESeqDataSetFromMatrix(countData = RNAcountdat, colData = RNAmetadata, design = ~ biological_rep+ Person)
RNAddstest <- DESeq(RNAddstest)
#sizeFactors(RNAdds_no21dds)
sizeFactors(RNAddstest) <- sizeFactors(RNAdds_no21dds)
#RNAddstest <- estimateSizeFactors(RNAddstest)
RNAddstest <- estimateDispersions(RNAddstest)

#hmm <- replace_na(noEthanDispersions, .0001)
#dispersions(RNAddstest) <- hmm
#noEthanDispersions
RNAddstest <- nbinomWaldTest(RNAddstest)
person1="Ethan"
person2="Eric"

RNAddsres<-results(RNAddstest,  contrast=c("Person", person1, person2))
RNAddsresOrdered <- RNAddsres[order(RNAddsres$pvalue),]
RNAddsresSig <- subset(RNAddsresOrdered, padj < 0.1)
RNAddsresSig_no21 <- as.data.frame(RNAddsresSig) %>%   
  rownames_to_column('gene') %>%
  filter(gene %in% annotationmergeno21$name) %>%
  arrange(gene) %>%
  column_to_rownames('gene')

#View(resdata)


resdata <- as.data.frame(RNAddsres)
#resdata <- na.omit(resdata)

# Remove really low genes maybe?
#resdata <- resdata[resdata$baseMean > 20,]
#nrow(resdata)
signifGenes <- resdata[resdata$padj < 0.1,]
#signifGenes <- na.omit(signifGenes)

ggplot(resdata, aes(log(baseMean,2), log2FoldChange)) + 
  geom_point(color="gray20", alpha = 0.5) + ##change opacity of points
  theme_classic() + ##use the classic theme template 
  xlab("log2 Average Normalized Counts") + ##label the x-axis
  ylab("log2 Fold Change") + ##label the y-axis
  ggtitle("MA plot") + ##to add title
  theme(plot.title = element_text(hjust = 0.5),
        title = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 18,face = "bold"),
        axis.text = element_text(size = 16)) + ##to center title
  geom_hline(aes(yintercept=0), colour="blue", linetype="dashed") + ##plot a horizontal line
  geom_point(data=signifGenes,aes(log(baseMean,2), log2FoldChange), color="red",size=1, alpha=0.75)

rownames(annotationnew) <- annotationnew$name
nrow(annotationnew)
lesschrs <- c("chr21","chr22")
fullresdata <- merge(resdata,annotationnew,by.x=0,by.y=4)
nrow(fullresdata)
fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]
quantile(fullresdata$baseMean)
fullresdata <- fullresdata[fullresdata$baseMean > 5,]
sigresdata <- fullresdata[fullresdata$padj<0.1,]
sigres21 <- sigresdata[sigresdata$chr == "chr21",]
fullresdata$coloring <- "red"
fullresdata$coloring[fullresdata$padj<0.1] <- "red"
fullresdata$coloring[fullresdata$padj>=0.1] <- "black"

ggplot(fullresdata, aes(x=chr, y=log2FoldChange)) + 
  geom_violin(trim=TRUE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_hline(yintercept=log2(1.5),linetype="dashed",color="blue") +
  #  geom_hline(yintercept=log2(0.66667),linetype="dashed",color="red") +
  geom_hline(yintercept=0) +
  ylab(paste0("log2FC(T21/D21)")) +
  xlab("Chromosome") +
  theme_classic(base_size = 24) +
  geom_jitter(aes(color=coloring,fill=coloring),shape=23, position=position_jitter(0.1,seed = 1),size=0.9, alpha=0.5) +
  scale_colour_manual(values=c("black","black"), guide=FALSE) +
  scale_fill_manual(values=c("black","black"), guide = FALSE) +
  stat_summary(fun.y=median, geom="crossbar", size=0.25, color="orange") +
  #  stat_summary(fun.y=mean, geom="crossbar", size=0.25, color="green") +
  theme(legend.title = element_blank())
#View(fullresdata)
medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]
medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians

ggsave("rnaseq_deseq2_nocorrection_noshrinkage_med100med139.svg",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "svg")
ggsave("rnaseq_deseq2_nocorrection_noshrinkage_med100med139.png",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "png")




RNAtestcountdat <- RNAcountdat
# Inform DESeq2 of copy number

count_for_dispersions_df <- as.data.frame(counts(RNAddstest, normalized=TRUE))
#count_for_dispersions_df <- count_for_dispersions_df[,colnames(count_for_dispersions_df) %like% "Ethan" | colnames(count_for_dispersions_df) %like% "Eric"]
count_for_dispersions_df$mean <- rowMeans(count_for_dispersions_df)

dispersion_dataframe <- as.data.frame(cbind(rownames(RNAddstest),dispersions(RNAddstest)))
head(count_for_dispersions_df)
head(dispersion_dataframe)
head(annotationmerge)
dispersion_final <- merge(dispersion_dataframe,annotationmerge,by.x=1,by.y=1)
dispersion_final2 <- merge(dispersion_final,count_for_dispersions_df,by.x=1,by.y=0)
head(dispersion_final2)

dispersion_final2$color <- ifelse(dispersion_final2$chr=="chr21","red",ifelse(dispersion_final2$chr=="chr20","green","black"))
dispersion_final2$V2 <- as.numeric(as.character(dispersion_final2$V2)) 

chr22dispersion <- dispersion_final2[dispersion_final2$chr=="chr22",]
chr21dispersion <- dispersion_final2[dispersion_final2$chr=="chr21",]

RNAddstestres<-results(RNAddstest,  contrast=c("Person", person1, person2))

RNAddsres<-results(RNAddstest,  contrast=c("Person", person1, person2))
RNAddsresOrdered <- RNAddsres[order(RNAddsres$pvalue),]
RNAddsresSig <- subset(RNAddsresOrdered, padj < 0.1)
RNAddsresSig_no21 <- as.data.frame(RNAddsresSig) %>%   
  rownames_to_column('gene') %>%
  filter(gene %in% annotationmergeno21$name) %>%
  arrange(gene) %>%
  column_to_rownames('gene')


resdata <- as.data.frame(RNAddsres)
resdata <- na.omit(resdata)
signifGenes <- resdata[resdata$padj < 0.1,]
signifGenes <- na.omit(signifGenes)

ggplot(resdata, aes(log(baseMean,2), log2FoldChange)) + 
  geom_point(color="gray20", alpha = 0.5) + ##change opacity of points
  theme_classic() + ##use the classic theme template 
  xlab("log2 Average Normalized Counts") + ##label the x-axis
  ylab("log2 Fold Change") + ##label the y-axis
  ggtitle("MA plot") + ##to add title
  theme(plot.title = element_text(hjust = 0.5),
        title = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 18,face = "bold"),
        axis.text = element_text(size = 16)) + ##to center title
  geom_hline(aes(yintercept=0), colour="blue", linetype="dashed") + ##plot a horizontal line
  geom_point(data=signifGenes,aes(log(baseMean,2), log2FoldChange), color="red",size=1, alpha=0.75)

rownames(annotationnew) <- annotationnew$name

lesschrs <- c("chr21","chr22")
fullresdata <- merge(resdata,annotationnew,by.x=0,by.y=4)
fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]
fullresdata <- fullresdata[fullresdata$baseMean > 1,]
sigresdata <- fullresdata[fullresdata$padj<0.1,]
sigres21 <- sigresdata[sigresdata$chr == "chr21",]
fullresdata$coloring <- "red"
fullresdata$coloring[fullresdata$padj<0.1] <- "red"
fullresdata$coloring[fullresdata$padj>=0.1] <- "black"

ggplot(fullresdata, aes(x=chr, y=log2FoldChange)) + 
  geom_violin(trim=TRUE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_hline(yintercept=log2(1.5),linetype="dashed",color="blue") +
  #  geom_hline(yintercept=log2(0.66667),linetype="dashed",color="red") +
  geom_hline(yintercept=0) +
  ylab(paste0("log2FC(FD21/D21)")) +
  xlab("Chromosome") +
  ggtitle(paste0("RNAseq DESeq2 Results")) +
  theme_classic(base_size = 18) +
  geom_jitter(aes(color=coloring,fill=coloring),shape=23, position=position_jitter(0.1,seed = 1),size=0.9, alpha=0.5) +
  scale_colour_manual(values=c("black","black"), guide=FALSE) +
  scale_fill_manual(values=c("black","black"), guide = FALSE) +
  stat_summary(fun.y=median, geom="crossbar", size=0.25, color="orange") +
  #  stat_summary(fun.y=mean, geom="crossbar", size=0.25, color="green") +
  theme(legend.title = element_blank())

medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]
medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians
ggsave("rnaseq_deseq2_nocorrection.svg",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "svg")
ggsave("rnaseq_deseq2_nocorrection.png",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "png")


nrow(fullresdata[fullresdata$log2FoldChange <= 0.25 & fullresdata$chr == "chr21",])
nrow(fullresdata)



###### GRO SEQ ######## 
#### UNCORRECTED REAL DATA GRO-SEQ RUN####
#traditional Deseq2 on on GRO count matrix's

GROseq_indir<-"/Shares/down/mixed/RNAGRO01132021/data/GROfiles/"
GROmetadata=read.table("/Shares/down/mixed/RNAGRO01132021/scripts/Deseq2_analysis/GROinfo.txt", sep="\t", header=TRUE)
GROmetadata$samplegroup <-paste(GROmetadata$person, "_", GROmetadata$libprep, sep="")
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
annotationmerge <-merge(annotationRNA, annotationGRO, all = FALSE,by="row.names")
row.names(annotationmerge)<-annotationmerge$GeneID

#keep only the genes on the main chromosomes
minichrs <-c("chr1", "chr2", "chr3","chr4", "chr5", "chr6", "chr7", "chr8", "chr9","chr10","chr11", "chr12", "chr13","chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")

library(dplyr)

annotationmerge <-annotationmerge %>%
  arrange(GeneID) %>%
  filter(chr %in% minichrs)
#add a column that has TRUE if the gene is not on chr21 for use as a control gene
annotationmerge
chr21genes <- annotationmerge[annotationmerge$chr=="chr21",]
head(chr21genes)
annotationmerge$controlgenes <- TRUE
annotationmerge$controlgenes18 <- annotationmerge$chr!="chr18"
annotationmerge$controlgenes19 <- annotationmerge$chr!="chr19"
annotationmerge$controlgenes20 <- annotationmerge$chr!="chr20"
annotationmerge$controlgenes21 <- annotationmerge$chr!="chr21"
annotationmerge$controlgenes21shuffle <- sample(annotationmerge$chr!="chr21")


annotationmergeno21 <- annotationmerge %>% filter(controlgenes21)
annotationmergeno21shuffle <- annotationmerge %>% filter(controlgenes21shuffle)
GROseq_indir
#read in data for GRO-seq
GRObodycountdat<- read.csv(paste0(GROseq_indir, "body.csv", sep=""))
GROmetadata
GRObodyddsFull <- DESeqDataSetFromMatrix(countData = GRObodycountdat, colData = GROmetadata, design = ~libprep+ person)
GRObodydds <- collapseReplicates( GRObodyddsFull,groupby = GRObodyddsFull$samplegroup,run = GRObodyddsFull$samplegroup)
GRObodydds <-DESeq(GRObodydds)

person1 = "Ethan"
person2 = "Eric"

GROddstestres<-results(GRObodydds,  contrast=c("person", person1, person2))
resdata <- as.data.frame(GROddstestres)
signifGenes <- resdata[resdata$padj < 0.01,]
signifGenes <- na.omit(signifGenes)

ggplot(resdata, aes(log(baseMean,2), log2FoldChange)) + 
  geom_point(color="gray20", alpha = 0.5) + ##change opacity of points
  theme_classic() + ##use the classic theme template 
  xlab("log2 Average Normalized Counts") + ##label the x-axis
  ylab("log2 Fold Change") + ##label the y-axis
  ggtitle("MA plot") + ##to add title
  theme(plot.title = element_text(hjust = 0.5),
        title = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 18,face = "bold"),
        axis.text = element_text(size = 16)) + ##to center title
  geom_hline(aes(yintercept=0), colour="blue", linetype="dashed") + ##plot a horizontal line
  geom_point(data=signifGenes,aes(log(baseMean,2), log2FoldChange), color="red",size=1, alpha=0.75)


fullresdata <- merge(resdata,annotationnew,by.x=0,by.y=4)
fullresdata <- merge(fullresdata,common_ids,by.x=1,by.y=1)
fullresdata <- na.omit(fullresdata)
GRO_uncorrected_fullresdata <- fullresdata

# All Chr21 Genes
nrow(fullresdata[!(duplicated(fullresdata$V2)) & fullresdata$chr=="chr21",])

# Not Transcribed
nrow(fullresdata[(is.na(fullresdata$padj)) & !(duplicated(fullresdata$V2)) & fullresdata$chr=="chr21",])

# all Significant
nrow(fullresdata[fullresdata$padj<.01 & !(is.na(fullresdata$padj)) & !(duplicated(fullresdata$V2)),])

# Significant 21
nrow(fullresdata[fullresdata$padj<.01 & fullresdata$chr=="chr21" & !(is.na(fullresdata$padj)) & !(duplicated(fullresdata$V2)),])

# Below 1.5 on 21
nrow(fullresdata[fullresdata$chr=="chr21" & !(duplicated(fullresdata$V2)) & !(is.na(fullresdata$padj)) & fullresdata$log2FoldChange<0.5849,])


common_fullresdata_test <- merge(fullresdata,common_ids,by.x=1,by.y=1)

nrow(common_fullresdata_test[!(duplicated(common_fullresdata_test$V2)) & common_fullresdata_test$padj<.01,])
nrow(common_fullresdata_test[!(duplicated(common_fullresdata_test$V2)) & common_fullresdata_test$chr=="chr21",])
nrow(common_fullresdata_test[!(duplicated(common_fullresdata_test$V2)) & common_fullresdata_test$padj<.01 & common_fullresdata_test$chr=="chr21",])
nrow(common_fullresdata_test[!(duplicated(common_fullresdata_test$V2)) & common_fullresdata_test$log2FoldChange<.5849 & common_fullresdata_test$chr=="chr21",])
nrow(common_fullresdata_test[!(duplicated(common_fullresdata_test$V2)) & common_fullresdata_test$log2FoldChange>=.5849 & common_fullresdata_test$chr=="chr21",])

fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]

fullresdata <- merge(resdata,annotationnew,by.x=0,by.y=4)
fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]
#quantile(fullresdata$baseMean, probs = seq(0,1,0.2))
#fullresdata <- fullresdata[fullresdata$baseMean > 60,] # Remove genes with very low expression

medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]

#dosagelowgenes <- medresdata[medresdata$log2FoldChange < 0.13 & medresdata$chr=="chr21",]
common_fullresdata_test <- merge(fullresdata,common_ids,by.x=1,by.y=1)

nrow(common_fullresdata_test[!(duplicated(common_fullresdata_test$V2)) & common_fullresdata_test$padj<.01,])
nrow(common_fullresdata_test[!(duplicated(common_fullresdata_test$V2)) & common_fullresdata_test$chr=="chr21",])
nrow(common_fullresdata_test[!(duplicated(common_fullresdata_test$V2)) & common_fullresdata_test$padj<.01 & common_fullresdata_test$chr=="chr21",])
nrow(common_fullresdata_test[!(duplicated(common_fullresdata_test$V2)) & common_fullresdata_test$log2FoldChange<.48 & common_fullresdata_test$chr=="chr21",])


#nrow(common_fullresdata_rna[common_fullresdata_rna$padj<.01,])
medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
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

ggsave("groseq_unchanged_ethan_eliz_108onchr21_3935total.svg",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "svg")
ggsave("groseq_unchanged_ethan_eliz_108onchr21_3935total.png",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "png")


#fullresdata <- fullresdata[fullresdata$baseMean > 100,]

sigresdata <- fullresdata[fullresdata$padj<0.01,]
sigres21 <- sigresdata[sigresdata$chr == "chr21",]
fullresdata$coloring <- "red"
fullresdata$coloring[fullresdata$padj<0.01] <- "red"
fullresdata$coloring[fullresdata$padj>=0.01] <- "black"

ggplot(fullresdata, aes(x=chr, y=log2FoldChange)) + 
  geom_violin(trim=TRUE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_hline(yintercept=log2(1.5),linetype="dashed",color="blue") +
  geom_hline(yintercept=log2(0.66667),linetype="dashed",color="blue") +
  geom_hline(yintercept=0) +
  ylab("log2FC(Eddie/Eli)") +
  xlab("Chromosome") +
  ggtitle("GROseq DESeq2 Results of simulated dataset (Eddie vs Eli)") +
  theme_classic(base_size = 18) +
  geom_jitter(aes(color=coloring,fill=coloring),shape=23, position=position_jitter(0.1,seed = 1),size=0.9, alpha=0.75) +
  scale_colour_manual(values=c("black","black"), guide=FALSE) +
  scale_fill_manual(values=c("black","black"), guide = FALSE) +
  stat_summary(fun.y=median, geom="crossbar", size=0.25, color="orange") +
  theme(legend.title = element_blank())

#Volcano Plot
fullresdata <- merge(resdata,annotationnew,by.x=0,by.y=4)
signif_fullresdata <- fullresdata[fullresdata$padj<.01,]
notsignif_fullresdata <- fullresdata[fullresdata$padj>=.01,]
fullresdata_chr21 <- fullresdata[fullresdata$chr=="chr21",]

fullresdata <- na.omit(fullresdata)
head(fullresdata)
ggplot() + 
  geom_point(data=fullresdata,aes(x=log2FoldChange, y=-log10(pvalue)),size=1,alpha=0.5) +
  geom_point(data=signif_fullresdata,aes(x=log2FoldChange, y=-log10(pvalue)),size=1,alpha=0.5,color="red") +
  geom_point(data=fullresdata_chr21,aes(x=log2FoldChange, y=-log10(pvalue)),size=1,alpha=0.5,color="green") +
  ylab(paste0("-log10(pvalue)")) +
  geom_vline(xintercept = log2(1.5), linetype="dotted", color="blue" ) +
  geom_vline(xintercept = log2(1), linetype="dotted", color="black" ) +
  xlab("log2FC") +
  theme_classic(base_size = 24) +
  theme(legend.title = element_blank())

ggsave("groseq_ethaneric_volcano.svg",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "svg")
ggsave("groseq_ethaneric_volcano.png",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "png")

#  Quintile of FC violin plots
fullresdata <- merge(resdata,annotationmerge,by.x=0,by.y=3)
fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]
fullresdata[is.na(fullresdata)] = 0

as.numeric(quantile(fullresdata$baseMean,probs=seq(0,1,.2))[2])

quantile1 <- fullresdata[fullresdata$baseMean<=as.numeric(quantile(fullresdata$baseMean,probs=seq(0,1,.2))[2]),]
quantile1$quantile <- 1
medians1 <- ddply(quantile1, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians1$quantile <- 1
medians1
quantile1
quantile2 <- fullresdata[fullresdata$baseMean <= as.numeric(quantile(fullresdata$baseMean,probs=seq(0,1,.2))[3]) & fullresdata$baseMean> as.numeric(quantile(fullresdata$baseMean,probs=seq(0,1,.2))[2]),]
quantile2$quantile <- 2
medians2 <- ddply(quantile2, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians2$quantile <- 2

quantile3 <- fullresdata[fullresdata$baseMean<= as.numeric(quantile(fullresdata$baseMean,probs=seq(0,1,.2))[4]) & fullresdata$baseMean> as.numeric(quantile(fullresdata$baseMean,probs=seq(0,1,.2))[3]),]
quantile3$quantile <- 3
medians3 <- ddply(quantile3, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians3$quantile <- 3

quantile4 <- fullresdata[fullresdata$baseMean<= as.numeric(quantile(fullresdata$baseMean,probs=seq(0,1,.2))[5]) & fullresdata$baseMean> as.numeric(quantile(fullresdata$baseMean,probs=seq(0,1,.2))[4]),]
quantile4$quantile <- 4
medians4 <- ddply(quantile4, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians4$quantile <- 4

quantile5 <- fullresdata[fullresdata$baseMean<= as.numeric(quantile(fullresdata$baseMean,probs=seq(0,1,.2))[6]) & fullresdata$baseMean> as.numeric(quantile(fullresdata$baseMean,probs=seq(0,1,.2))[5]),]
quantile5$quantile <- 5
medians5 <- ddply(quantile5, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians5$quantile <- 5

quantileall <- fullresdata
quantileall$quantile <- 6
medians6 <- ddply(quantileall, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians6$quantile <- 6

newfullres <- rbind(quantile1,quantile2,quantile3,quantile4,quantile5,quantileall)
newmedians <- rbind(medians1,medians2,medians3,medians4,medians5,medians6)
newmedians

# Facet Grid of Violins
colnames(newfullres)
newfullres <- newfullres[newfullres$chr %in% lesschrs,]

ggplot() + 
  geom_violin(data=newfullres,aes(x=chr, y=log2FoldChange)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_hline(yintercept=0) +
  ylab(paste0("Log2FC")) +
  geom_hline(yintercept = log2(1.5), linetype="dotted", color="blue" ) +
  xlab("Chromosome") +
  theme_classic(base_size = 24) +
  #geom_jitter(aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=2.5, alpha=0.85,color="black",fill="black") +
  #stat_summary(data=newfullres,aes(x=chr, y=log2FoldChange),fun.y=median, geom="crossbar", size=0.25, color="orange") +
  theme(legend.title = element_blank()) + 
  facet_wrap(~quantile, ncol = 5) +
  geom_hline(data = newmedians, aes(yintercept = log2(med)),
             colour = "orange")

ggsave("groseq_ethaneric_facets.svg",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "svg")
ggsave("groseq_ethaneric_facets.png",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "png")


nrow(fullresdata[fullresdata$log2FoldChange <= 0.25 & fullresdata$chr == "chr21",])

### GRO-seq nopromoters, no multimaps, ploidy correction ####
GRObodycountdat<- read.csv("/scratch/Users/sahu0957/ds_normalization/GRO/no_repeats/counts/featurecounts_gro_noproms_nomultis.sense.txt",
                           skip=1, sep="\t")

# colnames are different due to my featurecounts script. fixing that here by removing the path. Kinda hacky, go remove in script!
names <- colnames(GRObodycountdat)
names <- gsub(".*E","E",names)
names
nrow(GRObodycountdat[GRObodycountdat$Chr=="chr21",])
colnames(GRObodycountdat) <- names
colnames(GRObodycountdat)
nomultis_col <- colnames(GRObodycountdat)[7:22]

GROmetadata=read.table("/Shares/down/mixed/RNAGRO01132021/scripts/Deseq2_analysis/GROinfo.txt", sep="\t", header=TRUE)
GROmetadata$samplegroup <-paste(GROmetadata$person, "_", GROmetadata$libprep, sep="")
GROmetadata$nomultis <- nomultis_col

head(GRObodycountdat)
minichrs <-c("chr1", "chr2", "chr3","chr4", "chr5", "chr6", "chr7", "chr8", "chr9","chr10","chr11", "chr12", "chr13","chr14", "chr15", "chr16", "chr17", "chr18", "chr19","chr20","chr21","chr22")

masterannotationdf <- annotationmerge
masterannotationdf
masterannotationdf_only21 <- masterannotationdf %>% filter(chr=="chr21")
nrow(masterannotationdf_only21)
anuplodygenes <-as.vector(masterannotationdf_only21[["name"]]) #Only chr21 genes
anuplodygenes
baseploidy <- 2
alt_ploidy <-3
masterannotationdf
head(GRObodycountdat)
colnames(GRObodycountdat)
GRObodycountdat <- GRObodycountdat[GRObodycountdat$Geneid %in% masterannotationdf$GeneID,]
GRObodycountdat <- GRObodycountdat[!(duplicated(GRObodycountdat$Geneid)),]
rownames(GRObodycountdat) <- GRObodycountdat$Geneid


GRObodycountdat <- GRObodycountdat[,-c(1:6)]

nrow(GROmetadata)
ncol(GRObodycountdat)
GRObodycountdat
#run Deseq on RNA-seq with 21
ddsFull <- DESeqDataSetFromMatrix(countData = GRObodycountdat, colData = GROmetadata, design = ~ libprep + person) #use this for all samples
#ddsFull$type <- relevel(ddsFull$type, "WT")
dds <- collapseReplicates( ddsFull,groupby = ddsFull$samplegroup,run = ddsFull$samplegroup )


dds <- ddsFull

ddsFull <- DESeq(ddsFull)

person1 = "Ethan"
person2 = "Eric"
GROddsres<-results(ddsFull,  contrast=c("person", person1, person2))
resdata <- as.data.frame(GROddsres)
resdata
fullresdata <- merge(resdata,annotationnew,by.x=0,by.y=4)
fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]
medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]

medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
controlgenes <-ifelse(rownames(dds) %in% anuplodygenes, FALSE, TRUE) # Only normalize on non-aneuploidy genes
samplepinfo<-as.data.frame(colData(dds))
dds$ploidy
ploidy_trisomy21 = ifelse(rownames(dds) %in% anuplodygenes, alt_ploidy, baseploidy)
ploidy_trisomy21
ploidy_typical <- rep(baseploidy, nrow(dds))
ploidy_trisomy21
normFactors <- matrix(do.call(cbind, mget(paste0(dds$ploidy))),ncol=ncol(dds),nrow=nrow(dds),dimnames=list(1:nrow(dds),1:ncol(dds)))
normFactors <- normFactors/baseploidy


# we've already collapsed replicates for this object, so we can just make a soft copy
nrow(normFactors)
nrow(dds)
ddsCollapsed<-DESeq(dds)
ddsCollapsed_normfactor <-estimateSizeFactors(dds,normMatrix=normFactors, controlGenes=controlgenes) #adds a column to colData(ddsCollapsed) called sizeFactor if you use normmatrix then normalizationFactors(ddsCollapsed) instead.
ddsCollapsed_normfactor <- estimateDispersionsGeneEst(ddsCollapsed_normfactor) #adds mu to assays(ddsCollapsed) mu-has one column per sample and fills in 
#elementMetadata(ddsCollapsed) with baseMean, baseVar, allZero, dispGeneEst, dispGeneIter
ddsCollapsed_normfactor<-estimateDispersionsFit(ddsCollapsed_normfactor) # and fills in dispersionFunction(ddsCollapsed) and adds a column in  elementMetadata(ddsCollapsed) called dispFit
ddsCollapsed_normfactor <- estimateDispersionsMAP(ddsCollapsed_normfactor) #adds a attr(,"dispPriorVar") to dispersionFunction(ddsCollapsed) 
#and adds 4 columns to elementMetadata(ddsCollapsed): dispersion, dispIter, dispOutlier, dispMAP
ddsCollapsed_normfactor <- nbinomWaldTest(ddsCollapsed_normfactor) #adds both H and cooks to assays(ddsCollapsed)  both are for each sample you have.
#Adds a ton of columns to elementMetadata(ddsCollapsed) 
#including Intercept and all values they expect you to compare 
#SE_Intercept (SE stands for standard error)  and SE for every thing they expect you to compare. 
#Wald static for what they expect you to compare. 
#betaConv  betaIter         deviance  maxCooks

person1 = "Ethan"
person2 = "Eric"
RNAddsres<-results(ddsCollapsed_normfactor,  contrast=c("person", person1, person2))
resdata <- as.data.frame(RNAddsres)
fullresdata <- merge(resdata,annotationnew,by.x=0,by.y=4)
fullresdata <- merge(fullresdata,common_ids,by.x=1,by.y=1)
#fullresdata <- fullresdata[!(duplicated(fullresdata$V2)),]
#fullresdata <- na.omit(fullresdata)
nrow(fullresdata)
nrow(fullresdata[fullresdata$padj<.01,])
nrow(fullresdata[fullresdata$padj>=.01 & fullresdata$chr=="chr21",])
nrow(fullresdata[fullresdata$padj<.01 & fullresdata$chr=="chr21",])

fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]
quantile(fullresdata$baseMean, probs = seq(0,1,0.2))

nrow(fullresdata[!(duplicated(fullresdata$V2)) & fullresdata$chr=="chr21" & fullresdata$log2FoldChange<0 & !(is.na(fullresdata$log2FoldChange)),])
nrow(fullresdata[!(duplicated(fullresdata$V2)) & fullresdata$chr=="chr21" & fullresdata$log2FoldChange<0 & !(is.na(fullresdata$log2FoldChange)) & fullresdata$padj<.01 & fullresdata$baseMean>15,])
nrow(fullresdata[!(duplicated(fullresdata$V2)) & fullresdata$chr=="chr21" & fullresdata$log2FoldChange<0 & !(is.na(fullresdata$log2FoldChange)) & fullresdata$baseMean<15,])


nrow(fullresdata[fullresdata$baseMean < 15 & fullresdata$chr=="chr21" & fullresdata$log2FoldChange< 0,])
nrow(fullresdata[fullresdata$baseMean >= 15 & fullresdata$chr=="chr21" & fullresdata$padj < .01 & fullresdata$log2FoldChange<0,])
elieliz_ploidyresilientdiffs <- (fullresdata[fullresdata$chr=="chr21" & fullresdata$padj < .01,])
ethaneric_ploidyresilientdiffs$V2 %in% elieliz_ploidyresilientdiffs$V2

fullresdata <- fullresdata[fullresdata$baseMean > 15,] # Remove genes with very low expression

nrow(fullresdata)
nrow(fullresdata[fullresdata$padj<.01,])
nrow(fullresdata[fullresdata$padj<.01 & fullresdata$chr=="chr21",])
nrow(fullresdata[fullresdata$chr=="chr21",])

medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]
medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians
nrow(medresdata[medresdata$chr=="chr21" & medresdata$log2FoldChange<0,])
signif_medresdata <- medresdata[medresdata$padj<.01,]
notsignif_medresdata <- medresdata[medresdata$padj>=.01,]
nrow(medresdata[medresdata$chr=="chr21",])
nrow(signif_medresdata[signif_medresdata$chr=="chr21",])
nrow(notsignif_medresdata)
ggplot() + 
  geom_violin(data=fullresdata,trim=TRUE,aes(x=chr, y=log2FoldChange)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  #  geom_hline(yintercept=log2(0.66667),linetype="dashed",color="red") +
  geom_hline(yintercept=0) +
  ylab(paste0("Log2FC")) +
  xlab("Chromosome") +
  theme_classic(base_size = 24) +
  geom_jitter(data = signif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=1.5, alpha=0.85,color="red",fill="red") +
  geom_jitter(data = notsignif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=1.5, alpha=0.85,color="black",fill="black") +
  stat_summary(data=fullresdata,aes(x=chr, y=log2FoldChange),fun.y=median, geom="crossbar", size=0.25, color="orange") +
  theme(legend.title = element_blank())

ggsave("groseq_changed_ethan_eric_51onchr21_4311total.svg",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "svg")
ggsave("groseq_changed_ethan_eric_51onchr21_4311total.png",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "png")


ggplot(fullresdata, aes(x=chr, y=log2FoldChange)) + 
  geom_violin(trim=TRUE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_hline(yintercept=log2(1.5),linetype="dashed",color="blue") +
  geom_hline(yintercept=log2(0.66667),linetype="dashed",color="blue") +
  geom_hline(yintercept=0) +
  ylab("log2FC") +
  xlab("Chromosome") +
  theme_classic(base_size = 18) +
  geom_jitter(shape=23, position=position_jitter(0.1,seed = 1),size=0.9, alpha=0.75) +
  scale_colour_manual(values=c("black","black"), guide=FALSE) +
  scale_fill_manual(values=c("black","black"), guide = FALSE) +
  stat_summary(fun.y=median, geom="crossbar", size=0.25, color="orange") +
  theme(legend.title = element_blank())

####GRO-seq, alternative hypothesis adjustment ####
GRObodycountdat<- read.csv("/scratch/Users/sahu0957/ds_normalization/GRO/no_repeats/counts/featurecounts_gro_noproms_nomultis.sense.txt",
                           skip=1, sep="\t")

# colnames are different due to my featurecounts script. fixing that here by removing the path. Kinda hacky, go remove in script!
names <- colnames(GRObodycountdat)
names <- gsub(".*E","E",names)
names
nrow(GRObodycountdat[GRObodycountdat$Chr=="chr21",])
colnames(GRObodycountdat) <- names
colnames(GRObodycountdat)
nomultis_col <- colnames(GRObodycountdat)[7:22]

GROmetadata=read.table("/Shares/down/mixed/RNAGRO01132021/scripts/Deseq2_analysis/GROinfo.txt", sep="\t", header=TRUE)
GROmetadata$samplegroup <-paste(GROmetadata$person, "_", GROmetadata$libprep, sep="")
GROmetadata$nomultis <- nomultis_col

minichrs <-c("chr1", "chr2", "chr3","chr4", "chr5", "chr6", "chr7", "chr8", "chr9","chr10","chr11", "chr12", "chr13","chr14", "chr15", "chr16", "chr17", "chr18", "chr19","chr20","chr21","chr22")

GRObodycountdat <- GRObodycountdat[GRObodycountdat$Geneid %in% masterannotationdf$GeneID,]
GRObodycountdat <- GRObodycountdat[!(duplicated(GRObodycountdat$Geneid)),]
rownames(GRObodycountdat) <- GRObodycountdat$Geneid


GRObodycountdat <- GRObodycountdat[,-c(1:6)]
GRObodycountdat <- GRObodycountdat[,!(colnames(GRObodycountdat) %like% "Eric")]
GROmetadata <- GROmetadata[!(GROmetadata$person %like% "Eric"),]

nrow(GROmetadata)
ncol(GRObodycountdat)
GRObodycountdat
#run Deseq on GRO-seq with 21
ddsFull <- DESeqDataSetFromMatrix(countData = GRObodycountdat, colData = GROmetadata, design = ~ libprep + person) #use this for all samples
#ddsFull$type <- relevel(ddsFull$type, "WT")
ddsFull <- collapseReplicates( ddsFull,groupby = ddsFull$samplegroup,run = ddsFull$samplegroup )
ddsFull <- DESeq(ddsFull)

metrics_info_unchanged <- as.data.frame(mcols(ddsFull))
metrics_info_unchanged
metrics_info_unchanged <- metrics_info_unchanged[!(rownames(metrics_info_unchanged) %in% rownames(annotationmergeno21)),]
nrow(metrics_info_unchanged)

# Fetch all gene dispersions leaving out Ethan
GRObodycountdat<- read.csv("/scratch/Users/sahu0957/ds_normalization/GRO/no_repeats/counts/featurecounts_gro_noproms_nomultis.sense.txt",
                           skip=1, sep="\t")

# colnames are different due to my featurecounts script. fixing that here by removing the path. Kinda hacky, go remove in script!
names <- colnames(GRObodycountdat)
names <- gsub(".*E","E",names)
names
nrow(GRObodycountdat[GRObodycountdat$Chr=="chr21",])
colnames(GRObodycountdat) <- names
colnames(GRObodycountdat)
nomultis_col <- colnames(GRObodycountdat)[7:22]

GROmetadata=read.table("/Shares/down/mixed/RNAGRO01132021/scripts/Deseq2_analysis/GROinfo.txt", sep="\t", header=TRUE)
GROmetadata$samplegroup <-paste(GROmetadata$person, "_", GROmetadata$libprep, sep="")
GROmetadata$nomultis <- nomultis_col

minichrs <-c("chr1", "chr2", "chr3","chr4", "chr5", "chr6", "chr7", "chr8", "chr9","chr10","chr11", "chr12", "chr13","chr14", "chr15", "chr16", "chr17", "chr18", "chr19","chr20","chr21","chr22")

GRObodycountdat <- GRObodycountdat[GRObodycountdat$Geneid %in% masterannotationdf$GeneID,]
GRObodycountdat <- GRObodycountdat[!(duplicated(GRObodycountdat$Geneid)),]
rownames(GRObodycountdat) <- GRObodycountdat$Geneid


GRObodycountdat <- GRObodycountdat[,-c(1:6)]

nrow(GROmetadata)

#run Deseq on GRO-seq with Ethan removed
GRObodycountdat <- GRObodycountdat[,!(colnames(GRObodycountdat) %like% "Ethan")]
GROmetadata <- GROmetadata[!(GROmetadata$person %like% "Ethan"),]
GROmetadata
ddsFull <- DESeqDataSetFromMatrix(countData = GRObodycountdat, colData = GROmetadata, design = ~ libprep + person) #use this for all samples
#ddsFull$type <- relevel(ddsFull$type, "WT")
ddsFull <- collapseReplicates( ddsFull,groupby = ddsFull$samplegroup,run = ddsFull$samplegroup )

ddsFull <- DESeq(ddsFull)

metrics_info_noethan <- mcols(ddsFull)
metrics_info_noethan <- metrics_info_noethan[!(rownames(metrics_info_noethan) %in% rownames(annotationmergeno21)),]

nrow(metrics_info_noethan)
nrow(metrics_info_unchanged)
dispersion_nochr21_noethan <- merge(as.data.frame(metrics_info_unchanged),as.data.frame(metrics_info_noethan),by.x=0,by.y=0)
colnames(dispersion_nochr21_noethan)

eq <- function(x,y) {
  m <- lm(y ~ x)
  as.character(
    as.expression(
      substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                 list(a = format(coef(m)[1], digits = 4),
                      b = format(coef(m)[2], digits = 4),
                      r2 = format(summary(m)$r.squared, digits = 3)))
    )
  )
}

eq(dispersion_nochr21_noethan$dispMAP.x,dispersion_nochr21_noethan$dispMAP.y)

ggplot(df,aes(x = wt, y = hp)) + 
  geom_point() + 
  geom_smooth(method = "lm", se=FALSE) 
View(dispersion_nochr21_noethan)
colnames(dispersion_nochr21_noethan)
ggplot() + 
  geom_point(data=dispersion_nochr21_noethan,aes(x=as.numeric(as.character(dispersion.x)), y=as.numeric(as.character(dispersion.y)))) +
  geom_smooth(method="lm",color="blue") +
  geom_abline(color="red") +
  ylab(paste0("T21 Removed")) +
  xlab("D21 Brother Removed") +
  scale_x_log10() +
  scale_y_log10() +
  geom_text(x = 0, y = 0, label = eq(dispersion_nochr21_noethan$dispMAP.y,dispersion_nochr21_noethan$dispMAP.y), parse = TRUE) +
  ggtitle("Gene-wise Dispersion Fitted Estimates (GRO-seq)") +
  theme_classic(base_size = 24)
#geom_jitter(aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=2.5, alpha=0.85,color="black",fill="black") +
#stat_summary(data=newfullres,aes(x=chr, y=log2FoldChange),fun.y=median, geom="crossbar", size=0.25, color="orange") +
# theme(legend.title = element_blank()) 

ggsave("groseq_genewise_dispersion_brothersremoved.svg",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "svg")
ggsave("groseq_genewise_dispersion_brothersremoved.png",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "png")











GROddsres<-results(ddsFull,  contrast=c("person", person1, person2),lfcThreshold = log2(1.5),altHypothesis = "lessAbs")
resdata <- as.data.frame(GROddsres)
resdata
fullresdata <- merge(resdata,annotationnew,by.x=0,by.y=4)
fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]
medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]
medresdata
medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians

signif_medresdata <- medresdata[medresdata$padj<.1,]
notsignif_medresdata <- medresdata[medresdata$padj>=.1,]
View(medresdata)
nrow(signif_medresdata)
nrow(notsignif_medresdata)
ggplot() + 
  geom_violin(data=fullresdata,trim=TRUE,aes(x=chr, y=log2FoldChange)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  #  geom_hline(yintercept=log2(0.66667),linetype="dashed",color="red") +
  geom_hline(yintercept=0) +
  ylab(paste0("Log2FC")) +
  geom_hline(yintercept = log2(1.5), linetype="dotted", color="blue" ) +
  xlab("Chromosome") +
  theme_classic(base_size = 24) +
  geom_jitter(data = notsignif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=1.5, alpha=0.85,color="black",fill="black") +
  geom_jitter(data = signif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=1.5, alpha=0.85,color="red",fill="red") +
  stat_summary(data=fullresdata,aes(x=chr, y=log2FoldChange),fun.y=median, geom="crossbar", size=0.25, color="orange") +
  theme(legend.title = element_blank())

ggsave("groseq_adjustalthypothesis_56signif.svg",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "svg")
ggsave("groseq_adjustalthypothesis_56signif.png",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "png")
###########

#### GRO-sims alternative hypothesis####
  RNAcountdat <- read.csv("/Users/sahu0957/trisomy_normalization/RNAGRO01132021/data/RNAseqfiles/GROSIMS_varydepth0.08_p8_r3_s1_n3.csv", 
                          sep=",", row.names=1)
head(RNAcountdat)  
#Remove genes that are not in annotation (which was filtered to remove genes not in both and not on main chromosome above) 
  #put the data frames in the same order
  annotationmerge
  library(dplyr)
  RNAcountdat<- RNAcountdat %>%
    filter(Gene_id %in% annotationmerge$name) %>%
    arrange(Gene_id)
  RNAcountdat <- RNAcountdat[,-c(1)]
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
  RNAdds <-DESeq(RNAdds)

  person1 = "Eddie"
  person2 = "Elvis"
  RNAddsres<-results(RNAdds,  contrast=c("Person", person1, person2))
  resdata <- as.data.frame(RNAddsres)
  fullresdata <- merge(resdata,annotationnew,by.x=0,by.y=4)
  fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]
  medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]

  medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]
  medresdata
  medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
  medians
  
  signif_medresdata <- medresdata[medresdata$padj<.1,]
  notsignif_medresdata <- medresdata[medresdata$padj>=.1,]
  View(medresdata)
  nrow(signif_medresdata)
  nrow(notsignif_medresdata)
  ggplot() + 
    geom_violin(data=fullresdata,trim=TRUE,aes(x=chr, y=log2FoldChange)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    #  geom_hline(yintercept=log2(0.66667),linetype="dashed",color="red") +
    geom_hline(yintercept=0) +
    ylim(-5,5) +
    ylab(paste0("Log2FC")) +
    geom_hline(yintercept = log2(1.5), linetype="dotted", color="blue" ) +
    xlab("Chromosome") +
    theme_classic(base_size = 24) +
    geom_jitter(data = notsignif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=1.5, alpha=0.85,color="black",fill="black") +
    geom_jitter(data = signif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=1.5, alpha=0.85,color="red",fill="red") +
    stat_summary(data=fullresdata,aes(x=chr, y=log2FoldChange),fun.y=median, geom="crossbar", size=0.25, color="orange") +
    theme(legend.title = element_blank())
  
  ggsave("grosims_08disp_8poiss_37_1scale_3reps_0signif.svg",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "svg")
  ggsave("grosims_08disp_8poiss_37_1scale_3reps_0signif.png",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "png")
###########

# colnames are different due to my featurecounts script. fixing that here by removing the path. Kinda hacky, go remove in script!
names <- colnames(GRObodycountdat)
names <- gsub(".*E","E",names)
names
nrow(GRObodycountdat[GRObodycountdat$Chr=="chr21",])
colnames(GRObodycountdat) <- names
colnames(GRObodycountdat)
nomultis_col <- colnames(GRObodycountdat)[7:22]

GROmetadata=read.table("/Shares/down/mixed/RNAGRO01132021/scripts/Deseq2_analysis/GROinfo.txt", sep="\t", header=TRUE)
GROmetadata$samplegroup <-paste(GROmetadata$person, "_", GROmetadata$libprep, sep="")
GROmetadata$nomultis <- nomultis_col

minichrs <-c("chr1", "chr2", "chr3","chr4", "chr5", "chr6", "chr7", "chr8", "chr9","chr10","chr11", "chr12", "chr13","chr14", "chr15", "chr16", "chr17", "chr18", "chr19","chr20","chr21","chr22")

GRObodycountdat <- GRObodycountdat[GRObodycountdat$Geneid %in% masterannotationdf$GeneID,]
GRObodycountdat <- GRObodycountdat[!(duplicated(GRObodycountdat$Geneid)),]
rownames(GRObodycountdat) <- GRObodycountdat$Geneid


GRObodycountdat <- GRObodycountdat[,-c(1:6)]

nrow(GROmetadata)
ncol(GRObodycountdat)
GRObodycountdat
#run Deseq on RNA-seq with 21
ddsFull <- DESeqDataSetFromMatrix(countData = GRObodycountdat, colData = GROmetadata, design = ~ libprep + person) #use this for all samples
#ddsFull$type <- relevel(ddsFull$type, "WT")
dds <- collapseReplicates( ddsFull,groupby = ddsFull$samplegroup,run = ddsFull$samplegroup )


dds <- ddsFull

ddsFull <- DESeq(ddsFull)

person1 = "Ethan"
person2 = "Eric"
GROddsres<-results(ddsFull,  contrast=c("person", person1, person2),lfcThreshold = log2(1.5),altHypothesis = "lessAbs")
resdata <- as.data.frame(GROddsres)
resdata
fullresdata <- merge(resdata,annotationnew,by.x=0,by.y=4)
fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]
medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]
medresdata
medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians

signif_medresdata <- medresdata[medresdata$padj<.1,]
notsignif_medresdata <- medresdata[medresdata$padj>=.1,]
View(medresdata)
nrow(signif_medresdata)
nrow(notsignif_medresdata)
ggplot() + 
  geom_violin(data=fullresdata,trim=TRUE,aes(x=chr, y=log2FoldChange)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  #  geom_hline(yintercept=log2(0.66667),linetype="dashed",color="red") +
  geom_hline(yintercept=0) +
  ylab(paste0("Log2FC")) +
  geom_hline(yintercept = log2(1.5), linetype="dotted", color="blue" ) +
  xlab("Chromosome") +
  theme_classic(base_size = 24) +
  geom_jitter(data = notsignif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=1.5, alpha=0.85,color="black",fill="black") +
  geom_jitter(data = signif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=1.5, alpha=0.85,color="red",fill="red") +
  stat_summary(data=fullresdata,aes(x=chr, y=log2FoldChange),fun.y=median, geom="crossbar", size=0.25, color="orange") +
  theme(legend.title = element_blank())

ggsave("groseq_adjustalthypothesis_56signif.svg",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "svg")
ggsave("groseq_adjustalthypothesis_56signif.png",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "png")


###########



# Some leftover comparisons
common_fullresdata_rna <- merge(fullresdata,common_ids,by.x=1,by.y=1)
common_fullresdata_rna_gro <- merge(common_fullresdata_rna,common_fullresdata_gro,by.x=1,by.y=1)
colnames(common_fullresdata_rna_gro)

eq <- function(x,y) {
  m <- lm(y ~ x)
  as.character(
    as.expression(
      substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                 list(a = format(coef(m)[1], digits = 4),
                      b = format(coef(m)[2], digits = 4),
                      r2 = format(summary(m)$r.squared, digits = 3)))
    )
  )
}

eq(common_fullresdata_rna_gro$log2FoldChange.x,common_fullresdata_rna_gro$log2FoldChange.y)

common_fullresdata_rna_gro_sig <- common_fullresdata_rna_gro[common_fullresdata_rna_gro$padj.x < .01 | common_fullresdata_rna_gro$padj.y < .01,]

common_fullresdata_rna_gro_chr21 <- common_fullresdata_rna_gro[common_fullresdata_rna_gro$chr.x == "chr21",]
rna_21_hist <- common_fullresdata_rna[!duplicated(common_fullresdata_rna$V2)& common_fullresdata_rna$chr=="chr21",]

nrow(common_fullresdata_rna[!duplicated(common_fullresdata_rna$V2)&common_fullresdata_rna$log2FoldChange<.48 & common_fullresdata_rna$chr=="chr21",])
nrow(common_fullresdata_rna[!duplicated(common_fullresdata_rna$V2)&common_fullresdata_rna$chr=="chr21",])

ggplot(rna_21_hist, aes(x=log2FoldChange)) + 
  geom_histogram(aes(y=..density..),color="black", fill="skyblue") +
  geom_vline(xintercept = 0.58, color = "red") +
  geom_vline(xintercept = 0,color="black") +
  ggtitle("RNA-seq Chr21 Fold Change") + 
  ylab("Density") +
  xlab("log2(T21/D21)") +
  theme_classic(base_size=18)

ggplot(data=common_fullresdata_rna_gro_chr21,aes(x=log2FoldChange.x,y=log2FoldChange.y)) +
  geom_point(alpha=.25,size=.5) +
  geom_abline(color="red") +
  geom_hline(yintercept = 0) +
  theme_classic(base_size=18) +
  xlab("Log2FC(RNA-Seq)") +
  ylab("Log2FC(GRO-Seq") +
  xlim(-5,5) +
  ylim(-5,5) +
  geom_smooth(color="blue",method = "lm") +
  #  geom_text(x = 2, y = 2, label = eq(common_fullresdata_rna_gro$log2FoldChange.x,common_fullresdata_rna_gro$log2FoldChange.y), parse = TRUE) +
  geom_vline(xintercept = 0)

ggplot(data=common_fullresdata_rna_gro_sig,aes(x=log(baseMean.x),y=log(baseMean.y))) +
  geom_point(alpha=.25,size=.5) +
  geom_abline(color="red") +
  geom_hline(yintercept = 0) +
  theme_classic(base_size=18) +
  xlab("Log2FC(RNA-Seq)") +
  ylab("Log2FC(GRO-Seq") +
  #  xlim(-5,5) +
  #  ylim(-5,5) +
  geom_smooth(color="blue",method = "lm") +
  #  geom_text(x = 2, y = 2, label = eq(common_fullresdata_rna_gro$log2FoldChange.x,common_fullresdata_rna_gro$log2FoldChange.y), parse = TRUE) +
  geom_vline(xintercept = 0)

nrow(common_fullresdata_rna[common_fullresdata_rna$padj<.01,])
nrow(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$padj<.01 & common_fullresdata_rna$chr=="chr21",])


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
fullresdata <- merge(resdata,annotationnew,by.x=0,by.y=4)
common_fullresdata_rna <- merge(fullresdata,common_ids,by.x=1,by.y=1)

fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]
medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]

log2(1.5)
log2(1.3)
log2(1.7)
nrow(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$chr=="chr21",])
nrow(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$chr=="chr21" & common_fullresdata_rna$log2FoldChange<0.58,])
nrow(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$chr=="chr21" & common_fullresdata_rna$log2FoldChange>=0.58,])


nrow(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$chr=="chr21" & common_fullresdata_rna$log2FoldChange>0.3785 & common_fullresdata_rna$log2FoldChange<0.7655,])
nrow(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$chr=="chr21" & is.na(common_fullresdata_rna$log2FoldChange),])
#View(common_fullresdata_rna)
nrow(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$chr=="chr21" & !(is.na(common_fullresdata_rna$padj)),])
nrow(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$padj<.01 & common_fullresdata_rna$chr=="chr21" & !(is.na(common_fullresdata_rna$padj)) & common_fullresdata_rna$baseMean>60,])
nrow(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$chr=="chr21" & common_fullresdata_rna$baseMean>60,])
nrow(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$chr=="chr21" & common_fullresdata_rna$log2FoldChange>0 & !(is.na(common_fullresdata_rna$log2FoldChange)),]) 


#View(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$chr=="chr21" & common_fullresdata_rna$log2FoldChange<0 & common_fullresdata_rna$baseMean>60 & common_fullresdata_rna$padj>=.01 & !(is.na(common_fullresdata_rna$padj)),])

nrow(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$log2FoldChange<.48 & common_fullresdata_rna$chr=="chr21",])


nrow(common_fullresdata_rna[common_fullresdata_rna$padj<.01,])




library(plyr)

medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians

nrow(medresdata[medresdata$chr=="chr21" & medresdata$log2FoldChange<0,])
signif_medresdata <- medresdata[medresdata$padj<.01,]
notsignif_medresdata <- medresdata[medresdata$padj>=.01,]
nrow(medresdata[medresdata$chr=="chr21",])
nrow(signif_medresdata[signif_medresdata$chr=="chr21",])
nrow(notsignif_medresdata)
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

ggsave("simseq_highreps_08disp_8poiss_uncorrected.svg",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "svg")
ggsave("simseq_highreps_08disp_8poiss_uncorrected.png",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "png")

controlgenes <-ifelse(rownames(dds) %in% anuplodygenes, FALSE, TRUE) # Only normalize on non-aneuploidy genes

samplepinfo<-as.data.frame(colData(dds))
ploidy_trisomy21 = ifelse(rownames(dds) %in% anuplodygenes, alt_ploidy, baseploidy)
ploidy_typical <- rep(baseploidy, nrow(dds))
normFactors <- matrix(do.call(cbind, mget(paste0(dds$ploidy))),ncol=ncol(dds),nrow=nrow(dds),dimnames=list(1:nrow(dds),1:ncol(dds)))
normFactors <- normFactors/baseploidy
unique(normFactors)


ddsCollapsed<-DESeq(dds)
ddsCollapsed_normfactor <-estimateSizeFactors(dds,normMatrix=normFactors, controlGenes=controlgenes) #adds a column to colData(ddsCollapsed) called sizeFactor if you use normmatrix then normalizationFactors(ddsCollapsed) instead.
ddsCollapsed_normfactor <- estimateDispersionsGeneEst(ddsCollapsed_normfactor) #adds mu to assays(ddsCollapsed) mu-has one column per sample and fills in 
#elementMetadata(ddsCollapsed) with baseMean, baseVar, allZero, dispGeneEst, dispGeneIter
ddsCollapsed_normfactor<-estimateDispersionsFit(ddsCollapsed_normfactor) # and fills in dispersionFunction(ddsCollapsed) and adds a column in  elementMetadata(ddsCollapsed) called dispFit
ddsCollapsed_normfactor <- estimateDispersionsMAP(ddsCollapsed_normfactor) #adds a attr(,"dispPriorVar") to dispersionFunction(ddsCollapsed) 
#and adds 4 columns to elementMetadata(ddsCollapsed): dispersion, dispIter, dispOutlier, dispMAP
ddsCollapsed_normfactor <- nbinomWaldTest(ddsCollapsed_normfactor) #adds both H and cooks to assays(ddsCollapsed)  both are for each sample you have.

person1 = "Eddie"
person2 = "Elvis"
RNAddsres<-results(ddsCollapsed_normfactor,  contrast=c("Person", person1, person2))
resdata <- as.data.frame(RNAddsres)
fullresdata <- merge(resdata,annotationnew,by.x=0,by.y=4)
fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]
quantile(fullresdata$baseMean, probs = seq(0,1,0.1))
fullresdata <- fullresdata[fullresdata$baseMean > 6,] # Remove genes with very low expression

medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]

medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians

nrow(fullresdata)
nrow(fullresdata[fullresdata$padj<.01,])
nrow(fullresdata[fullresdata$chr=="chr21",])

nrow(fullresdata[fullresdata$padj<.01 & fullresdata$chr=="chr21",])
nrow(fullresdata[fullresdata$padj<.01 & fullresdata$chr=="chr21",])

nrow(fullresdata[fullresdata$chr=="chr21",])



medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]
medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians
nrow(medresdata[medresdata$chr=="chr21" & medresdata$log2FoldChange<0,])
signif_medresdata <- medresdata[medresdata$padj<.01,]
notsignif_medresdata <- medresdata[medresdata$padj>=.01,]
nrow(medresdata[medresdata$chr=="chr21",])
nrow(signif_medresdata[signif_medresdata$chr=="chr21",])
nrow(notsignif_medresdata)
ggplot() + 
  geom_violin(data=fullresdata,trim=TRUE,aes(x=chr, y=log2FoldChange)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  #  geom_hline(yintercept=log2(0.66667),linetype="dashed",color="red") +
  geom_hline(yintercept=0) +
  ylab(paste0("Log2FC")) +
  xlab("Chromosome") +
  theme_classic(base_size = 24) +
 # geom_hline(yintercept=log2(1.5),color="blue",linetype="dashed") +
  geom_jitter(data = signif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=1.5, alpha=0.85,color="red",fill="red") +
  geom_jitter(data = notsignif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=1.5, alpha=0.85,color="black",fill="black") +
  stat_summary(data=fullresdata,aes(x=chr, y=log2FoldChange),fun.y=median, geom="crossbar", size=0.25, color="orange") +
  theme(legend.title = element_blank())

ggsave("simseq_highdisp_lowreps_corrected.svg",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "svg")
ggsave("simseq_highdisp_lowreps_corrected.png",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "png")

common_fullresdata_rna <- merge(fullresdata,common_ids,by.x=1,by.y=1)
log2(1.5)
nrow(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$chr=="chr21" & common_fullresdata_rna$log2FoldChange<0.58,])
nrow(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$chr=="chr21" & is.na(common_fullresdata_rna$log2FoldChange),])

nrow(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$chr=="chr21" & !(is.na(common_fullresdata_rna$padj)),])
(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$padj<.01 & common_fullresdata_rna$chr=="chr21" & !(is.na(common_fullresdata_rna$padj)) & common_fullresdata_rna$baseMean>60,])
nrow(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$chr=="chr21" & common_fullresdata_rna$baseMean>60,])
nrow(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$chr=="chr21" & common_fullresdata_rna$log2FoldChange>0 & !(is.na(common_fullresdata_rna$log2FoldChange)),]) 


View(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$chr=="chr21" & common_fullresdata_rna$log2FoldChange<0 & common_fullresdata_rna$baseMean>60 & common_fullresdata_rna$padj>=.01 & !(is.na(common_fullresdata_rna$padj)),])

nrow(common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$log2FoldChange<.48 & common_fullresdata_rna$chr=="chr21",])


nrow(common_fullresdata_rna[common_fullresdata_rna$padj<.01,])















