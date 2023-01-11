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
lesschrs <- c("chr21")

RNAcountdat <- read.csv("/Users/sahu0957/trisomy_normalization/RNAGRO01132021/data/RNAseqfiles/NEWSCALEreplicatefinalsims_a0.001_p1_r2_s1_n3.csv", 
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
RNAdds <-DESeq(RNAdds, betaPrior = TRUE) # Earlier DESeq2 versions include a betaPrior for the intercept term... wild!

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
library(plyr)
medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians

common_fullresdata_rna <- merge(medresdata,common_ids,by.x=1,by.y=1)

signif_medresdata <- medresdata[medresdata$padj<.01,]
notsignif_medresdata <- medresdata[medresdata$padj>=.01,]
signif_medresdata
nrow(notsignif_medresdata)
ggplot() + 
  geom_violin(data=fullresdata,trim=TRUE,aes(x=chr, y=log2FoldChange)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  #  geom_hline(yintercept=log2(0.66667),linetype="dashed",color="red") +
  geom_hline(yintercept=0) +
  ylab(paste0("Log2FC")) +
  ylim(-3,3) +
  geom_hline(yintercept = log2(1.5), linetype="dotted", color="blue" ) +
  xlab("Chromosome") +
  theme_classic(base_size = 24) +
  geom_jitter(data = signif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=1.5, alpha=0.85,color="red",fill="red") +
  geom_jitter(data = notsignif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=1.5, alpha=0.85,color="black",fill="black") +
  stat_summary(data=fullresdata,aes(x=chr, y=log2FoldChange),fun.y=median, geom="crossbar", size=0.25, color="orange") +
  theme(legend.title = element_blank())

ggsave("rnaseq_deseq2_realsims_new_violins_map_med124.svg",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "svg")
ggsave("rnaseq_deseq2_realsims_new_violins_map_med124.png",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "png")

### Real Data ###
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
RNAmetadata$samplegroup <-paste(RNAmetadata$Person, "_", RNAmetadata$biological_rep, sep="")
head(RNAcountdat)


#this is where you would remove samples if they are bad
#RNAcountdat <- read.csv(RNAcoveragedat,sep="\t",skip=1)
RNAcountdat <- read.csv("/scratch/Users/sahu0957/ds_normalization/RNA/nextflowhg38_20220214/counts/featurecounts_rna_withmultis.sense.txt",
                        sep="\t", skip=1)

colnames(RNAcountdat)
RNAcountdat <- RNAcountdat[!(duplicated(RNAcountdat$Geneid)),]
rownames(RNAcountdat) <- RNAcountdat$Geneid
library(dplyr)
RNAcountdat<- RNAcountdat %>%
  filter(Geneid %in% masterannotationdf$name) %>%
  arrange(Geneid)
nrow(RNAcountdat)

#run Deseq on RNA-seq with 21
RNAcountdat <- RNAcountdat[,-c(1:6)]
ddsFull <- DESeqDataSetFromMatrix(countData = RNAcountdat, colData = RNAmetadata, design = ~biological_rep + Person)

# Can run DESeq2 unaltered here to compare results
ddsFull <- DESeq(ddsFull, betaPrior = TRUE)
RNAddsres<-results(ddsFull,contrast=c("Person","Ethan","Eric"))
resdata <- as.data.frame(RNAddsres)

fullresdata <- merge(resdata,annotationnew,by.x=0,by.y=4)
#fullresdata <- merge(fullresdata,common_ids,by.x=1,by.y=1)

common_fullresdata_test <- merge(fullresdata,common_ids,by.x=1,by.y=1)


fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]

medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]

library(plyr)
medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians

signif_medresdata <- medresdata[medresdata$padj<.01,]
notsignif_medresdata <- medresdata[medresdata$padj>=.01,]

fullresdata$chr <- factor(fullresdata$chr, levels=minichrs)

ggplot() + 
  geom_violin(data=fullresdata,trim=TRUE,aes(x=chr, y=log2FoldChange)) +
  #  geom_hline(yintercept=log2(0.66667),linetype="dashed",color="red") +
  geom_hline(yintercept=0) +
  ylab(paste0("Log2FC")) +
  geom_hline(yintercept = log2(1.5), linetype="dotted", color="blue" ) +
  theme_classic(base_size = 30) +
  ylim(-3,3) +
  geom_jitter(data = signif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=2.5, alpha=0.85,color="red",fill="red") +
  geom_jitter(data = notsignif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=2.5, alpha=0.85,color="black",fill="black") +
  stat_summary(data=fullresdata,aes(x=chr, y=log2FoldChange),fun.y=median, geom="crossbar", size=0.25, color="orange") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.title = element_blank()) +
  theme(axis.title.x = element_blank())

ggsave("rnaseq_map_med124.svg",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "svg")
ggsave("rnaseq_map_med124.png",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "png")


