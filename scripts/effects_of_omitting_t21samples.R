library("reshape2")
library("plyr")
library("ggplot2")
library("miscTools")
library("stringr")
library("DESeq2")
library("tidyverse")
library("UpSetR")
theme_set(
  theme_classic(base_size = 18)
)

RNAindir<-"DS_Normalization/counts/rna"
GROseq_indir<-"DS_Normalization/counts/gro"
RNAcoveragedat<-"res_featureCounts_gene_idfull_143138.coverage.csv"

#Below is the bed that got used for RNA-seq. 
#This is not the bed file that got used for GRO-seq exactly becuase to be acuurate in the GROS-seq I had to remove any gene bodys/tss that were in the data more than once. 
ori_worldbed <- "DS_Normalization/annotation/hg38_refseq.bed"
#Below are the bed file I used for GRO-seq anaysis
#in brief, I keep the longest isoform with uniq body start and stop, body(+1000, -500) start and stop, tss(-500, +500) start and stop
#I also removed gene <3000bp and genes whose body or tss <0 in coordinates
#below I will only keep genes that were analyzed by both GR0-seq and RNA-seq
beddir<-"DS_Normalization/annotation"
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
annotationmerge <- annotationmerge %>%
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

# Assign Common IDs to genes, for maximal isoform filtering
common_ids <- read.table("DS_Normalization/annotation/refseq_to_common_id.txt",sep="\t")
common_ids <- common_ids[common_ids$V1 %in% annotationmerge$GeneID,]
####Uncorrected Real RNA-seq data####
RNAbed <-ori_worldbed
filetable <- RNAmetadata
masterannotationdf <- annotationmerge

anuplodygenes <-masterannotationdf_only21[["name"]] #Only chr21 genes
baseploidy <- 2
alt_ploidy <-3

#this is loading the metadata
### RNA-seq ####
RNAmetadata=read.table("DS_Normalization/metadata/RNAinfo.txt", sep="\t", header=TRUE)
RNAmetadata$samplegroup <-paste(RNAmetadata$Person, "_", RNAmetadata$biological_rep, sep="")
RNAcountdat <- read.csv(RNAcoveragedat,sep="\t",skip=1)
RNAcountdat$Gene_id <- rownames(RNAcountdat)
RNAcountdat<- RNAcountdat %>%
  filter(Gene_id %in% annotationmerge$name) %>%
  arrange(Gene_id)

rownames(RNAcountdat) <- RNAcountdat$Gene_id
RNAcountdat <- RNAcountdat[,-c(ncol(RNAcountdat))]

RNAcountdat <- RNAcountdat[!(duplicated(RNAcountdat$Geneid)),]
rownames(RNAcountdat) <- RNAcountdat$Geneid
RNAcountdat<- RNAcountdat %>%
  filter(Geneid %in% masterannotationdf$name) %>%
  arrange(Geneid)

rownames(RNAcountdat) <- RNAcountdat$Geneid
RNAcountdat <- RNAcountdat[-c(1:6)]
#run Deseq on RNA-seq with 21

RNAcountdat <- RNAcountdat[,!(colnames(RNAcountdat) %like% "Ethan")]
RNAmetadata <- RNAmetadata[!(RNAmetadata$Person %like% "Ethan"),]
ddsFull <- DESeqDataSetFromMatrix(countData = RNAcountdat, colData = RNAmetadata, design = ~biological_rep + Person)

# Can run DESeq2 unaltered here to compare results
ddsFull <- DESeq(ddsFull)
RNAddsres<-results(ddsFull,contrast=c("Person","Eli","Eric"))
noethan_resdata <- as.data.frame(RNAddsres)

#this is loading the metadata
RNAmetadata=read.table("DS_Normalization/metadata/RNAinfo.txt", sep="\t", header=TRUE)
RNAmetadata$samplegroup <-paste(RNAmetadata$Person, "_", RNAmetadata$biological_rep, sep="")
RNAcountdat <- read.csv(RNAcoveragedat,sep="\t",skip=1)
RNAcountdat$Gene_id <- rownames(RNAcountdat)
RNAcountdat<- RNAcountdat %>%
  filter(Gene_id %in% annotationmerge$name) %>%
  arrange(Gene_id)

rownames(RNAcountdat) <- RNAcountdat$Gene_id
RNAcountdat <- RNAcountdat[!(duplicated(RNAcountdat$Geneid)),]
RNAcountdat<- RNAcountdat %>%
  filter(Geneid %in% masterannotationdf$name) %>%
  arrange(Geneid)
RNAcountdat <- RNAcountdat[,-c(ncol(RNAcountdat))]
RNAcountdat <- RNAcountdat[-c(1:6)]

RNAcountdat <- RNAcountdat[,!(colnames(RNAcountdat) %like% "Elizabeth")]
RNAmetadata <- RNAmetadata[!(RNAmetadata$Person %like% "Elizabeth"),]
ddsFull <- DESeqDataSetFromMatrix(countData = RNAcountdat, colData = RNAmetadata, design = ~biological_rep + Person)

# Can run DESeq2 unaltered here to compare results
ddsFull <- DESeq(ddsFull)
RNAddsres<-results(ddsFull,contrast=c("Person","Eli","Eric"))
noeliz_resdata <- as.data.frame(RNAddsres)
noethan_fullresdata <- merge(noethan_resdata,annotationnew,by.x=0,by.y=4)
noethan_fullresdata <- merge(noethan_fullresdata,common_ids,by.x=1,by.y=1)
noethan_fullresdata <- noethan_fullresdata[!(duplicated(noethan_fullresdata$V2)),]

noeliz_fullresdata <- merge(noeliz_resdata,annotationnew,by.x=0,by.y=4)
noeliz_fullresdata <- merge(noeliz_fullresdata,common_ids,by.x=1,by.y=1)
noeliz_fullresdata <- noeliz_fullresdata[!(duplicated(noeliz_fullresdata$V2)),]

bothrna_fullresdata <- merge(both_resdata,annotationnew,by.x=0,by.y=4)

fullresdata <- noethan_fullresdata[noethan_fullresdata$chr %in% lesschrs,]
medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]

medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))

signif_medresdata <- medresdata[medresdata$padj<.01,]
notsignif_medresdata <- medresdata[medresdata$padj>=.01,]
fullresdata$chr <- factor(fullresdata$chr, levels=minichrs)

ggplot() + 
  geom_violin(data=fullresdata,trim=TRUE,aes(x=chr, y=log2FoldChange)) +
  geom_hline(yintercept=0) +
  ylab(paste0("Log2FC")) +
  geom_hline(yintercept = log2(1.5), linetype="dotted", color="blue" ) +
  theme_classic(base_size = 30) +
  ylim(-5,5) +
  geom_jitter(data = signif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=2.5, alpha=0.85,color="red",fill="red") +
  geom_jitter(data = notsignif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=2.5, alpha=0.85,color="black",fill="black") +
  stat_summary(data=fullresdata,aes(x=chr, y=log2FoldChange),fun.y=median, geom="crossbar", size=0.25, color="orange") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.title = element_blank()) +
  theme(axis.title.x = element_blank())

listInput <- list(no_ethan = na.omit(noethan_fullresdata[noethan_fullresdata$chr=="chr21" & noethan_fullresdata$padj<.01,]$Row.names), 
                  no_elizabeth = na.omit(noeliz_fullresdata[noeliz_fullresdata$chr=="chr21" & noeliz_fullresdata$padj<.01,]$Row.names), 
                  both = na.omit(bothrna_fullresdata[bothrna_fullresdata$chr=="chr21" & bothrna_fullresdata$padj<.01,]$Row.names))

upset(fromList(listInput), order.by = "freq")

### GRO-seq ####
GROseq_indir<-"DS_Normalization/counts/gro/"
GROmetadata=read.table("DS_Normalization/metadata/GROinfo.txt", sep="\t", header=TRUE)
GROmetadata$samplegroup <-paste(GROmetadata$person, "_", GROmetadata$libprep, sep="")
#read in data for GRO-seq
GRObodycountdat<- read.csv(paste0(GROseq_indir, "body.csv", sep=""))

# colnames are different due to my featurecounts script. fixing that here by removing the path.
GRObodycountdat <- GRObodycountdat[GRObodycountdat$Geneid %in% masterannotationdf$GeneID,]
GRObodycountdat <- GRObodycountdat[!(duplicated(GRObodycountdat$Geneid)),]
rownames(GRObodycountdat) <- GRObodycountdat$Geneid
GRObodycountdat <- GRObodycountdat[,-c(1:6)]
#run Deseq on RNA-seq with 21


GRObodycountdat <- GRObodycountdat[,!(colnames(GRObodycountdat) %like% "Ethan")]
GROmetadata <- GROmetadata[!(GROmetadata$person %like% "Ethan"),]

ddsFull <- DESeqDataSetFromMatrix(countData = GRObodycountdat, colData = GROmetadata, design = ~ libprep + person) #use this for all samples
dds <- collapseReplicates( ddsFull,groupby = ddsFull$samplegroup,run = ddsFull$samplegroup )
dds <- ddsFull

ddsFull <- DESeq(ddsFull)

person1 = "Eli"
person2 = "Eric"
GROddsres<-results(ddsFull,  contrast=c("person", person1, person2))
resdata <- as.data.frame(GROddsres)
fullresdata <- merge(resdata,annotationnew,by.x=0,by.y=4)
both_fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]

GRObodycountdat<- read.csv(paste0(GROseq_indir, "body.csv", sep=""))

GROmetadata=read.table("DS_Normalization/metadata/GROinfo.txt", sep="\t", header=TRUE)
GROmetadata$samplegroup <-paste(GROmetadata$person, "_", GROmetadata$libprep, sep="")
GRObodycountdat <- GRObodycountdat[GRObodycountdat$Geneid %in% masterannotationdf$GeneID,]
GRObodycountdat <- GRObodycountdat[!(duplicated(GRObodycountdat$Geneid)),]
rownames(GRObodycountdat) <- GRObodycountdat$Geneid

GRObodycountdat <- GRObodycountdat[,-c(1:6)]
#run Deseq on RNA-seq with 21

GRObodycountdat <- GRObodycountdat[,!(colnames(GRObodycountdat) %like% "Elizabeth")]
GROmetadata <- GROmetadata[!(GROmetadata$person %like% "Elizabeth"),]

ddsFull <- DESeqDataSetFromMatrix(countData = GRObodycountdat, colData = GROmetadata, design = ~ libprep + person) #use this for all samples
dds <- collapseReplicates( ddsFull,groupby = ddsFull$samplegroup,run = ddsFull$samplegroup )
dds <- ddsFull

ddsFull <- DESeq(ddsFull)

person1 = "Eli"
person2 = "Eric"
GROddsres<-results(ddsFull,  contrast=c("person", person1, person2))
resdata <- as.data.frame(GROddsres)

fullresdata <- merge(resdata,annotationnew,by.x=0,by.y=4)
noeliz_fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]

medresdata <- noeliz_fullresdata[!(is.na(noeliz_fullresdata$log2FoldChange)),]
medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
signif_medresdata <- medresdata[medresdata$padj<.01,]
notsignif_medresdata <- medresdata[medresdata$padj>=.01,]

ggplot() + 
  geom_violin(data=noeliz_fullresdata,trim=TRUE,aes(x=chr, y=log2FoldChange)) +
  #  geom_hline(yintercept=log2(0.66667),linetype="dashed",color="red") +
  geom_hline(yintercept=0) +
  ylab(paste0("Log2FC")) +
  geom_hline(yintercept = log2(1.5), linetype="dotted", color="blue" ) +
  theme_classic(base_size = 30) +
  ylim(-5,5) +
  geom_jitter(data = signif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=2.5, alpha=0.85,color="red",fill="red") +
  geom_jitter(data = notsignif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=2.5, alpha=0.85,color="black",fill="black") +
  stat_summary(data=fullresdata,aes(x=chr, y=log2FoldChange),fun.y=median, geom="crossbar", size=0.25, color="orange") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.title = element_blank()) +
  theme(axis.title.x = element_blank())

listInput <- list(no_ethan = na.omit(noethan_fullresdata[noethan_fullresdata$chr=="chr21" & noethan_fullresdata$padj<.01,]$Row.names), 
                  no_elizabeth = na.omit(noeliz_fullresdata[noeliz_fullresdata$chr=="chr21" & noeliz_fullresdata$padj<.01,]$Row.names), 
                  both = na.omit(both_fullresdata[both_fullresdata$chr=="chr21" & both_fullresdata$padj<.01,]$Row.names))

upset(fromList(listInput), order.by = "freq")


