# The following script was used to generate CDF distributions of the fold-changes observed in differential analysis

library("ggplot2")
library("tidyverse")
theme_set(
  theme_classic(base_size = 18)
)
########RNA-SEQ#######
RNAindir<-"DS_Normalization/counts/rna"
GROseq_indir<-"DS_Normalization/counts/gro"
RNAcoveragedat<-"res_featureCounts_gene_idfull_143138.coverage.csv"

#Below is the bed that got used for RNA-seq. 
# To standardize RNA-seq and GRO-seq comparisons, any gene bodys/tss that were in the data more than once were removed. 
ori_worldbed <- "DS_Normalization/annotation/hg38_refseq.bed"
# Below are the bed file used for GRO-seq anaysis
# Filtered to Longest isoform with uniq body start and stop, body(+1000, -500) start and stop, tss(-500, +500) start and stop
# Also removed gene <3000bp and genes whose body or tss <0 in coordinates
# only keep genes that were analyzed by both GR0-seq and RNA-seq
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
annotationmerge

annotationmerge <- annotationmerge %>%
  arrange(GeneID) %>%
  filter(chr %in% minichrs)


########RNA-SEQ#######
####Uncorrected Real RNA-seq data####
lesschrs <- c("chr21","chr22")
masterannotationdf <- annotationmerge

masterannotationdf_only21 <- masterannotationdf %>% filter(chr=="chr21")

#this is loading the metadata
RNAmetadata=read.table("/DS_Normalization/annotation/RNAinfo.txt", sep="\t", header=TRUE)
RNAmetadata$samplegroup <-paste(RNAmetadata$Person, "_", RNAmetadata$biological_rep, sep="")

RNAcountdat <- read.csv("/DS_Normalization/counts/featurecounts_rna_withmultis.sense.txt",
                        sep="\t", skip=1)

RNAcountdat<- RNAcountdat %>%
    filter(Geneid %in% annotationmerge$name) %>%
  arrange(Geneid)

# Filter to maximal isoform
RNAcountdat$sum <- rowSums(RNAcountdat[,7:ncol(RNAcountdat)])/RNAcountdat$Length
RNAcountdat <- merge(RNAcountdat,common_ids,by.x=1,by.y=1)
RNAcountdat <- (RNAcountdat[order(RNAcountdat$V2,-RNAcountdat$sum),])
RNAcountdat <- (RNAcountdat[!(duplicated(RNAcountdat$V2)),])
RNAcountdat <- subset(RNAcountdat, select=-c(V2,sum))
nrow(RNAcountdat)

# Further filter to those common between RNA and GRO
RNAcountdat<- RNAcountdat %>%
  filter(Geneid %in% masterannotationdf$name) %>%
  arrange(Geneid)
rownames(RNAcountdat) <- RNAcountdat$Geneid

# If there are some duplicated entries
RNAcountdat <- RNAcountdat[!(duplicated(RNAcountdat$Geneid)),]
rownames(RNAcountdat) <- RNAcountdat$Geneid


#run Deseq on RNA-seq with 21
RNAcountdat <- RNAcountdat[,-c(1:6)]
ddsFull <- DESeqDataSetFromMatrix(countData = RNAcountdat, colData = RNAmetadata, design = ~biological_rep + Person)

# Can run DESeq2 unaltered here to compare results
ddsFull <- DESeq(ddsFull)
RNAddsres<-results(ddsFull,contrast=c("Person","Ethan","Eric"))
resdata <- as.data.frame(RNAddsres)
resdata$geneid <- rownames(resdata)
ncol(resdata)
rownames(annotationmerge) <- annotationmerge$GeneID

fullresdata <- merge(resdata,annotationmerge,by.x=0,by.y=0)
RNA_uncorrected_fullresdata <- fullresdata
fullresdata <- fullresdata[fullresdata$chr %in% minichrs,]

# fetching medians data

medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]
medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
fullresdata$chr <- factor(fullresdata$chr, levels=minichrs)

# Draw ECDF Plot
CDF <- ecdf(fullresdata[fullresdata$chr=="chr21",]$log2FoldChange)
CDF
# draw CDF Plot
fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]

ggplot(fullresdata, aes(x = log2FoldChange,color=chr)) +
  stat_ecdf() +
  xlim(-2,2) +
  geom_vline(xintercept = log2(1.5),color="red",linetype="dotted") +
  geom_vline(xintercept=log2(1.0),color="blue",linetype="dotted") +
  geom_vline(xintercept=log2(1.5)-(2*sd(na.omit(fullresdata[fullresdata$chr!="chr21",]$log2FoldChange))),
             color="purple",linetype="dotted") +
  theme_classic()

#### Uncorrected GRO-seq ####
#### UNCORRECTED REAL DATA GRO-SEQ RUN####
#traditional Deseq2 on on GRO count matrix's

GROseq_indir<-"DS_Normalization/counts/gro"
GROmetadata=read.table("DS_Normalization/metadata/GROinfo.txt", sep="\t", header=TRUE)
GROmetadata$samplegroup <-paste(GROmetadata$person, "_", GROmetadata$libprep, sep="")

########RNA-SEQ#######
RNAindir<-"DS_Normalization/counts/rna"
GROseq_indir<-"DS_Normalization/counts/gro"
RNAcoveragedat<-"res_featureCounts_gene_idfull_143138.coverage.csv"

#Below is the bed that got used for RNA-seq. 
# To standardize RNA-seq and GRO-seq comparisons, any gene bodys/tss that were in the data more than once were removed. 
ori_worldbed <- "DS_Normalization/annotation/hg38_refseq.bed"
# Below are the bed file used for GRO-seq anaysis
# Filtered to Longest isoform with uniq body start and stop, body(+1000, -500) start and stop, tss(-500, +500) start and stop
# Also removed gene <3000bp and genes whose body or tss <0 in coordinates
# only keep genes that were analyzed by both GR0-seq and RNA-seq
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
annotationmerge <-merge(annotationRNA, annotationGRO, all = FALSE,by="row.names")
row.names(annotationmerge)<-annotationmerge$GeneID

#keep only the genes on the main chromosomes
minichrs <-c("chr1", "chr2", "chr3","chr4", "chr5", "chr6", "chr7", "chr8", "chr9","chr10","chr11", "chr12", "chr13","chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")

library(dplyr)

annotationmerge <-annotationmerge %>%
  arrange(GeneID) %>%
  filter(chr %in% minichrs)
#read in data for GRO-seq
GRObodycountdat<- read.csv(paste0(GROseq_indir, "body.csv", sep=""))
GRObodycountdat$Gene_id <- rownames(GRObodycountdat)
GRObodycountdat<- GRObodycountdat %>%
  filter(Gene_id %in% annotationmerge$name) %>%
  arrange(Gene_id)


rownames(GRObodycountdat) <- GRObodycountdat$Gene_id
GRObodycountdat <- merge(GRObodycountdat,annotationmerge,by.x=0,by.y=3)
colnames(GRObodycountdat)
GRObodycountdat$sum <- rowSums(GRObodycountdat[,2:17])/GRObodycountdat$Length
GRObodycountdat <- merge(GRObodycountdat,common_ids,by.x=1,by.y=1)
GRObodycountdat <- (GRObodycountdat[order(GRObodycountdat$V2,-GRObodycountdat$sum),])
GRObodycountdat <- (GRObodycountdat[!(duplicated(GRObodycountdat$V2)),])
GRObodycountdat <- subset(GRObodycountdat, select=-c(V2,sum))
rownames(GRObodycountdat) <- GRObodycountdat$Row.names
GRObodycountdat <- GRObodycountdat[,2:17]

GRObodyddsFull <- DESeqDataSetFromMatrix(countData = GRObodycountdat, colData = GROmetadata, design = ~libprep+ person)
GRObodydds <- collapseReplicates( GRObodyddsFull,groupby = GRObodyddsFull$samplegroup,run = GRObodyddsFull$samplegroup)
GRObodydds <-DESeq(GRObodydds)

person1 = "T21"
person2 = "D21"

GROddstestres<-results(GRObodydds,  contrast=c("person", person1, person2))
resdata <- as.data.frame(GROddstestres)
fullresdata <- merge(resdata,annotationmerge,by.x=0,by.y=1)
fullresdata <- fullresdata[fullresdata$chr %in% minichrs,]

# draw CDF Plot
fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]

ggplot(fullresdata, aes(x = log2FoldChange,color=chr)) +
  stat_ecdf() +
  geom_vline(xintercept = log2(1.50),color="red",linetype="dotted") +
  geom_vline(xintercept=log2(1.0),color="blue",linetype="dotted") +
  xlim(-2,2) +
  geom_vline(xintercept=log2(1.5)-(2*sd(fullresdata[fullresdata$chr!="chr21" & !(is.na(fullresdata$log2FoldChange)),]$log2FoldChange))
             ,color="purple",linetype="dotted") +
  theme_classic()


