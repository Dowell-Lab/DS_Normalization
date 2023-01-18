# Asymptotic Dispersion and Extra Poisson Noise extraction
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

### RNA ####
RNAmetadata=read.table("DS_Normalization/metadata/RNAinfo.txt", sep="\t", header=TRUE)
RNAmetadata$samplegroup <-paste(RNAmetadata$Person, "_", RNAmetadata$biological_rep, sep="")

#this is where you would remove samples if they are bad
RNAcountdat <- read.csv(RNAcoveragedat,sep="\t",skip=1)

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

RNAcountdat$Gene_id <- rownames(RNAcountdat)
RNAcountdat<- RNAcountdat %>%
  filter(Gene_id %in% annotationmerge$name) %>%
  arrange(Gene_id) %>%
  column_to_rownames('Gene_id')

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
colnames(disp_mat)=c("asymptDisp", "extraPois", "varLogDispEsts", "dispPriorVar")
for (i in 1:25)
{rdf=as.data.frame(rbind(RNAdds_shuffleremove21info(i)))
colnames(rdf)=c("asymptDisp", "extraPois", "varLogDispEsts", "dispPriorVar")
disp_mat = rbind(disp_mat, rdf)}

disp_mat$sample <-"shuffled removal chr21 genes"
realdispinfo = c(attributes(RNAdds@dispersionFunction)$coefficients[1], attributes(RNAdds@dispersionFunction)$coefficients[2], attributes(RNAdds@dispersionFunction)$varLogDispEsts, attributes(RNAdds@dispersionFunction)$dispPriorVar, "real")
disp_nochr21info = c(attributes(RNAdds_no21dds@dispersionFunction)$coefficients[1], attributes(RNAdds_no21dds@dispersionFunction)$coefficients[2], attributes(RNAdds_no21dds@dispersionFunction)$varLogDispEsts, attributes(RNAdds_no21dds@dispersionFunction)$dispPriorVar, "no chr21")
rdf = as.data.frame(rbind(realdispinfo, disp_nochr21info))
colnames(rdf)=c("asymptDisp", "extraPois", "varLogDispEsts", "dispPriorVar", "sample")
disp_mat <-rbind(disp_mat, rdf)
disp_mat$asymptDisp <- as.numeric(disp_mat$asymptDisp)
disp_mat$extraPois <- as.numeric(disp_mat$extraPois)
disp_mat$varLogDispEsts <- as.numeric(disp_mat$varLogDispEsts)
disp_mat$dispPriorVar <- as.numeric(disp_mat$dispPriorVar)

minidisp_mat <- disp_mat %>% filter(sample=="shuffled removal chr21 genes")
realdisp_mat <-disp_mat %>% filter(sample=="real")
nochr21disp_mat <-disp_mat %>% filter(sample=="no chr21")

ggplot(disp_mat, aes(x=sample, y=asymptDisp)) + 
  geom_boxplot() +
  theme_classic(base_size=18)

ggplot(disp_mat, aes(x=sample, y=extraPois)) + 
  geom_boxplot() +
  theme_classic(base_size=18)

ggplot(disp_mat, aes(x=sample, y=varLogDispEsts)) + 
  geom_boxplot() +
  theme_classic(base_size=18)

ggplot(disp_mat, aes(x=sample, y=dispPriorVar)) + 
  geom_boxplot() +
  theme_classic(base_size=18)

### GRO ###
RNAcountdat<- read.csv(paste0(GROseq_indir, "body.csv", sep=""))

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

RNAcountdat$Gene_id <- rownames(RNAcountdat)
RNAcountdat<- RNAcountdat %>%
  filter(Gene_id %in% annotationmerge$name) %>%
  arrange(Gene_id) %>%
  column_to_rownames('Gene_id')

RNAmetadata=read.table("DS_Normalization/metadata/GROinfo.txt", sep="\t", header=TRUE)
RNAmetadata$samplegroup <-paste(RNAmetadata$person, "_", RNAmetadata$libprep, sep="")

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
colnames(disp_mat)=c("asymptDisp", "extraPois", "varLogDispEsts", "dispPriorVar")

for (i in 1:25)
{rdf=as.data.frame(rbind(RNAdds_shuffleremove21info(i)))
colnames(rdf)=c("asymptDisp", "extraPois", "varLogDispEsts", "dispPriorVar")
disp_mat = rbind(disp_mat, rdf)}
disp_mat$sample <-"shuffled removal chr21 genes"
realdispinfo = c(attributes(RNAdds@dispersionFunction)$coefficients[1], attributes(RNAdds@dispersionFunction)$coefficients[2], attributes(RNAdds@dispersionFunction)$varLogDispEsts, attributes(RNAdds@dispersionFunction)$dispPriorVar, "real")
disp_nochr21info = c(attributes(RNAdds_no21dds@dispersionFunction)$coefficients[1], attributes(RNAdds_no21dds@dispersionFunction)$coefficients[2], attributes(RNAdds_no21dds@dispersionFunction)$varLogDispEsts, attributes(RNAdds_no21dds@dispersionFunction)$dispPriorVar, "no chr21")
rdf = as.data.frame(rbind(realdispinfo, disp_nochr21info))
colnames(rdf)=c("asymptDisp", "extraPois", "varLogDispEsts", "dispPriorVar", "sample")
disp_mat <-rbind(disp_mat, rdf)
disp_mat$asymptDisp <- as.numeric(disp_mat$asymptDisp)
disp_mat$extraPois <- as.numeric(disp_mat$extraPois)
disp_mat$varLogDispEsts <- as.numeric(disp_mat$varLogDispEsts)
disp_mat$dispPriorVar <- as.numeric(disp_mat$dispPriorVar)

minidisp_mat <- disp_mat %>% filter(sample=="shuffled removal chr21 genes")
realdisp_mat <-disp_mat %>% filter(sample=="real")
nochr21disp_mat <-disp_mat %>% filter(sample=="no chr21")

ggplot(disp_mat, aes(x=sample, y=asymptDisp)) + 
  geom_boxplot() +
  theme_classic(base_size=18)

ggplot(disp_mat, aes(x=sample, y=extraPois)) + 
  geom_boxplot() +
  theme_classic(base_size=18)

ggplot(disp_mat, aes(x=sample, y=varLogDispEsts)) + 
  geom_boxplot() +
  theme_classic(base_size=18)

ggplot(disp_mat, aes(x=sample, y=dispPriorVar)) + 
  geom_boxplot() +
  theme_classic(base_size=18)
