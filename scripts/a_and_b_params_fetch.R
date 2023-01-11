# Asymptotic Dispersion and Extra Poisson Noise
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


### RNA ####
### Various trial things here... looking at hyperparameters of the data

RNAcountdat <- read.csv("/Users/sahu0957/trisomy_normalization/RNAGRO01132021/data/RNAseqfiles/res_featureCounts_gene_idfull_143138.coverage.csv", 
                        sep=",", row.names=1)

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

### GRO ###

RNAcountdat <- read.csv("/Users/sahu0957/trisomy_normalization/RNAGRO01132021/data/GROfiles/body.csv", 
                        sep=",", row.names=1)

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

head(RNAcountdat)
RNAcountdat$Gene_id <- rownames(RNAcountdat)
RNAcountdat<- RNAcountdat %>%
  filter(Gene_id %in% annotationmerge$name) %>%
  arrange(Gene_id) %>%
  column_to_rownames('Gene_id')
nrow(RNAcountdat)
nrow(chr21genes)
#View(chr21genes)
RNAmetadata=read.table("/Shares/down/mixed/RNAGRO01132021/scripts/Deseq2_analysis/GROinfo.txt", sep="\t", header=TRUE)
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

ggsave("groseq_asymptdisp.svg",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "svg")
ggsave("groseq_asymptdisp.png",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "png")


ggplot(disp_mat, aes(x=sample, y=extraPois)) + 
  geom_boxplot() +
  theme_classic(base_size=18)

ggsave("groseq_extrapois.svg",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "svg")
ggsave("groseq_extrapois.png",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "png")

ggplot(disp_mat, aes(x=sample, y=varLogDispEsts)) + 
  geom_boxplot() +
  theme_classic(base_size=18)

ggplot(disp_mat, aes(x=sample, y=dispPriorVar)) + 
  geom_boxplot() +
  theme_classic(base_size=18)
