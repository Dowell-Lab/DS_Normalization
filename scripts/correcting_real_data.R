####Correcting Real Data ####
# Which genes are tested is very important- right now, filtered to those with signal in RNA and GRO (314 chr21 genes)
# In the paper, I've also filtered isoforms (262 chr21 genes). At the end of the day, the percentages remain the same


########RNA-SEQ#######
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
annotationmerge

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

common_ids <- read.table("/Users/sahu0957/backup/genome_files/refseq_to_common_id.txt",sep="\t")
common_ids <- common_ids[common_ids$V1 %in% annotationmerge$GeneID,]
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
RNAmetadata$samplegroup <-paste(RNAmetadata$Person, "_", RNAmetadata$biological_rep, sep="")

#this is where you would remove samples if they are bad
#RNAcountdat <- read.csv(RNAcoveragedat,sep="\t",skip=1)
RNAcountdat <- read.csv("/scratch/Users/sahu0957/ds_normalization/RNA/nextflowhg38_20220214/counts/featurecounts_rna_withmultis.sense.txt",
                        sep="\t", skip=1)
# Maximal isoform filter, if you wanna
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

#run Deseq on RNA-seq with 21
RNAcountdat <- RNAcountdat[,-c(1:6)]
head(RNAcountdat)
ddsFull <- DESeqDataSetFromMatrix(countData = RNAcountdat, colData = RNAmetadata, design = ~biological_rep + Person)

# Can run DESeq2 unaltered here to compare results
ddsFull <- DESeq(ddsFull)
RNAddsres<-results(ddsFull,contrast=c("Person","Ethan","Eric"))
resdata <- as.data.frame(RNAddsres)
fullresdata <- merge(resdata,annotationmerge,by.x=0,by.y=3)
write.csv(fullresdata,"/scratch/Shares/dowell/for_marya/ds_normalization/deseq2/uncorrected_RNA_results.csv")

fullresdata <- merge(fullresdata,common_ids,by.x=1,by.y=1)
# If we can't get a FC estimate on it, just remove it (genes with no reads, for example, shouldn't be called 0 LFC)
fullresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]

RNA_uncorrected_fullresdata <- fullresdata
nrow(RNA_uncorrected_fullresdata)

lesschrs
fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]

medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]

library(plyr)

# Depending on which genes you pick, RNA-seq MFC ranges from ~1.34-1.41.
medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians

signif_medresdata <- medresdata[medresdata$padj<.01,]
notsignif_medresdata <- medresdata[medresdata$padj>=.01,]
fullresdata$chr <- factor(fullresdata$chr, levels=minichrs)

fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]

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

# Correcting things

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


# Maximal isoform filter, if you wanna
RNAcountdat$sum <- rowSums(RNAcountdat[,7:ncol(RNAcountdat)])/RNAcountdat$Length
RNAcountdat <- merge(RNAcountdat,common_ids,by.x=1,by.y=1)
RNAcountdat <- (RNAcountdat[order(RNAcountdat$V2,-RNAcountdat$sum),])
RNAcountdat <- (RNAcountdat[!(duplicated(RNAcountdat$V2)),])
RNAcountdat <- subset(RNAcountdat, select=-c(V2,sum))
nrow(RNAcountdat)

# We've got 1300 duplicate entries... How was this filtered out before?
#RNAcountdat <- RNAcountdat[!(duplicated(RNAcountdat$Geneid)),]
# Filter to genes common in both RNA and GRO seq

library(tibble)

# Depending on which genes you use, you'll get different FC results. Can either keep all or only those in both RNA/GRO files, for example

RNAcountdat<- RNAcountdat %>%
  filter(Geneid %in% annotationmerge$GeneID) %>%
  arrange(Geneid) %>%
  column_to_rownames('Geneid')
colnames(RNAcountdat)
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
fullresdata <- merge(nomultisresdata,annotationmerge,by.x=0,by.y=3)
medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]

medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians

signif_medresdata <- medresdata[medresdata$padj<.01,]
notsignif_medresdata <- medresdata[medresdata$padj>=.01,]
nrow(signif_medresdata)
nrow(notsignif_medresdata)

fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]

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
norm_data <- as.data.frame(normalizationFactors(ddsCollapsed_normfactor))
norm_dataframe_withannos <- merge(confused,masterannotationdf,by.x=0,by.y=0)


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
fullresdata <- merge(resdata,annotationmerge,by.x=0,by.y=3)
fullresdata <- fullresdata[fullresdata$chr %in% minichrs,]
nrow(fullresdata)
quantile(fullresdata$baseMean, probs = seq(0,1,0.2))
#View(fullresdata)
fullresdata <- fullresdata[fullresdata$baseMean > 12,] # Remove genes with very low expression
#fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]
medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]
common_fullresdata_rna <- merge(fullresdata,common_ids,by.x=1,by.y=1)
write.csv(fullresdata,"/scratch/Shares/dowell/for_marya/ds_normalization/deseq2/normalization_matrix_basemeanfilter_RNA_results.csv")

common_ids


medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians
colnames(fullresdata)

signif_medresdata <- (common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$padj<.01,])

notsignif_medresdata <- (common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$padj>=.01,])
View(signif_medresdata[signif_medresdata$chr=="chr21" & signif_medresdata$log2FoldChange< 0,])
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

fullresdata
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
RNAddsres<-results(ddsFull,  contrast=c("Person", person1, person2),altHypothesis="lessAbs",alpha=.1,lfcThreshold=log2(1.5))
resdata <- as.data.frame(RNAddsres)
fullresdata <- merge(resdata,annotationmerge,by.x=0,by.y=3)

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

#### GRO-SEQ####
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
nrow(GRObodycountdat)
View(annotationmerge)
# Annoyingly, this file doesn't include annotations -_-
GRObodycountdat <- merge(GRObodycountdat,annotationmerge,by.x=0,by.y=3)
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

person1 = "Ethan"
person2 = "Eric"

GROddstestres<-results(GRObodydds,  contrast=c("person", person1, person2))
resdata <- as.data.frame(GROddstestres)
fullresdata <- merge(resdata,annotationmerge,by.x=0,by.y=3)
#write.csv(fullresdata,"/scratch/Shares/dowell/for_marya/ds_normalization/deseq2/uncorrected_GRO_results.csv")

#fullresdata <- merge(fullresdata,common_ids,by.x=1,by.y=1)
GRO_uncorrected_fullresdata <- fullresdata

fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]

medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]
fullresdata
common_fullresdata <- merge(fullresdata,common_ids,by.x=1,by.y=1)
common_fullresdata
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

ggsave("groseq_deseq2_violins_uncorrected_chr2122_med138_med104.svg",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "svg")
ggsave("groseq_deseq2_violins_uncorrected_chr2122_med138_med104.png",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "png")


### GRO-seq nopromoters, no multimaps, ploidy correction ####
GRObodycountdat<- read.csv("/scratch/Users/sahu0957/ds_normalization/GRO/no_repeats/counts/featurecounts_gro_noproms_nomultis.sense.txt",
                           skip=1, sep="\t")

# colnames are different due to my featurecounts script. fixing that here by removing the path. Kinda hacky, go remove in script!
names <- colnames(GRObodycountdat)
names <- gsub(".*E","E",names)
nrow(GRObodycountdat[GRObodycountdat$Chr=="chr21",])
colnames(GRObodycountdat) <- names
colnames(GRObodycountdat)
GRObodycountdat$sum <- rowSums(GRObodycountdat[,7:ncol(GRObodycountdat)])/GRObodycountdat$Length
GRObodycountdat <- merge(GRObodycountdat,common_ids,by.x=1,by.y=1)
GRObodycountdat <- (GRObodycountdat[order(GRObodycountdat$V2,-GRObodycountdat$sum),])
GRObodycountdat <- (GRObodycountdat[!(duplicated(GRObodycountdat$V2)),])
GRObodycountdat <- subset(GRObodycountdat, select=-c(V2,sum))
head(GRObodycountdat)
rownames(GRObodycountdat) <- GRObodycountdat$Geneid

GRObodycountdat <- GRObodycountdat[,-c(1:6)]
nomultis_col
GROmetadata=read.table("/Shares/down/mixed/RNAGRO01132021/scripts/Deseq2_analysis/GROinfo.txt", sep="\t", header=TRUE)
GROmetadata$samplegroup <-paste(GROmetadata$person, "_", GROmetadata$libprep, sep="")
GROmetadata$nomultis <- nomultis_col
minichrs <-c("chr1", "chr2", "chr3","chr4", "chr5", "chr6", "chr7", "chr8", "chr9","chr10","chr11", "chr12", "chr13","chr14", "chr15", "chr16", "chr17", "chr18", "chr19","chr20","chr21","chr22")

masterannotationdf <- annotationmerge
masterannotationdf_only21 <- masterannotationdf %>% filter(chr=="chr21")
anuplodygenes <-as.vector(masterannotationdf_only21[["name"]]) #Only chr21 genes
anuplodygenes
baseploidy <- 2
alt_ploidy <-3

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
fullresdata <- merge(resdata,annotationmerge,by.x=0,by.y=3)
medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]

medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians
controlgenes <-ifelse(rownames(dds) %in% anuplodygenes, FALSE, TRUE) # Only normalize on non-aneuploidy genes
samplepinfo<-as.data.frame(colData(dds))
ploidy_trisomy21 = ifelse(rownames(dds) %in% anuplodygenes, alt_ploidy, baseploidy)
ploidy_typical <- rep(baseploidy, nrow(dds))
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

person1 = "Eli"
person2 = "Elizabeth"
RNAddsres<-results(ddsCollapsed_normfactor,  contrast=c("person", person1, person2))
resdata <- as.data.frame(RNAddsres)
fullresdata <- merge(resdata,annotationmerge,by.x=0,by.y=3)
commonfullresdata <- merge(fullresdata,common_ids,by.x=1,by.y=1)
#fullresdata <- fullresdata[!(duplicated(fullresdata$V2)),]
#fullresdata <- na.omit(fullresdata)

nrow(fullresdata)
quantile(fullresdata$baseMean, probs = seq(0,1,0.2))

commonfullresdata <- commonfullresdata[commonfullresdata$baseMean > 12,] # Remove genes with very low expression
View(commonfullresdata[commonfullresdata$padj<.1 & commonfullresdata$chr=="chr21",])

fullresdata <- fullresdata[fullresdata$baseMean > 12,] # Remove genes with very low expression
#write.csv(fullresdata,"/scratch/Shares/dowell/for_marya/ds_normalization/deseq2/normmatrix_basemean12_GRO_results.csv")
fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]

medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]
medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians
signif_medresdata <- medresdata[medresdata$padj<.01,]
notsignif_medresdata <- medresdata[medresdata$padj>=.01,]
View(signif_medresdata[signif_medresdata$log2FoldChange<0 & signif_medresdata$chr=="chr21",])

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


ggsave("groseq_deseq2_violins_corrected_chr2122_med101_med104.svg",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "svg")
ggsave("groseq_deseq2_violins_corrected_chr2122_med101_med104.png",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "png")


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
#run Deseq on GRO-seq with 21
ddsFull <- DESeqDataSetFromMatrix(countData = GRObodycountdat, colData = GROmetadata, design = ~ libprep + person) #use this for all samples
#ddsFull$type <- relevel(ddsFull$type, "WT")
ddsFull <- collapseReplicates( ddsFull,groupby = ddsFull$samplegroup,run = ddsFull$samplegroup )
ddsFull <- DESeq(ddsFull)

GROddsres<-results(ddsFull,  contrast=c("person", person1, person2),lfcThreshold = log2(1.5),altHypothesis = "lessAbs")
resdata <- as.data.frame(GROddsres)
fullresdata <- merge(resdata,annotationmerge,by.x=0,by.y=3)
fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]
medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]
medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))

signif_medresdata <- medresdata[medresdata$padj<.1,]
notsignif_medresdata <- medresdata[medresdata$padj>=.1,]

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

