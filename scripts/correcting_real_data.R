####Correcting Real Data ####
# This script accounts for trisomy using two different strategies, such that differential analysis is more applicable for these samples
# First strategy: read counts over genes are divided by their ploidy number, such that FCs over all chromosomes is expected to
# be a normal distribution around LFC=0. Significantly low genes are then interpretable as below the expected value (i.e., dosage compensated)
# Second, a new alternative hypothesis is utilized. LFC shrinkage is inapplicable in this instance

library(dplyr)
library(plyr)
library(ggplot2)
library(tibble)
########RNA-SEQ#######
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
RNAmetadata=read.table("DS_Normalization/metadata/RNAinfo.txt", sep="\t", header=TRUE)
RNAmetadata$samplegroup <-paste(RNAmetadata$Person, "_", RNAmetadata$biological_rep, sep="")

#this is where you would remove samples if they are bad
RNAcountdat <- read.csv(RNAcoveragedat,sep="\t",skip=1)

# Filter by isoform based on Common IDs
RNAcountdat$sum <- rowSums(RNAcountdat[,7:ncol(RNAcountdat)])/RNAcountdat$Length
RNAcountdat <- merge(RNAcountdat,common_ids,by.x=1,by.y=1)
RNAcountdat <- (RNAcountdat[order(RNAcountdat$V2,-RNAcountdat$sum),])
RNAcountdat <- (RNAcountdat[!(duplicated(RNAcountdat$V2)),])
RNAcountdat <- subset(RNAcountdat, select=-c(V2,sum))
# Further filter to those common between RNA and GRO
RNAcountdat<- RNAcountdat %>%
    filter(Geneid %in% masterannotationdf$name) %>%
  arrange(Geneid)
rownames(RNAcountdat) <- RNAcountdat$Geneid

#run Deseq on RNA-seq with 21
RNAcountdat <- RNAcountdat[,-c(1:6)]
ddsFull <- DESeqDataSetFromMatrix(countData = RNAcountdat, colData = RNAmetadata, design = ~biological_rep + Person)

# Can run DESeq2 unaltered here to compare results
# We use pseudonyms to identify each sample:
# Eli: Father
# Elizabeth: Mother
# Eric: Son with D21
# Ethan: Son with T21
ddsFull <- DESeq(ddsFull)
RNAddsres<-results(ddsFull,contrast=c("Person","Ethan","Eric"))
resdata <- as.data.frame(RNAddsres)
fullresdata <- merge(resdata,annotationmerge,by.x=0,by.y=3)
fullresdata <- merge(fullresdata,common_ids,by.x=1,by.y=1)
# If we can't get a FC estimate on it, just remove it (genes with no reads, for example, shouldn't be called 0 LFC)
fullresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]

RNA_uncorrected_fullresdata <- fullresdata
fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]
medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]

# Calculate median FC by chromosome. Depending on your filter, RNA-seq results in a FC of 1.34-1.41
medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))

#Padj cutoff. Adjust as desired
padjusted_cutoff <- 0.01
signif_medresdata <- medresdata[medresdata$padj < padjusted_cutoff,]
notsignif_medresdata <- medresdata[medresdata$padj >= padjusted_cutoff,]


fullresdata$chr <- factor(fullresdata$chr, levels=minichrs)

# For easier visualization, filter to only two chromosomes (21 and 22)
fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]

# Violin plot visualization
ggplot() + 
  geom_violin(data=fullresdata,trim=TRUE,aes(x=chr, y=log2FoldChange)) +
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

#this is where you would remove samples if they are bad
RNAcountdat <- read.csv("DS_Normalization/counts/rna/featurecounts_rna_rmsked.sense.txt",
                        sep="\t", skip=1)

# filter to only main chrs
RNAcountdat <- RNAcountdat[RNAcountdat$Chr %in% minichrs,]


# Remove long colnames to be more readable
names <- colnames(RNAcountdat)
names <- gsub(".*removed..","",names)
colnames(RNAcountdat) <- names

# Maximal isoform filter, based on Common ID
RNAcountdat$sum <- rowSums(RNAcountdat[,7:ncol(RNAcountdat)])/RNAcountdat$Length
RNAcountdat <- merge(RNAcountdat,common_ids,by.x=1,by.y=1)
RNAcountdat <- (RNAcountdat[order(RNAcountdat$V2,-RNAcountdat$sum),])
RNAcountdat <- (RNAcountdat[!(duplicated(RNAcountdat$V2)),])
RNAcountdat <- subset(RNAcountdat, select=-c(V2,sum))

RNAcountdat<- RNAcountdat %>%
  filter(Geneid %in% annotationmerge$GeneID) %>%
  arrange(Geneid) %>%
  column_to_rownames('Geneid')

# Remove annotation info from count file
RNAcountdat <- RNAcountdat[,-c(1:5)]
#run Deseq on RNA-seq with 21
colnames(RNAcountdat)
ddsFull <- DESeqDataSetFromMatrix(countData = RNAcountdat, colData = RNAmetadata, design = ~ biological_rep + Person) #use this for all samples
dds <- ddsFull
ddsFull <- DESeq(ddsFull)

person1 = "Ethan"
person2 = "Eric"
RNAddsres<-results(ddsFull,  contrast=c("Person", person1, person2))
nomultisresdata <- as.data.frame(RNAddsres)
fullresdata <- merge(nomultisresdata,annotationmerge,by.x=0,by.y=3)
medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]

medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians

padjusted_cutoff <- 0.01
signif_medresdata <- medresdata[medresdata$padj<padjusted_cutoff,]
notsignif_medresdata <- medresdata[medresdata$padj>=padjusted_cutoff,]

fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]

ggplot() + 
  geom_violin(data=fullresdata,trim=TRUE,aes(x=chr, y=log2FoldChange)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_hline(yintercept=0) +
  ylab(paste0("Log2FC")) +
  geom_hline(yintercept = log2(1.5), linetype="dotted", color="blue" ) +
  xlab("Chromosome") +
  theme_classic(base_size = 24) +
  geom_jitter(data = signif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=1.5, alpha=0.85,color="red",fill="red") +
  geom_jitter(data = notsignif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=1.5, alpha=0.85,color="black",fill="black") +
  stat_summary(data=fullresdata,aes(x=chr, y=log2FoldChange),fun.y=median, geom="crossbar", size=0.25, color="orange") +
  theme(legend.title = element_blank())


# Normalizing by ploidy. This will make chr21 FC distribution comparable to disomic chromosomes, to assess if there is truly
# dosage compensation or if low FC genes are just expected by chance
controlgenes <-ifelse(rownames(dds) %in% anuplodygenes, FALSE, TRUE) # Only normalize on non-aneuploidy genes
samplepinfo<-as.data.frame(colData(dds))

ploidy_trisomy21 = ifelse(rownames(dds) %in% anuplodygenes, alt_ploidy, baseploidy)
ploidy_typical <- rep(baseploidy, nrow(dds))

normFactors <- matrix(do.call(cbind, mget(paste0(dds$ploidy))),ncol=ncol(dds),nrow=nrow(dds),dimnames=list(1:nrow(dds),1:ncol(dds)))
normFactors <- normFactors/baseploidy

ddsCollapsed<-DESeq(dds,betaPrior = FALSE)
ddsCollapsed_normfactor <-estimateSizeFactors(dds,normMatrix=normFactors, controlGenes=controlgenes) #adds a column to colData(ddsCollapsed) called sizeFactor if you use normmatrix then normalizationFactors(ddsCollapsed) instead.
norm_data <- as.data.frame(normalizationFactors(ddsCollapsed_normfactor))
norm_dataframe_withannos <- merge(confused,masterannotationdf,by.x=0,by.y=0)


### Correcting for ploidy should yield a distribution tightly around 1.0
ddsCollapsed_normfactor <- estimateDispersionsGeneEst(ddsCollapsed_normfactor) #adds mu to assays(ddsCollapsed) mu-has one column per sample and fills in 
ddsCollapsed_normfactor<-estimateDispersionsFit(ddsCollapsed_normfactor) # and fills in dispersionFunction(ddsCollapsed) and adds a column in  elementMetadata(ddsCollapsed) called dispFit
ddsCollapsed_normfactor <- estimateDispersionsMAP(ddsCollapsed_normfactor) #adds a attr(,"dispPriorVar") to dispersionFunction(ddsCollapsed) 
ddsCollapsed_normfactor <- nbinomWaldTest(ddsCollapsed_normfactor,betaPrior = FALSE) #adds both H and cooks to assays(ddsCollapsed)  both are for each sample you have.

person1 = "Ethan"
person2 = "Eric"
RNAddsres<-results(ddsCollapsed_normfactor,  contrast=c("Person", person1, person2))
resdata <- as.data.frame(RNAddsres)
fullresdata <- merge(resdata,annotationmerge,by.x=0,by.y=3)
fullresdata <- fullresdata[fullresdata$chr %in% minichrs,]

# Adding in a minimum baseMean filter can be helpful for removing false positives (noisy genes). Can be done pre- or post-differential analysis
fullresdata <- fullresdata[fullresdata$baseMean > 30,] # Remove genes with very low expression
medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]
common_fullresdata_rna <- merge(fullresdata,common_ids,by.x=1,by.y=1)

medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
signif_medresdata <- (common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$padj<.01,])
notsignif_medresdata <- (common_fullresdata_rna[!(duplicated(common_fullresdata_rna$V2)) & common_fullresdata_rna$padj>=.01,])

ggplot() + 
  geom_violin(data=fullresdata,trim=TRUE,aes(x=chr, y=log2FoldChange)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_hline(yintercept=0) +
  ylab(paste0("Log2FC")) +
  xlab("Chromosome") +
  theme_classic(base_size = 24) +
  geom_jitter(data = signif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=1.5, alpha=0.85,color="red",fill="red") +
  geom_jitter(data = notsignif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=1.5, alpha=0.85,color="black",fill="black") +
  stat_summary(data=fullresdata,aes(x=chr, y=log2FoldChange),fun.y=median, geom="crossbar", size=0.25, color="orange") +
  theme(legend.title = element_blank())

#### Change the null hypothesis- Real RNA data ####

#this is where you would remove samples if they are bad
RNAcountdat <- read.csv("DS_Normalization/counts/rna/featurecounts_rna_rmsked.sense.txt",
                        sep="\t", skip=1)

# filter to only canonical chrs
RNAcountdat <- RNAcountdat[RNAcountdat$Chr %in% minichrs,]

names <- colnames(RNAcountdat)
names <- gsub(".*removed..","",names)
colnames(RNAcountdat) <- names

RNAcountdat<- RNAcountdat %>%
  filter(Geneid %in% annotationmerge$GeneID) %>%
  arrange(Geneid) %>%
  column_to_rownames('Geneid')

RNAcountdat <- RNAcountdat[,-c(1:5)]
#run Deseq on RNA-seq with 21
ddsFull <- DESeqDataSetFromMatrix(countData = RNAcountdat, colData = RNAmetadata, design = ~ biological_rep + Person) #use this for all samples
dds <- ddsFull
ddsFull <- DESeq(ddsFull,)

person1 = "Ethan"
person2 = "Eric"

# Adjusting alternative hypothesis. We identify significant genes using altHypothesis LFC < |log2(1.5)|
RNAddsres<-results(ddsFull,  contrast=c("Person", person1, person2),altHypothesis="lessAbs",alpha=.01,lfcThreshold=log2(1.5))
resdata <- as.data.frame(RNAddsres)
fullresdata <- merge(resdata,annotationmerge,by.x=0,by.y=3)

# Make a violin
fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]

medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]

common_fullresdata_rna <- merge(medresdata,common_ids,by.x=1,by.y=1)
medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))

signif_medresdata <- medresdata[medresdata$padj < padjusted_cutoff,]
notsignif_medresdata <- medresdata[medresdata$padj >= padjusted_cutoff,]

ggplot() + 
  geom_violin(data=fullresdata,trim=TRUE,aes(x=chr, y=log2FoldChange)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
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

GROseq_indir<-"DS_Normalization/counts/gro/"
GROmetadata=read.table("DS_Normalization/metadata/GROinfo.txt", sep="\t", header=TRUE)
GROmetadata$samplegroup <-paste(GROmetadata$person, "_", GROmetadata$libprep, sep="")


#read in data for GRO-seq
GRObodycountdat<- read.csv(paste0(GROseq_indir, "body.csv", sep=""))

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
GRO_uncorrected_fullresdata <- fullresdata

fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]
medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]
common_fullresdata <- merge(fullresdata,common_ids,by.x=1,by.y=1)
medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
signif_medresdata <- medresdata[medresdata$padj < padjusted_cutoff,]
notsignif_medresdata <- medresdata[medresdata$padj >= padjusted_cutoff,]

ggplot() + 
  geom_violin(data=fullresdata,trim=TRUE,aes(x=chr, y=log2FoldChange)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_hline(yintercept=0) +
  ylab(paste0("Log2FC")) +
  geom_hline(yintercept = log2(1.5), linetype="dotted", color="blue" ) +
  xlab("Chromosome") +
  theme_classic(base_size = 24) +
  geom_jitter(data = signif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=1.5, alpha=0.85,color="red",fill="red") +
  geom_jitter(data = notsignif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=1.5, alpha=0.85,color="black",fill="black") +
  stat_summary(data=fullresdata,aes(x=chr, y=log2FoldChange),fun.y=median, geom="crossbar", size=0.25, color="orange") +
  theme(legend.title = element_blank())


### GRO-seq nopromoters, no multimaps, ploidy correction ####
GRObodycountdat<- read.csv("DS_Normalization/counts/gro/featurecounts_gro_noproms_nomultis.sense.txt",
                           skip=1, sep="\t")

# colnames are different due to my featurecounts script. fixing that here by removing the path. Kinda hacky, go remove in script!
names <- colnames(GRObodycountdat)
names <- gsub(".*E","E",names)
colnames(GRObodycountdat) <- names
GRObodycountdat$sum <- rowSums(GRObodycountdat[,7:ncol(GRObodycountdat)])/GRObodycountdat$Length
GRObodycountdat <- merge(GRObodycountdat,common_ids,by.x=1,by.y=1)
GRObodycountdat <- (GRObodycountdat[order(GRObodycountdat$V2,-GRObodycountdat$sum),])
GRObodycountdat <- (GRObodycountdat[!(duplicated(GRObodycountdat$V2)),])
GRObodycountdat <- subset(GRObodycountdat, select=-c(V2,sum))
rownames(GRObodycountdat) <- GRObodycountdat$Geneid

GRObodycountdat <- GRObodycountdat[,-c(1:6)]
GROmetadata=read.table("DS_Normalization/metadata/GROinfo.txt", sep="\t", header=TRUE)
GROmetadata$samplegroup <-paste(GROmetadata$person, "_", GROmetadata$libprep, sep="")
GROmetadata$nomultis <- nomultis_col

#run Deseq on GRO-seq with 21
ddsFull <- DESeqDataSetFromMatrix(countData = GRObodycountdat, colData = GROmetadata, design = ~ libprep + person) #use this for all samples
dds <- collapseReplicates( ddsFull,groupby = ddsFull$samplegroup,run = ddsFull$samplegroup )

dds <- ddsFull
ddsFull <- DESeq(ddsFull)

person1 = "Ethan"
person2 = "Eric"
GROddsres<-results(ddsFull,  contrast=c("person", person1, person2))
resdata <- as.data.frame(GROddsres)
fullresdata <- merge(resdata,annotationmerge,by.x=0,by.y=3)
medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]

medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
controlgenes <-ifelse(rownames(dds) %in% anuplodygenes, FALSE, TRUE) # Only normalize on non-aneuploidy genes
samplepinfo<-as.data.frame(colData(dds))
ploidy_trisomy21 = ifelse(rownames(dds) %in% anuplodygenes, alt_ploidy, baseploidy)
ploidy_typical <- rep(baseploidy, nrow(dds))
normFactors <- matrix(do.call(cbind, mget(paste0(dds$ploidy))),ncol=ncol(dds),nrow=nrow(dds),dimnames=list(1:nrow(dds),1:ncol(dds)))
normFactors <- normFactors/baseploidy

ddsCollapsed<-DESeq(dds)
ddsCollapsed_normfactor <-estimateSizeFactors(dds,normMatrix=normFactors, controlGenes=controlgenes) #adds a column to colData(ddsCollapsed) called sizeFactor if you use normmatrix then normalizationFactors(ddsCollapsed) instead.
ddsCollapsed_normfactor <- estimateDispersionsGeneEst(ddsCollapsed_normfactor) #adds mu to assays(ddsCollapsed) mu-has one column per sample and fills in 
ddsCollapsed_normfactor<-estimateDispersionsFit(ddsCollapsed_normfactor) # and fills in dispersionFunction(ddsCollapsed) and adds a column in  elementMetadata(ddsCollapsed) called dispFit
ddsCollapsed_normfactor <- estimateDispersionsMAP(ddsCollapsed_normfactor) #adds a attr(,"dispPriorVar") to dispersionFunction(ddsCollapsed) 
ddsCollapsed_normfactor <- nbinomWaldTest(ddsCollapsed_normfactor) #adds both H and cooks to assays(ddsCollapsed)  both are for each sample you have.


person1 = "Ethan"
person2 = "Eric"
RNAddsres<-results(ddsCollapsed_normfactor,  contrast=c("person", person1, person2))
resdata <- as.data.frame(RNAddsres)
fullresdata <- merge(resdata,annotationmerge,by.x=0,by.y=3)
commonfullresdata <- merge(fullresdata,common_ids,by.x=1,by.y=1)

commonfullresdata <- commonfullresdata[commonfullresdata$baseMean > 12,] # Remove genes with very low expression

fullresdata <- fullresdata[fullresdata$baseMean > 30,] # Remove genes with very low expression
fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]

medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]
medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))

signif_medresdata <- medresdata[medresdata$padj < padjusted_cutoff,]
notsignif_medresdata <- medresdata[medresdata$padj >= padjusted_cutoff,]

ggplot() + 
  geom_violin(data=fullresdata,trim=TRUE,aes(x=chr, y=log2FoldChange)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_hline(yintercept=0) +
  ylab(paste0("Log2FC")) +
  xlab("Chromosome") +
  theme_classic(base_size = 24) +
  geom_jitter(data = signif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=1.5, alpha=0.85,color="red",fill="red") +
  geom_jitter(data = notsignif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=1.5, alpha=0.85,color="black",fill="black") +
  stat_summary(data=fullresdata,aes(x=chr, y=log2FoldChange),fun.y=median, geom="crossbar", size=0.25, color="orange") +
  theme(legend.title = element_blank())

####GRO-seq, alternative hypothesis adjustment ####
GRObodycountdat<- read.csv("DS_Normalization/counts/gro/featurecounts_gro_noproms_nomultis.sense.txt",
                           skip=1, sep="\t")


# Cleaning up column names so that they're more legible.
names <- colnames(GRObodycountdat)
names <- gsub(".*E","E",names)
colnames(GRObodycountdat) <- names

nomultis_col <- colnames(GRObodycountdat)[7:22]

GROmetadata=read.table("DS_Normalization/metadata/GROinfo.txt", sep="\t", header=TRUE)
GROmetadata$samplegroup <-paste(GROmetadata$person, "_", GROmetadata$libprep, sep="")
GROmetadata$nomultis <- nomultis_col

GRObodycountdat <- GRObodycountdat[GRObodycountdat$Geneid %in% masterannotationdf$GeneID,]
GRObodycountdat <- GRObodycountdat[!(duplicated(GRObodycountdat$Geneid)),]
rownames(GRObodycountdat) <- GRObodycountdat$Geneid

# Remove annotation info columns
GRObodycountdat <- GRObodycountdat[,-c(1:6)]
#run Deseq on GRO-seq with 21
ddsFull <- DESeqDataSetFromMatrix(countData = GRObodycountdat, colData = GROmetadata, design = ~ libprep + person) #use this for all samples
#ddsFull$type <- relevel(ddsFull$type, "WT")
ddsFull <- collapseReplicates( ddsFull,groupby = ddsFull$samplegroup,run = ddsFull$samplegroup )
ddsFull <- DESeq(ddsFull)

# Run with alternative hypothesis instead. While genes on other chrs aren't easily interpretable, anything significant on 21 is 
# potentially "dosage compensated
GROddsres<-results(ddsFull,  contrast=c("person", person1, person2),lfcThreshold = log2(1.5),altHypothesis = "lessAbs")
resdata <- as.data.frame(GROddsres)
fullresdata <- merge(resdata,annotationmerge,by.x=0,by.y=3)
fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]
medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]
medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))

signif_medresdata <- medresdata[medresdata$padj < padjust_cutoff,]
notsignif_medresdata <- medresdata[medresdata$padj >= padjust_cutoff,]

ggplot() + 
  geom_violin(data=fullresdata,trim=TRUE,aes(x=chr, y=log2FoldChange)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_hline(yintercept=0) +
  ylab(paste0("Log2FC")) +
  geom_hline(yintercept = log2(1.5), linetype="dotted", color="blue" ) +
  xlab("Chromosome") +
  theme_classic(base_size = 24) +
  geom_jitter(data = notsignif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=1.5, alpha=0.85,color="black",fill="black") +
  geom_jitter(data = signif_medresdata,aes(x=chr, y=log2FoldChange),shape=1, position=position_jitter(0.15,seed = 1),size=1.5, alpha=0.85,color="red",fill="red") +
  stat_summary(data=fullresdata,aes(x=chr, y=log2FoldChange),fun.y=median, geom="crossbar", size=0.25, color="orange") +
  theme(legend.title = element_blank())

