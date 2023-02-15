####Facet Grids ####
# Identifying quantiles of expression levels and their fold changes for each gene
library(plyr)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(tibble)


########RNA-SEQ#######
#pointers to original files Mary made
RNAindir<-"DS_Normalization/counts/rna"
GROseq_indir<-"DS_Normalization/counts/gro"
RNAcoveragedat<-"res_featureCounts_gene_idfull_143138.coverage.csv"

#Below is the bed that got used for RNA-seq. 
#This is not the bed file that got used for GRO-seq exactly becuase to be acuurate in the GROS-seq I had to remove any gene bodys/tss that were in the data more than once. 
ori_worldbed <- "DS_Normalization/annotation/hg38_refseq.bed"
#Below are the bed file I used for GRO-seq anaysis
#genes.bed and tss.bed came from /Users/allenma/humangtf.ipynb
#in breif i keep the longest isoform with uniq body start and stop, body(+1000, -500) start and stop, tss(-500, +500) start and stop
#I also removed gene <3000bp and genes whose body or tss <0 in coordinates
#below I will only keep genes that were anayized by both GR0-seq and RNA-seq
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
lesschrs <- c("chr20","chr21")

common_ids <- read.table("DS_Normalization/annotation/refseq_to_common_id.txt",sep="\t")
common_ids <- common_ids[common_ids$V1 %in% annotationmerge$GeneID,]
####Uncorrected Real RNA-seq data####
RNAbed <-ori_worldbed
filetable <- RNAmetadata
minichrs <-c("chr1", "chr2", "chr3","chr4", "chr5", "chr6", "chr7", "chr8", "chr9","chr10","chr11", "chr12", "chr13","chr14", "chr15", "chr16", "chr17", "chr18", "chr19","chr20","chr21","chr22")
masterannotationdf <- annotationmerge

masterannotationdf_only21 <- masterannotationdf %>% filter(chr=="chr21")
anuplodygenes <-masterannotationdf_only21[["name"]] #Only chr21 genes
baseploidy <- 2
alt_ploidy <-3

#this is loading the metadata
RNAmetadata=read.table("DS_Normalization/metadata/RNAinfo.txt", sep="\t", header=TRUE)
RNAmetadata$samplegroup <-paste(RNAmetadata$Person, "_", RNAmetadata$biological_rep, sep="")

#this is where you would remove samples if they are bad
RNAcountdat <- read.csv(RNAcoveragedat,sep="\t",skip=1)

# Maximal isoform filter
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
ddsFull <- DESeq(ddsFull)
RNAddsres<-results(ddsFull,contrast=c("Person","Ethan","Eric"))
resdata <- as.data.frame(RNAddsres)
fullresdata <- merge(resdata,annotationmerge,by.x=0,by.y=3)
fullresdata <- merge(fullresdata,common_ids,by.x=1,by.y=1)

# Depending on which genes you pick, RNA-seq MFC ranges from ~1.34-1.41.
medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
#  Quintile of FC violin plots
fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]
fullresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]

# Separate into quintiles
# FIXME: change to iterable
quantile(fullresdata$baseMean,probs=seq(0,1,.2))
quantile1 <- fullresdata[fullresdata$baseMean<=(as.list(quantile(fullresdata$baseMean,probs=seq(0,1,.2)))[[2]][1]),]
quantile1$quantile <- 1
medians1 <- ddply(quantile1, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians1$quantile <- 1
quantile2 <- fullresdata[fullresdata$baseMean <= (as.list(quantile(fullresdata$baseMean,probs=seq(0,1,.2)))[[3]][1]) & fullresdata$baseMean> (as.list(quantile(fullresdata$baseMean,probs=seq(0,1,.2)))[[2]][1]),]
quantile2$quantile <- 2
medians2 <- ddply(quantile2, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians2$quantile <- 2
quantile3 <- fullresdata[fullresdata$baseMean<= (as.list(quantile(fullresdata$baseMean,probs=seq(0,1,.2)))[[4]][1]) & fullresdata$baseMean> (as.list(quantile(fullresdata$baseMean,probs=seq(0,1,.2)))[[3]][1]),]
quantile3$quantile <- 3
medians3 <- ddply(quantile3, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians3$quantile <- 3
quantile4 <- fullresdata[fullresdata$baseMean<= (as.list(quantile(fullresdata$baseMean,probs=seq(0,1,.2)))[[5]][1]) & fullresdata$baseMean> (as.list(quantile(fullresdata$baseMean,probs=seq(0,1,.2)))[[4]][1]),]
quantile4$quantile <- 4
medians4 <- ddply(quantile4, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians4$quantile <- 4
quantile5 <- fullresdata[fullresdata$baseMean<= (as.list(quantile(fullresdata$baseMean,probs=seq(0,1,.2)))[[6]][1]) & fullresdata$baseMean> (as.list(quantile(fullresdata$baseMean,probs=seq(0,1,.2)))[[5]][1]),]
quantile5$quantile <- 5
medians5 <- ddply(quantile5, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians5$quantile <- 5

# And a total to compare to:
quantileall <- fullresdata
quantileall$quantile <- 6
medians6 <- ddply(quantileall, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians6$quantile <- 6

newfullres <- rbind(quantile1,quantile2,quantile3,quantile4,quantile5,quantileall)
newmedians <- rbind(medians1,medians2,medians3,medians4,medians5,medians6)

# Facet Grid of Violins
newfullres <- newfullres[newfullres$chr %in% lesschrs,]

ggplot() + 
  geom_violin(data=newfullres,aes(x=chr, y=log2FoldChange)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_hline(yintercept=0) +
  ylab(paste0("Log2FC")) +
  geom_hline(yintercept = log2(1.5), linetype="dotted", color="blue" ) +
  xlab("Chromosome") +
  theme_classic(base_size = 24) +
  theme(legend.title = element_blank()) + 
  facet_wrap(~quantile, ncol = 6) +
  geom_hline(data = newmedians, aes(yintercept = log2(med)),
             colour = "orange")

### Corrected RNA-seq ###
# Correcting things

####REMOVING REPEATS####
# We noticed that a lot of samples had specific genes that always appeared near 1.0, especially rRNA. Could it be
# that repeats are causing multi-mapping issues for these genes? To test, we removed reads over RNA repeats in genome

RNAbed <-ori_worldbed
filetable <- RNAmetadata

masterannotationdf <- annotationmerge
masterannotationdf_only21 <- masterannotationdf %>% filter(chr=="chr21")
head(masterannotationdf_only21)
anuplodygenes <-masterannotationdf_only21[["name"]] #Only chr21 genes
baseploidy <- 2
alt_ploidy <-3

#this is loading the metadata
RNAmetadata=read.table("DS_Normalization/metadata/RNAinfo.txt", sep="\t", header=TRUE)
RNAcountdat <- read.csv("DS_Normalization/counts/rna/featurecounts_rna_rmsked.sense.txt",
                        sep="\t", skip=1)

# filter to only canonical chrs
RNAcountdat <- RNAcountdat[RNAcountdat$Chr %in% minichrs,]
names <- colnames(RNAcountdat)
names <- gsub(".*removed..","",names)
colnames(RNAcountdat) <- names


# Maximal isoform filter, if you wanna
RNAcountdat$sum <- rowSums(RNAcountdat[,7:ncol(RNAcountdat)])/RNAcountdat$Length
RNAcountdat <- merge(RNAcountdat,common_ids,by.x=1,by.y=1)
RNAcountdat <- (RNAcountdat[order(RNAcountdat$V2,-RNAcountdat$sum),])
RNAcountdat <- (RNAcountdat[!(duplicated(RNAcountdat$V2)),])
RNAcountdat <- subset(RNAcountdat, select=-c(V2,sum))

# Depending on which genes you use, you'll get different FC results. Can either keep all or only those in both RNA/GRO files, for example

RNAcountdat<- RNAcountdat %>%
  filter(Geneid %in% annotationmerge$GeneID) %>%
  arrange(Geneid) %>%
  column_to_rownames('Geneid')

RNAcountdat <- RNAcountdat[,-c(1:5)]
#run Deseq on RNA-seq with 21
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

fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]
fullresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]

# FIXME: Change to iterable
quantile(fullresdata$baseMean,probs=seq(0,1,.2))
quantile1 <- fullresdata[fullresdata$baseMean<=(as.list(quantile(fullresdata$baseMean,probs=seq(0,1,.2)))[[2]][1]),]
quantile1$quantile <- 1
medians1 <- ddply(quantile1, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians1$quantile <- 1
quantile2 <- fullresdata[fullresdata$baseMean <= (as.list(quantile(fullresdata$baseMean,probs=seq(0,1,.2)))[[3]][1]) & fullresdata$baseMean> (as.list(quantile(fullresdata$baseMean,probs=seq(0,1,.2)))[[2]][1]),]
quantile2$quantile <- 2
medians2 <- ddply(quantile2, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians2$quantile <- 2

quantile3 <- fullresdata[fullresdata$baseMean<= (as.list(quantile(fullresdata$baseMean,probs=seq(0,1,.2)))[[4]][1]) & fullresdata$baseMean> (as.list(quantile(fullresdata$baseMean,probs=seq(0,1,.2)))[[3]][1]),]
quantile3$quantile <- 3
medians3 <- ddply(quantile3, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians3$quantile <- 3

quantile4 <- fullresdata[fullresdata$baseMean<= (as.list(quantile(fullresdata$baseMean,probs=seq(0,1,.2)))[[5]][1]) & fullresdata$baseMean> (as.list(quantile(fullresdata$baseMean,probs=seq(0,1,.2)))[[4]][1]),]
quantile4$quantile <- 4
medians4 <- ddply(quantile4, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians4$quantile <- 4

quantile5 <- fullresdata[fullresdata$baseMean<= (as.list(quantile(fullresdata$baseMean,probs=seq(0,1,.2)))[[6]][1]) & fullresdata$baseMean> (as.list(quantile(fullresdata$baseMean,probs=seq(0,1,.2)))[[5]][1]),]
quantile5$quantile <- 5
medians5 <- ddply(quantile5, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians5$quantile <- 5

quantileall <- fullresdata
quantileall$quantile <- 6
medians6 <- ddply(quantileall, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians6$quantile <- 6

newfullres <- rbind(quantile1,quantile2,quantile3,quantile4,quantile5,quantileall)
newmedians <- rbind(medians1,medians2,medians3,medians4,medians5,medians6)

# Facet Grid of Violins
newfullres <- newfullres[newfullres$chr %in% lesschrs,]

ggplot() + 
  geom_violin(data=newfullres,aes(x=chr, y=log2FoldChange)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_hline(yintercept=0) +
  ylab(paste0("Log2FC")) +
  geom_hline(yintercept = log2(1.5), linetype="dotted", color="blue" ) +
  xlab("Chromosome") +
  theme_classic(base_size = 24) +
  theme(legend.title = element_blank()) + 
  facet_wrap(~quantile, ncol = 6) +
  geom_hline(data = newmedians, aes(yintercept = log2(med)),
             colour = "orange")

### Uncorrected GRO ###
#### GRO-SEQ####
#### UNCORRECTED REAL DATA GRO-SEQ RUN####
#traditional Deseq2 on on GRO count matrix's

GROseq_indir<-"DS_Normalization/counts/gro/"
GROmetadata=read.table("DS_Normalization/metadata/GROinfo.txt", sep="\t", header=TRUE)
GROmetadata$samplegroup <-paste(GROmetadata$person, "_", GROmetadata$libprep, sep="")

library(dplyr)

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

fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]
fullresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]

quantile(fullresdata$baseMean,probs=seq(0,1,.2))

quantile1 <- fullresdata[fullresdata$baseMean<=(as.list(quantile(fullresdata$baseMean,probs=seq(0,1,.2)))[[2]][1]),]
quantile1$quantile <- 1
medians1 <- ddply(quantile1, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians1$quantile <- 1
quantile2 <- fullresdata[fullresdata$baseMean <= (as.list(quantile(fullresdata$baseMean,probs=seq(0,1,.2)))[[3]][1]) & fullresdata$baseMean> (as.list(quantile(fullresdata$baseMean,probs=seq(0,1,.2)))[[2]][1]),]
quantile2$quantile <- 2
medians2 <- ddply(quantile2, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians2$quantile <- 2

quantile3 <- fullresdata[fullresdata$baseMean<= (as.list(quantile(fullresdata$baseMean,probs=seq(0,1,.2)))[[4]][1]) & fullresdata$baseMean> (as.list(quantile(fullresdata$baseMean,probs=seq(0,1,.2)))[[3]][1]),]
quantile3$quantile <- 3
medians3 <- ddply(quantile3, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians3$quantile <- 3

quantile4 <- fullresdata[fullresdata$baseMean<= (as.list(quantile(fullresdata$baseMean,probs=seq(0,1,.2)))[[5]][1]) & fullresdata$baseMean> (as.list(quantile(fullresdata$baseMean,probs=seq(0,1,.2)))[[4]][1]),]
quantile4$quantile <- 4
medians4 <- ddply(quantile4, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians4$quantile <- 4

quantile5 <- fullresdata[fullresdata$baseMean<= (as.list(quantile(fullresdata$baseMean,probs=seq(0,1,.2)))[[6]][1]) & fullresdata$baseMean> (as.list(quantile(fullresdata$baseMean,probs=seq(0,1,.2)))[[5]][1]),]
quantile5$quantile <- 5
medians5 <- ddply(quantile5, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians5$quantile <- 5

quantileall <- fullresdata
quantileall$quantile <- 6
medians6 <- ddply(quantileall, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians6$quantile <- 6

newfullres <- rbind(quantile1,quantile2,quantile3,quantile4,quantile5,quantileall)
newmedians <- rbind(medians1,medians2,medians3,medians4,medians5,medians6)

# Facet Grid of Violins
newfullres <- newfullres[newfullres$chr %in% lesschrs,]

ggplot() + 
  geom_violin(data=newfullres,aes(x=chr, y=log2FoldChange)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_hline(yintercept=0) +
  ylab(paste0("Log2FC")) +
  geom_hline(yintercept = log2(1.5), linetype="dotted", color="blue" ) +
  xlab("Chromosome") +
  theme_classic(base_size = 24) +
  theme(legend.title = element_blank()) + 
  facet_wrap(~quantile, ncol = 6) +
  geom_hline(data = newmedians, aes(yintercept = log2(med)),
             colour = "orange")

### Corrected GRO ###
GRObodycountdat<- read.csv("DS_Normalization/counts/gro/featurecounts_gro_noproms_nomultis.sense.txt",
                           skip=1, sep="\t")

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

#run Deseq on RNA-seq with 21
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

fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]
fullresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]
quantile(fullresdata$baseMean,probs=seq(0,1,.2))

quantile1 <- fullresdata[fullresdata$baseMean<=(as.list(quantile(fullresdata$baseMean,probs=seq(0,1,.2)))[[2]][1]),]
quantile1$quantile <- 1
medians1 <- ddply(quantile1, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians1$quantile <- 1
quantile2 <- fullresdata[fullresdata$baseMean <= (as.list(quantile(fullresdata$baseMean,probs=seq(0,1,.2)))[[3]][1]) & fullresdata$baseMean> (as.list(quantile(fullresdata$baseMean,probs=seq(0,1,.2)))[[2]][1]),]
quantile2$quantile <- 2
medians2 <- ddply(quantile2, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians2$quantile <- 2

quantile3 <- fullresdata[fullresdata$baseMean<= (as.list(quantile(fullresdata$baseMean,probs=seq(0,1,.2)))[[4]][1]) & fullresdata$baseMean> (as.list(quantile(fullresdata$baseMean,probs=seq(0,1,.2)))[[3]][1]),]
quantile3$quantile <- 3
medians3 <- ddply(quantile3, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians3$quantile <- 3

quantile4 <- fullresdata[fullresdata$baseMean<= (as.list(quantile(fullresdata$baseMean,probs=seq(0,1,.2)))[[5]][1]) & fullresdata$baseMean> (as.list(quantile(fullresdata$baseMean,probs=seq(0,1,.2)))[[4]][1]),]
quantile4$quantile <- 4
medians4 <- ddply(quantile4, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians4$quantile <- 4

quantile5 <- fullresdata[fullresdata$baseMean<= (as.list(quantile(fullresdata$baseMean,probs=seq(0,1,.2)))[[6]][1]) & fullresdata$baseMean> (as.list(quantile(fullresdata$baseMean,probs=seq(0,1,.2)))[[5]][1]),]
quantile5$quantile <- 5
medians5 <- ddply(quantile5, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians5$quantile <- 5

quantileall <- fullresdata
quantileall$quantile <- 6
medians6 <- ddply(quantileall, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians6$quantile <- 6

newfullres <- rbind(quantile1,quantile2,quantile3,quantile4,quantile5,quantileall)
newmedians <- rbind(medians1,medians2,medians3,medians4,medians5,medians6)

# Facet Grid of Violins
newfullres <- newfullres[newfullres$chr %in% lesschrs,]

ggplot() + 
  geom_violin(data=newfullres,aes(x=chr, y=log2FoldChange)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_hline(yintercept=0) +
  ylab(paste0("Log2FC")) +
  geom_hline(yintercept = log2(1.5), linetype="dotted", color="blue" ) +
  xlab("Chromosome") +
  ylim(-4,4) +
  theme_classic(base_size = 24) +
  theme(legend.title = element_blank()) + 
  facet_wrap(~quantile, ncol = 6) +
  geom_hline(data = newmedians, aes(yintercept = log2(med)),
             colour = "orange")