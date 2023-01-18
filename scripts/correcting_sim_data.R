### Correcting Sim Data ###
library(tidyverse)

RNAcountdat <- read.csv("/path/to/simcounts",
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

# We have to make fake metadata here for the design matrix. Note that sims incorporate NO batch effects
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
anuplodygenes <-masterannotationdf_only21[["name"]] #Only chr21 genes
baseploidy <- 2
alt_ploidy <-3

#run Deseq on RNA-seq with 21
ddsFull <- DESeqDataSetFromMatrix(countData = RNAcountdat, colData = RNAmetadata, design = ~ Person) #use this for all samples
dds <- ddsFull

ddsnocorr <- DESeq(ddsFull)

person1 = "Eddie" #T21
person2 = "Elvis" #D21
RNAddsres<-results(ddsnocorr,  contrast=c("Person", person1, person2))
resdata <- as.data.frame(RNAddsres)
fullresdata <- merge(resdata,annotationnew,by.x=0,by.y=4)

common_fullresdata_rna <- merge(fullresdata,common_ids,by.x=1,by.y=1)

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

person1 = "Eddie"
person2 = "Elvis"
RNAddsres<-results(ddsCollapsed_normfactor,  contrast=c("Person", person1, person2))
resdata <- as.data.frame(RNAddsres)
fullresdata <- merge(resdata,annotationnew,by.x=0,by.y=4)
fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]
fullresdata <- fullresdata[fullresdata$baseMean > 30,] # Remove genes with very low expression

medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]
medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
signif_medresdata <- medresdata[medresdata$padj<.01,]
notsignif_medresdata <- medresdata[medresdata$padj>=.01,]

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
