library("reshape2")
library("plyr")
library("ggplot2")
library("miscTools")
library("stringr")
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
lesschrs <- c("chr20","chr21","chr22")

###Initial DESeq2 Run with NO adjustments####
###UNCORRECTED SIM COMPARISONS####
# I'm purposely picking a high dispersion/noise file because we can see just how well our methods work to correct it
# Note that since this is simulated data, my gut feeling is that the corrections won't be quite as good as with real data
RNAcountdat <- read.csv("/Users/sahu0957/trisomy_normalization/RNAGRO01132021/data/RNAseqfiles/NEWSCALEreplicatefinalsims_a0.982_p1_r2_s1_n3.csv", 
                        sep=",", row.names=1)
head(RNAcountdat)
#RNAcountdat$Gene_id <- rownames(RNAcountdat)
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
RNAdds <- DESeqDataSetFromMatrix(countData = RNAcountdat, colData = RNAmetadata, design = ~Person)
RNAdds <-DESeq(RNAdds, betaPrior = FALSE) # Earlier DESeq2 versions include a betaPrior for the intercept term... wild!

person1 = "Eddie"
person2 = "Elvis"
RNAddsres<-results(RNAdds,  contrast=c("Person", person1, person2))
resdata <- as.data.frame(RNAddsres)
fullresdata <- merge(resdata,annotationnew,by.x=0,by.y=4)
fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]
medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]
medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians

ggplot(fullresdata, aes(x=chr, y=log2FoldChange)) + 
  geom_violin(trim=TRUE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_hline(yintercept=log2(1.5),linetype="dashed",color="blue") +
  #  geom_hline(yintercept=log2(0.66667),linetype="dashed",color="red") +
  geom_hline(yintercept=0) +
  ylab(paste0("log2FC(ST21/SD21)")) +
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


###looping runs with NO correction#### 
# We can see whether there's a pattern associated with more/less dispersion by looping DESeq2 runs through all of our simulated files
# each count file is in the same format, where the sim info is in the name of the file and the counts in the file itself
# so we can parse out the parameters of the simulation first and then run DESeq2
# a = asymptotic dispersion
# p = extra-poisson noise
# n = number of replicates in a single experiment
# r = number of experiments
# s = scalar for depth
countfiles <- list.files(path = "/Users/sahu0957/trisomy_normalization/RNAGRO01132021/data/RNAseqfiles/",pattern = "NEWSCALE")
countfiles
median_colnames <- c("chr","med","asympt_Disp","extraPoisson","rep","scale","n")
median_scaleddata <- read.table(text = "",
                                col.names = median_colnames,
                                stringsAsFactors = FALSE)

disp_info <- str_split(k, "_")

disp_info

for (k in countfiles){
  disp_info <- str_split(k, "_")
  asympt_Dispersion = substring(disp_info[[1]][2], 2)
  extraPoisson_noise= substring(disp_info[[1]][3], 2)
  rep_number= substring(disp_info[[1]][4], 2)
  scale_number=substring(disp_info[[1]][5], 2)
  n_number <- substring(disp_info[[1]][6], 2, nchar(disp_info[[1]][6])-4)
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
  median_scaleddata <- rbind(median_scaleddata,medians)
}

median_scaleddata

median_scaleddata$asympt_Disp <- as.numeric(as.character(median_scaleddata$asympt_Disp))
median_scaleddata$extraPoisson <- as.numeric(as.character(median_scaleddata$extraPoisson))
median_scaleddata$rep <- as.numeric(as.character(median_scaleddata$rep))
median_scaleddata$scale <- as.numeric(as.character(median_scaleddata$scale))
median_scaleddata$n <- as.numeric(as.character(median_scaleddata$n))

# To visualize, we first "marginalize" on one of the parameters so we're comparing like with like on all samples
test_mediandata <- median_scaleddata[median_scaleddata$asympt_Disp==.01,]
test_mediandata
ggplot(test_mediandata, aes(x=extraPoisson,y=med,color=chr)) +
  geom_point() +
  geom_smooth(method="loess",span=1) +
  labs(colour = "Chromosome",title = "No LFC Shrinkage") +
  geom_hline(yintercept = 1.5,linetype="dashed") +
  geom_hline(yintercept = 1.0, linetype="dashed") +
  ylab("Median Fold Change") +
  xlab("Extra-Poisson Noise") +
  #  annotate("segment", x = 0.001, xend = 0.001, y = 1.2, yend = 1.28, colour = "red", size=0.5, alpha=1, arrow=arrow()) +
  #  annotate("segment", x = 0.16, xend = 0.16, y = 1.12, yend = 1.2, colour = "red", size=0.5, alpha=1, arrow=arrow()) +
  #xlab(expression(paste(alpha))) +
  theme_classic(base_size=18)

ggsave("rnaseq_deseq2_realsims_varypoisson.svg",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "svg")
ggsave("rnaseq_deseq2_realsims_varypoinsson.png",path = "/scratch/Users/sahu0957/ds_normalization/figs",device = "png")

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
RNAbed <-ori_worldbed
filetable <- RNAmetadata
minichrs <-c("chr1", "chr2", "chr3","chr4", "chr5", "chr6", "chr7", "chr8", "chr9","chr10","chr11", "chr12", "chr13","chr14", "chr15", "chr16", "chr17", "chr18", "chr19","chr20","chr21","chr22")

masterannotationdf <- annotationmerge
masterannotationdf_only21 <- masterannotationdf %>% filter(chr=="chr21")
anuplodygenes <-masterannotationdf_only21[["name"]] #Only chr21 genes
baseploidy <- 2
alt_ploidy <-3

#this is loading the metadata
RNAmetadata=read.table("/Shares/down/mixed/RNAGRO01132021/scripts/Deseq2_analysis/RNAinfo.txt", sep="\t", header=TRUE)

#this is where you would remove samples if they are bad
RNAcountdat <- read.csv(paste0(RNAindir, RNAcoveragedat),
                        sep=",", row.names=1)

RNAcountdat$Gene_id <- rownames(RNAcountdat)
RNAcountdat<- RNAcountdat %>%
  filter(Gene_id %in% masterannotationdf$name) %>%
  arrange(Gene_id) %>%
  column_to_rownames('Gene_id')


#run Deseq on RNA-seq with 21
ddsFull <- DESeqDataSetFromMatrix(countData = RNAcountdat, colData = RNAmetadata, design = ~ biological_rep + Person) #use this for all samples
# We don't need to collapse replicates with RNA-seq but we will with GRO
#ddsFull$type <- relevel(ddsFull$type, "WT")
#dds <- collapseReplicates( ddsFull,groupby = ddsFull$samplegroup,run = ddsFull$samplegroup )
dds <- ddsFull

# Can run DESeq2 unaltered here to compare results
#ddsFull <- DESeq(ddsFull)

#person1 = "Ethan"
#person2 = "Eric"
#RNAddsres<-results(ddsFull,  contrast=c("Person", person1, person2))
#resdata <- as.data.frame(RNAddsres)
#fullresdata <- merge(resdata,annotationnew,by.x=0,by.y=4)
#fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]
#medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]

#medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
#medians

controlgenes <-ifelse(rownames(dds) %in% anuplodygenes, FALSE, TRUE) # Only normalize on non-aneuploidy genes

samplepinfo<-as.data.frame(colData(dds))

ploidy_fortrisomy = ifelse(rownames(dds) %in% anuplodygenes, alt_ploidy, baseploidy)
ploidy_typical <- rep(baseploidy, nrow(dds))
ploidy_fortrisomy
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
fullresdata <- fullresdata[fullresdata$baseMean > 12,] # Remove genes with very low expression
medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]

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

# filter to only canonical chrs
head(RNAcountdat)
RNAcountdat <- RNAcountdat[RNAcountdat$Chr %in% minichrs,]

# colnames are different due to my featurecounts script. fixing that here by removing the path. Kinda hacky, go remove in script!
names <- colnames(RNAcountdat)
names <- gsub(".*removed..","",names)
colnames(RNAcountdat) <- names

# We've got 1300 duplicate entries... How was this filtered out before?
RNAcountdat <- RNAcountdat[!(duplicated(RNAcountdat$Geneid)),]
# Filter to genes common in both RNA and GRO seq
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

ploidy_fortrisomy = ifelse(rownames(dds) %in% anuplodygenes, alt_ploidy, baseploidy)
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
fullresdata <- fullresdata[fullresdata$baseMean > 30,] # Remove genes with very low expression

medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]
dosagelowgenes <- medresdata[medresdata$log2FoldChange < 0.13 & medresdata$chr=="chr21",]
nrow(dosagelowgenes)
medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians


###REPEAT-REMOVED DATA DONE!####
#### START CORRECTING SIM DATA####

# Again, deliberately picking a real noisy file to try this out on.
# These sims have NOT had any repeats removed. I'll need to remake these sims to reflect that later
RNAcountdat <- read.csv("/Users/sahu0957/trisomy_normalization/RNAGRO01132021/data/RNAseqfiles/NEWSCALEreplicatefinalsims_a1.2_p10_r3_s1_n3.csv",
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
fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]
medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]

medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians


controlgenes <-ifelse(rownames(dds) %in% anuplodygenes, FALSE, TRUE) # Only normalize on non-aneuploidy genes

samplepinfo<-as.data.frame(colData(dds))
ploidy_fortrisomy = ifelse(rownames(dds) %in% anuplodygenes, alt_ploidy, baseploidy)
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
lesschrs <- c("chr20","chr21")
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

lesschrs <- c("chr20","chr21")
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
#traditional Deseq2 on on GRO count matrix's

GRObodyddsFull <- DESeqDataSetFromMatrix(countData = GRObodycountdat, colData = GROmetadata, design = ~libprep+ person)
GRObodydds <- collapseReplicates( GRObodyddsFull,groupby = GRObodyddsFull$samplegroup,run = GRObodyddsFull$samplegroup)
GRObodydds <-DESeq(GRObodydds)

GROtssddsFull <- DESeqDataSetFromMatrix(countData = GROtsscountdat, colData = GROmetadata, design = ~libprep+person)
GROtssdds <- collapseReplicates( GROtssddsFull,groupby = GROtssddsFull$samplegroup,run = GROtssddsFull$samplegroup )
GROtssdds <-DESeq(GROtssdds)

GROgeneddsFull <- DESeqDataSetFromMatrix(countData = GROgenecountdat, colData = GROmetadata, design = ~libprep+person)
GROgenedds <- collapseReplicates( GROgeneddsFull,groupby = GROgeneddsFull$samplegroup,run = GROgeneddsFull$samplegroup )
GROgenedds <-DESeq(GROgenedds)

removegenesfromnorm_GRObody<- function(dropchr){
  GRObodyddsFull <- DESeqDataSetFromMatrix(countData = GRObodycountdat, colData = GROmetadata, design = ~libprep+ person)
  GRObodydds <- collapseReplicates( GRObodyddsFull,groupby = GRObodyddsFull$samplegroup,run = GRObodyddsFull$samplegroup)
  GRObodydds_nochrincontroldds <-estimateSizeFactors(GRObodydds, controlGenes=dropchr) 
  GRObodydds_nochrincontroldds <- estimateDispersionsGeneEst(GRObodydds_nochrincontroldds) 
  GRObodydds_nochrincontroldds<-estimateDispersionsFit(GRObodydds_nochrincontroldds) 
  GRObodydds_nochrincontroldds <- estimateDispersionsMAP(GRObodydds_nochrincontroldds) 
  GRObodydds_nochrincontroldds <- nbinomWaldTest(GRObodydds_nochrincontroldds)
  sf <-as.data.frame(sizeFactors(GRObodydds_nochrincontroldds))} 

person1 = "Ethan"
person2 = "Eddie"

RNAseq_alteredgenes <- function(person1, person2){
  RNAdds <- DESeqDataSetFromMatrix(countData = GROgenecountdat, colData = GROmetadata, design = ~ libprep+ person)
  RNAdds <-DESeq(RNAdds)
  RNAddsres<-results(RNAdds,  contrast=c("Person", person1, person2))
  RNAdds_no21dds <- DESeqDataSetFromMatrix(countData = GROgenecountdatdropchr21, colData = GROmetadata, design = ~ libprep+ person)
  RNAdds_no21dds <-DESeq(RNAdds_no21dds)
  RNAdds_no21ddsres<-results(RNAdds_no21dds,contrast=c("Person", person1, person2))
  RNAdds_no21ddsshuffle <- DESeqDataSetFromMatrix(countData = GROgenecountdatdropchr21shuffle, 
                                                  colData = GROmetadata, design = ~ libprep+ person)
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

diffexpdf <- rbind(RNAseq_alteredgenes("Eli", "Eric"), 
                   RNAseq_alteredgenes("Elizabeth", "Eric"), 
                   RNAseq_alteredgenes("Eli", "Ethan"),
                   RNAseq_alteredgenes("Elizabeth", "Ethan"), 
                   RNAseq_alteredgenes("Elizabeth", "Eli"), 
                   RNAseq_alteredgenes("Eric", "Ethan"), 
                   RNAseq_alteredgenes("Eddie", "Ethan")) 

diffexpdf$compairison <-rownames(diffexpdf)
diffexpdf$diff <- diffexpdf$RNAseq - diffexpdf$RNAseq_no_21
diffexpdf_long =  gather(diffexpdf, sample, ngenes, -compairison, -diff) %>%separate(compairison, into=c("comparison", "direction"), sep="__")
diffexpdf_long$comparison <- factor(diffexpdf_long$comparison, 
                                    levels=c("Eli_vs_Eric", "Elizabeth_vs_Eric", "Elizabeth_vs_Eli", "Eli_vs_Ethan", "Elizabeth_vs_Ethan", "Eric_vs_Ethan","Eddie_vs_Ethan"))

ggplot(diffexpdf_long, aes(x=comparison, y=ngenes, color=sample)) + geom_point()+ theme(axis.text.x = element_text(angle = 90))
ggplot(diffexpdf_long, aes(x=comparison, y=diff, color=direction)) + geom_point()+ theme(axis.text.x = element_text(angle = 90))

person1 = "Eddie"
person2 = "Eli"

RNAtestcountdat <- GRObodycountdat[,colnames(GRObodycountdat) %like% person1 | colnames(GRObodycountdat) %like% person2] 
RNAtestmetadata <- GROmetadata[GROmetadata$file %like% person1 | GROmetadata$file %like% person2,]
RNAtestcountdatdropchr21 <- GRObodycountdatdropchr21[,colnames(GRObodycountdatdropchr21) %like% person1 | colnames(GRObodycountdatdropchr21) %like% person2]
RNAtestcountdatdropchr21shuffle <- GRObodycountdatdropchr21shuffle[,colnames(GRObodycountdatdropchr21shuffle) %like% person1 | colnames(GRObodycountdatdropchr21shuffle) %like% person2]

RNAddstest_no21dds <- DESeqDataSetFromMatrix(countData = RNAtestcountdatdropchr21, colData = RNAtestmetadata, design = ~ libprep+ person)
RNAddstest_no21dds <- collapseReplicates( RNAddstest_no21dds,groupby = RNAddstest_no21dds$samplegroup,run = RNAddstest_no21dds$samplegroup)
RNAddstest_no21dds <-DESeq(RNAddstest_no21dds)

RNAddstest <- DESeqDataSetFromMatrix(countData = RNAtestcountdat, colData = RNAtestmetadata, design = ~ libprep+ person)
RNAddstest <- collapseReplicates( RNAddstest,groupby = RNAddstest$samplegroup,run = RNAddstest$samplegroup)
sizeFactors(RNAddstest) <- sizeFactors(RNAddstest_no21dds)
RNAddstest <-DESeq(RNAddstest)

RNAddstestres<-results(RNAddstest,  contrast=c("person", person1, person2))

RNAddstest_no21ddsres<-results(RNAddstest_no21dds,contrast=c("person", person1, person2))

RNAddstest_no21ddsshuffle <- DESeqDataSetFromMatrix(countData = RNAtestcountdatdropchr21shuffle, colData = RNAtestmetadata, design = ~ libprep+ person)
RNAddstest_no21ddsshuffle <- collapseReplicates( RNAddstest_no21ddsshuffle,groupby = RNAddstest_no21ddsshuffle$samplegroup,
                                                 run = RNAddstest_no21ddsshuffle$samplegroup)

RNAddstest_no21ddsshuffle <-DESeq(RNAddstest_no21dds)
RNAddstest_no21ddsshuffleres<-results(RNAddstest_no21ddsshuffle,contrast=c("person", person1, person2))

RNAddsres<-results(RNAddstest,  contrast=c("person", person1, person2))
RNAddsresOrdered <- RNAddsres[order(RNAddsres$pvalue),]
RNAddsresSig <- subset(RNAddsresOrdered, padj < 0.01)
RNAddsresSig_no21 <- as.data.frame(RNAddsresSig) %>%   
  rownames_to_column('gene') %>%
  filter(gene %in% annotationmergeno21$name) %>%
  arrange(gene) %>%
  column_to_rownames('gene')

RNAdds_no21ddsres<-results(RNAddstest_no21dds,contrast=c("person", person1, person2))


RNAdds_no21ddsresOrdered <- RNAdds_no21ddsres[order(RNAdds_no21ddsres$pvalue),]
RNAdds_no21ddsresSig <- subset(RNAdds_no21ddsresOrdered, padj < 0.01)
RNAdds_no21ddsresSig <- as.data.frame(RNAdds_no21ddsresSig)
RNAdds_no21ddsresSig_up <- RNAdds_no21ddsresSig %>% filter(log2FoldChange>0)
RNAdds_no21ddsresSig_down <- RNAdds_no21ddsresSig %>% filter(log2FoldChange<=0)
RNAddsresSig_no21_up <- RNAddsresSig_no21 %>% filter(log2FoldChange>0)
RNAddsresSig_no21_down <- RNAddsresSig_no21 %>% filter(log2FoldChange<=0)

resdata <- as.data.frame(RNAddsres)
resdata <- na.omit(resdata)
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

rownames(annotationnew) <- annotationnew$name

lesschrs <- c("chr19","chr20","chr21","chr22")

fullresdata <- merge(resdata,annotationnew,by.x=0,by.y=4)
fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]
fullresdata <- fullresdata[fullresdata$baseMean > 100,]

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

nrow(fullresdata[fullresdata$log2FoldChange <= 0.25 & fullresdata$chr == "chr21",])

### With Real Data:
RNAbed <-ori_worldbed
filetable <- RNAmetadata
RNAcoveragedat
minichrs <-c("chr1", "chr2", "chr3","chr4", "chr5", "chr6", "chr7", "chr8", "chr9","chr10","chr11", "chr12", "chr13","chr14", "chr15", "chr16", "chr17", "chr18", "chr19","chr20","chr21","chr22")

masterannotationdf <- annotationmerge
masterannotationdf_only21 <- masterannotationdf %>% filter(chr=="chr21")
head(masterannotationdf_only21)
anuplodygenes <-masterannotationdf_only21[["name"]] #Only chr21 genes
baseploidy <- 2
alt_ploidy <-3

#this is loading the metadata
RNAmetadata=read.table("/Shares/down/mixed/RNAGRO01132021/scripts/Deseq2_analysis/RNAinfo.txt", sep="\t", header=TRUE)
RNAmetadata

#this is where you would remove samples if they are bad
RNAcountdat <- read.csv(paste0(RNAindir, RNAcoveragedat),
                        sep=",", row.names=1)

RNAcountdat$Gene_id <- rownames(RNAcountdat)
head(RNAcountdat)
RNAcountdat<- RNAcountdat %>%
  filter(Gene_id %in% masterannotationdf$name) %>%
  arrange(Gene_id) %>%
  column_to_rownames('Gene_id')
dim(RNAcountdat)


#run Deseq on RNA-seq with 21
ddsFull <- DESeqDataSetFromMatrix(countData = RNAcountdat, colData = RNAmetadata, design = ~ biological_rep + Person) #use this for all samples
#ddsFull$type <- relevel(ddsFull$type, "WT")
#dds <- collapseReplicates( ddsFull,groupby = ddsFull$samplegroup,run = ddsFull$samplegroup )
dds <- ddsFull

#ddsFull <- DESeq(ddsFull)

#person1 = "Ethan"
#person2 = "Eric"
#RNAddsres<-results(ddsFull,  contrast=c("Person", person1, person2))
#resdata <- as.data.frame(RNAddsres)
#fullresdata <- merge(resdata,annotationnew,by.x=0,by.y=4)
#fullresdata <- fullresdata[fullresdata$chr %in% lesschrs,]
#medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]

#medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
#medians


controlgenes <-ifelse(rownames(dds) %in% anuplodygenes, FALSE, TRUE) # Only normalize on non-aneuploidy genes

samplepinfo<-as.data.frame(colData(dds))

ploidy_fortrisomy = ifelse(rownames(dds) %in% anuplodygenes, alt_ploidy, baseploidy)
ploidy_typical <- rep(baseploidy, nrow(dds))
ploidy_fortrisomy
normFactors <- matrix(do.call(cbind, mget(paste0(dds$ploidy))),ncol=ncol(dds),nrow=nrow(dds),dimnames=list(1:nrow(dds),1:ncol(dds)))
normFactors <- normFactors/baseploidy
normFactors
head(RNAcountdat)
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
fullresdata <- fullresdata[fullresdata$baseMean > 12,] # Remove genes with very low expression
medresdata <- fullresdata[!(is.na(fullresdata$log2FoldChange)),]
medians <- ddply(medresdata, .(chr), summarise, med = 2^(median(log2FoldChange)))
medians$asympt_Disp <- asympt_Dispersion
medians$extraPoisson <- extraPoisson_noise
medians$rep <- rep_number
medians$scale <- scale_number
#  medians$n <- n_number
median_scaleddata <- rbind(median_scaleddata,medians)
medians





