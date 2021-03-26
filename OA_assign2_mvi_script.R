#Set wd
setwd("/Users/mariaingersoll/Desktop/BU2020-2021/Ecol_and_Enviro_Genomics/BI586-git/OA_assign2")

#Load packages
library(DESeq2) 
library(affycoretools)
library(arrayQualityMetrics)
library(genefilter)
library(Biobase)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(vegan)
library(ggrepel)
library(tidyverse)
library(RColorBrewer)
library(gplots)
library(VennDiagram)

#Check for outliers (Open the html index file and analyze for outliers)
countData <- read.table("1.Assign2_GE_Sym_OA.csv", sep=",", header=T)
head(countData)
row.names(countData)=countData$X
countData$X=NULL
length(countData[,1])
#41850
treat=c( "Control", "P2553", "P2553", "Control")
g=data.frame(treat)
g
colData<- g
dds<-DESeqDataSetFromMatrix(countData=countData, colData=g, design=~treat)
vsd.ge=assay(vst(dds))
rl=vst(dds)
e=ExpressionSet(assay(rl), AnnotatedDataFrame(as.data.frame(colData(rl))))
arrayQualityMetrics(e,intgroup=c("treat"),force=T)

#Re-read in counts and Plot counts.
countData <- read.table("1.Assign2_GE_Sym_OA.csv", sep=",", header=T)
head(countData)
row.names(countData)=countData$X
countData$X=NULL
length(countData[,1])
totalCounts=colSums(countData)
totalCounts
barplot(totalCounts, col=c("coral", "blue", "blue", "coral"), ylab="raw counts")

#DESeq analysis to see how our data varies by treatment. Re-read in counts then perform DESeq function.
countData <- read.table("1.Assign2_GE_Sym_OA.csv", sep=",", header=T)
head(countData)
row.names(countData)=countData$X
countData$X=NULL
length(countData[,1])
treat=c( "Control", "P2553", "P2553", "Control")
g=data.frame(treat)
g
colData<- g


dds<-DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~treat) 

dds<-DESeq(dds)
head(dds)
res<- results(dds)

#Look at dispersion plot, which tells you the spread of your data.
plotDispEsts(dds, main="Dispersion plot")

#P2553 vs. Control pairwise comparison 
#(we're looking at gene expression RELATIVE to the control). MA plot to look at the depth of convergence and general trends of log fold change.
resP2553 <- results(dds, contrast=c("treat","P2553","Control"))
table(resP2553$padj<0.01)
summary(resP2553)

nrow(resP2553[resP2553$padj<0.05 & !is.na(resP2553$padj),])

plotMA(resP2553, main="Control vs P2553")
plotMA(resP2553, main="Control vs P2553", ylim=c(-2,2))

#Use pairwise comparison to determine the amount of up and downregulated genes in the P2553 treatment.
results <- as.data.frame(resP2553)
head(results)

nrow(resP2553[resP2553$padj<0.1 & resP2553$log2FoldChange > 0 & !is.na(resP2553$padj),])
nrow(resP2553[resP2553$padj<0.1 & resP2553$log2FoldChange < 0 & !is.na(resP2553$padj),])

write.table(resP2553, file="P2553.txt", quote=F, sep="\t")

cd <- read.table("P2553.txt")
head(cd)

#Make the GO table.
go_input_P2553 = cd %>%
  tibble::rownames_to_column(var = "iso") %>%
  mutate(mutated_p = -log(pvalue)) %>%
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  na.omit() %>%
  select(iso, mutated_p_updown)
head(go_input_P2553)
colnames(go_input_P2553) <- c("gene", "pval")
head(go_input_P2553)
write.csv(go_input_P2553, file="P2553_GO.csv", quote=F, row.names=FALSE)

#Get p values.
val2553=cbind(resP2553$pvalue, resP2553$padj)
head(val2553)
colnames(val2553)=c("pval.2553", "padj.2553")
length(val2553[,1])
table(complete.cases(val2553))

#Make rlogdata and pvals table.
rlog=rlogTransformation(dds, blind=TRUE) 
rld=assay(rlog)
head(rld)
colnames(rld)=paste(colData$treat)
head(rld)
length(rld[,1])

rldpvals=cbind(rld,val2553)
head(rldpvals)
dim(rldpvals)
table(complete.cases(rldpvals))

write.csv(rldpvals, "OA_assign2_RLDandPVALS.csv", quote=F)

colnames(rld)=paste(colData$treat)
head(rld)

#Create a sample distance heat map.
sampleDists <- as.matrix(dist(t(rld)))
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          margin=c(10, 10), main="Sample Distance Matrix")

#heat map of sample distances for pco2, overall similarity between samples based on the rlog values
rldpvals <- read.csv(file="OA_assign2_RLDandPVALS.csv", row.names=1)
head(rldpvals)
rld=rldpvals[,1:4]
head(rld)

sampleDists <- dist(t(rld))
sampleDistMatrix <- as.matrix( sampleDists )
treat=c( "Control_1", "P2553_1", "P2553_2", "Control_2")
colnames(sampleDistMatrix)=paste(treat)
rownames(sampleDistMatrix)=paste(treat)

heat.colors = colorRampPalette(rev(c("blue","yellow","red")),bias=0.3)(100)
pheatmap(sampleDistMatrix,color = heat.colors,cex=0.9,border_color=NA,cluster_rows=T,cluster_cols=T)

#PCA time
rld_t=t(rld)
pca <- prcomp(rld_t,center = TRUE)
head(pca)
li <- pca$sdev^2 / sum(pca$sdev^2)
pc1v <- round(li[1] * 100, 1)
pc2v <- round(li[2] * 100, 1)
pca_s <- as.data.frame(pca$x)
head(pca_s)
pca_s <- pca_s[,c(1,2)]
pca_s$Samples = row.names(pca_s)
pca_s$treat=colData$treat
head(pca_s)

cbPalette <- c("darkgoldenrod2",  "darkolivegreen3", "dodgerblue3")
ggplot(pca_s, aes(PC1, PC2, color = treat, pch = treat)) +
  geom_point(size=3) +
  #  geom_text_repel(aes(label=Samples)) +
  scale_colour_manual(values=cbPalette)+
  theme_bw() +
  # geom_density2d(alpha=.5)+
  geom_polygon(alpha=.2)+
  xlab(paste0("PC1: ",pc1v,"% variance")) +
  ylab(paste0("PC2: ",pc2v,"% variance")) 
head(pca)

#Creating a list of all the genes that have a pval <=0.1
head(rldpvals)
rldpvals_10 = filter(rldpvals, pval.2553 <= 0.1, preserve = TRUE)
length(rldpvals_10[,1])
#903
length(rldpvals[,1])
rld_10 = rldpvals_10[,1:4]
rld_10 = rownames_to_column(rld_10, var = "Gene_ID")
head(rld_10)

#Make the variable "gene" with the gene symbol called description
head(gene)
gene = read.delim("davies_cladeC_iso2gene.tab.txt", sep = "\t")%>%
  mutate(Description = gsub(".* GN=", "", Description)) %>%
  mutate(Description = gsub(" .*", "", Description))
head(gene)

#Join your two variables rld_10 and gene by "Gene_ID"
sigDEG = rld_10 %>%
  left_join(gene) %>%
  mutate(Description = make.names(Description, unique = TRUE)) %>%
  column_to_rownames(var = "Description") %>%
  dplyr::select(-Gene_ID) %>%
  drop_na() %>%
  dplyr::select(sort(current_vars()))
head(sigDEG)
nrow(sigDEG)
#903, nice!

#Now make a heatmap of all the DEGs
col0=colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
heatmap_sigDEG = heatmap.2(as.matrix(sigDEG), col = col0, Rowv = TRUE, Colv = TRUE, scale = "row",
                           dendrogram = "both",
                           trace = "none",
                           main = "Significant DEG",
                           margin = c(5,15))



