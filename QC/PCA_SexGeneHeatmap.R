## @Author Jennifer Beane, Kahkeshan Hijazi 
## Section of Computational Biomedicine, Department of Medicine, Boston University School of Medicine, Boston

## randomized Principal Component Analysis
## Clustering by sex-gene expression

rm(list=ls())

if (!require("heatmap3",character.only = TRUE)) {
  install.packages("heatmap3",dep=TRUE)
}
  # library(heatmap3)
if (!require("ggplot2",character.only = TRUE)) {
  install.packages("ggplot2",dep=TRUE)
}

load("inputData.RData")
# d - expression data on probe level
# s - clinical information
# annot - microarray annotation file

data<-d[,2:123]
rownames(data)<-as.character(d[,1])
d.gene<-d[,124:127]
rownames(d.gene)<-as.character(d[,1])
annot.gene<-annot[match(rownames(data),annot[,1]),]

#define phenotypic variables
pat<-s$patient
pat<-as.factor(pat)
stage<-s$phenotypical.stage
stage<-gsub("normal normofluorescent",1,stage)
stage<-gsub("normal hypofluorescent",2,stage)
stage<-gsub("hyperplasia",3,stage)
stage<-gsub("metaplasia",4,stage)
stage<-gsub("mild dysplasia",5,stage)
stage<-gsub("moderate dysplasia",6,stage)
stage<-gsub("severe dysplasia",7,stage)
stage<-gsub("carcinoma in situ",8,stage)
stage<-gsub("squamous cell carcinoma",9,stage)
stage<-as.factor(stage)
gender<-as.factor(s$gender)
age<-as.numeric(s$age)
pat<-as.factor(s$patient)
smoke<-as.factor(s$smoking.status)
p.lc<-as.factor(s$previous.lung.cancer)
batch<-as.factor(s$Hybridization.day.Run)

pdf(file="PCA_plot.pdf")
z.data<-t(scale(t(data),center=T,scale=T))
pca<-prcomp(z.data,scale=F,center=F)
mean.pc1<-mean(pca$rotation[,1])
sd.pc1<-sd(pca$rotation[,1])
mean.pc2<-mean(pca$rotation[,2])
sd.pc2<-sd(pca$rotation[,2])
outliers<-c(which(pca$rotation[,1]<(mean.pc1-sd.pc1*2)),
which(pca$rotation[,1]>(mean.pc1+sd.pc1*2)),
which(pca$rotation[,2]<(mean.pc2-sd.pc2*2)),
which(pca$rotation[,2]>(mean.pc2+sd.pc2*2)))
out.color<-rep("no",length=ncol(data))
out.color[outliers]<-"yes"
perc.var<-(pca$sdev)^2 / sum(pca$sdev^2)
percentVar <- round(100 * perc.var)
print(ggplot(as.data.frame(pca$rotation[,c(1,2)]), aes(PC1, PC2, color=stage, shape=smoke)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")))
print(ggplot(as.data.frame(pca$rotation[,c(1,2)]), aes(PC1, PC2, color=out.color)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")))
print(ggplot(as.data.frame(pca$rotation[,c(1,2)]), aes(PC1, PC2,color=batch)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")))
dev.off()

sex.genes <- read.table("gender_genes.txt",sep="\t",header=T)
sex.genes<-sex.genes[,1]
pdf(file="heatmap_sex_genes.pdf")
com.genes<-intersect(sex.genes,annot.gene$GeneSymbol)
sex.data<-data[match(com.genes,annot.gene$GeneSymbol),]
sex.col<-rep("pink",length(gender))
sex.col[which(gender=="Male")]<-"blue"
heatmap3(as.matrix(sex.data),margins = c(10, 10),labRow=sex.genes,method="ward.D2",cexRow=0.3,ColSideColors=sex.col,col=colorRampPalette(c("blue","white","red"), space="rgb")(255), balanceColor = T)
dev.off()
