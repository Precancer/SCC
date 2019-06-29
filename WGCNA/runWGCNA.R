## @Author Jennifer Beane, Kahkeshan Hijazi
## Section of Computational Biomedicine, Department of Medicine, Boston University School of Medicine, Boston

## 1. Runs mixed-effects model on probes on 9 stages (factor variable)
## 2. Weighted-gene correlation network analysis on the significant genes from step 1
rm(list=ls())

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if(!require("GO.db",character.only = TRUE)) {
  BiocManager::install("GO.db")
}
if(!require("impute")) {
  BiocManager::install("impute")
}
if(!require("minet")) {
  BiocManager::install("minet")
}
if(!require("org.Hs.eg.db")) {
  BiocManager::install("org.Hs.eg.db")
}
if(!require(nlme)) {
  install.packages("nlme")
}
if(!require(WGCNA)) {
  install.packages("WGCNA", dep = T)
}
if (!require("heatmap3",character.only = TRUE)) {
  install.packages("heatmap3",dep=TRUE)
}

load("inputDataExpression.RData")

data<-d[,2:123]
rownames(data)<-as.character(d[,1])
d.gene<-d[,124:127]
rownames(d.gene)<-as.character(d[,1])

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

#Run mixed effect model
res<-c();
good.ind<-c();
for(i in 1:nrow(data))  {
  # y is the expression levels of gene i
  y <- as.numeric(data[i,])
  # runs the linear model with patient as the random effect
  if(class(try( lme(y ~ gender+age+smoke+p.lc+stage, random = ~1|pat, method="ML"))) != "try-error") {
    if(class(try( lme(y ~ gender+age+smoke+p.lc, random = ~1|pat, method="ML"))) != "try-error") {
      
      model1 <- lme(y~ gender+age+smoke+p.lc+stage, random = ~1|pat, method="ML")
      # runs the baseline model, removing the cancer term that you're interested in
      model2 <- lme(y ~ gender+age+smoke+p.lc, random = ~1|pat, method="ML")
      good.ind<-c(good.ind,i);
      # save the coefficient (betas) from the full model for the site and cancer terms
      coefficients <- c(summary(model1)$coefficients$fixed)
      # save the T Values from the full model for the site and cancer terms
      tvals <- c(summary(model1)$tTable[,4:5])
      # calculate the ANOVA model comparison for the cancer term
      pvals <- anova(model1, model2)$p[2]
      next.row <- c(coefficients,tvals,pvals)
      res <- rbind(res, next.row)
    }
  }
}
colnames(res)<-c(rep(names(summary(model1)$coefficients$fixed),3),"ANOVAp")
res.fdr <- p.adjust(res[,dim(res)[2]], method="fdr")
res <- cbind(res, res.fdr)
rownames(res)<-rownames(data)[good.ind]
write.table(res, "9stage_factor.txt", sep="\t")

sig.ind<-which(res[,41]<0.001)
write.table(res[sig.ind,], "9stage_factor_fdr001.txt", sep="\t")


#conduct WGCNA analysis
allowWGCNAThreads()
#define input data matrix
set<-t(data[match(rownames(res)[sig.ind],rownames(data)),])

#Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
#Call the network topology analysis function
sft = pickSoftThreshold(set, powerVector = powers, verbose = 5)
#Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
pdf(file="9stage_factor_results_fdr001_4734genes_PowerSelecting_graphs.pdf")
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
#this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
#Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

net = WGCNA::blockwiseModules(set, power = 12, minModuleSize = 50,
                       reassignThreshold = 0, mergeCutHeight = 0.3, networkType="signed",
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "9stage_results_fdr001_4734genes_min50GenesPerCluster_mergecut_0.3_KH",
                       verbose = 3)
MEs <- net$MEs
#"ME7" "ME1" "ME5" "ME6" "ME2" "ME3" "ME4" "ME0"

###############Make dendogram#####################
# open a graphics window
#sizeGrWindow(12, 9)
pdf(file="9stage_results_fdr001_4734genes_50GenesPerCluster_mergecut_0.3_7modules_dendogram.pdf")
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
#########################
moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs

#moduleColors
#    black      blue     brown     green      grey       red turquoise    yellow 
#       64      1256       280       183        18        79      2662       192 
#table(net$colors)
#   0    1    2    3    4    5    6    7 
#  18 2662 1256  280  192  183   79   64 
#grey module is the junk module and contains genes that cannot be placed in co-expression modules

###see how similar MEs are
dissimME=(1-t(cor(MEs, method="p")))/2
hclustdatME=hclust(as.dist(dissimME), method="average" )
# Plot the eigengene dendrogram
pdf(file="dendogram_MEs.pdf")
par(mfrow=c(1,1))
plot(hclustdatME, main="Clustering tree based of the module eigengenes")
dev.off()

