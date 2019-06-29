## @Author: Mihaela Angelova
## Email: mihaela.angelova@crick.ac.uk

options(stringsAsFactors = F)
if(!require("WGCNA")) {
  install.packages("WGCNA", dep = T)
}
if(!require("flashClust")) {
  install.packages("flashClust")
}

replace = FALSE # Sampling with (TRUE) or without (FALSE) replacement; logical value
fraction = if (replace) 1 else 0.63

cat("Importing data...\n")
load("inputDataModule.RData")
## data.clin - clinical information with developmental stage
## data <- Expression data for the genes in the modules
stages = data.clin$stage[match(rownames(data), data.clin$file)]
colors <- read.csv( "wgcna_clusters.txt", sep = "\t", header = T )$module
useSamples = unlist(lapply(split(1:nrow(data), stages), function(x) {
  sample(x, size = length(x) * fraction, replace = replace)
}))
adj1 <- adjacency(data[useSamples, ], power = power, type = "signed" )
adj2 <- adjacency(data[-useSamples, ], power = power, type = "signed" )
multiData <- multiData(Reference = adj1, Test = adj2)
multiColor <- list(Reference = colors)
test <- modulePreservation(multiData = multiData, testNetworks = 2, verbose=5,
                          calculateClusterCoeff = FALSE, nPermutations = 1,
                          maxModuleSize = max(table(colors)),
                          multiColor = multiColor, networkType = "signed", quickCor = 1)
save(test, file = "modulePreservation.RData")

##Plotting results

ref = 1 # reference set
test = 2 # test set

plotData = as.data.frame(
  do.call(rbind, lapply(1:length(res), function(x) {
    cbind(medianRank.pres= res[[x]]$preservation$observed[[ref]][[test]]$medianRank.pres,
          Zsummary.pres = res[[x]]$preservation$Z[[ref]][[test]]$Zsummary.pres,
          moduleSize = res[[x]]$preservation$Z[[ref]][[test]]$moduleSize,
          ßcluster = rownames(res[[x]]$preservation$Z[[ref]][[test]]))
  })))
plotData$moduleSize = as.numeric(as.character(plotData$moduleSize))
plotData$Zsummary.pres = as.numeric(as.character(plotData$Zsummary.pres))
plotData$medianRank.pres = as.numeric(as.character(plotData$medianRank.pres))
plotData = ddply(plotData, .(cluster), summarize, 
                 medianRank.pres=mean(medianRank.pres),
                 Zsummary.pres=mean(Zsummary.pres),
                 moduleSize = unique(moduleSize))
plotData = plotData[, c(2:ncol(plotData), 1)]
modColors = cols = c( "darkgrey", "gold", "black", "red", "blue", 
                      "#00B0F6",  "yellow", "green","magenta") 
names(modColors)=  rownames(res[[1]]$preservation$observed[[ref]][[test]])
plotMods = !(modColors %in% c("grey", "gold"));
ext = modColors[plotMods];
mains = c("Preservation Median rank", "Preservation Zsummary")
par(mfrow = c(1,2))
pdf("Module_preservation.pdf", width = 5, height = 5,
    useDingbats = F)
for (p in 1:2) {
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  if (p==2) {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(plotData[plotMods, 3], plotData[plotMods, p], col = 1, 
       bg = modColors[plotMods], pch = 21,
       main = "", #mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, max(plotData$moduleSize) ), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  legend("topleft", legend = plotData[plotMods, 4], 
         col = sapply( plotData[plotMods, 4], 
                       function(x) modColors[[x]]), bty = "n", pch = 20)
  if (p==2)  {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  }
}
dev.off()
