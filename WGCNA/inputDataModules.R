data.clin <- read.csv("/Volumes/My Passport/work/projects/precancer/data/metadata/lung/metadata.lung.txt", sep = "\t", header = T )
file.exp <- file.path("/Volumes/My Passport/work/projects/precancer/analyses/lung/WGCNA/module_gene_expression.txt")

data.exp <- read.csv( file.exp, sep = "\t", header = T, row.names = 1)
data.exp = data.exp[, !colnames(data.exp) %in% c("Common", "Genbank", "Gene.Symbol", "Map")]
data = t(data.exp)

save(data.clin, data, file = "inputModuleData.RData")
