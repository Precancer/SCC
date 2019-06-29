# prepare input data and save to RData file
clinFile = "../Celine_Paper/mascaux_annot.txt" #"/Volumes/My Passport/work/projects/precancer/data/metadata/lung/metadata.lung.txt"
expProbeFile = "/Volumes/My Passport/work/projects/precancer/analyses/lung/expression/expression.lung.probe.txt"
d<-read.delim(file=expProbeFile,sep="\t",header=T)
s<-read.table(file=clinFile,sep="\t",header=T)
save(d, s, file = "inputDataExpression.RData")
