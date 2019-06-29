## @Author: Mihaela Angelova
## Email: mihaela.angelova@crick.ac.uk

# 'signature' accepts the value of LM22  or HD (in-house defined HD signature, Table S4)
signature = "LM22"  # HD or LM22
absolute = TRUE
analysisID <- paste("CIBERSORT", signature, absolute, sep = ".")
groupLevels <-  c( "normal", "low grade", "high grade", "SCC" )
# Load color settings
source("color_settings.R")

if(!require(reshape2)) {
  install.packages("reshape2")
}

if(!require(ggplot2)) {
  install.packages("ggplot2")
}

if(!require(plyr)) {
  install.packages("ggplot2")
}

# Cells excluded from Fig 2a, but included in Fig. S4a
cellsExclude = c("Tregs", "Monocytes", "Plasma.cells", "Eosinophils")
## Cell types renamed
# T cells follicular helper  "Tfh"
# T cells regulatory (Tregs) "Tregs"
# T cells gamma delta        "Tgd"
# Dendritic cells resting    "DC resting"
# Dendritic cells activated  "DC activated"

fixed_random_factors <- c("group", "stage", "patient", "previous.lung.cancer", "smoking.status")
outputDir = "."

##Load clinical data and CIBERSORT output for the selected signature (LM22 or HD)
load(file.path(outputDir, paste0("CIBERSORT.Output_Abs_", signature, ".RData")))

dataMelt = data.frame( melt(data.ciber, id = "Input Sample") )
match.melt = match( dataMelt$Input.Sample, data.clin$file )
dataMelt = cbind(dataMelt, 
                  data.clin[match.melt, fixed_random_factors])

dataMelt = dataMelt[ , grep( "Correlation|P-value|RMSE|Absolute", colnames(dataMelt), value = T, invert = T ) ]

dataMelt$group = factor(dataMelt$group, levels = groupLevels)
dataMelt$Input.Sample = factor(dataMelt$Input.Sample, levels = data.ciber$`Input Sample`)
dataMelt = dataMelt[! dataMelt$variable %in% cellsExclude, ]
dataMelt$CellType = dataMelt$variable = 
  factor(dataMelt$variable, levels = rev(names(ciberRainbow)))
dataMelt = dataMelt[ !is.na(dataMelt$variable), ]

dataMelt$matched = dataMelt$patient %in% 
  data.clin$patient[ duplicated(data.clin$patient )]
dataMelt = droplevels(dataMelt)

print("Plotting average stacked bars per stage")
mrg = ddply(dataMelt, .( CellType, group, stage ), summarize, median = mean(value))
mrg = mrg[ !mrg$CellType %in% cellsExclude, ]
pdf(file.path(outputDir, paste0( "barplot.", analysisID, ".pdf")),
    width = 6.5, height = 5 ) #7.5 for LM22
g <- ggplot( mrg, aes( as.factor(stage), median, fill = CellType ) )
plot = g + geom_bar( width = 1, stat = "identity", position = "fill", color = "white" )+
  theme_classic(base_size = 14) + theme(legend.position = "left") +
  theme( axis.text.x = element_text( size = 14), 
         legend.key.size = unit( 2, "mm" ), 
         legend.text     = element_text( size = 4),
         panel.spacing   = unit(0, "lines"),
         legend.position = "right",
         legend.title    = element_blank()) +
  scale_fill_manual(values = ciberRainbow ) +
  xlab( "Developmental stage" ) + ylab( "Absolute abundance")
print(plot)
dev.off()

pdf( file.path( outputDir, paste0( "Intralesional_heterogeneity_Fig5a_", analysisID, ".pdf" ) ),
     width = 12, height = 5 )
g <- ggplot(dataMelt[ dataMelt$matched, ], aes(as.factor(stage), value, fill = CellType))
plot = g + geom_bar(stat = "identity", position = "fill", color = "transparent") +
  # guides(fill = FALSE) +
  theme_bw() + theme(legend.position = "top") +
  facet_grid( . ~ as.factor(patient), scales = "free", space = "free" ) +
  theme(axis.text.x     = element_text( angle = 90, size = 12 ),
        legend.key.size = unit( 0.8, "mm" ),
        legend.text     = element_text( size= 6),
        legend.title    = element_blank()) +
  scale_fill_manual(values = ciberRainbow) +
  xlab("Developmental stage")
print(plot)
dev.off()
