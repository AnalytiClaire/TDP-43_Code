library(WGCNA)
### C9orf72 ###
# Display the current working directory
setwd ("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GeneExpressionAnalysis/TopGenes_2016-02-15/")

#Read in desired genes
C9Results <- read.csv ("C9rankeduniqueresult.csv", header=TRUE) #Taking only the genes we deemed acceptable through criteria. See details on
#gene expression analysis to find criteria

C9ID <- C9Results[,1]
C9ID <- data.matrix(C9ID)

#Read in raw expression values
setwd ("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/C9orf72_LCM/")
C9RawExp <- read.csv("eset_NineP_150612_exprs.csv") 

C9Exp <-merge(C9ID, C9RawExp, by.x="V1", by.y="X") #merge raw expression data with accepted genes
rownames(C9Exp) <- C9Exp[,1] #make probeset IDs row names
colnames(C9Exp) <- colnames(C9RawExp) #make file names column names
C9Exp <- cbind(C9Exp[,2:12]) #remove ID column

C9Pat <- C9Exp[,1:8]
C9Con <- C9Exp[,9:11]

C9Pat <- t(C9Exp)
C9Con <- t(C9Exp)

###PATIENT ANALYSIS###

# ###Choosing soft threshold
# # Choose a set of soft-thresholding powers
# powers = c(c(1:10), seq(from = 12, to=20, by=2))
# # Call the network topology analysis function
# sft = pickSoftThreshold(C9Pat, powerVector = powers, verbose = 5)
# # Plot the results:
# sizeGrWindow(9, 5)
# par(mfrow = c(1,2));
# cex1 = 0.9;
# # Scale-free topology fit index as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
#      main = paste("Scale independence"));
# text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      labels=powers,cex=cex1,col="red");
# # this line corresponds to using an R^2 cut-off of h
# abline(h=0.70,col="red")
# # Mean connectivity as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], sft$fitIndices[,5],
#      xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
#      main = paste("Mean connectivity"))
# text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

##SOFT THRESHOLD VALUE OF 6 SELECTED##

# C9Exp <- data.matrix(C9Exp) #csv files contain character matrices, the following code requires numeric

##One-step network construction and module detection
setwd ("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/WGCNA/C9orf72/")
net = blockwiseModules(C9Pat, power = 6,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "C9TOM",
                       verbose = 3)
table(net$colors)

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(C9Pat, power = 6);
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;

# Call the plot function
sizeGrWindow(9,9)
TOMplot(plotTOM, moduleColors, main = "Network heatmap plot, all genes")

# Recalculate module eigengenes
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
# Isolate weight from the clinical traits
weight = as.data.frame(datTraits$weight_g);
names(weight) = "weight"
# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, weight))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)