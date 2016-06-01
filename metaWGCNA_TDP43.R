setwd(dir="/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/DEG_Test2/")

library(WGCNA)
library(impute)
library(dynamicTreeCut)
library(qvalue)
library(flashClust)
library(lattice)
library(survival)
library(Formula)
library(ggplot2)
library(Hmisc)

C9Expr <- read.csv(file = "C9rankeduniqueresult.csv")
rownames(C9Expr) <- C9Expr$Probe.Set.ID
datExprA1 <- C9Expr[,52:59] #C9 patients
datExprA2 <- C9Expr[,49:51] #C9 controls

CHExpr <- read.csv(file = "CHrankeduniqueresult.csv")
rownames(CHExpr) <- CHExpr$Probe.Set.ID
datExprB1 <- CHExpr[,55:57] #C9 patients
datExprB2 <- CHExpr[,49:54] #C9 controls

IDs <- read.csv("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/WGCNA/C9orf72/C9ID geneIDs.csv")
IDs <- C9Expr[,c("Probe.Set.ID", "Gene.Symbol")]
genes <- IDs$Gene.Symbol
probeID <- IDs$Probe.Set.ID

#For matching probes if from different platforms
#  (Section will take ~5-10 minutes to run)
# source("collapseRows_NEW.R") # ONLY uncomment this line if you get an error with it commented
datExprB1g = (collapseRows(datExprB1,genes,probeID))[[1]]
datExprB2g = (collapseRows(datExprB2,genes,probeID))[[1]]

#Limit analysis to common probes
commonProbesA = intersect (rownames(datExprA1),rownames(datExprA2))
datExprA1p = datExprA1[commonProbesA,]
datExprA2p = datExprA2[commonProbesA,]

commonGenesB = intersect (rownames(datExprB1g),rownames(datExprB2g))
datExprB1g = datExprB1g[commonGenesB,]
datExprB2g = datExprB2g[commonGenesB,]

Exp <- t(datExprA2p)

###Choosing soft threshold
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(Exp, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.50,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/WGCNA/")
#Correlating general network properties
softPower = 6 # (Read WGCNA tutorial to learn how to pick your power)
rankExprA1= rank(rowMeans(datExprA1p))
rankExprA2= rank(rowMeans(datExprA2p))
random5000= sample(commonProbesA,5000)
rankConnA1= rank(softConnectivity(t(datExprA1p[random5000,]),type="signed",power=softPower, minNSamples = 3))
rankConnA2= rank(softConnectivity(t(datExprA2p[random5000,]),type="signed",power=softPower, minNSamples = 3))

rankExprB1= rank(rowMeans(datExprB1g))
rankExprB2= rank(rowMeans(datExprB2g))
random5000= sample(commonGenesB,5000)
rankConnB1= rank(softConnectivity(t(datExprB1g[random5000,]),type="signed",power=softPower, minNSamples = 3))
rankConnB2= rank(softConnectivity(t(datExprB2g[random5000,]),type="signed",power=softPower, minNSamples = 3))

pdf("generalNetworkProperties.pdf", height=10, width=9)
par(mfrow=c(2,2))
verboseScatterplot(rankExprA1,rankExprA2, xlab="Ranked Expression (A1)", 
                   ylab="Ranked Expression (A2)")
verboseScatterplot(rankConnA1,rankConnA2, xlab="Ranked Connectivity (A1)", 
                   ylab="Ranked Connectivity (A2)")
verboseScatterplot(rankExprB1,rankExprB2, xlab="Ranked Expression (B1)", 
                   ylab="Ranked Expression (B2)")
verboseScatterplot(rankConnB1,rankConnB2, xlab="Ranked Connectivity (B1)", 
                   ylab="Ranked Connectivity (B2)")
dev.off()

#calculate all of the necessary values to run WGCNA
#(this will take around 10 minutes)
adjacencyA1 = adjacency(t(datExprA1p),power=softPower,type="signed");
diag(adjacencyA1)=0
dissTOMA1   = 1-TOMsimilarity(adjacencyA1, TOMType="signed")
geneTreeA1  = flashClust(as.dist(dissTOMA1), method="average")

adjacencyA2 = adjacency(t(datExprA2p),power=softPower,type="signed");
diag(adjacencyA2)=0
dissTOMA2   = 1-TOMsimilarity(adjacencyA2, TOMType="signed")
geneTreeA2  = flashClust(as.dist(dissTOMA2), method="average")

adjacencyB1 = adjacency(t(datExprB1g),power=softPower,type="signed");
diag(adjacencyB1)=0
dissTOMB1   = 1-TOMsimilarity(adjacencyB1, TOMType="signed")
geneTreeB1  = flashClust(as.dist(dissTOMB1), method="average")

adjacencyB2 = adjacency(t(datExprB2g),power=softPower,type="signed");
diag(adjacencyB2)=0
dissTOMB2   = 1-TOMsimilarity(adjacencyB2, TOMType="signed")
geneTreeB2  = flashClust(as.dist(dissTOMB2), method="average")

# save.image("tutorial.RData")  #  (Section will take ~5-15 minutes to run)

pdf("dendrogram.pdf",height=6,width=16)
par(mfrow=c(1,2))
plot(geneTreeA1,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity (A1)", labels=FALSE,hang=0.04);
plot(geneTreeA2,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity (A2)", labels=FALSE,hang=0.04); 
plot(geneTreeB1,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity (B1)", labels=FALSE,hang=0.04);
plot(geneTreeB2,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity (B2)", labels=FALSE,hang=0.04); 
dev.off()

#determine modules based on control data set

mColorh=NULL
for (ds in 0:3){
  tree = cutreeHybrid(dendro = geneTreeB2, pamStage=FALSE,
                      minClusterSize = (30-3*ds), cutHeight = 0.99, 
                      deepSplit = ds, distM = dissTOMB2)
  mColorh=cbind(mColorh,labels2colors(tree$labels));
}
pdf("Module_choices.pdf", height=10,width=25); 
plotDendroAndColors(geneTreeB2, mColorh, paste("dpSplt =",0:3), main = "",dendroLabels=FALSE);
dev.off()

modulesA1 =  mColorh[,4] #choose based on deepslit values in plot

#calculate the principle components for visualizations 
PCs1A    = moduleEigengenes(t(datExprB2g),  colors=modulesA1) 
ME_1A    = PCs1A$eigengenes
distPC1A = 1-abs(cor(ME_1A,use="p"))
distPC1A = ifelse(is.na(distPC1A), 0, distPC1A)
pcTree1A = hclust(as.dist(distPC1A),method="a") 
MDS_1A   = cmdscale(as.dist(distPC1A),2)
colorsA1 = names(table(modulesA1))
pdf("ModuleEigengeneVisualizations.pdf",height=6,width=8)
par(mfrow=c(1,1), mar=c(0, 3, 1, 1) + 0.1, cex=1)

plot(pcTree1A, xlab="",ylab="",main="",sub="")

plot(MDS_1A, col= colorsA1,  main="MDS plot", cex=2, pch=19)

ordergenes = geneTreeB2$order
plotMat(scale(log(datExprB2g[ordergenes,])) , rlabels= modulesA1[ordergenes], clabels= colnames(datExprB2g), rcols=modulesA1[ordergenes])

for (which.module in names(table(modulesA1))){
  ME = ME_1A[, paste("ME",which.module, sep="")] 
  barplot(ME, col=which.module, main="", cex.main=2, 
          ylab="eigengene expression",xlab="array sample") 
} 

dev.off()

##Qualitatively and quantitatively measure network preservation at the module level##

#assess how well modules in network 1 are preserved in network 2
pdf("Final_modules.pdf",height=8,width=12)
plotDendroAndColors(geneTreeA1, modulesA1, "Modules", dendroLabels=FALSE, hang=0.03, addGuide=TRUE, 
                    guideHang=0.05, main="Gene dendrogram and module colors (A1)") 
plotDendroAndColors(geneTreeA2, modulesA1, "Modules", dendroLabels=FALSE, hang=0.03, addGuide=TRUE, 
                    guideHang=0.05, main="Gene dendrogram and module colors (A2)") 
plotDendroAndColors(geneTreeB1, modulesA1, "Modules", dendroLabels=FALSE, hang=0.03, addGuide=TRUE, 
                    guideHang=0.05, main="Gene dendrogram and module colors (A1)") 
plotDendroAndColors(geneTreeB2, modulesA1, "Modules", dendroLabels=FALSE, hang=0.03, addGuide=TRUE, 
                    guideHang=0.05, main="Gene dendrogram and module colors (A2)") 
dev.off()

