setwd(dir="/Users/clairegreen/Downloads/metaAnalysisFiles")

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

load("metaAnalysisData.RData")

#  (Section will take ~5-10 minutes to run)
# source("collapseRows_NEW.R") # ONLY uncomment this line if you get an error with it commented
datExprB1g = (collapseRows(datExprB1,genesI,probesI))[[1]]
datExprB2g = (collapseRows(datExprB2,genesA,probesA))[[1]]

commonProbesA = intersect (rownames(datExprA1),rownames(datExprA2))
datExprA1p = datExprA1[commonProbesA,]
datExprA2p = datExprA2[commonProbesA,]

commonGenesB = intersect (rownames(datExprB1g),rownames(datExprB2g))
datExprB1g = datExprB1g[commonGenesB,]
datExprB2g = datExprB2g[commonGenesB,]


Exp <- t(datExprA1)
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


softPower = 10 # (Read WGCNA tutorial to learn how to pick your power)
rankExprA1= rank(rowMeans(datExprA1p))
rankExprA2= rank(rowMeans(datExprA2p))
random5000= sample(commonProbesA,5000)
rankConnA1= rank(softConnectivity(t(datExprA1p[random5000,]),type="signed",power=softPower))
rankConnA2= rank(softConnectivity(t(datExprA2p[random5000,]),type="signed",power=softPower))

rankExprB1= rank(rowMeans(datExprB1g))
rankExprB2= rank(rowMeans(datExprB2g))
random5000= sample(commonGenesB,5000)
rankConnB1= rank(softConnectivity(t(datExprB1g[random5000,]),type="signed",power=softPower))
rankConnB2= rank(softConnectivity(t(datExprB2g[random5000,]),type="signed",power=softPower))

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

