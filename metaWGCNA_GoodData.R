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


PETExpr <- read.csv(file = "PETEXPRSrankeduniqueresult.csv")
rownames(PETExpr) <- PETExpr$hgnc_symbol
datExprPET <- PETExpr[,9:26]

RAVExpr <- read.csv(file = "RAVEXPRSrankeduniqueresult.csv")
rownames(RAVExpr) <- RAVExpr$hgnc_symbol
datExprRAV <- RAVExpr[,17:28]

FTLDExpr <- read.csv(file = "FTLDrankeduniqueresult.csv")
rownames(FTLDExpr) <- FTLDExpr$Probe.Set.ID
datExprFTLD <- FTLDExpr[,57:72]

IDs <- read.csv("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/WGCNA/GoodData/FTLDSymbols.csv")
genes <- IDs$Gene.Symbol
probeID <- IDs$Probe.Set.ID

keepGenesDups = (collapseRows(datExprFTLD,genes,probeID))[[2]]
rownames(datExprFTLD)<-keepGenesDups[,1]

commonProbesA = intersect (rownames(datExprA1),rownames(datExprA2))
commonProbesA = intersect (commonProbesA, rownames(datExprA3))