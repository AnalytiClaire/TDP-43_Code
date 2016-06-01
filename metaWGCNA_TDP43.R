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
