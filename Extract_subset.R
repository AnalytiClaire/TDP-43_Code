#### Method for extracting genes from expression tables and taking the average across samples #####

library(matrixStats)

#Read in files
setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/noMedian/")

C9 <- read.csv("C9_unique.csv")
C9 <- C9[order(C9$P.Value),]
CH <- read.csv("CH_unique.csv")
CH <- CH[order(CH$P.Value),]
sals <- read.csv("sals_unique.csv")
sals <- sals[order(sals$P.Value),]
ftld <- read.csv("ftld_unique.csv")
ftld <- ftld[order(ftld$P.Value),]
vcp <- read.csv("vcp_unique.csv")
vcp <- vcp[order(vcp$P.Value),]

setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/TDP-43_DEseq2/")

pet <- read.csv("PET_results_keepfiltering.csv")
pet <- pet[!duplicated(pet$hgnc_symbol),]
rav <- read.csv("RAV_results_keepfiltering.csv")
rav <- rav[!duplicated(rav$hgnc_symbol),]

#Read in genes of interest
setwd("/Users/clairegreen/Desktop/")
names <- read.table("List5996.txt")
names <- names$V1


##### MICROARRAY #####
mat <- vcp

#Take the subset of genes from the tables
subset <- subset(mat, mat$Gene.Symbol %in% names, drop = TRUE)
rownames(subset) <- subset$Gene.Symbol
subset[,1] = NULL

#Take one condition (patients or controls)
setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Expression_DEGonly/5996/")
write.csv(subset, file = "vcp_5996_all.csv")
subset[,1:7] <- NULL
subset[,length(subset)] <- NULL
write.csv(subset, file = "vcp_5996_expr.csv")


#### RNA-SEQ ######
mat <- rav

#Take the subset of genes from the tables
subset <- subset(mat, mat$hgnc_symbol %in% names, drop = TRUE)
rownames(subset) <- subset$hgnc_symbol
subset[,length(mat)] = NULL

#Take one condition (patients or controls)
setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Expression_DEGonly/5996/")
write.csv(subset, file = "rav_5996_all.csv")
subset[,1:8] <- NULL
write.csv(subset, file = "rav_5996_expr.csv")
