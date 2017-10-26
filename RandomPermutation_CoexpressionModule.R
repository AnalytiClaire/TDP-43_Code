###Random permutations for Differential Expression##
options(scipen=999)

DEGPPI <- readLines("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/DEG_PPI_Genes_nofib.txt")

DEGs <- readLines("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/FoldChangeResults/Filtered_upanddown.txt")

#Load file with all genes
library(hgu133plus2.db)
sym <- hgu133plus2SYMBOL
sym1 <- mappedkeys(sym)
sym2 <- as.list (sym[c(sym1)]) 
sym3 <- data.frame (sym2)
sym.probes <- names (sym2)
sym.genes <- sym3[1,]

sym.genes <- t(sym.genes)

allgenes <- sym.genes[!duplicated(sym.genes),]


setwd("~/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/PPI_Coexpression/CorrelationValue/")
C9_cor.5 <- read.csv("C9_cor.5_cytoscape.csv")
sALS_cor.5 <-read.csv("sALS_cor.5_cytoscape.csv")
FTLD_cor.5 <-read.csv("FTLD_cor.5_cytoscape.csv")
VCP_cor.5 <-read.csv("VCP_cor.5_cytoscape.csv")
PET_cor.5 <-read.csv("PET_cor.5_cytoscape.csv")
RAV_cor.5 <-read.csv("RAV_cor.5_cytoscape.csv")

#indicate the number of overlapping genes identified by DE analysis
test <- 7
samplenum <- 104
samplelist <- DEGPPI

m=10000 #number of repetitions 
r <- c(1:m) #store repetition numbers in vector "r"


for (j in 1:m){
  sample <- sample(samplelist, size=samplenum, replace=FALSE)
  random <- Reduce(intersect, list(sample, DEGs))
  r[j] <- length(random)
}


# ## Gene numbers found by running "Foldchangeupdown.R"
# for (j in 1:m){
#   random1 <- sample (samplelist, size=samplenum, replace=FALSE)
#   random2 <- sample (samplelist, size=samplenum, replace=FALSE)
#   random3 <- sample (samplelist, size=samplenum, replace=FALSE)
#   random4 <- sample (samplelist, size=samplenum, replace=FALSE)
#   random5 <- sample (samplelist, size=samplenum, replace=FALSE)
#   random6 <- sample (samplelist, size=samplenum, replace=FALSE)
#   random <- Reduce(intersect, list(random1, random2, random3, random4, random5, random6))

test1 <- which(r > test)  # count number of times r is larger than test value
result <- (length(test1)/m) # calculate P value
result
mean(r)
range(r)
