setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/DEG_Test2/")

A1 <- read.csv("C9rankeduniqueresult.csv")
A2 <- read.csv("CHrankeduniqueresult.csv")
A3 <- read.csv("sALSrankeduniqueresult.csv")
A4 <- read.csv("FTLDrankeduniqueresult.csv")
A5 <- read.csv("VCPrankeduniqueresult.csv")

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/")
names <- read.table("5500genes.txt")
names <- names$V1

subsetC9 <- subset(A5, A5$Gene.Symbol %in% names, drop = TRUE)
subsetC9 <-subsetC9[!duplicated(subsetC9[,15]),]
rownames(subsetC9) <- subsetC9$Gene.Symbol
subsetC9 <- subsetC9[,52:58]
subsetC9[,(ncol(subsetC9)+1)] <- rowMedians(subsetC9) 
subsetC9[,(ncol(subsetC9)-1)] <- rownames(subsetC9)
subsetC9 <- subsetC9[,(ncol(subsetC9)-1):ncol(subsetC9)]

setwd("/Users/clairegreen/Desktop/")
write.csv(subsetC9, file = "subset.csv")
