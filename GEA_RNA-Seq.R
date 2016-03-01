##RNA-Seq Gene Expression Analysis using Limma##

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GeneExpressionAnalysis/RNA-seq/")
# Counts <- read.table(file = 'GSE67196_Petrucelli2015_ALS_genes.rawcount.txt', header = TRUE)
# 
# write.csv(x = Counts, file = "counts_petrucelli.csv")

Counts <- read.csv(file = "counts_petrucelli.csv", header = TRUE)
Counts[,1] <- NULL #do twice
Counts[,1] <- NULL #do twice

Counts[Counts == 0] <- NA
Counts<-subset( Counts, subset=(GeneID !="NA")) #if no gene symbol, discount

# Countszero <-subset(Counts, subset=(row !=0))
# Countszero <- apply(Counts, 1, function(row) all(row !="NA"))
# Counts <- Counts[Countszero,]

library()


