
#Load database of associations
setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GWAS/")

A <- read.table(file = "signif.snp.GWAScentral.txt")
a <- A$V1

B <- read.table(file = "signif.snp.GWAScentral.p0.0001.1.txt")
b <- B$V1

C <- read.table(file = "GWAScentral_<=0.0001.txt")
c <- c$x

D <- read.table(file = "GWAScentral_<=0.00001.txt")
d <- D$x

E <- read.table(file = "signif.snp.NeuroX.txt")
e <- E$V1

F <- read.table(file = "signif.snp.NeuroX.p5E08.txt")
f <- F$V1


#load test file

setwd (dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GSEA/PCxN Example/")
Z <- read.csv(file = "C9pcxngenes.csv")
y <- Z$NEGATIVE_REGULATION_OF

#remove any duplicates
y <- y[!duplicated(y)]

#Intersect
overlap <- Reduce(intersect, list(o2, o3))
print(overlap)



#Load file with all genes
library(hgu133plus2.db)
sym <- hgu133plus2SYMBOL
sym1 <- mappedkeys(sym)
sym2 <- as.list (sym[c(sym1)])
sym3 <- data.frame (sym2)
sym.probes <- names (sym2)
sym.genes <- sym3[1,]

#calculate counts

# n <- max(length(a), length(b), length(c), length(d))
# length(a) <- n
# length(b) <- n
# length(c) <- n
# length(d) <- n
# 
# x <- cbind(a,b,c,d)
# 
# for (i in x[,1:4]) {

#How many PCxN genes contain snps
x.in <- length (which(y %in% c)) 
#how many do not
x.out <- length(y) - x.in
#total number of snps
tot.in <- length (GC.04$HGNC.Gene.Symbol)
#total number of all genes
tot.out <- length (sym.genes)

#create count matrix
counts <- matrix (nrow=2, ncol=2)
counts [1,] <- c(x.in, tot.in)
counts [2,] <- c(x.out, tot.out)

#Conduct fisher's exact test for count data
a5 <-fisher.test (counts)
enrich <- a5$p
print(enrich)

