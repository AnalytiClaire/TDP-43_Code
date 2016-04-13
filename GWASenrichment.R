
#Load database of associations
setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GWAS/")

a <- read.table(file = "signif.snp.GWAScentral.txt")
a <- a$V1

b <- read.table(file = "signif.snp.GWAScentral.p0.0001.1.txt")
b <- b$V1

c <- read.table(file = "signif.snp.NeuroX.txt")
c <- c$V1

d <- read.table(file = "signif.snp.NeuroX.p5E08.txt")
d <- d$V1

#load test file

setwd (dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GSEA/PCxN Example/")
Z <- read.csv(file = "allgenes2.csv")
y <- Z$RHO_PROTEIN_SIGNAL_TRANSDUCTION

#remove any duplicates
y <- y[!duplicated(y)]

#Intersect
overlap <- Reduce(intersect, list(y, a))
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

x.in <- length (which(y %in% d))
x.out <- length(y) - x.in
tot.in <- length (d)
tot.out <- length (sym.genes)

counts <- matrix (nrow=2, ncol=2)
counts [1,] <- c(x.in, tot.in)
counts [2,] <- c(x.out, tot.out)

#Conduct fisher's exact test for count data

a5 <-fisher.test (counts)
enrich <- a5$p
print(enrich)

