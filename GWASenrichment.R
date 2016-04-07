setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GWAS/")

GWAScen <- read.table(file = "signif.snp.GWAScentral.txt")
GWAScen <- GWAScen$V1

neuroX <- read.table(file = "signif.snp.NeuroX.p5E08.txt")
neuroX <- neuroX$V1

pthresh <- 1.000e-05

gwasthresh <- subset(GWAScen, p.value <= pthresh)

overlap <- Reduce(intersect, list(x, GWAScen))
print(overlap)






x <- read.table(file = "DEG+50.txt")
x <- x$V1

library(hgu133plus2.db)
sym <- hgu133plus2SYMBOL
sym1 <- mappedkeys(sym)
sym2 <- as.list (sym[c(sym1)])
sym3 <- data.frame (sym2)
sym.probes <- names (sym2)
sym.genes <- sym3[1,]

x.in <- length (which(x %in% GWAScen))
x.out <- length(x) - x.in
tot.in <- length (GWAScen)
tot.out <- length (sym.genes)

counts <- matrix (nrow=2, ncol=2)
counts [1,] <- c(x.in, tot.in)
counts [2,] <- c(x.out, tot.out)

a5 <-fisher.test (counts)
lists.enrich[4,i] <- a5$p
