setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/")

GWAScen <- read.table(file = "signif.snp.GWAScentral.p0.0001.1.txt")
GWAScen <- GWAScen$V1

neuroX <- read.table(file = "signif.snp.NeuroX.p5E08.txt")
neuroX <- neuroX$V1

GWAScatALS <- read.table(file = "ALSGWAScatalog.txt")
GWAScatALS <- GWAScatALS$V1

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/Pathprint/")
pathprint <- read.table(file = "Pathprintgenes.txt")
pathprint <- pathprint$V1

pathprintunique <- pathprint[!duplicated(pathprint)]



# pthresh <- 1.000e-05
# 
# gwasthresh <- subset(GWAScen, p.value <= pthresh)
# 
overlap <- Reduce(intersect, list(pathprint, GWAScatALS))
print(overlap)




x <- read.table(file = "DEGs.txt")
x <- x$V1

library(hgu133plus2.db)
sym <- hgu133plus2SYMBOL
sym1 <- mappedkeys(sym)
sym2 <- as.list (sym[c(sym1)])
sym3 <- data.frame (sym2)
sym.probes <- names (sym2)
sym.genes <- sym3[1,]

x.in <- length (which(pathprintunique %in% neuroX))
x.out <- length(pathprintunique) - x.in
tot.in <- length (neuroX)
tot.out <- length (sym.genes)

counts <- matrix (nrow=2, ncol=2)
counts [1,] <- c(x.in, tot.in)
counts [2,] <- c(x.out, tot.out)

a5 <-fisher.test (counts)
enrich <- a5$p
