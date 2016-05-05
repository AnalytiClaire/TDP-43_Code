# set working directory
# set working directory
setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GSEA/PCxN Example/probesets/")

#Load individual gene names for each significance threshold
A <- read.table(file = "threegenes.txt")
a <- A$V1

B <- read.table(file = "fourgenes.txt")
b <- B$V1

C <- read.table(file = "fivegenes.txt")
c <- C$V1

D <- read.table(file = "sixgenes")
d <- D$V1

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GWAS")

E <- read.table(file = "signif.snp.NeuroX.txt")
e <- E$V1

F <- read.table(file = "signif.snp.NeuroX.p5E08.txt")
f <- F$V1

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GeneExpressionAnalysis")

G <- read.table(file = "allgenes.txt")
g <- G$V1

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GWAS")
H <- read.table(file = "signif.snp.AD.GWASCentralp5E08.txt")
h <- H$V1

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression")
K <- read.table(file = "subnet.28.GM.txt")
k <- K$V1

####Load full gwas datasets ####
setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GWAS/")

GCEN <- read.csv("ALS.gwascentral.martquery_0301121419_683.csv")

three <- GCEN[(GCEN$p.value <= 0.001),]

four <- GCEN[(GCEN$p.value <= 0.0001),]

five <- GCEN[(GCEN$p.value <= 0.00001),]

six <- GCEN[(GCEN$p.value <= 0.00001),]




setwd (dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GSEA/PCxN Example")

Z <- read.csv(file = "Pathprint/PP+20.csv", na.strings = c("", "NA)"))
Z <- as.list(Z)


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
sym.genes <- t(sym.genes)

allgenes <- sym.genes[!duplicated(sym.genes),]
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

y <- Z$ABC.transporters..KEGG.

#How many PCxN genes contain snps
x.in <- length (which(k %in% y)) 
#how many do not
x.out <- length(y) - x.in
#total number of snps
tot.in <- length(k)
#total number of all genes
tot.out <- length(allgenes)-length(tot.in)


#create count matrix
counts <- matrix (nrow=2, ncol=2)
counts [1,] <- c(x.in, tot.in)
counts [2,] <- c(x.out, tot.out)

#Conduct fisher's exact test for count data
a5 <-fisher.test (counts)
enrich <- a5$p
print(enrich)

