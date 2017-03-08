# Example script for hyperPathway

library(pathprint)

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

setwd(dir = "/Users/clairegreen/Desktop/")

G <- read.table(file = "TDP-43DEGs.txt")
g <- G$V1

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GWAS")
H <- read.table(file = "signif.snp.AD.GWASCentralp5E08.txt")
h <- H$V1

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/M&R/")
K <- read.table(file = "subnet.28.GM.txt")
k <- K$V1

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/M&R/ExAC/")
L <- read.table(file = "exac.pli.0.95.txt")
l <- L$V1

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/M&R/")
M <- read.table(file = "Cirulli.txt")
m <- M$V1

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/M&R/")
N <- read.table(file = "GeneCardsAD.txt")
n <- N$V1

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/M&R/")
O <- read.table(file = "GeneCardsALS.txt")
o <- O$V1

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/QQ/Test3/")
P <- read.table(file = "genemania-genes.txt")
p <- P$V1

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/M&R/")
Q <- read.table(file = "Pasterkamp_TDP43.txt")
q <- Q$V1

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/M&R/")
R <- read.table(file = "Taylor_TDP43.txt")
r <- R$V1

setwd(dir = "/Users/clairegreen/Desktop/")
S <- read.table(file = "ProteintargetingtoER.txt")
s <- S$V1

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/M&R/6500")
U <- read.table("6500 + 200.txt")
u <- U$V1

setwd(dir = "/Users/clairegreen/Desktop/Desktop Folder/")
V <- read.table("MicrogliaGenes_doi:10.1038:nn.3599.txt")
v <- V$V1

# 
# ####Load full gwas datasets ####
# setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GWAS/")
# 
# GCEN <- read.csv("ALS.gwascentral.martquery_0301121419_683.csv")
# 
# three <- GCEN[(GCEN$p.value <= 0.001),]
# 
# four <- GCEN[(GCEN$p.value <= 0.0001),]
# 
# five <- GCEN[(GCEN$p.value <= 0.00001),]
# 
# six <- GCEN[(GCEN$p.value <= 0.00001),]

# load script and geneset file
# The geneset is a list of pathways including
## Reactome high level processes
## Pathways from KEGG and Wikipathways
## Netpath transcriptionally regulated gene sets
## Lincoln Stein's static modules
# format is Entrez Gene

setwd (dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/Pathprint/Pathprint 25.04.16/hyperPathway/All.Pathways/PathprintPathways(29)/")

Z <- read.csv(file = "pathprintgenes.csv", na.strings = c("", "NA)"))
Z <- as.list(Z)
Z <- lapply(Z, function(x) x[!is.na(x)])

z <- read.csv(file = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/Pathprint/PP+20.csv", na.strings = c("", "NA)"))
z <- as.list(z)
z <- lapply(z, function(x) x[!is.na(x)])

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

#Read in geneset
setwd (dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/DEG_Test2/")
Y <- read.csv(file = "AllgenesNO.csv", na.strings = c("", "NA)"))
Y <- as.list(Y)
Y<- lapply(Y, function(x) x[!is.na(x)])

setwd (dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/DEG_Test2/")
X <- read.csv(file = "Allgenes+GM.csv", na.strings = c("", "NA)"))
X <- as.list(X)
X<- lapply(X, function(x) x[!is.na(x)])

setwd (dir = "/Users/clairegreen/Desktop/")
y <- read.table(file = "GM+58.txt")
y <- y$V1

setwd (dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/")
W <- read.csv(file = "BenchmarkGenes.csv", na.strings = c("", "NA)"))
W <- as.list(W)
W<- lapply(W, function(x) x[!is.na(x)])


# run script
pathwayEnrichment2 <- hyperPathway(
                genelist = x,
								geneset = Z,
								Nchip = length(allgenes)
							 )
setwd (dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/TDP-43_DEseq2/")
write.csv(pathwayEnrichment2, file = "pathprint_enrich.csv")


o2 <- Reduce(intersect, list(x, W$ALSOD))
print(o2)
write.csv(o2, file = "intersect.csv")


