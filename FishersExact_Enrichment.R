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

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/M&R/")
K <- read.table(file = "subnet.28.GM.txt")
k <- K$V1

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/M&R/ExAC/")
L <- read.table(file = "exac.pli.0.95.txt")
l <- L$V1


####Load full gwas datasets ####
setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/M&R/")

GCEN <- read.csv("ALS.gwascentral.martquery_0301121419_683.csv")

three <- GCEN[(GCEN$p.value <= 0.001),]

four <- GCEN[(GCEN$p.value <= 0.0001),]

five <- GCEN[(GCEN$p.value <= 0.00001),]

six <- GCEN[(GCEN$p.value <= 0.00001),]








setwd (dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/Pathprint/Pathprint 25.04.16/FishersExact/FE.All.Pathways/FE.PathprintPathways(29)/")

#Read in geneset
Z <- read.csv(file = "FE.pathprintgenes.csv", na.strings = c("", "NA)"))
Z <- as.list(Z)
Z<- lapply(Z, function(x) x[!is.na(x)])

geneset <- Z
genelist <- k
genelist.length <- length(h)

hyper <- as.data.frame(matrix(nrow = length(geneset), ncol = 1))
rownames(hyper) <- names(geneset)
colnames(hyper) <- c("p-value")

#Conduct fisher's exact test

for (i in 1:length(geneset)) {
  if (length(intersect(genelist, unlist(geneset[i]))) < 
      1) {
    hyper[i, 1] <- 1
  }
  else if (length(intersect(genelist, unlist(geneset[i]))) > 
           0) {
    
    #Pathway gene/GWAS genes intersection
    x.in <- length(intersect(genelist, unlist(geneset[i])))
    #Remaining pathway genes
    x.out <- length(unlist(geneset[i])) - x.in
    #total number of snps
    tot.in <- genelist.length
    #total number of all genes
    tot.out <- length(allgenes)-length(tot.in)
    
    #create count matrix
    counts <- matrix (nrow=2, ncol=2)
    counts [1,] <- c(x.in, tot.in)
    counts [2,] <- c(x.out, tot.out)
    
    #Conduct fisher's exact test for count data
    a5 <-fisher.test (counts)
    hyper[i, 1] <- a5$p
  }
}

hyper[, 2] <- p.adjust(hyper[, 1], method = "BH")
overlap <- vector("list", 0)
for (i in 1:length(geneset)) {
  temp.overlap <- list(intersect(genelist, unlist(geneset[[i]])))
  overlap <- append(overlap, temp.overlap)
}
names(overlap) <- rownames(hyper)
for (i in 1:length(geneset)) {
  hyper[i, 3] <- length(overlap[[i]])
  hyper[i, 4] <- length(geneset[[i]])
}
hyper[, 5] <- rownames(hyper)
hyper <- cbind((1:length(hyper[, 1])), hyper)
colnames(hyper) <- c("ID", "P-value", "BHadjP-value", "nGenes", 
                     "nPathway", "Name")

write.csv(hyper, file = "FE.EP20.subnetwork28.csv")


o1 <- Reduce(intersect, list(k, Z$Prion.diseases..KEGG.))
print(o1)
