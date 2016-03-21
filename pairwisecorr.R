setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GeneExpressionAnalysis/Microarray/TopGenes_2016-02-15")

#load dataset
exprsC9 <- read.csv("C9rankeduniqueresult.csv")

#Load list of interesting genes
setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GeneExpressionAnalysis/Microarray/")
Genelist <- read.csv("MADEGs.csv")

#Make gene symbol row names
rownames(exprsC9) <- exprsC9$Gene.Symbol
exprsC9pat <- exprsC9[,52:59]

#Make gene symbol a column
exprsC9pat <- cbind(exprsC9pat, exprsC9$Gene.Symbol)
colnames(exprsC9pat)[9] <- "Gene.Symbol"

#Merge by interesting gene names with expression to form matrix
C9patgene <- merge(Genelist, exprsC9pat, by.x = "Gene", by.y = "Gene.Symbol")
rownames(C9patgene) <- C9patgene$Gene
C9patgene[,1] <- NULL


##For loop for generating regression values and p values
C9patgene <- t(C9patgene)

reg <- matrix(0, ncol(C9patgene), ncol(C9patgene))
p.value <- matrix(0, ncol(C9patgene), ncol(C9patgene))

for (i in 1:ncol(C9patgene)){
  for (j in 1:ncol(C9patgene)){
    reg[i,j] <- cor.test(C9patgene[,i], C9patgene[,j], method = "kendall")$estimate
    rownames(reg) <- colnames(reg) <- colnames(C9patgene)
  }}

for (i in 1:ncol(C9patgene)){
  for (j in 1:ncol(C9patgene)){
    p.value[i,j] <- cor.test(C9patgene[,i], C9patgene[,j], method = "kendall")$p.value
  rownames(p.value) <- colnames(p.value) <- colnames(C9patgene)
}}

##Only take upper triangle without diagonal (all comparisons are currently doubled)
ptri <- p.value
ptri[lower.tri(ptri, diag = TRUE)] <- NA

#Turn into vector
library(gdata)
p.vec <- unmatrix(ptri)
#Remove NA values
p.vec <- na.omit(p.vec)
#Multiple hypothesis testing correction 
p.adj <- p.adjust(p.vec, method = "fdr", n = length(p.vec))

#Create results table
reg.mat <- unmatrix(reg)
reg.mat <- as.data.frame(reg.mat)
p.adj <- as.data.frame(p.adj)
p.mat <- as.data.frame(p.vec)

pvals <- merge(p.adj, p.mat, by.x = "row.names", by.y = "row.names")
rownames(pvals)<- pvals$Row.names
pvals[,1] <- NULL
results <- merge(pvals, reg.mat, by.x = "row.names", by.y = "row.names")
rownames(results)<- results$Row.names
results[,1] <- NULL
results <- results[order(results$p.vec),]

##PSYCH METHOD### 
library(psych)

#Generate matrix of correlation and p values
cortest <- corr.test(t(C9patgene), use = "pairwise", method = "spearman", adjust = "fdr") #if genes are rownames, the t is important

#Extract results table
# cortestoutput <- cortest$ci
cortestpadjust <- cortest$p

coradj <- cortestpadjust
coradj[lower.tri(coradj, diag = TRUE)] <- NA

#Turn into vector
library(gdata)
coradj <- as.matrix(coradj)
corvec <- unmatrix(coradj)
#Remove NA values
corvec <- na.omit(corvec)
corvec <- as.data.frame(corvec)

#Select significant results
sigoutput <- subset(cortestoutput, cortestoutput$p < 0.05)
write.csv(sigoutput, file = "C9allsigcoexpr1.csv")


