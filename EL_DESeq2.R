#### Analysis of Eddie Lee's TDP-43 data ####

# TDPneg = nuclei absent of TDP-43
# TDPpos = nuclei containing TDP-43

library(DESeq2)
library(biomaRt)

setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/EL_TDP-43/")
load("Count_Matrix_Genes.rda")


#Assign condition
condition <- c("Neg", "Pos","Neg", "Pos","Neg","Pos","Neg", "Pos","Neg", "Pos",
                  "Neg", "Pos","Neg", "Pos")


Counts <- counts.g
#Rows must have at least 3 samples with scores of 10 or higher 
keep <- rowSums(Counts>=10) >= 3
Counts<-Counts[keep,]

#Create a coldata frame
coldata <- data.frame(row.names=colnames(counts.g), condition)
coldata

#Create a DESeqDataSet object 
dds <- DESeqDataSetFromMatrix(countData=counts.g, colData=coldata, design=~condition)
dds

#Run DEseq2 pipleline
dds <- DESeq(dds)

res <- results(dds)
table(res$padj<0.05)
## Order by p-value
res <- res[order(res$pvalue), ]
## Merge with normalized count data
resdata_raw <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)

names(resdata_raw)[1] <- "Gene"
head(resdata_raw)
genes <- as.vector(resdata_raw$Gene)
mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
mart_back <- getBM(attributes =c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values=genes,  mart=mart)
result <- merge(resdata_raw, mart_back, by.x = "Gene", by.y = "ensembl_gene_id")
genesort <- result[order(result$pvalue),]

write.csv(genesort, "EL_results.csv")

Sig.padj <- subset(genesort, subset=(padj < 0.05))
Sig.padj <- Sig.padj$hgnc_symbol
Sig.padj <- Sig.padj[!duplicated(Sig.padj)]

write.table(Sig.padj, "sig.padj.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
