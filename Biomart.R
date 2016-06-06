##Convert information using biomart

#Important terms:
# Affymetrix Human Genome U133 Plus 2.0 array probe set = affy_hg_u133_plus_2
# Affymetrix Human Genome U133A 2 array probe set = affy_hg_u133a_2
# Ensembl ID = ensembl_gene_ensembl
# HGNC symbol = hgnc_symbol
# Entrez ID = entrezgene

genes <- read.csv("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/WGCNA/GoodData/FTLDSymbols.csv", header = TRUE)

library (biomaRt)
mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")

analysis.name <- "FTLDID"
genes <- genes$Probe.Set.ID
mart_attribute <- listAttributes(mart)
mart_filter <- listFilters(mart)
annotation <- getBM(attributes=c("affy_hg_u133_plus_2","ensembl_gene_id", "hgnc_symbol"),
                   filters = "affy_hg_u133_plus_2", values = genes, mart = mart)

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/WGCNA/GoodData/")
write.csv(annotation, paste(analysis.name,"geneIDs.csv"), sep = "")




C9file <- read.csv("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/DEG Test2/C9rankeduniqueresult.csv")
C9 <- merge(C9file, annotation, by.x = "Ensembl", by.y = "ensembl_gene_id")
