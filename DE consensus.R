#####DIFFERENTIAL GENE EXPRESSION INTERSECT
#takes csv files of top X DE genes and identifies any consensus genes 

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GeneExpressionAnalysis/TopGenes_2016-02-12/")
C9orf72 <- read.csv(file = "C9_ap_4000")
 
CHMP2B <- read.csv(file = "CH_ap_4000")
 
sALS<- read.csv(file = "sALS_ap_4000")

FTLD<- read.csv(file = "FTLD_ap_4000")

VCP<- read.csv(file = "VCP_ap_4000")


C9_DE<- C9orf72$Gene.Symbol
CHMP2B_DE <- CHMP2B$Gene.Symbol
sALS_DE <- sALS$Gene.Symbol
FTLD_DE <- FTLD$Gene.Symbol
VCP_DE <- VCP$Gene.Symbol
overlap <- Reduce(intersect, list(C9_DE, CHMP2B_DE, sALS_DE, FTLD_DE, VCP_DE))
print(overlap)


# c9setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43 Data Sets/Pathprint/")
# PP_gene_list <- read.csv(file = "GeneList_2.csv")
# PP_gene_list <- PP_gene_list$GENE
# setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43 Data Sets/DE Genes")
# DE_gene_list <- read.csv(file = "DE_genelist_5000.csv")
# DE_gene_list <- DE_gene_list$GENE
# overlap <- Reduce(intersect, list(PP_gene_list, DE_gene_list))
# print(overlap)
