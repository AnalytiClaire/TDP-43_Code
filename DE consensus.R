#####DIFFERENTIAL GENE EXPRESSION INTERSECT
#takes csv files of top X DE genes and identifies any consensus genes 

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GeneExpressionAnalysis/Microarray/TopGenes_2016-02-15/")
A <- read.csv(file =   "C9_ap_5000")

B <- read.csv(file =   "CH_ap_5000")
 
C <- read.csv(file = "sALS_ap_5000")
 
D <- read.csv(file = "FTLD_ap_5000")

E <- read.csv(file =  "VCP_ap_5000")

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GeneExpressionAnalysis/RNA-seq/16_3_17/")

F <- read.csv(file =  "PET_ap_5000")

G <- read.csv(file =  "RAV_ap_5000")


A_DE <- A$Gene.Symbol
B_DE <- B$Gene.Symbol
C_DE <- C$Gene.Symbol
D_DE <- D$Gene.Symbol
E_DE <- E$Gene.Symbol
F_DE <- F$X
G_DE <- G$X

overlap <- Reduce(intersect, list(A_DE, B_DE, C_DE, D_DE, E_DE, F_DE, G_DE))
print(overlap)

write.csv(overlap, file = "overlap.csv")

# c9setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43 Data Sets/Pathprint/")
# PP_gene_list <- read.csv(file = "GeneList_2.csv")
# PP_gene_list <- PP_gene_list$GENE
# setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43 Data Sets/DE Genes")
# DE_gene_list <- read.csv(file = "DE_genelist_5000.csv")
# DE_gene_list <- DE_gene_list$GENE
# overlap <- Reduce(intersect, list(PP_gene_list, DE_gene_list))
# print(overlap)

RNA <- overlap
Micro <- overlap

snpenrich <- Reduce(intersect, list(overlap, gwas3))
print(snpenrich)
