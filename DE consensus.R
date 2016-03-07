#####DIFFERENTIAL GENE EXPRESSION INTERSECT
#takes csv files of top X DE genes and identifies any consensus genes 

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GeneExpressionAnalysis/RNA-seq/16_4_1/")
A <- read.csv(file = "PetMar1_ap_2000")
 
B <- read.csv(file = "RavMar1_ap_2000")
 
C <- read.csv(file = "sALS_ap_5000")

D <- read.csv(file = "FTLD_ap_5000")

E <- read.csv(file = "VCP_ap_5000")


A_DE <- A$X
B_DE <- B$X
C_DE <- C$Gene.Symbol
D_DE <- D$Gene.Symbol
E_DE <- E$Gene.Symbol
overlap <- Reduce(intersect, list(A_DE, B_DE, C_DE, D_DE, E_DE))
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

overlap <- Reduce(intersect, list(RNA, Micro))
