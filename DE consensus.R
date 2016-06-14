#####DIFFERENTIAL GENE EXPRESSION INTERSECT
#takes csv files of top X DE genes and identifies any consensus genes 

setwd("/Users/clairegreen/Documents/PhD/TDP-43/non-TDP-43 Data Sets/JK2011_SOD1/")
A <- read.csv(file =   "JKSOD1 _ap_500.csv")

setwd("/Users/clairegreen/Documents/PhD/TDP-43/non-TDP-43 Data Sets/J-Y_FUS/")
B <- read.csv(file =   "FUS _ap_500.csv")
 
# C <- read.csv(file =   "sALS _ap_6500.csv")
#  
# D <- read.csv(file =   "FTLD _ap_6500.csv")
# 
# E <- read.csv(file =   "VCP _ap_6500.csv")
# 
# # setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GeneExpressionAnalysis/RNA-seq/16_3_17/")
# 
# F <- read.csv(file =   "PET _ap_6500.csv")
# 
# G <- read.csv(file =   "RAV _ap_6500.csv")
# 
# H <- read.csv(file =   "FSHD _ap_6500.csv")

A_DE <- A$Gene.Symbol
B_DE <- B$Gene.Symbol
C_DE <- C$Gene.Symbol
D_DE <- D$Gene.Symbol
E_DE <- E$Gene.Symbol
F_DE <- F$hgnc_symbol
G_DE <- G$hgnc_symbol
H_DE <- H$Gene.Symbol

overlap <- Reduce(intersect, list(A_DE, B_DE))
print(overlap)

setwd("/Users/clairegreen/Desktop/")
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
