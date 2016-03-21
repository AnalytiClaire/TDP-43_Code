#####DIFFERENTIAL GENE EXPRESSION INTERSECT
#takes csv files of top X DE genes and identifies any consensus genes 

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GeneExpressionAnalysis/Microarray/TopGenes_ensemblID/")
A <- read.csv(file =   "C9b_ap_5000.csv")

B <- read.csv(file =   "CHb_ap_5000.csv")
 
C <- read.csv(file = "sALSb_ap_5000.csv")
 
D <- read.csv(file = "FTLDb_ap_5000.csv")

E <- read.csv(file =  "VCPb_ap_5000.csv")

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GeneExpressionAnalysis/RNA-seq/16_3_21/")

F <- read.csv(file =  "PETb_ap_5000")

G <- read.csv(file =  "RAVb_ap_5000")


A_DE <- A$Ensembl
B_DE <- B$Ensembl
C_DE <- C$Ensembl
D_DE <- D$Ensembl
E_DE <- E$Ensembl
F_DE <- F$Row.names
G_DE <- G$Row.names

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

overlap <- Reduce(intersect, list(RNA, Micro))
