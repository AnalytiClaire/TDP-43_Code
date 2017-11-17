library(biomaRt)

setwd("~/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/PPI_Coexpression/CorrelationValue/RMETHOD/")

iref14 <- read.table("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/iref14_Human_UP_noDup_table_nodash.txt", header = T)
braingenes <- read.csv("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Zhang_BrainCelltype_Markers_braingenes.csv", header = T)
genelist <- readLines("~/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/PPI_Coexpression/CorrelationValue/RMETHOD/NuggetGenes.txt")

mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
attributes <- listAttributes(mart)
mart_back <- getBM(attributes =c("hgnc_symbol", "uniprotswissprot"), filters="hgnc_symbol", values=genelist,  mart=mart)

genelist_Uniprot <- subset(mart_back, !(mart_back$uniprotswissprot == ""))
swiss <- subset(genelist_Uniprot, genelist_Uniprot$hgnc_symbol %in% genelist)
write.csv(swiss, "martback.csv", row.names = F)

###### IDENTIFY MISSING GENES AND FIND UNIPROT CODES FOR THEM. NB SOME GENES MAY NOT BE PROTEIN CODING#####

disgenes <- readLines("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/TDP-43genes.txt")

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/PPI_Coexpression/CorrelationValue/RMETHOD/")
mart_table <- read.csv("martback.csv", header = T)

# for (i in 1:length(disgenes)){
#   mutgene <- disgenes[i]
#   mutgene_back <- getBM(attributes =c("hgnc_symbol", "uniprotswissprot"), filters="hgnc_symbol", values=mutgene,  mart=mart)
#   martplusgene <- rbind(mart_table, mutgene_back[1,])
#   uniprot_gene <- mart_table$uniprotswissprot
#   PPI2 <- subset(iref14, iref14$V1 %in% uniprot_gene & iref14$V2 %in% uniprot_gene)
#   mutgeneprot <- mutgene_back$uniprotswissprot
#   interactions <- subset(PPI, PPI$V1 %in% mutgeneprot | PPI$V2 %in% mutgeneprot)
#   write.csv(interactions, file = paste(mutgene, "_PPI.csv", sep = ""), row.names = F, quote = F)
# }

mutgene <- disgenes
mutgene_back <- getBM(attributes =c("hgnc_symbol", "uniprotswissprot"), filters="hgnc_symbol", values=mutgene,  mart=mart)
martplusgene <- rbind(mart_table, mutgene_back)
genelist_Uniprot <- subset(martplusgene, !(martplusgene$uniprotswissprot == ""))
genelist_Uniprot <- subset(genelist_Uniprot,!(duplicated(genelist_Uniprot$hgnc_symbol)))
uniprot_gene <- genelist_Uniprot$uniprotswissprot
PPI_All <- subset(iref14, iref14$V1 %in% uniprot_gene & iref14$V2 %in% uniprot_gene)
write.csv(PPI_All, "ALL_PPI.csv", row.names = F, quote = F)

#symbol conversion https://biodbnet-abcc.ncifcrf.gov/db/db2db.php

### Group1.1 only
setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/Coexpression/PPI_Coexpression/CorrelationValue/PPI_mutgenes")
mart_table <- read.csv("martback.csv", header = T)
group1.1 <- mart_table$uniprotswissprot
PPI <- subset(iref14, iref14$V1 %in% group1.1 & iref14$V2 %in% group1.1)
write.csv(PPI, "Group1.1_PPI.csv", row.names = F, quote = F)



