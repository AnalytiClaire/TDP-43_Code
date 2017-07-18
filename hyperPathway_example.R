# Example script for hyperPathway

library(pathprint)

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

setwd(dir = "/Users/clairegreen/Desktop/")

G <- read.table(file = "TDP-43DEGs.txt")
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

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/M&R/")
M <- read.table(file = "Cirulli.txt")
m <- M$V1

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/M&R/")
N <- read.table(file = "GeneCardsAD.txt")
n <- N$V1

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/M&R/")
O <- read.table(file = "GeneCardsALS.txt")
o <- O$V1

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/QQ/Test3/")
P <- read.table(file = "genemania-genes.txt")
p <- P$V1

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/M&R/")
Q <- read.table(file = "Pasterkamp_TDP43.txt")
q <- Q$V1

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/M&R/")
R <- read.table(file = "Taylor_TDP43.txt")
r <- R$V1

setwd(dir = "/Users/clairegreen/Desktop/")
S <- read.table(file = "ProteintargetingtoER.txt")
s <- S$V1

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/M&R/6500")
U <- read.table("6500 + 200.txt")
u <- U$V1

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression")
V <- read.table("OneBenchmarkList.txt")
v <- V$V1

setwd(dir = "/Users/clairegreen/Desktop/")
List5996 <- read.table("List5996.txt")
List5996 <- List5996$V1

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression")
List6732 <- read.table("LIST_6732.txt")
List6732 <- List6732$V1

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/FoldChangeResults/")
up <- read.table("intersect_up_1.txt")
up <- up$V1

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/FoldChangeResults/")
down <- read.table("intersect_down_1.txt")
down <- down$V1

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/FoldChangeResults/")
up_filter <- read.table("Filtered_up_genes.txt")
up_filter <- up_filter$V1

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/FoldChangeResults/")
down_filter <-read.table("Filtered_down_genes.txt")
down_filter <- down_filter$V1


# ####Load full gwas datasets ####
# setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GWAS/")
# 
# GCEN <- read.csv("ALS.gwascentral.martquery_0301121419_683.csv")
# 
# three <- GCEN[(GCEN$p.value <= 0.001),]
# 
# four <- GCEN[(GCEN$p.value <= 0.0001),]
# 
# five <- GCEN[(GCEN$p.value <= 0.00001),]
# 
# six <- GCEN[(GCEN$p.value <= 0.00001),]

# load script and geneset file
# The geneset is a list of pathways including
## Reactome high level processes
## Pathways from KEGG and Wikipathways
## Netpath transcriptionally regulated gene sets
## Lincoln Stein's static modules
# format is Entrez Gene

setwd (dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/Pathprint/Pathprint 25.04.16/hyperPathway/All.Pathways/PathprintPathways(29)/")

Z <- read.csv(file = "pathprintgenes.csv", na.strings = c("", "NA)"))
Z <- as.list(Z)
Z <- lapply(Z, function(x) x[!is.na(x)])

z <- read.csv(file = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/Pathprint/PathprintThreshold/GeneList_200.csv", na.strings = c("", "NA)"))
z <- as.list(z)
z <- lapply(z, function(x) x[!is.na(x)])


#Pathprint gene lists
load("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/Pathprint/Pathprint_biomart.RData")

# #pathprint result dataframe
# 
# IL2_down_reg._targets <- pathprint_list$`IL-3 down reg. targets (Netpath)`
# Mitochondrial_LCFatty_Acid_BetaOxidation <- pathprint_list$`Mitochondrial LC-Fatty Acid Beta-Oxidation (Wikipathways)`
# Vitamin_B12_Metabolism <- pathprint_list$`Vitamin B12 Metabolism (Wikipathways)`
# Vitamin_D_synthesis <- pathprint_list$`Vitamin D synthesis (Wikipathways)`
# BRCA1_28 <- pathprint_list$`{BRCA1,28} (Static Module)`
# CD4_14 <- pathprint_list$`{CD4,14} (Static Module)`
# SEPT2_21 <- pathprint_list$`{SEPT2,21} (Static Module)`
# Fatty_Acid_Biosynthesis <- pathprint_list$`Fatty Acid Biosynthesis (Wikipathways)`
# FCGR2B_50 <- pathprint_list$`{FCGR2B,50} (Static Module)`
# RB1_11 <- pathprint_list$`{RB1,11} (Static Module)`
# Muscle_contraction <- pathprint_list$`Muscle contraction (Reactome)`
# Steroid_Biosynthesis <- pathprint_list$`Steroid Biosynthesis (Wikipathways)`
# Glutathione_metabolism <- pathprint_list$`Glutathione metabolism (KEGG)`
# Fat_digestion_and_absorption <- pathprint_list$`Fat digestion and absorption (KEGG)`
# AminoacyltRNA_biosynthesis <- pathprint_list$`Aminoacyl-tRNA biosynthesis (KEGG)`
# VCP_17 <- pathprint_list$`{VCP,17} (Static Module)`
# Osteoclast_differentiation <- pathprint_list$`Osteoclast differentiation (KEGG)`
# Fluoropyrimidine_Activity <- pathprint_list$`Fluoropyrimidine Activity (Wikipathways)`
# IL7_down_reg._targets <- pathprint_list$`IL-7 down reg. targets (Netpath)`
# Heme_Biosynthesis <- pathprint_list$`Heme Biosynthesis (Wikipathways)`
# SP1_88 <- pathprint_list$`{SP1,88} (Static Module)`
# Signaling_by_Notch <- pathprint_list$`Signaling by Notch (Reactome)`
# Ubiquinone_and_other_terpenoidquinone_biosynthesis <- pathprint_list$`Ubiquinone and other terpenoid-quinone biosynthesis (KEGG)`
# Sulfur_metabolism <- pathprint_list$`Sulfur metabolism (KEGG)`
# Caffeine_metabolism <- pathprint_list$`Caffeine metabolism (KEGG)`
# Phototransduction <- pathprint_list$`Phototransduction (KEGG)`
# Parkinsons_disease <- pathprint_list$`Parkinson's disease (KEGG)`
# PPP2CA_20 <- pathprint_list$`{PPP2CA,20} (Static Module)`
# ABL1_15 <- pathprint_list$`{ABL1,15} (Static Module)`
# Type_II_diabetes_mellitus <- pathprint_list$`Type II diabetes mellitus (KEGG)`
# Signaling_by_EGFR <- pathprint_list$`Signaling by EGFR (Reactome)`
# Glucocorticoid_Mineralcorticoid_Metabolism <- pathprint_list$`Glucocorticoid &amp; Mineralcorticoid Metabolism (Wikipathways)`
# Alzheimers_disease <- pathprint_list$`Alzheimer's disease (KEGG)`
# Blood_Clotting_Cascade <- pathprint_list$`Blood Clotting Cascade (Wikipathways)`
# Terpenoid_backbone_biosynthesis <- pathprint_list$`Terpenoid backbone biosynthesis (KEGG)`
# Signaling_by_GPCR <- pathprint_list$`Signaling by GPCR (Reactome)`
# 
# PP_up <- list(IL2_down_reg._targets,Mitochondrial_LCFatty_Acid_BetaOxidation,Vitamin_B12_Metabolism,Vitamin_D_synthesis,BRCA1_28,CD4_14,
#       SEPT2_21,Fatty_Acid_Biosynthesis,FCGR2B_50,RB1_11,Muscle_contraction,Steroid_Biosynthesis,Glutathione_metabolism,Fat_digestion_and_absorption,
#       AminoacyltRNA_biosynthesis,VCP_17,Osteoclast_differentiation,Fluoropyrimidine_Activity,IL7_down_reg._targets,Heme_Biosynthesis,SP1_88,
#       Signaling_by_Notch)
# PP_down <- list(Ubiquinone_and_other_terpenoidquinone_biosynthesis,Sulfur_metabolism,Caffeine_metabolism,Phototransduction,Parkinsons_disease,
#       PPP2CA_20,ABL1_15,Type_II_diabetes_mellitus,Signaling_by_EGFR,Glucocorticoid_Mineralcorticoid_Metabolism,Alzheimers_disease,
#       Blood_Clotting_Cascade,Terpenoid_backbone_biosynthesis,Signaling_by_GPCR)
# 
# names(PP_up) <- names


#Load file with all genes
library(hgu133plus2.db)
sym <- hgu133plus2SYMBOL
sym1 <- mappedkeys(sym)
sym2 <- as.list (sym[c(sym1)]) 
sym3 <- data.frame (sym2)
sym.probes <- names (sym2)
sym.genes <- sym3[1,]

sym.genes <- t(sym.genes)

allgenes <- sym.genes[!duplicated(sym.genes),]

#Read in geneset
setwd (dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/DEG_Test2/")
Y <- read.csv(file = "AllgenesNO.csv", na.strings = c("", "NA)"))
Y <- as.list(Y)
Y<- lapply(Y, function(x) x[!is.na(x)])

setwd (dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/LF_miRNA/")
X <- read.csv(file = "miRNAtargetGenes.csv", na.strings = c("", "NA)"))
X <- as.list(X)
X<- lapply(X, function(x) x[!is.na(x)])

setwd (dir = "/Users/clairegreen/Desktop/")
y <- read.table(file = "GM+58.txt")
y <- y$V1

setwd (dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/")
W <- read.csv(file = "BenchmarkGenes2.csv", na.strings = c("", "NA)"))
W <- as.list(W)
W<- lapply(W, function(x) x[!is.na(x)])


# run script
pathwayEnrichment <- hyperPathway(
                genelist = up_filter,
								geneset = X,
								Nchip = length(allgenes)
							 )
setwd (dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/FoldChangeResults/")
write.csv(pathwayEnrichment, file = "pathprint_enrich_up.csv")


o2 <- Reduce(intersect, list(up_filter, Pathprint_results$VCP_17))
print(o2)
write.csv(o2, file = "intersect.csv")


