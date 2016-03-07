##RNA-Seq Gene Expression Analysis using Limma##

analysis.name<-"PET" #Label analysis
setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GeneExpressionAnalysis/RNA-seq/Petrucelli/")
# Counts <- read.table(file = 'GSE67196_Petrucelli2015_ALS_genes.rawcount.txt', header = TRUE)
# 
# write.csv(x = Counts, file = "counts_petrucelli.csv")

Counts <- read.csv(file = "FCX_Petrucelli.csv", header = TRUE)

Counts[Counts == 0] <- NA
# Counts[Counts<30] <- NA
Counts <- na.omit(Counts)
# Counts<-subset(Counts, subset=(GeneID !="NA")) #if no gene symbol, discount

# Countszero <-subset(Counts, subset=(row !=0))
# Countszero <- apply(Counts, 1, function(row) all(row !="NA"))
# Counts <- Counts[Countszero,]

library(limma)
library(edgeR)

# rownames(Counts)<-Counts[,1]

# Counts[,1] <- NULL
Countnum <- Counts[,2:28]
# Counts <- data.matrix(Counts)

#DGElist
dge <- DGEList(counts=Countnum)
dge <- calcNormFactors(dge)

#Design
Treat<-factor(rep(c("Control", "Patient"),c(9,18)), levels=c("Control", "Patient"))
design<-model.matrix(~Treat)
rownames(design)<-colnames(Countnum)
design

#Voom transformation
v <- voom(dge,design,plot=TRUE)

#Limma fitting
fit <- lmFit(v,design)
fit <- eBayes(fit)
result<-topTable(fit, coef="TreatPatient", adjust="BH", number=nrow(Countnum)) #"BH" adjust for multiple hypothesis testing
result <- merge(result, Counts, by="row.names", all=TRUE)
result <- result[,2:8]

uniqueresult <- result[!duplicated(result[,7]),]
genesort <- uniqueresult[order(uniqueresult$adj.P.Val),]

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GeneExpressionAnalysis/RNA-seq/16_3_2a/")
write.csv(genesort, file=paste(analysis.name, "rankeduniqueresult.csv", sep=""), sep="\t", row.names=TRUE, quote = FALSE)

topgene <- genesort[1:1000,]
write.csv(x = topgene, file = paste(analysis.name,"_ap_1000", sep = ""))
topgene <- genesort[1:2000,]
write.csv(x = topgene, file = paste(analysis.name,"_ap_2000", sep = ""))
topgene <- genesort[1:3000,]
write.csv(x = topgene, file = paste(analysis.name,"_ap_3000", sep = ""))
topgene <- genesort[1:4000,]
write.csv(x = topgene, file = paste(analysis.name,"_ap_4000", sep = ""))
topgene <- genesort[1:5000,]
write.csv(x = topgene, file = paste(analysis.name,"_ap_5000", sep = ""))
