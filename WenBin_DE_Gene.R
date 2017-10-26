##Differential Expression of Genes##

library(affy)
library(tcltk)
library(widgetTools)
library(DynDoc)
library(tools)
library(Biobase)
library(tkWidgets)
library(plyr)

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GeneExpressionAnalysis/Microarray/FTLD/")

#run program to choose .CEL files from directory
celfiles <- fileBrowser(textToShow = "Choose CEL files", testFun = hasSuffix("[cC][eE][lL]"))
#celfiles<-basename(celfiles)
Data<-ReadAffy(filenames=celfiles) #read in files
rmaEset<-rma(Data) #normalise using RMA
analysis.name<-"FTLD_SP" #Label analysis
dataMatrixAll<-exprs(rmaEset) #takes expression from normalised expression set

# setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/NormalisedExpressionMatrices/")
# write.csv(dataMatrixAll, file = paste(analysis.name,"eset.csv", sep = ""))

#mas5call generates presence/absence calls for each probeset
mas5call<-mas5calls(Data)
callMatrixAll<-exprs(mas5call)
colnames(callMatrixAll)<-sub(".CEL", ".mas5-Detection", colnames(callMatrixAll),fixed=TRUE)
colnames(callMatrixAll)<-sub(".cel", ".mas5-Detection", colnames(callMatrixAll),fixed=TRUE)
callMatrixAll<-as.data.frame(callMatrixAll)
callMatrixAll$ProbeSetID<-rownames(callMatrixAll)
countPf<-function(x){
  sum(x=="P")
}

#count how many samples have presence calls
countPl<-apply(callMatrixAll, 1, countPf)
callMatrixAll$ProbeSetID<-rownames(callMatrixAll)
countPdf<-data.frame(ProbeSetID=names(countPl), countP=countPl) 

#read annotation file

###USING BIOMART
# library (biomaRt)
# mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org") 
# x <- rownames(dataMatrixAll) #create vector containing probe IDs
# mart_attribute <- listAttributes(mart)
# annotation <- getBM(attributes=c("affy_hg_u133a_2", "hgnc_symbol", "description"), 
#                    filters = "affy_hg_u133a_2", values = x, mart = mart)
# annotation<-subset(annotation, subset=(hgnc_symbol !="")) #if no gene symbol, discount

# USING ANNOTATION FILE (if .csv, convert to .txt using excel)
# annotation.file<-"/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/HG-U133_Plus_2.na35.annot.csv/HG-U133_Plus_2.na35_SHORT.annot.txt"
annotation.file<-"/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/HG-U133A_2.na35.annot.csv/HG-U133A_2.na35.annot_SHORT.txt"
#annotation.file<-"/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/HG-U133A.na35.annot.csv/HG-U133A.na35.annot.txt"
annotation<-read.table(annotation.file, header = TRUE, row.names=NULL, sep="\t", skip=0, stringsAsFactors=F, quote = "", comment.char="!", fill = TRUE, as.is = TRUE)
dim(annotation)
nrow(annotation)
annotation<-subset( annotation, subset=(Gene.Symbol !="---")) #if no gene symbol, discount

# Remove rows in which genes are noted to have negative strand matching probes
idxNegativeStrand<-grep("Negative Strand Matching Probes", annotation$Annotation.Notes)
if(length(idxNegativeStrand)>0)
{
  annotation<-annotation[-idxNegativeStrand,]
}


expressionMatrix<-exprs(rmaEset)
colnames(expressionMatrix)

#this is for matched samples
Treat<-factor(rep(c("Control", "Patient"),c(8,10)), levels=c("Control", "Patient"))
design<-model.matrix(~Treat)
rownames(design)<-colnames(expressionMatrix)
design

#Conduct statistical analysis of expression
library(limma)
fit<-lmFit(expressionMatrix, design) #linear model fit
fit<-eBayes(fit) 
result<-topTable(fit, coef="TreatPatient", adjust="BH", number=nrow(expressionMatrix)) #"BH" adjust for multiple hypothesis testing
#toptable normally takes top number but this takes all


result$"ProbeSetID"<-rownames(result) #make probeset IDs the row names
head(result$"ProbeSetID") 
result$"Fold Change"<-2^result$logFC 
result$"Fold Change"[result$"Fold Change"<1]<-(-1)/result$"Fold Change"[result$"Fold Change"<1] #converts log fold change into a linear value above or below 0
expressionLinear<-as.data.frame(2^expressionMatrix)
expressionLinear$ProbeSetID<-rownames(expressionLinear)
result<-merge(result, expressionLinear, by.x="ProbeSetID", by.y="ProbeSetID") #merge values into one array
result<-merge(annotation, result, by.x="Probe.Set.ID", by.y="ProbeSetID")
result<-merge(result, countPdf, by.x="Probe.Set.ID", by.y="ProbeSetID")

# setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/PPI_Network/SimilarityNetworkFusion/Microarray_RNAseq/")
# write.csv(result, file=paste(analysis.name, "result.csv", sep=""), sep="\t", row.names=FALSE, quote = FALSE)

# result<-subset(result, Gene.Symbol!="") #removes any probes for which there are no gene symbols
# result<-subset(result, subset=(countP>2)) #only takes results that have at least 2 samples with a presence call for a probe


genesort <- result[order(result$P.Value),]

genesort$Gene.Symbol <- sapply(strsplit(genesort$Gene.Symbol,"///"), `[`, 1)

genesort$Ensembl <- sapply(strsplit(genesort$Ensembl,"///"), `[`, 1)

# genesort <- result[order(result$adj.P.Val),]
uniqueresult <- genesort[!duplicated(genesort[,5]),]

unique_patient <- uniqueresult[,c(5,20:37)]

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/Microarray_AllGenesExpression/")
write.csv(genesort, file=paste(analysis.name, "rankeduniqueresult.csv", sep=""), sep="\t", row.names=FALSE, quote = FALSE)
write.csv(unique_patient, file=paste(analysis.name, "uniquegene_samples.csv", sep=""), sep="\t", row.names=FALSE, quote = FALSE)




# foldchange<-1.1
# pvalue<-0.05
# #adj_P_Val<-0.05
# siggenes<-subset(uniqueresult, subset=(P.Value < pvalue) & abs(logFC) > log2(foldchange))
# #siggenes<-subset(result, subset=(adj.P.Val < 0.05))
# nrow(siggenes)
# siggenesup<-subset(siggenes, subset= logFC > 0)
# siggenesdown<-subset(siggenes, subset=logFC < 0)
# colnames(siggenesup)
# nrow(siggenesup)
# nrow(siggenesdown)
# UpandDown<-intersect(siggenesup$"Gene.Symbol", siggenesdown$"Gene.Symbol")
# length(UpandDown)
# 
# UporDown<-subset(siggenes, subset=(!siggenes$"Gene.Symbol"%in% UpandDown))
# upsiggenes<-subset(siggenesup, subset=(!siggenesup$"Gene.Symbol"%in% UpandDown))
# downsiggenes<-subset(siggenesdown, subset=(!siggenesdown$"Gene.Symbol"%in% UpandDown))
# length(unique(siggenes$"Gene.Symbol"))
# uniquesiggenes <- unique(siggenes$Gene.Symbol)
# length(unique(upsiggenes$"Gene.Symbol"))
# length(unique(downsiggenes$"Gene.Symbol"))
# length(unique(UporDown$"Gene.Symbol"))


###Write results to CSV files for consensus analysis
setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/noMedian/")
#dir.create(paste("TopGenes", Sys.Date(), sep = "_")) #create directory using the day's date
#Take results, remove duplicate rows for genes, order by adjusted p value and take top X number of genes


#For ordering by adjusted p value
genesort <- uniqueresult[order(uniqueresult$adj.P.Val),]
# uniqueresult <- result[!duplicated(result[,15]),]
write.csv(uniqueresult, file=paste(analysis.name, "rankeduniqueresult.csv", sep=""), sep="\t", row.names=FALSE, quote = FALSE)

