###Random permutations for Differential Expression##

#load the total number of unique genes for each data set
A <- read.csv("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GeneExpressionAnalysis/Microarray/TopGenes_2016-02-15/C9rankeduniqueresult.csv")
B <- read.csv("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GeneExpressionAnalysis/Microarray/TopGenes_2016-02-15/CHrankeduniqueresult.csv")
C <- read.csv("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GeneExpressionAnalysis/Microarray/TopGenes_2016-02-15/sALSrankeduniqueresult.csv")
D <- read.csv("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GeneExpressionAnalysis/Microarray/TopGenes_2016-02-15/FTLDrankeduniqueresult.csv")
E <- read.csv("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GeneExpressionAnalysis/Microarray/TopGenes_2016-02-15/VCPrankeduniqueresult.csv")
F <- read.csv("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GeneExpressionAnalysis/RNA-seq/16_3_17/PET16_3_17rankeduniqueresult.csv")
G <- read.csv("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/GeneExpressionAnalysis/RNA-seq/16_3_17/RAV16_3_17rankeduniqueresult.csv")

# ## Save annotation file locations to variable
# annotation.U133plus2<-"/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/HG-U133_Plus_2.na35.annot.csv/HG-U133_Plus_2.na35.annot.txt"
# annotation.U133A2<-"/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/HG-U133A_2.na35.annot.csv/HG-U133A_2.na35.annot.txt"
# 
# # Read in annotation files
# annotationU133plus2<-read.table(annotation.U133plus2, header = TRUE, row.names=NULL, sep="\t", skip=0, stringsAsFactors=F, quote = "", comment.char="!", fill = TRUE )
# annotationU133A2<-read.table(annotation.U133A2, header = TRUE, row.names=NULL, sep="\t", skip=0, stringsAsFactors=F, quote = "", comment.char="!", fill = TRUE )
# 
# # Remove rows in which gene symbol is absent
# annotationU133plus2<-subset(annotationU133plus2, subset=(Gene.Symbol !="---")) #if no gene symbol, discount
# annotationU133A2<-subset(annotationU133A2, subset=(Gene.Symbol !="---")) #if no gene symbol, discount
# 
# # Remove rows in which gene symbol is duplicated
# annotationU133plus2 <- annotationU133plus2[!duplicated(annotationU133plus2[,15]),]
# annotationU133A2 <- annotationU133A2[!duplicated(annotationU133A2[,15]),]
# 
# 
# # Remove rows in which genes are noted to have negative strand matching probes
# idxNegativeStrand<-grep("Negative Strand Matching Probes", annotationU133plus2$Annotation.Notes)
# if(length(idxNegativeStrand)>0)
# {
#   annotationU133plus2<-annotationU133plus2[-idxNegativeStrand,]
# }
# nrow(annotationU133plus2)
# 
# 
# 
# idxNegativeStrand<-grep("Negative Strand Matching Probes", annotationU133A2$Annotation.Notes)
# if(length(idxNegativeStrand)>0)
# {
#   annotationU133A2<-annotationU133A2[-idxNegativeStrand,]
# }
# nrow(annotationU133A2)


#indicate the number of overlapping genes identified by DE analysis
test <- 3

m=10000 #number of repetitions 
r <- c(1:m) #store repetition numbers in vector "r"

for (j in 1:m)
{
  random1 <- sample (A$Gene.Symbol, size=5000, replace=FALSE)
  random2 <- sample (B$Gene.Symbol, size=5000, replace=FALSE)
  random3 <- sample (C$Gene.Symbol, size=5000, replace=FALSE)
  random4 <- sample (D$Gene.Symbol, size=5000, replace=FALSE)
  random5 <- sample (E$Gene.Symbol, size=5000, replace=FALSE)
  random6 <- sample (F$X,           size=5000, replace=FALSE)
  random7 <- sample (G$X,           size=5000, replace=FALSE)
  random <- Reduce(intersect, list(random1, random2, random3, random4, random5, random6, random7))
  r[j] <- length(random)
}

test1 <- which(r > test)  # count number of times r is larger than test value
result <- (length(test1)/m) # calculate P value
mean(r)
