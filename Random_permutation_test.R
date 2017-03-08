###Random permutations for Differential Expression##

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/DEG_Test2/")

#load the total number of unique genes for each data set
Group1 <- read.csv("C9rankeduniqueresult.csv")
Group2 <- read.csv("CHrankeduniqueresult.csv")
Group3 <- read.csv("sALSrankeduniqueresult.csv")
Group4 <- read.csv("FTLDrankeduniqueresult.csv")
Group5 <- read.csv("VCPrankeduniqueresult.csv")
Group6 <- read.csv("PETNOrankeduniqueresult.csv")
Group7 <- read.csv("RAVrankeduniqueresult.csv")

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
test <- 178
samplenum <- 6500

m=10000 #number of repetitions 
r <- c(1:m) #store repetition numbers in vector "r"

for (j in 1:m)
{
  random1 <- sample (Group1$Gene.Symbol, size=samplenum, replace=FALSE)
  random2 <- sample (Group2$Gene.Symbol, size=samplenum, replace=FALSE)
  random3 <- sample (Group3$Gene.Symbol, size=samplenum, replace=FALSE)
  random4 <- sample (Group4$Gene.Symbol, size=samplenum, replace=FALSE)
  random5 <- sample (Group5$Gene.Symbol, size=samplenum, replace=FALSE)
  random6 <- sample (Group6$symbol, size=samplenum, replace=FALSE)
  random7 <- sample (Group7$hgnc_symbol, size=samplenum, replace=FALSE)
  random <- Reduce(intersect, list(random1, random2, random3, random4, random5, random6, random7))
  r[j] <- length(random)
}

test1 <- which(r > test)  # count number of times r is larger than test value
result <- (length(test1)/m) # calculate P value
mean(r)
