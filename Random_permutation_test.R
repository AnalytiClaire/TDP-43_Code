###Random permutations for Differential Expression##

# #load the total number of unique genes for each data set
# load("/Users/clairegreen/Documents/PhD/TDP-43 Sig/TDP-43 Data Sets/DE Genes/C9orf72/DE_Workspace.RData")
# unique_C9 <- uniqueresult
# load("/Users/clairegreen/Documents/PhD/TDP-43 Sig/TDP-43 Data Sets/DE Genes/CHMP2B/DE_Workspace.RData")
# unique_CHMP2B <- uniqueresult
# load("/Users/clairegreen/Documents/PhD/TDP-43 Sig/TDP-43 Data Sets/DE Genes/sALS/DE_Workspace.RData")
# unique_sALS <- uniqueresult
# load("/Users/clairegreen/Documents/PhD/TDP-43 Sig/TDP-43 Data Sets/DE Genes/FTLD/DE_Workspace.RData")
# unique_FTLD <- uniqueresult
# load("/Users/clairegreen/Documents/PhD/TDP-43 Sig/TDP-43 Data Sets/DE Genes/VCP/DE_Workspace.RData")
# unique_VCP <- uniqueresult

annotation.file<-"/Users/clairegreen/Documents/PhD/TDP-43 Sig/TDP-43 Data Sets/HG-U133_Plus_2.na35.annot.csv/HG-U133_Plus_2.na35.annot.txt"
annotation<-read.table(annotation.file, header = TRUE, row.names=NULL, sep="\t", skip=0, stringsAsFactors=F, quote = "", comment.char="!", fill = TRUE )
annotation<-subset( annotation, subset=(Gene.Symbol !="---")) #if no gene symbol, discount

#indicate the number of overlapping genes identified by DE analysis
test <- 181

m=10000 #number of repetitions 
r <- c(1:m) #store repetition numbers in vector "r"

for (j in 1:m)
{
  random1 <- sample (annotation$Gene.Symbol, size=5000, replace=F)
  random2 <- sample (annotation$Gene.Symbol, size=5000, replace=F)
  random3 <- sample (annotation$Gene.Symbol, size=5000, replace=F)
  random4 <- sample (annotation$Gene.Symbol, size=5000, replace=F)
  random5 <- sample (annotation$Gene.Symbol, size=5000, replace=F)
  random <- Reduce(intersect, list(random1, random2, random3, random4, random5))
  r[j] <- length(random)
}

test1 <- which(r > test) 
result <- (length(test1)/m)

