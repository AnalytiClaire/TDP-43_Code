#QQplot for p Values#

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/DEG Test2/")
C9DEgene <- read.csv("C9rankeduniqueresult.csv")
CHDEgene <- read.csv("CHrankeduniqueresult.csv")
sALSDEgene <- read.csv("sALSrankeduniqueresult.csv")
FTLDDEgene <- read.csv("FTLDrankeduniqueresult.csv")
VCPDEgene <- read.csv("VCPrankeduniqueresult.csv")

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/DEG Test2/")
PETDEgene <- read.csv("PETrankeduniqueresult.csv")
RAVDEgene <- read.csv("RAVrankeduniqueresult.csv")

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/HG-U133_Plus_2.na35.annot.csv/")
Annotation <- read.csv("HG-U133_Plus_2.na35.annot.csv")

#Create list for C9orf72#
genelist <- C9DEgene$Ensembl
thresh <- c(500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000, 7500, 8000)
resultC9 <- list(a = 1:500, b = 1:1000, c = 1:1500, d = 1:2000, e = 1:2500, f = 1:3000, g = 1:3500, 
               h = 1:4000, i = 1:4500, j = 1:5000, k = 1:5500, l = 1:6000, m = 1:6500, n = 1:7000, 
               o = 1:7500, p = 1:8000)

for (i in 1:length(thresh)) {
  resultC9[[i]] <- genelist[1:thresh[[i]]]
}

#Create list for CHMP2B#
genelist <- CHDEgene$Ensembl
thresh <- c(500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000, 7500, 8000)
resultCH <- list(a = 1:500, b = 1:1000, c = 1:1500, d = 1:2000, e = 1:2500, f = 1:3000, g = 1:3500, 
                h = 1:4000, i = 1:4500, j = 1:5000, k = 1:5500, l = 1:6000, m = 1:6500, n = 1:7000, 
                o = 1:7500, p = 1:8000)

for (i in 1:length(thresh)) {
  resultCH[[i]] <- genelist[1:thresh[[i]]]
}

#Create list for sALS#
genelist <- sALSDEgene$Ensembl
thresh <- c(500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000, 7500, 8000)
resultsALS <- list(a = 1:500, b = 1:1000, c = 1:1500, d = 1:2000, e = 1:2500, f = 1:3000, g = 1:3500, 
                h = 1:4000, i = 1:4500, j = 1:5000, k = 1:5500, l = 1:6000, m = 1:6500, n = 1:7000, 
                o = 1:7500, p = 1:8000)

for (i in 1:length(thresh)) {
  resultsALS[[i]] <- genelist[1:thresh[[i]]]
}

#Create list for FTLD#
genelist <- FTLDDEgene$Ensembl
thresh <- c(500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000, 7500, 8000)
resultFTLD <- list(a = 1:500, b = 1:1000, c = 1:1500, d = 1:2000, e = 1:2500, f = 1:3000, g = 1:3500, 
                h = 1:4000, i = 1:4500, j = 1:5000, k = 1:5500, l = 1:6000, m = 1:6500, n = 1:7000, 
                o = 1:7500, p = 1:8000)

for (i in 1:length(thresh)) {
  resultFTLD[[i]] <- genelist[1:thresh[[i]]]
}

#Create list for VCP#
genelist <- VCPDEgene$Ensembl
thresh <- c(500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000, 7500, 8000)
resultVCP <- list(a = 1:500, b = 1:1000, c = 1:1500, d = 1:2000, e = 1:2500, f = 1:3000, g = 1:3500, 
                h = 1:4000, i = 1:4500, j = 1:5000, k = 1:5500, l = 1:6000, m = 1:6500, n = 1:7000, 
                o = 1:7500, p = 1:8000)

for (i in 1:length(thresh)) {
  resultVCP[[i]] <- genelist[1:thresh[[i]]]
}

#Create list for Petrucelli#
genelist <- PETDEgene$Row.names
thresh <- c(500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000, 7500, 8000)
resultPET <- list(a = 1:500, b = 1:1000, c = 1:1500, d = 1:2000, e = 1:2500, f = 1:3000, g = 1:3500, 
                h = 1:4000, i = 1:4500, j = 1:5000, k = 1:5500, l = 1:6000, m = 1:6500, n = 1:7000, 
                o = 1:7500, p = 1:8000)

for (i in 1:length(thresh)) {
  resultPET[[i]] <- genelist[1:thresh[[i]]]
}

#Create list for Ravits#
genelist <- RAVDEgene$Row.names.1
thresh <- c(500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000, 7500, 8000)
resultRAV <- list(a = 1:500, b = 1:1000, c = 1:1500, d = 1:2000, e = 1:2500, f = 1:3000, g = 1:3500, 
                h = 1:4000, i = 1:4500, j = 1:5000, k = 1:5500, l = 1:6000, m = 1:6500, n = 1:7000, 
                o = 1:7500, p = 1:8000)

for (i in 1:length(thresh)) {
  resultRAV[[i]] <- genelist[1:thresh[[i]]]
}

##Intersect lists of same length##
result <- list(list(a = 1:5000, b = 1:5000, c = 1:5000, d = 1:5000, e = 1:5000, f = 1:5000, g = 1:5000, 
                    h = 1:5000, i = 1:5000, j = 1:5000, k = 1:5000, l = 1:5000, m = 1:5000, n = 1:5000, 
                    o = 1:5000, p = 1:5000))

for (i in 1:16) {
  result[[i]] <- Reduce(intersect, list(resultC9[[i]],resultCH[[i]],resultsALS[[i]],resultFTLD[[i]],
                                        resultVCP[[i]],resultPET[[i]],resultRAV[[i]]))
}

##Make list values equal, convert to data frame and save##
for (i in 1:16) {
  length(result[[i]]) <- 665
}

setwd(dir = "/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/DEG Test2/")
resultdf <- as.data.frame(result)
write.csv(result, file = "Ensconsensus.csv")


# result2000 <- as.data.frame(result[[4]])
# result2500 <- as.data.frame(result[[5]])
result3000 <- as.data.frame(result[[6]])
result3500 <- as.data.frame(result[[7]])
result4000 <- as.data.frame(result[[8]])
result4500 <- as.data.frame(result[[9]])
result5000 <- as.data.frame(result[[10]])
result5500 <- as.data.frame(result[[11]])
result6000 <- as.data.frame(result[[12]])
result6500 <- as.data.frame(result[[13]])
result7000 <- as.data.frame(result[[14]])
result7500 <- as.data.frame(result[[15]])
result8000 <- as.data.frame(result[[16]])

Z <- read.csv(file = "Allgenes.csv", header = TRUE, na.strings = c("", "NA)"))
Z <- as.list(Z)
Z<- lapply(Z, function(x) x[!is.na(x)])
Z <- lapply(Z, function(x) x[!duplicated(x)])
