setwd("/Users/clairegreen/Desktop/")

source("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/JiantaoGRAIL/GRAIL.R")

load("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/JiantaoGRAIL/text_2006_12_Symbol.RData")
load("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/JiantaoGRAIL/Hsa_v12_08_CoExpression_Symbol.RData")

# read seed genes
sTable = read.table(file = "KnownTDP43.txt", sep = "\t") #table with sig.DE gene output (fc, pval etc)
SeedUP = as.character(sTable[,"Symbol"])[as.numeric(sTable$logFC) > 0] #upregulated genes
SeedDN = as.character(sTable[,"Symbol"])[as.numeric(sTable$logFC) < 0] #downregulated genes
SeedUPDN = c(SeedUP, SeedDN) #gene list ordered by foldchange
SList  = list(SeedUP = unique(SeedUP), SeedDN = unique(SeedDN), SeedUPDN = unique(SeedUPDN)) #take unique values of each

# query genes
queryVector = as.character(read.table(file = "QueryTDP43.txt", row.names = NULL, header = F, sep = "\t")[,1])

# filter by Whole network

for(i in 1:length(sTable)){

	SeedVector = as.character(sTable[[i]])
	SeedList   = as.list(SeedVector)
	names(SeedList) = SeedVector

	tmList = GrailWorkflow(SeedList, queryVector, TmDB,   method = "coexpression", Nthreads = 4, Pfcutoff = 0.1, minN = 20, ABS = FALSE)
	coList = GrailWorkflow(SeedList, queryVector, CoexDB, method = "coexpression", Nthreads = 4, Pfcutoff = 0.1, minN = 20, ABS = FALSE)

	# merge results
	tmTable = tmList$qTable
	coTable = coList$qTable
	rownames(tmTable) = as.character(tmTable[,1])
	rownames(coTable) = as.character(coTable[,1])

	ID = unique(c(SeedVector, queryVector))
	NN = length(ID)
	tmBestGene = coBestGene = rep("", NN)
	tmBestPv   = coBestPv   = rep(NA,  NN)
	names(tmBestGene) = names(tmBestPv) = names(coBestGene) = names(coBestPv) = ID

	tmBestGene[ rownames(tmTable)] = as.character(tmTable[,"queryVector"])
	tmBestPv[ rownames(tmTable)]   = as.numeric(tmTable[,"pvalueVector"])

	coBestGene[ rownames(coTable)] = as.character(coTable[,"queryVector"])
	coBestPv[ rownames(coTable)]   = as.numeric(coTable[,"pvalueVector"])
	
	Flag   <- ID %in% queryVector
	mTable <- data.frame(ID, tmBestGene, tmBestPv, coBestGene, coBestPv, Flag)	

	write.table(mTable, file = paste0( names(sTable)[i], ".txt"), row.names = F, col.names = T, sep = "\t")
}
