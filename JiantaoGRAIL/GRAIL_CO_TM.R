setwd("/n/hsphS10/hsphfs1/hide_lab/jshi/Development/GRAIL/WM")

source("/n/hsphS10/hsphfs1/hide_lab/jshi/Script/GRAIL.R")

load("/n/hsphS10/hsphfs1/hide_lab/jshi/software/GRAIL_ranking_v3/NetworkDB/TextMing/text_2006_12_Symbol.RData")
load("/n/hsphS10/hsphfs1/hide_lab/jshi/software/GRAIL_ranking_v3/NetworkDB/CoExpression/Hsa_v12_08_CoExpression_Symbol.RData")

# read seed genes
sTable = read.table(file = "RAW/WM_BL_vs_HD_BL.txt", row.names = 1, header = T, sep = "\t")
SeedUP = as.character(sTable[,"Symbol"])[as.numeric(sTable$logFC) > 0]
SeedDN = as.character(sTable[,"Symbol"])[as.numeric(sTable$logFC) < 0]
SeedUPDN = c(SeedUP, SeedDN)
SList  = list(SeedUP = unique(SeedUP), SeedDN = unique(SeedDN), SeedUPDN = unique(SeedUPDN))

# query genes
queryVector = as.character(read.table(file = "RAW/Query_Genes_New.txt", row.names = NULL, header = F, sep = "\t")[,1])

# filter by Whole network

for(i in 1:length(SList)){

	SeedVector = as.character(SList[[i]])
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

	write.table(mTable, file = paste0("Results/", names(SList)[i], ".txt"), row.names = F, col.names = T, sep = "\t")
}
