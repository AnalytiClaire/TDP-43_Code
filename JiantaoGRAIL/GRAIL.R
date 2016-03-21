GeneSignificance <- function(Nseed, Cg, Pfcutoff = 0.1) {
# Nseed is the number of seed region used, integer
# Cg is the weight for gene g
	
	if(Nseed < 5)
		stop("Less than 5 regions found.\n ")
  
	Pg <- 0
	lambda <- Nseed*Pfcutoff

	for(i in 0:Nseed)
		Pg <- Pg + dpois(i, lambda)* pgamma(Cg, shape = i, rate = 1, lower.tail = FALSE, log.p = FALSE)

	return(Pg)
}

GeneWeightByRelatedness <- function(relatednessVector, GeneID, seedList, Pfcutoff = 0.1, ABS = FALSE) {
# relatednessVector is the relatedness vector, the names are entrez gene ID
# GeneID is the entrez gene needs to be assessed
# seedList is genes in each seed region, returned from geneFiler. 
# seedList must be filtered to make sure it only contains the genes in database

  sumNa <- sum(is.na(relatednessVector))
  if(sumNa > 100)
    stop("Too many missing values in the relatednessVector.\n")
  relatednessVector[is.na(relatednessVector)] <- 0
   
	# if abosolute value
	if(ABS)
		relatednessVector <- abs(relatednessVector)
  
	relatednessVector   <- sort(relatednessVector, decreasing = TRUE)
	
	# get the ranked gene list
	gID <- names(relatednessVector)
	N   <- length(gID)
	Ng  <- 1: N

	# process region by region
	M <- length(seedList)
	GeneOverlap <- rep(TRUE, M)
	adjustP     <- rep(NaN,  M)
	
	for(i in 1:M) {
		
		g <- seedList[[i]]
		if( GeneID %in% g)
			next
		GeneOverlap[i] <- FALSE
		
		# check the seed genes in the rank
		Index <- gID %in% g
		unAdjustP  <- Ng[Index][1]/N
		adjustP[i] <- 1- (1- unAdjustP)^sum(Index)
	}

	# calculate the weight
	adjustP <- adjustP[ ! GeneOverlap]
	Nseed   <- length(adjustP)
	Cg <- 0
	
  for(i in 1:Nseed) {
    
		if( adjustP[i] > Pfcutoff)
			next
		Cg <- Cg - log(adjustP[i]/Pfcutoff)
	}
  
	# return the value
	rList <- list(Nseed = Nseed, Cg = Cg)
	return(rList)
}

GetComputedRelatednessVector <- function(GeneID, MX){
# GeneID is the ID for query Gene
# M is a pre-computed relatedness matrix
  
  if( ! GeneID %in% colnames(MX))
    stop(paste(GeneID,"is not found in matrix M.\n"))
  
  rValue <- MX[,GeneID]
  names(rValue) <- rownames(MX)
  
  return(rValue)
}

GetCorrelationFromExpression <- function(GeneID, MX, minN = 20){
# GeneID is the ID for query Gene
# M is a gene expression matrix with genes in the row and samples in the column
  
  N     <- nrow(MX)
  Index <- rownames(MX) %in% GeneID
  
  if(sum(Index) != 1)
    stop(paste(sum(Index),"genes found for gene", GeneID,"\n"))
  
  TagEx    <- MX[GeneID,]
  TagIndex <- ! is.na(TagEx)
  
  # test all other genes
  rValue <- rep(0,N)
  for(i in 1:N){
    
    TestEx <- MX[i,]
    TestIndex <- ! is.na(TestEx)
    
    PairIndex <- TagIndex & TestIndex
    
    if(sum(PairIndex) < minN)
      next
    
    rValue[i] <- cor(TagEx[PairIndex], TestEx[PairIndex])
  }
  names(rValue) <- rownames(MX)
  
  return(rValue)
}

GeneFilter <- function(seedList, PG){
  
  rList <- list()
  M <- length(seedList)
  
  for(i in 1:M)
    rList[[i]] <- intersect(seedList[[i]], PG)
  names(rList) <- names(seedList)
  
  flag <- sapply(rList,function(x) length(x)) > 0
  rList<- rList[flag]
  
  return(rList)
}

ChooseBestGene <- function(seedList, queryVector, pvalueVector){
  
  if(length(queryVector) != length(pvalueVector))
    stop("queryVector and pvalueVector must be of the same length.\n")
  
  M <- length(seedList)
  chosenG <- rep("", M)
  chosenP <- rep(NaN,M)
  
  for(i in 1:M){
    
    temp <- seedList[[i]]
    if(length(temp) == 0)
      next
    
    Index <- queryVector %in% temp
    if(sum(Index) == 0)
      next
    
    tempG <- queryVector[Index]
    tempP <- pvalueVector[Index]
    
    chosenG[i] <- tempG[tempP == min(tempP)][1]
    chosenP[i] <- tempP[tempP == min(tempP)][1]
  }
  
  rTable <- data.frame(chosenG,chosenP)
  rownames(rTable) <- names(seedList)
  
  return(rTable)
}

GrailWorkflow <- function(seedList, queryVector, MX, method = "expression", Nthreads = 4, Pfcutoff = 0.1, minN = 20, ABS = FALSE){
# element in seedList, queryVector and row names of M must be of the same type, for example gene symbol/entrez ID
  
  # check required packges
  check1 <- ! require(foreach)
  check2 <- ! require(doMC)
  if(check1 | check2)
    stop("foreach and doMC are required for multi-threading.\n")
  registerDoMC(Nthreads)
  
  # clearn the seedList
  PG <- rownames(MX)
  seedList <- GeneFilter(seedList, PG)
  
  # clean queryVector
  seedVector <- unlist(seedList)
  queryVector<- intersect( c(seedVector,queryVector), PG )
  N <- length(queryVector)
  if(N < 5)
    stop("Nothing to query.\n")
  
  # process eqch query gene
  if(! method %in% c("expression","coexpression"))
    stop("Method must be either expression or coexpression.\n")
  
  pb<-txtProgressBar(min = 0, max = N, style = 3)
  xx <- foreach(i = 1:N) %dopar% {
    
    GeneID <- queryVector[i]
    
    if(method == "expression")
      relatednessVector <- GetCorrelationFromExpression(GeneID, MX, minN)
    else
      relatednessVector <- GetComputedRelatednessVector(GeneID, MX)
    
    rList <- GeneWeightByRelatedness(relatednessVector, GeneID, seedList, Pfcutoff, ABS)
    
    pvalue<- GeneSignificance(rList$Nseed, rList$Cg, Pfcutoff)
    
    # update progress
    setTxtProgressBar(pb, i)
    
    # return
    pvalue
  }
  cat("\n")
  
  pvalueVector <- unlist(xx)
  
  # get best gene
  bTable <- ChooseBestGene(seedList, queryVector, pvalueVector)
  rList  <- list(seedList = seedList, bTable = bTable, qTable = data.frame(queryVector,pvalueVector))
  
  return(rList)
}

