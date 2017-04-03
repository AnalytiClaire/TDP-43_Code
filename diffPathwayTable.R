diffPathwayTable <- function (fingerprints, fac, threshold) 
{
  fac <- as.factor(fac)
  levels <- levels(fac)
  if (length(levels) != 2) 
    stop("fac must contain 2 levels")
  if (length(fac) != ncol(fingerprints)) 
    stop("fac length must equal nubmer of fingerprints")
  if ((threshold < 0 | threshold > 2)) 
    stop("threshold value should be between 0 and 2")
  if (sum(fac == levels[1]) == 1) {
    print("N.B. only 1 sample for group 1")
    mean.1 <- fingerprints[, fac == levels[1]]
  }
  else {
    mean.1 <- apply(fingerprints[, fac == levels[1]], 1, 
                    mean)
  }
  if (sum(fac == levels[2]) == 1) {
    print("N.B. only 1 sample for group 2")
    mean.2 <- fingerprints[, fac == levels[2]]
  }
  else {
    mean.2 <- apply(fingerprints[, fac == levels[2]], 1, 
                    mean)
  }
  mean.diff <- mean.1 - mean.2
  meandiff <- abs(mean.diff)
  signifPathawys <- na.omit(names(mean.diff)[abs(mean.diff) > 
                                               threshold])
  return(as.data.frame(meandiff))
  return(as.vector(signifPathawys))
}


c9.lcm <- diffPathwayTable(C9.LCM_pathprint, vec.c9, thres)
CHMP2B.lcm <- diffPathwayTable(CHMP2B.LCM_pathprint, vec.ch, thres)
SALS.lcm <- diffPathwayTable(SALS.LCM_pathprint, vec.sals, thres)
FTLD_FCx <- diffPathwayTable(FTLD_pathprint, vec.FTLD, thres)
VCP.m <- diffPathwayTable(VCP_pathprint, vec.vcp, thres)

# setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/Pathprint/PathprintThreshold")
# write.csv(c9.lcm, "C9.csv")
# write.csv(SALS.lcm, "sALS.csv")
# write.csv(CHMP2B.lcm, "CHMP2B.csv")
# write.csv(FTLD_FCx, "FTLD.csv")
# write.csv(VCP.m, "VCP.csv")


library(data.table)
setDT(c9.lcm, keep.rownames = TRUE)
setDT(CHMP2B.lcm, keep.rownames = TRUE)
setDT(SALS.lcm, keep.rownames = TRUE)
setDT(FTLD_FCx, keep.rownames = TRUE)
setDT(VCP.m, keep.rownames = TRUE)

c9.lcm <- c9.lcm[order(c9.lcm[,2], decreasing = TRUE),]
CHMP2B.lcm <- CHMP2B.lcm[order(CHMP2B.lcm$meandiff, decreasing = TRUE),]
SALS.lcm <- SALS.lcm[order(SALS.lcm$meandiff, decreasing = TRUE),]
FTLD_FCx <- FTLD_FCx[order(FTLD_FCx$meandiff, decreasing = TRUE),]
VCP.m <- VCP.m[order(VCP.m$meandiff, decreasing = TRUE),]

thresh <- 200

c9thresh <- c9.lcm[1:thresh,]
chthresh <- CHMP2B.lcm[1:thresh,]
salsthresh <- SALS.lcm[1:thresh,]
ftldthresh <- FTLD_FCx[1:thresh,]
vcpthresh <- VCP.m[1:thresh,]



overlap <- Reduce(intersect, list(c9thresh$rn, chthresh$rn, salsthresh$rn, ftldthresh$rn, vcpthresh$rn)) #selects pathways that are present in all data sets listed
overlap

setwd("/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/Pathprint/PathprintThreshold")
write.table(overlap, "Pathprint_thresholdrevision.txt", quote = F, row.names = F, col.names = F)

# 
# AllResult <- merge(c9.lcm, CHMP2B.lcm, by = "row.names")
# AllResult <- merge(AllResult, SALS.lcm, by.x = "Row.names", by.y = "row.names")
# AllResult <- merge(AllResult, FTLD_FCx, by.x = "Row.names", by.y = "row.names")
# AllResult <- merge(AllResult, VCP.m, by.x = "Row.names", by.y = "row.names")
# 
# # rownames(AllResult) <- AllResult[,1]
# AllResult[,1] <- NULL
# colnames(AllResult) <- c("C9orf72", "CHMP2B", "sALS", "FTLD", "VCP")
# 
# 
# AllResult <- AllResult[order(AllResult$C9orf72, decreasing = TRUE),]
